#!/usr/bin/env python3
"""Load genomes, phenotypes, merged features and resistance calls into SQLite.

Creates the database from database/schema.sql, then loads:

    genomes, phenotypes   from data/samples.tsv and data/genomes_manifest.tsv
    sequences             from the FASTA files (accession and length)
    features              from the merged GFF3 files (merge_annotations.py)
    resistance_calls      from the portal (data/portal_calls.tsv, with its
                          per-antibiotic repeat rows collapsed), the AMRFinderPlus
                          and Kleborate reports, and the split-gene screen
    tool_versions         one row per call source

Loading is all-or-nothing: one transaction, so a failed file leaves no
half-loaded database. The schema's constraints (contig names that exist in the
FASTA, start >= 1, start <= end, a valid strand, unique IDs) reject bad rows.

Usage:
    python scripts/load_database.py --fasta data/genomes/*.fasta \\
        --gff3 results/merged/*.klebtrace.gff3 \\
        --amrfinderplus results/phase2_galaxy/amrfinderplus_report/*.tsv \\
        --kleborate results/kleborate/general__species_output.txt \\
                    results/kleborate_negatives/general__species_output.txt \\
        --out results/klebtrace.sqlite
"""

from __future__ import annotations

import argparse
import re
import sqlite3
import sys
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

from merge_annotations import read_fasta_lengths

ANTIBIOTIC = "meropenem"
CARBAPENEMASE_REFERENCE = "carbapenem-hydrolyzing"

TOOL_VERSIONS = [
    ("portal", "AMRFinderPlus (AMR Portal pipeline, protein input)",
     "portal annotation_tool_version 1.5.0 / 1.6.0", "not stated in release 2026-07", "AMR Portal"),
    ("amrfinderplus", "AMRFinderPlus (nucleotide input)", "4.2.7", "2026-05-15.1",
     "usegalaxy.org or the KlebTrace Nextflow pipeline (same version and database)"),
    ("kleborate", "Kleborate", "3.3.0b1", "Kleborate AMR database v3.3", "local (Bioconda)"),
    ("split_gene_screen", "screen_carbapenemases.py with minimap2", "minimap2 2.31",
     "Kleborate AMR database v3.3 (568 carbapenemase alleles)", "local"),
]


def create(db: Path, schema: Path) -> sqlite3.Connection:
    """Create a fresh database from the schema (an existing file is replaced)."""
    if db.exists():
        db.unlink()
    db.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db, isolation_level=None)
    con.executescript(schema.read_text())
    con.execute("BEGIN")
    return con


def _none(value):
    return None if pd.isna(value) or value in ("", "-", "NA") else value


def _flag(value) -> int | None:
    if pd.isna(value):
        return None
    return int(str(value).strip().lower() in ("true", "1", "yes"))


def is_carbapenemase(subtype, subclass, reference_name) -> int:
    """Acquired gene (not a point mutation) labelled CARBAPENEM, or an unnamed allele whose
    closest reference is carbapenem-hydrolyzing (the split OXA-232 in GCA_019928365.1)."""
    if str(subtype).upper().startswith("POINT"):
        return 0
    return int("CARBAPENEM" in str(subclass).upper()
               or CARBAPENEMASE_REFERENCE in str(reference_name).lower())


def load_genomes(con: sqlite3.Connection, samples: Path, manifest: Path) -> int:
    s = pd.read_csv(samples, sep="\t", dtype=str)
    m = pd.read_csv(manifest, sep="\t", dtype=str).set_index("assembly")
    missing = sorted(set(s["assembly"]) - set(m.index))
    if missing:
        raise ValueError(f"Genomes in {samples} but not in {manifest}: {missing}")
    for r in s.to_dict("records"):
        g = m.loc[r["assembly"]]
        con.execute(
            "INSERT INTO genomes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (r["assembly"], r["biosample"], _none(r.get("country")),
             int(float(r["collection_year"])) if _none(r.get("collection_year")) else None,
             _none(r.get("region")), r["selection_group"], _none(r.get("carbapenemase_family")),
             _flag(r.get("in_ensembl")), int(g["n_sequences"]), int(g["total_bp"]), int(g["n50"]),
             float(g["gc_percent"]), g["md5"]))
        con.execute("INSERT INTO phenotypes VALUES (?, ?, ?, ?, ?)",
                    (r["assembly"], ANTIBIOTIC, r["phenotype"], _none(r.get("mic_mg_l")),
                     _flag(r.get("strict_resistant"))))
    return len(s)


def load_sequences(con: sqlite3.Connection, fastas: list[Path]) -> int:
    known = {row[0] for row in con.execute("SELECT assembly FROM genomes")}
    n = 0
    for fasta in fastas:
        assembly = fasta.name.removesuffix(".fasta")
        if assembly not in known:
            continue
        rows = [(assembly, seqid, length) for seqid, length in read_fasta_lengths(fasta).items()]
        con.executemany("INSERT INTO sequences VALUES (?, ?, ?)", rows)
        n += len(rows)
    return n


def parse_attributes(column9: str) -> dict[str, str]:
    pairs = (item.split("=", 1) for item in column9.strip().split(";") if "=" in item)
    return {k: unquote(v) for k, v in pairs}


def load_features(con: sqlite3.Connection, gff3s: list[Path]) -> int:
    n = 0
    for path in gff3s:
        assembly = path.name.removesuffix(".klebtrace.gff3")
        rows = []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) != 9:
                    raise ValueError(f"{path}: expected 9 columns: {line[:80]}")
                a = parse_attributes(f[8])
                rows.append((assembly, f[0], f[1], f[2], int(f[3]), int(f[4]), f[6],
                             a.get("ID"), a.get("Name"), a.get("product"), f[8]))
        con.executemany(
            "INSERT INTO features (assembly, seqid, source, type, start, end, strand, gff_id, "
            "name, product, attributes) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", rows)
        n += len(rows)
    return n


def _insert_calls(con: sqlite3.Connection, rows: list[tuple]) -> int:
    con.executemany(
        "INSERT INTO resistance_calls (assembly, source, gene, element_subtype, subclass, "
        "is_carbapenemase, method, is_partial, seqid, start, end, coverage_pct, identity_pct, "
        "closest_reference) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", rows)
    return len(rows)


def _num(value, kind=float):
    value = _none(value)
    return None if value is None else kind(float(value))


def load_portal_calls(con: sqlite3.Connection, portal_calls: Path) -> int:
    df = pd.read_csv(portal_calls, sep="\t", dtype=str)
    # One row per antibiotic the gene affects (split_subclass): keep one per gene and position
    df = df.drop_duplicates(["assembly_ID", "amr_element_symbol", "region", "region_start", "region_end"])
    rows = [(r["assembly_ID"], "portal", r["amr_element_symbol"], r["element_subtype"], r["subclass"],
             is_carbapenemase(r["element_subtype"], r["subclass"], r["reference_name"]),
             r["amrfinderplus_method"], int(str(r["amrfinderplus_method"]).startswith("PARTIAL")),
             r["region"], _num(r["region_start"], int), _num(r["region_end"], int),
             _num(r["reference_sequence_coverage"]), _num(r["reference_sequence_identity"]),
             _none(r["reference_name"]))
            for r in df.to_dict("records")]
    return _insert_calls(con, rows)


def load_amrfinderplus(con: sqlite3.Connection, reports: list[Path]) -> int:
    known = {row[0] for row in con.execute("SELECT assembly FROM genomes")}
    n = 0
    for path in reports:
        assembly = re.sub(r"(\.amrfinderplus)?\.tsv$", "", path.name)
        if assembly not in known:
            continue
        df = pd.read_csv(path, sep="\t", dtype=str)
        rows = [(assembly, "amrfinderplus", r["Element symbol"], r["Subtype"], r["Subclass"],
                 is_carbapenemase(r["Subtype"], r["Subclass"], r.get("Closest reference name")),
                 r["Method"], int(str(r["Method"]).startswith("PARTIAL")),
                 str(r["Contig id"]).split()[0], int(r["Start"]), int(r["Stop"]),
                 _num(r["% Coverage of reference"]), _num(r["% Identity to reference"]),
                 _none(r.get("Closest reference name")))
                for r in df.to_dict("records") if r["Type"] == "AMR"]
        n += _insert_calls(con, rows)
    return n


def load_kleborate(con: sqlite3.Connection, reports: list[Path]) -> int:
    """Carbapenemases and porin changes from Kleborate (no coordinates).

    Kleborate marks an incomplete gene with "?" and an inexact match with "*";
    "?" is loaded as partial.
    """
    known = {row[0] for row in con.execute("SELECT assembly FROM genomes")}
    rows = []
    for path in reports:
        df = pd.read_csv(path, sep="\t", dtype=str)
        for r in df.to_dict("records"):
            if r["strain"] not in known:
                continue
            for gene in str(_none(r.get("kpsc__amr__Bla_Carb_acquired")) or "").split(";"):
                if gene:
                    rows.append((r["strain"], "kleborate", gene, "AMR", "CARBAPENEM", 1, None,
                                 int("?" in gene), None, None, None, None, None, None))
            for change in str(_none(r.get("kpsc__amr__Omp_mutations")) or "").split(";"):
                if change:
                    rows.append((r["strain"], "kleborate", change, "POINT", "CARBAPENEM", 0, None,
                                 0, None, None, None, None, None, None))
    return _insert_calls(con, rows)


def load_screen(con: sqlite3.Connection, screen: Path) -> int:
    known = {row[0] for row in con.execute("SELECT assembly FROM genomes")}
    df = pd.read_csv(screen, sep="\t", dtype=str)
    rows = [(r["assembly"], "split_gene_screen", f"bla{r['allele']}", "AMR", "CARBAPENEM", 1,
             r["call"], int(r["call"] in ("split", "partial", "interrupted")),
             _none(r["contigs"]), None, None, round(100 * float(r["coverage"]), 1), None, None)
            for r in df.to_dict("records") if r["assembly"] in known and r["call"] != "none"]
    return _insert_calls(con, rows)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--schema", type=Path, default=Path("database/schema.sql"))
    p.add_argument("--samples", type=Path, default=Path("data/samples.tsv"))
    p.add_argument("--manifest", type=Path, default=Path("data/genomes_manifest.tsv"))
    p.add_argument("--portal-calls", type=Path, default=Path("data/portal_calls.tsv"))
    p.add_argument("--screen", type=Path, default=Path("results/portal_negative_screen.tsv"))
    p.add_argument("--fasta", type=Path, nargs="+", required=True)
    p.add_argument("--gff3", type=Path, nargs="+", required=True)
    p.add_argument("--amrfinderplus", type=Path, nargs="+", required=True)
    p.add_argument("--kleborate", type=Path, nargs="+", required=True)
    p.add_argument("--out", type=Path, default=Path("results/klebtrace.sqlite"))
    args = p.parse_args(argv)

    con = create(args.out, args.schema)
    try:
        counts = {
            "genomes": load_genomes(con, args.samples, args.manifest),
            "sequences": load_sequences(con, args.fasta),
            "features": load_features(con, args.gff3),
            "portal calls": load_portal_calls(con, args.portal_calls),
            "AMRFinderPlus calls": load_amrfinderplus(con, args.amrfinderplus),
            "Kleborate calls": load_kleborate(con, args.kleborate),
            "split-gene screen calls": load_screen(con, args.screen) if args.screen.exists() else 0,
        }
        con.executemany("INSERT INTO tool_versions VALUES (?, ?, ?, ?, ?)", TOOL_VERSIONS)
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        con.close()
        args.out.unlink(missing_ok=True)
        raise
    con.close()
    for table, n in counts.items():
        print(f"{n:>7} {table}")
    print(f"Database written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
