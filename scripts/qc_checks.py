#!/usr/bin/env python3
"""Quality checks on the merged annotation and the loaded database.

Each check returns rows of (check, assembly, status, detail), with status PASS,
WARN or FAIL, and all rows go into one QC report. FAIL means the data are
wrong and the pipeline stops; WARN means the data are valid but a person
should look at them.

    gff3_valid             9 columns, integer coordinates, start <= end, strand
                           in + - . ?, CDS phase 0/1/2, ##gff-version header
    so_types               every type is a known Sequence Ontology term (WARN)
    seqids_match_fasta     every seqid is in the genome's FASTA and every
                           feature ends within its sequence
    unique_ids             no repeated ID attribute within a genome
    cds_density            Bakta CDS per kb within the range for
                           K. pneumoniae; scales with genome size, unlike a
                           fixed gene count (WARN)
    coordinates_agree      each AMRFinderPlus gene matches a Bakta CDS on the
                           same strand: 3' ends equal or 3 bp apart (Bakta's CDS
                           includes the stop codon), 5' ends equal or a whole
                           number of codons apart (a different start codon) (WARN)
    database_complete      every genome has sequences, features, a phenotype
                           and AMRFinderPlus calls
    duplicated_markers     more than one full copy of a single-copy chromosomal
                           gene, a sign of duplicated contigs (WARN)
    carbapenemase_review   a carbapenemase call that is partial or unnamed, or
                           that the portal does not report (WARN)

Usage:
    python scripts/qc_checks.py --database results/klebtrace.sqlite \\
        --fasta data/genomes/*.fasta --gff3 results/merged/*.klebtrace.gff3 --out results/qc_report.tsv
"""

from __future__ import annotations

import argparse
import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd

from merge_annotations import read_fasta_lengths

REPORT_COLUMNS = ["check", "assembly", "status", "detail"]

# Bakta CDS per kb. The seven study genomes range from 0.93 to 0.96 (including
# small ORFs); a value far outside this suggests a broken or contaminated
# assembly or annotation.
CDS_PER_KB = (0.85, 1.00)

# Sequence Ontology term names used by Bakta, AMRFinderPlus and PlasmidFinder
# in the merged files. crispr-spacer is Bakta's own name and is not in this list.
SO_TERMS = {"region", "gene", "CDS", "tRNA", "tmRNA", "rRNA", "ncRNA", "regulatory_region",
            "origin_of_replication", "origin_of_transfer", "CRISPR", "direct_repeat", "gap"}

SINGLE_COPY_MARKERS = ("fosA", "oqxA", "oqxB", "kdeA", "emrD")
STRANDS = {"+", "-", ".", "?"}


def assembly_of(path: Path) -> str:
    return re.sub(r"\.klebtrace\.gff3$|\.fasta$", "", path.name)


def read_gff3(path: Path) -> tuple[list[str], list[list[str]]]:
    header, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
            elif line.strip():
                rows.append(line.rstrip("\n").split("\t"))
    return header, rows


def gff3_frame(rows: list[list[str]]) -> pd.DataFrame:
    df = pd.DataFrame([r for r in rows if len(r) == 9],
                      columns=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["ID"] = df["attributes"].str.extract(r"(?:^|;)ID=([^;]+)")[0]
    df["Name"] = df["attributes"].str.extract(r"(?:^|;)Name=([^;]+)")[0]
    return df


def result(check: str, assembly: str, problems: list[str], ok: str, level: str = "FAIL") -> tuple:
    if not problems:
        return (check, assembly, "PASS", ok)
    shown = "; ".join(problems[:5]) + (f"; and {len(problems) - 5} more" if len(problems) > 5 else "")
    return (check, assembly, level, f"{len(problems)} problem(s): {shown}")


# --- GFF3 checks -------------------------------------------------------------

def check_gff3_valid(assembly: str, header: list[str], rows: list[list[str]]) -> tuple:
    problems = []
    if not header or header[0] != "##gff-version 3":
        problems.append("first line is not ##gff-version 3")
    for i, r in enumerate(rows, start=1):
        if len(r) != 9:
            problems.append(f"feature {i}: {len(r)} columns")
            continue
        if not (r[3].isdigit() and r[4].isdigit()):
            problems.append(f"feature {i}: non-integer coordinates {r[3]}..{r[4]}")
        elif not 1 <= int(r[3]) <= int(r[4]):
            problems.append(f"feature {i}: start {r[3]} after end {r[4]} or below 1")
        if r[6] not in STRANDS:
            problems.append(f"feature {i}: strand '{r[6]}'")
        if r[2] == "CDS" and r[7] not in ("0", "1", "2"):
            problems.append(f"feature {i}: CDS phase '{r[7]}'")
        if any("=" not in pair for pair in r[8].split(";") if pair):
            problems.append(f"feature {i}: attribute without '=' (unencoded ';'?)")
    return result("gff3_valid", assembly, problems, f"{len(rows)} features valid")


def check_so_types(assembly: str, gff: pd.DataFrame) -> tuple:
    unknown = gff.loc[~gff["type"].isin(SO_TERMS), "type"].value_counts()
    problems = [f"{t} ({n})" for t, n in unknown.items()]
    return result("so_types", assembly, problems, "all types are SO terms", "WARN")


def check_seqids_match_fasta(assembly: str, gff: pd.DataFrame, lengths: dict[str, int]) -> tuple:
    problems = [f"{s} not in FASTA" for s in sorted(set(gff["seqid"]) - set(lengths))]
    known = gff[gff["seqid"].isin(lengths)]
    beyond = known[known["end"] > known["seqid"].map(lengths)]
    problems += [f"{r.ID} ends at {r.end} on {r.seqid} ({lengths[r.seqid]} bp)" for r in beyond.itertuples()]
    return result("seqids_match_fasta", assembly, problems,
                  f"{gff['seqid'].nunique()} seqids, all in the FASTA ({len(lengths)} sequences)")


def check_unique_ids(assembly: str, gff: pd.DataFrame) -> tuple:
    missing = int(gff["ID"].isna().sum())
    dup = gff.loc[gff["ID"].duplicated() & gff["ID"].notna(), "ID"].unique().tolist()
    problems = ([f"{missing} features without ID"] if missing else []) + [f"duplicate ID {d}" for d in dup]
    return result("unique_ids", assembly, problems, f"{gff['ID'].nunique()} unique IDs")


def check_cds_density(assembly: str, gff: pd.DataFrame, genome_bp: int,
                      expected: tuple[float, float] = CDS_PER_KB) -> tuple:
    n_cds = int(((gff["source"] == "Bakta") & (gff["type"] == "CDS")).sum())
    density = n_cds / (genome_bp / 1000)
    detail = f"{n_cds} CDS on {genome_bp / 1e6:.2f} Mb = {density:.3f} per kb"
    problems = [] if expected[0] <= density <= expected[1] else [f"{detail}, outside {expected[0]}-{expected[1]}"]
    return result("cds_density", assembly, problems, detail, "WARN")


def compare_ends(strand: str, d_start: int, d_end: int) -> tuple[int, int]:
    """Convert Bakta-minus-AMRFinderPlus coordinate differences into differences at the
    gene's 5' (start codon) and 3' (stop codon) ends; positive means Bakta's CDS is longer."""
    return (-d_start, d_end) if strand == "+" else (d_end, -d_start)


def check_coordinates_agree(assembly: str, gff: pd.DataFrame) -> tuple:
    """Compare each full-length AMRFinderPlus gene with the overlapping Bakta CDS on the same strand.

    Expected: Bakta's CDS includes the stop codon and AMRFinderPlus's does not, so the 3' end
    differs by 3 bp. Bakta (Prodigal) may also choose a different start codon from the reference,
    which moves the 5' end by a multiple of 3 bp; that is common and recorded, not flagged.
    Flagged: a 3' end that differs otherwise, a 5' shift that breaks the reading frame, or an
    AMRFinderPlus gene with no Bakta CDS at all.
    """
    amr = gff[(gff["source"] == "AMRFinderPlus") & ~gff["attributes"].str.contains("element_subtype=POINT")
              & ~gff["attributes"].str.contains("method=PARTIAL")]
    cds = gff[(gff["source"] == "Bakta") & (gff["type"] == "CDS")]
    same_ends = other_start = 0
    problems = []
    for a in amr.itertuples():
        hits = cds[(cds["seqid"] == a.seqid) & (cds["strand"] == a.strand)
                   & (cds["start"] <= a.end) & (cds["end"] >= a.start)]
        if hits.empty:
            problems.append(f"{a.Name} {a.seqid}:{a.start}-{a.end} has no Bakta CDS")
            continue
        overlap = hits[["end"]].clip(upper=a.end)["end"] - hits[["start"]].clip(lower=a.start)["start"]
        b = hits.loc[overlap.idxmax()]
        d5, d3 = compare_ends(a.strand, int(b["start"] - a.start), int(b["end"] - a.end))
        if d3 in (0, 3) and d5 == 0:
            same_ends += 1
        elif d3 in (0, 3) and d5 % 3 == 0:
            other_start += 1
        else:
            problems.append(f"{a.Name} {a.seqid}: 5' end {d5:+d} bp, 3' end {d3:+d} bp")
    ok = (f"{same_ends + other_start} of {len(amr)} genes agree ({same_ends} same start, "
          f"{other_start} different start codon; 3' ends differ only by the stop codon)")
    return result("coordinates_agree", assembly, problems, ok, "WARN")


# --- Database checks ---------------------------------------------------------

def check_database_complete(con: sqlite3.Connection) -> list[tuple]:
    sql = """
        SELECT g.assembly,
               (SELECT count(*) FROM sequences s WHERE s.assembly = g.assembly),
               (SELECT count(*) FROM features f WHERE f.assembly = g.assembly),
               (SELECT count(*) FROM phenotypes p WHERE p.assembly = g.assembly),
               (SELECT count(*) FROM resistance_calls c WHERE c.assembly = g.assembly
                                                       AND c.source = 'amrfinderplus')
        FROM genomes g ORDER BY g.assembly"""
    out = []
    for assembly, n_seq, n_feat, n_pheno, n_amr in con.execute(sql):
        problems = [name for name, n in (("sequences", n_seq), ("features", n_feat), ("phenotype", n_pheno),
                                         ("AMRFinderPlus calls", n_amr)) if n == 0]
        out.append(result("database_complete", assembly, [f"no {p}" for p in problems],
                          f"{n_seq} sequences, {n_feat} features, {n_amr} AMRFinderPlus calls"))
    return out


def check_duplicated_markers(con: sqlite3.Connection) -> list[tuple]:
    marks = ",".join("?" * len(SINGLE_COPY_MARKERS))
    counts = pd.read_sql(
        f"SELECT g.assembly, c.gene, count(c.gene) AS n FROM genomes g "
        f"LEFT JOIN resistance_calls c ON c.assembly = g.assembly AND c.source = 'amrfinderplus' "
        f"AND c.is_partial = 0 AND c.gene IN ({marks}) GROUP BY g.assembly, c.gene",
        con, params=SINGLE_COPY_MARKERS)
    out = []
    for assembly, group in counts.groupby("assembly"):
        extra = group[group["n"] > 1]
        problems = [f"{r.gene} x{r.n}" for r in extra.itertuples()]
        out.append(result("duplicated_markers", assembly, problems,
                          "single-copy chromosomal genes found once", "WARN"))
    return out


def check_carbapenemase_review(con: sqlite3.Connection) -> list[tuple]:
    calls = pd.read_sql(
        "SELECT g.assembly, c.source, c.gene, c.is_partial FROM genomes g "
        "LEFT JOIN resistance_calls c ON c.assembly = g.assembly AND c.is_carbapenemase = 1", con)
    out = []
    for assembly, group in calls.groupby("assembly"):
        amr = group[group["source"] == "amrfinderplus"]
        portal = group[group["source"] == "portal"]
        problems = [f"{r.gene} partial" for r in amr.itertuples() if r.is_partial]
        problems += [f"{r.gene} allele not named" for r in amr.itertuples()
                     if not r.is_partial and not str(r.gene)[-1].isdigit()]
        if len(amr) and portal.empty:
            problems.append("carbapenemase not reported by the portal")
        n = len(amr)
        out.append(result("carbapenemase_review", assembly, problems,
                          f"{n} carbapenemase call(s), all named, complete and in the portal" if n
                          else "no carbapenemase", "WARN"))
    return out


def run_checks(database: Path, fastas: list[Path], gff3s: list[Path]) -> pd.DataFrame:
    fasta_by_assembly = {assembly_of(f): f for f in fastas}
    rows = []
    for path in sorted(gff3s):
        assembly = assembly_of(path)
        header, raw = read_gff3(path)
        gff = gff3_frame(raw)
        rows.append(check_gff3_valid(assembly, header, raw))
        rows.append(check_so_types(assembly, gff))
        if assembly not in fasta_by_assembly:
            rows.append(("seqids_match_fasta", assembly, "FAIL", "no FASTA given for this genome"))
            continue
        lengths = read_fasta_lengths(fasta_by_assembly[assembly])
        rows.append(check_seqids_match_fasta(assembly, gff, lengths))
        rows.append(check_unique_ids(assembly, gff))
        rows.append(check_cds_density(assembly, gff, sum(lengths.values())))
        rows.append(check_coordinates_agree(assembly, gff))
    with sqlite3.connect(database) as con:
        rows += check_database_complete(con)
        rows += check_duplicated_markers(con)
        rows += check_carbapenemase_review(con)
    return pd.DataFrame(rows, columns=REPORT_COLUMNS)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--database", type=Path, required=True)
    p.add_argument("--fasta", type=Path, nargs="+", required=True)
    p.add_argument("--gff3", type=Path, nargs="+", required=True)
    p.add_argument("--out", type=Path, default=Path("qc_report.tsv"))
    args = p.parse_args(argv)
    report = run_checks(args.database, args.fasta, args.gff3)
    report.to_csv(args.out, sep="\t", index=False)
    print(report["status"].value_counts().to_string())
    print(report[report["status"] != "PASS"].to_string(index=False))
    return 1 if (report["status"] == "FAIL").any() else 0


if __name__ == "__main__":
    sys.exit(main())
