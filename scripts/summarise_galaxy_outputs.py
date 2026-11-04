#!/usr/bin/env python3
"""Summarise the KlebTrace Phase 2 Galaxy outputs, one row per genome.

Reads the per-genome files downloaded from the Galaxy workflow v0.2
(AMRFinderPlus, PlasmidFinder, IntegronFinder and ISEScan) and writes:

    genome_summary.tsv     one row per genome: assembly statistics,
                           carbapenem resistance mechanisms, other
                           beta-lactamases, virulence loci, plasmid
                           replicons, integrons, IS elements and QC flags
    amrfinderplus_all.tsv  every AMRFinderPlus call for all genomes in one
                           table, with an `assembly` column

Assembly statistics come from data/genomes_manifest.tsv (written by
download_genomes.py from the FASTA files themselves). The IntegronFinder summary
is not used for this, because it leaves out some short sequences (for example,
294 of the 304 sequences in GCA_001182445.1).

Usage:
    python scripts/summarise_galaxy_outputs.py --indir results/phase2_galaxy
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Genes found once on the chromosome of almost every K. pneumoniae. More than
# one full-length copy in a draft assembly suggests duplicated contigs or a
# second strain in the sample.
CHROMOSOMAL_MARKERS = ("fosA", "oqxA", "oqxB", "kdeA", "emrD", "fieF")

FRAGMENTED_CONTIGS = 200

# When AMRFinderPlus cannot name the allele (for example, a gene split across
# contig ends), it reports a family symbol such as "blaOXA" with the subclass of
# the whole family (BETA-LACTAM), not CARBAPENEM. The closest reference name
# still says "carbapenem-hydrolyzing", so carbapenemases are found by either.
CARBAPENEMASE_REFERENCE = "carbapenem-hydrolyzing"


def accession(path: Path) -> str:
    return path.name.removesuffix(".tsv")


def read_amrfinder(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.insert(0, "assembly", accession(path))
    return df


def read_isescan_total(path: Path) -> int:
    df = pd.read_csv(
        path, sep=r"\s+", skiprows=1, header=None,
        names=["seqid", "family", "nIS", "pct_genome", "bps4IS", "dnaLen"],
    )
    total = df.loc[df["family"] == "total", "nIS"]
    return int(total.iloc[0]) if len(total) else int(df["nIS"].sum())


def join(values) -> str:
    return ";".join(values) if len(values) else "none"


def summarise_genome(acc: str, indir: Path, manifest: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    stats = manifest.loc[acc]
    sizes = pd.read_csv(indir / "integronfinder_summary" / f"{acc}.tsv", sep="\t", comment="#")
    amr = read_amrfinder(indir / "amrfinderplus_report" / f"{acc}.tsv")
    plasmids = pd.read_csv(indir / "plasmidfinder" / f"{acc}.tsv", sep="\t")
    n_is = read_isescan_total(indir / "isescan_summary" / f"{acc}.tsv")

    subclass = amr["Subclass"].fillna("")
    subtype = amr["Subtype"].fillna("")
    method = amr["Method"].fillna("")
    partial = method.str.startswith("PARTIAL")

    def label(rows: pd.DataFrame) -> list[str]:
        out = []
        for symbol, group in rows.groupby("Element symbol", sort=True):
            tag = symbol
            if len(group) > 1:
                tag += f" (x{len(group)})"
            if group["Method"].str.startswith("PARTIAL").any():
                tag += " [partial]"
            if symbol.startswith("bla") and not symbol[-1].isdigit() and "Closest reference name" in group:
                closest = group["Closest reference name"].iloc[0].split()[-1]
                tag += f" [closest {closest}]"
            out.append(tag)
        return out

    is_carb = subclass.str.contains("CARBAPENEM")
    names = amr.get("Closest reference name", pd.Series("", index=amr.index)).fillna("")
    carbapenemases = amr[
        (is_carb | names.str.contains(CARBAPENEMASE_REFERENCE)) & (subtype == "AMR")
    ]
    porins = amr[is_carb & subtype.str.startswith("POINT")]
    other_bla = amr[
        amr["Element symbol"].str.startswith("bla")
        & ~amr.index.isin(carbapenemases.index)
        & ~is_carb
        & subclass.str.contains("CEPHALOSPORIN")
    ]
    virulence = amr[amr["Type"] == "VIRULENCE"]

    flags = []
    if stats["n_sequences"] > FRAGMENTED_CONTIGS:
        flags.append(f"fragmented assembly ({stats['n_sequences']} sequences)")
    if (carbapenemases["Method"].str.startswith("PARTIAL")).any():
        flags.append("carbapenemase gene partial at contig end")
    if not carbapenemases["Element symbol"].str[-1].str.isdigit().all():
        flags.append("carbapenemase allele not named")
    full = amr[~partial]
    dup = [g for g in CHROMOSOMAL_MARKERS if (full["Element symbol"] == g).sum() > 1]
    if dup:
        flags.append("duplicated chromosomal genes: " + ",".join(dup))

    row = {
        "assembly": acc,
        "n_sequences": int(stats["n_sequences"]),
        "total_bp": int(stats["total_bp"]),
        "n50_bp": int(stats["n50"]),
        "largest_contig_bp": int(sizes["size"].max()),
        "carbapenemase_genes": join(label(carbapenemases)),
        "carbapenem_porin_changes": join(label(porins)),
        "other_esbl_ampc": join(label(other_bla)),
        "virulence_loci": join(sorted(virulence["Element symbol"].unique())),
        "n_amrfinder_calls": len(amr),
        "n_partial_calls": int(partial.sum()),
        "plasmid_replicons": join(sorted(plasmids["Plasmid"].unique())),
        "integrons_complete": int(sizes["complete"].sum()),
        "integrons_in0": int(sizes["In0"].sum()),
        "integrons_calin": int(sizes["CALIN"].sum()),
        "n_is_elements": n_is,
        "qc_flags": join(flags),
    }
    return row, amr


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--indir", type=Path, default=Path("results/phase2_galaxy"))
    p.add_argument("--outdir", type=Path, default=None, help="defaults to --indir")
    p.add_argument("--manifest", type=Path, default=Path("data/genomes_manifest.tsv"))
    args = p.parse_args(argv)
    outdir = args.outdir or args.indir
    outdir.mkdir(parents=True, exist_ok=True)

    accessions = sorted(accession(f) for f in (args.indir / "amrfinderplus_report").glob("*.tsv"))
    if not accessions:
        print(f"No AMRFinderPlus files found in {args.indir}", file=sys.stderr)
        return 1

    manifest = pd.read_csv(args.manifest, sep="\t").set_index("assembly")
    rows, tables = [], []
    for acc in accessions:
        row, amr = summarise_genome(acc, args.indir, manifest)
        rows.append(row)
        tables.append(amr)

    pd.DataFrame(rows).to_csv(outdir / "genome_summary.tsv", sep="\t", index=False)
    pd.concat(tables).to_csv(outdir / "amrfinderplus_all.tsv", sep="\t", index=False)
    print(f"Summarised {len(rows)} genomes into {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
