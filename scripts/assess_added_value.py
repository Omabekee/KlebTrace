#!/usr/bin/env python3
"""Compare KlebTrace's merged annotation with Ensembl's annotation of the same genome.

For each feature from a specialist tool (AMRFinderPlus, PlasmidFinder) and for
Bakta, find the Ensembl gene it overlaps on the same strand and classify it:

    same_name        Ensembl's gene has the same name (blaNDM-1 and blaNDM-1_2)
    other_name       Ensembl has a gene there, under another name or none
    no_ensembl_gene  Ensembl has no gene at this position
    not_a_gene       a feature Ensembl does not represent as a gene
                     (point mutation call, plasmid replicon)

Both files must use the same sequence names, so run to_ensembl_seqids.py first.

Usage:
    python scripts/assess_added_value.py --assembly GCA_000764615.1 \\
        --klebtrace results/ensembl_tracks/GCA_000764615.1.ensembl_names.gff3 \\
        --ensembl data/ensembl/Klebsiella_pneumoniae_gca_000764615.ASM76461v1.63.gff3.gz \\
        --outdir results/tool_assessment
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

COLUMNS = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]


def read_gff3(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.rstrip("\n").split("\t"))
    df = pd.DataFrame(rows, columns=COLUMNS)
    df[["start", "end"]] = df[["start", "end"]].astype(int)
    return df


def attr(df: pd.DataFrame, key: str) -> pd.Series:
    return df["attributes"].str.extract(rf"(?:^|;){key}=([^;]+)")[0].map(
        lambda v: unquote(v) if isinstance(v, str) else v)


def base_name(name) -> str:
    """Prokka numbers repeated gene names (blaNDM-1_2, dsbD_3); drop that suffix."""
    return re.sub(r"_\d+$", "", str(name)) if isinstance(name, str) else ""


def ensembl_genes(ensembl: pd.DataFrame) -> pd.DataFrame:
    genes = ensembl[ensembl["type"].isin(["gene", "ncRNA_gene"])].copy()
    genes["name"] = attr(genes, "Name").map(base_name)
    genes["description"] = attr(genes, "description").fillna("")
    return genes


def classify(features: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    out = []
    for f in features.itertuples():
        hits = genes[(genes["seqid"] == f.seqid) & (genes["strand"] == f.strand)
                     & (genes["start"] <= f.end) & (genes["end"] >= f.start)]
        subtype = f.element_subtype if isinstance(f.element_subtype, str) else ""
        if f.source == "PlasmidFinder" or subtype.startswith("POINT"):
            status, best = "not_a_gene", None
        elif hits.empty:
            status, best = "no_ensembl_gene", None
        else:
            overlap = hits["end"].clip(upper=f.end) - hits["start"].clip(lower=f.start)
            best = hits.loc[overlap.idxmax()]
            status = "same_name" if best["name"] and best["name"] == base_name(f.name) else "other_name"
        out.append({
            "source": f.source, "type": f.type, "seqid": f.seqid, "start": f.start, "end": f.end,
            "name": f.name, "subclass": f.subclass, "element_subtype": subtype, "status": status,
            "ensembl_name": None if best is None else best["name"],
            "ensembl_description": None if best is None else best["description"],
        })
    return pd.DataFrame(out)


def assess(klebtrace: Path, ensembl: Path) -> pd.DataFrame:
    kt = read_gff3(klebtrace)
    # Bakta keeps the gene symbol in gene= and the product in Name=
    kt["name"] = attr(kt, "gene").where(kt["source"] == "Bakta", attr(kt, "Name"))
    kt["subclass"] = attr(kt, "subclass")
    kt["element_subtype"] = attr(kt, "element_subtype")
    specialist = kt[kt["source"].isin(["AMRFinderPlus", "PlasmidFinder"])]
    bakta_cds = kt[(kt["source"] == "Bakta") & (kt["type"] == "CDS")]
    genes = ensembl_genes(read_gff3(ensembl))
    return pd.concat([classify(specialist, genes), classify(bakta_cds, genes)], ignore_index=True)


def summary(table: pd.DataFrame, assembly: str, ensembl: Path) -> pd.DataFrame:
    genes = ensembl_genes(read_gff3(ensembl))
    rows = []
    for source, group in table.groupby("source"):
        counts = group["status"].value_counts()
        rows.append({"assembly": assembly, "source": source, "features": len(group),
                     **{s: int(counts.get(s, 0)) for s in ("same_name", "other_name", "no_ensembl_gene", "not_a_gene")}})
    rows.append({"assembly": assembly, "source": "Ensembl (Prokka via MGnify)",
                 "features": len(genes), "same_name": None, "other_name": None,
                 "no_ensembl_gene": None, "not_a_gene": None})
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--assembly", required=True)
    p.add_argument("--klebtrace", type=Path, required=True)
    p.add_argument("--ensembl", type=Path, required=True)
    p.add_argument("--outdir", type=Path, default=Path("results/tool_assessment"))
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    table = assess(args.klebtrace, args.ensembl)
    table.to_csv(args.outdir / f"{args.assembly}.added_value.tsv", sep="\t", index=False)
    s = summary(table, args.assembly, args.ensembl)
    s.to_csv(args.outdir / f"{args.assembly}.added_value_summary.tsv", sep="\t", index=False)
    print(s.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
