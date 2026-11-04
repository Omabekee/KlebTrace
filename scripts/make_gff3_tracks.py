#!/usr/bin/env python3
"""Convert AMRFinderPlus and PlasmidFinder tables into GFF3 tracks for a genome browser.

Bakta already writes GFF3. AMRFinderPlus and PlasmidFinder write tables, so
this script turns one genome's tables into two GFF3 files that JBrowse (or an
Ensembl custom track) can display alongside the Bakta annotation:

    <assembly>.amrfinderplus.gff3
    <assembly>.plasmidfinder.gff3

Coordinates: both tools report 1-based, inclusive positions, which is also the
GFF3 convention, so positions are copied unchanged. (The Phase 0 conversion in
Galaxy shifted starts by +1; see the learning log.)

Contig names: PlasmidFinder reports the full FASTA header, for example
"CP008933.1 Klebsiella pneumoniae strain PMK1 plasmid ...". Only the first
word, the sequence accession, is kept so that it matches the FASTA and the
Bakta GFF3.

Usage:
    python scripts/make_gff3_tracks.py --assembly GCA_000764615.1 \
        --indir results/phase2_galaxy --outdir results/tracks
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import quote

import pandas as pd

GFF_COLUMNS = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]


def encode(value) -> str:
    """Percent-encode the characters GFF3 reserves in attribute values."""
    return quote(str(value), safe=" /()[]:.+-_'*")


def attributes(pairs: dict) -> str:
    return ";".join(f"{k}={encode(v)}" for k, v in pairs.items()
                    if v is not None and str(v) not in ("", "nan", "NA"))


def first_word(text: str) -> str:
    return str(text).strip().split()[0]


def amrfinder_to_gff(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for i, r in enumerate(df.to_dict("records"), start=1):
        subtype = str(r["Subtype"])
        rows.append({
            "seqid": first_word(r["Contig id"]),
            "source": "AMRFinderPlus",
            "type": "point_mutation" if subtype.startswith("POINT") else "gene",
            "start": int(r["Start"]),
            "end": int(r["Stop"]),
            "score": ".",
            "strand": r["Strand"] if r["Strand"] in ("+", "-") else ".",
            "phase": ".",
            "attributes": attributes({
                "ID": f"amrfinderplus_{i}",
                "Name": r["Element symbol"],
                "product": r["Element name"],
                "element_type": r["Type"],
                "element_subtype": subtype,
                "class": r["Class"],
                "subclass": r["Subclass"],
                "method": r["Method"],
                "coverage": r["% Coverage of reference"],
                "identity": r["% Identity to reference"],
            }),
        })
    return pd.DataFrame(rows, columns=GFF_COLUMNS)


def plasmidfinder_to_gff(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for i, r in enumerate(df.to_dict("records"), start=1):
        start, end = str(r["Position in contig"]).split("..")
        rows.append({
            "seqid": first_word(r["Contig"]),
            "source": "PlasmidFinder",
            "type": "origin_of_replication",
            "start": int(start),
            "end": int(end),
            "score": ".",
            "strand": ".",
            "phase": ".",
            "attributes": attributes({
                "ID": f"plasmidfinder_{i}",
                "Name": r["Plasmid"],
                "identity": r["Identity"],
                "length_query_template": r["Query / Template length"],
                "reference": r.get("Accession number"),
            }),
        })
    return pd.DataFrame(rows, columns=GFF_COLUMNS)


def write_gff(df: pd.DataFrame, path: Path) -> None:
    with path.open("w") as fh:
        fh.write("##gff-version 3\n")
        df.sort_values(["seqid", "start"]).to_csv(fh, sep="\t", header=False, index=False)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--assembly", required=True, help="e.g. GCA_000764615.1")
    p.add_argument("--indir", type=Path, default=Path("results/phase2_galaxy"))
    p.add_argument("--outdir", type=Path, default=Path("results/tracks"))
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    amr = pd.read_csv(args.indir / "amrfinderplus_report" / f"{args.assembly}.tsv", sep="\t", dtype=str)
    pf = pd.read_csv(args.indir / "plasmidfinder" / f"{args.assembly}.tsv", sep="\t", dtype=str)

    write_gff(amrfinder_to_gff(amr), args.outdir / f"{args.assembly}.amrfinderplus.gff3")
    write_gff(plasmidfinder_to_gff(pf), args.outdir / f"{args.assembly}.plasmidfinder.gff3")
    print(f"Wrote {len(amr)} AMRFinderPlus and {len(pf)} PlasmidFinder features to {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
