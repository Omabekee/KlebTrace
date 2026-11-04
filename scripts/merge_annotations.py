#!/usr/bin/env python3
"""Merge one genome's tool outputs into a single harmonised GFF3.

Inputs: the FASTA, Bakta's GFF3, and the AMRFinderPlus, Kleborate and
PlasmidFinder tables. Output: one GFF3 in which

- every seqid is the sequence accession from the FASTA header (for example
  CP008933.1), checked against the FASTA, with one ##sequence-region pragma per
  sequence and no embedded ##FASTA section;
- column 2 (source) names the tool (Bakta, AMRFinderPlus, PlasmidFinder); for
  Bakta, the program that predicted the feature (Prodigal, tRNAscan-SE, ...)
  is kept in a `predicted_by` attribute;
- column 3 (type) is a Sequence Ontology term (Bakta's oriC, oriV, oriT and
  crispr-repeat are renamed; crispr-spacer is kept and flagged by qc_checks.py);
- column 9 has ID, Name, product and Dbxref, with GFF3 reserved characters
  percent-encoded and IDs prefixed by tool so they are unique in the file.

Kleborate reports whole-genome results without coordinates (species, sequence
type, virulence and resistance scores), so these go in `# kleborate` comment
lines in the header.

AMRFinderPlus point mutations (for example ompK36_D135DGD) are reported with the
coordinates of the whole affected gene, not of the change itself, so they are
written as type `gene` with `element_subtype=POINT`, not as a 1-bp
point_mutation.

Usage:
    merge_annotations.py --assembly GCA_000764615.1 --fasta X.fasta --bakta X.gff3 \\
        --amrfinderplus X.tsv --kleborate X.tsv --plasmidfinder X.tsv --out X.klebtrace.gff3
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import quote

import pandas as pd

GFF_COLUMNS = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

# Bakta types that are not Sequence Ontology term names
SO_RENAME = {"oriC": "origin_of_replication", "oriV": "origin_of_replication",
             "oriT": "origin_of_transfer", "crispr-repeat": "direct_repeat"}

KLEBORATE_FIELDS = {
    "species": "general__species__species",
    "species_match": "general__species__species_match",
    "ST": "kpsc__mlst__ST",
    "virulence_score": "kpsc__virulence_score__virulence_score",
    "resistance_score": "kpsc__resistance_score__resistance_score",
    "carbapenemases": "kpsc__amr__Bla_Carb_acquired",
    "omp_mutations": "kpsc__amr__Omp_mutations",
    "K_locus": "kpsc__kaptive__K_locus",
    "O_type": "kpsc__kaptive__O_type",
}


def encode(value) -> str:
    """Percent-encode the characters GFF3 reserves in attribute values (; = & , tab, newline, %)."""
    return quote(str(value), safe=" /()[]:.+-_'*<>|!?#@^~")


def attributes(pairs: dict) -> str:
    return ";".join(f"{k}={encode(v)}" for k, v in pairs.items()
                    if v is not None and str(v) not in ("", "nan", "NA", "-"))


def read_fasta_lengths(fasta: Path) -> dict[str, int]:
    """Return {accession: length} for every sequence, in file order."""
    lengths: dict[str, int] = {}
    current = None
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in lengths:
                    raise ValueError(f"{fasta}: duplicate sequence {current}")
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len(line.strip())
    return lengths


def read_bakta_gff3(path: Path) -> pd.DataFrame:
    """Read Bakta's feature lines, stopping at the ##FASTA section."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}: expected 9 columns, got {len(fields)}: {line[:80]}")
            rows.append(fields)
    df = pd.DataFrame(rows, columns=GFF_COLUMNS)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def bakta_features(gff: pd.DataFrame) -> pd.DataFrame:
    """Harmonise Bakta features: source Bakta, SO type names, predictor kept as an attribute.

    Bakta's own attribute values are already percent-encoded, so they are kept as written.
    """
    out = gff.copy()
    predictor = out["source"]
    out["source"] = "Bakta"
    out["type"] = out["type"].replace(SO_RENAME)
    out["attributes"] = [
        attrs if pred == "Bakta" else f"{attrs};predicted_by={encode(pred)}"
        for attrs, pred in zip(out["attributes"], predictor)
    ]
    return out


def amrfinderplus_features(report: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for i, r in enumerate(report.to_dict("records"), start=1):
        accession = r.get("Closest reference accession")
        rows.append({
            "seqid": str(r["Contig id"]).split()[0],
            "source": "AMRFinderPlus",
            "type": "gene",
            "start": int(r["Start"]),
            "end": int(r["Stop"]),
            "score": ".",
            "strand": r["Strand"] if r["Strand"] in ("+", "-") else ".",
            "phase": ".",
            "attributes": attributes({
                "ID": f"amrfinderplus_{i}",
                "Name": r["Element symbol"],
                "product": r["Element name"],
                "Dbxref": f"RefSeq:{accession}" if accession not in (None, "", "NA") else None,
                "element_type": r["Type"],
                "element_subtype": r["Subtype"],
                "class": r["Class"],
                "subclass": r["Subclass"],
                "method": r["Method"],
                "coverage": r["% Coverage of reference"],
                "identity": r["% Identity to reference"],
                "closest_reference": r.get("Closest reference name"),
            }),
        })
    return pd.DataFrame(rows, columns=GFF_COLUMNS)


def plasmidfinder_features(report: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for i, r in enumerate(report.to_dict("records"), start=1):
        start, end = str(r["Position in contig"]).split("..")
        rows.append({
            "seqid": str(r["Contig"]).split()[0],
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
                "product": "plasmid replicon",
                "Dbxref": f"INSDC:{r['Accession number']}" if r.get("Accession number") else None,
                "identity": r["Identity"],
                "length_query_template": r["Query / Template length"],
            }),
        })
    return pd.DataFrame(rows, columns=GFF_COLUMNS)


def kleborate_header(report: pd.DataFrame, assembly: str) -> list[str]:
    """Genome-level Kleborate results as GFF3 comment lines."""
    rows = report[report["strain"].astype(str).isin({assembly, assembly.rsplit(".", 1)[0]})]
    if rows.empty:
        raise ValueError(f"No Kleborate row for {assembly}")
    row = rows.iloc[0]
    values = {k: row.get(col) for k, col in KLEBORATE_FIELDS.items()}
    return [f"# kleborate {k}={v}" for k, v in values.items() if pd.notna(v)]


def check_against_fasta(features: pd.DataFrame, lengths: dict[str, int]) -> None:
    unknown = sorted(set(features["seqid"]) - set(lengths))
    if unknown:
        raise ValueError(f"Features on sequences not in the FASTA: {unknown[:5]}")
    too_long = features[features["end"] > features["seqid"].map(lengths)]
    if not too_long.empty:
        raise ValueError(f"Features end beyond their sequence: {too_long.head(3).to_dict('records')}")


def merge(assembly: str, fasta: Path, bakta: Path, amrfinderplus: Path,
          kleborate: Path, plasmidfinder: Path) -> list[str]:
    """Return the merged GFF3 as lines."""
    lengths = read_fasta_lengths(fasta)
    read = lambda p: pd.read_csv(p, sep="\t", dtype=str)  # noqa: E731
    features = pd.concat([
        bakta_features(read_bakta_gff3(bakta)),
        amrfinderplus_features(read(amrfinderplus)),
        plasmidfinder_features(read(plasmidfinder)),
    ], ignore_index=True)
    check_against_fasta(features, lengths)

    order = {seqid: i for i, seqid in enumerate(lengths)}
    features = features.assign(_order=features["seqid"].map(order)).sort_values(
        ["_order", "start", "end", "source"], kind="stable").drop(columns="_order")

    ids = features["attributes"].str.extract(r"(?:^|;)ID=([^;]+)")[0]
    duplicated = ids[ids.duplicated() & ids.notna()]
    if not duplicated.empty:
        raise ValueError(f"Duplicate IDs: {sorted(set(duplicated))[:5]}")

    lines = ["##gff-version 3",
             "##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.1/so.obo",
             f"# KlebTrace merged annotation for {assembly}: Bakta, AMRFinderPlus, PlasmidFinder"]
    lines += kleborate_header(read(kleborate), assembly)
    lines += [f"##sequence-region {seqid} 1 {length}" for seqid, length in lengths.items()]
    lines += ["\t".join(map(str, row)) for row in features[GFF_COLUMNS].itertuples(index=False)]
    return lines


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--assembly", required=True)
    for name in ("fasta", "bakta", "amrfinderplus", "kleborate", "plasmidfinder", "out"):
        p.add_argument(f"--{name}", type=Path, required=True)
    args = p.parse_args(argv)
    lines = merge(args.assembly, args.fasta, args.bakta, args.amrfinderplus,
                  args.kleborate, args.plasmidfinder)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    n = sum(not line.startswith("#") for line in lines)
    print(f"{args.assembly}: {n} features written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
