#!/usr/bin/env python3
"""Rename the sequences in a KlebTrace GFF3 to the names Ensembl uses for the same assembly.

KlebTrace files use the ENA sequence accession (for example CP008933.1), while
Ensembl Bacteria uses other names for the same sequences: "Chromosome" and the
plasmid names for GCA_000764615.1 (pPMK1-NDM), and the submitter's contig names
for GCA_900774885.1 (ERS2179090SCcontig000001). A custom track only lines up in
the Ensembl browser if the names match.

Each Ensembl sequence is matched to the FASTA sequence of the same length. When
several sequences share a length, the one whose FASTA header mentions the
Ensembl name is used. Any sequence left unmatched stops the script.

Usage:
    python scripts/to_ensembl_seqids.py --fasta data/genomes/GCA_000764615.1.fasta \\
        --ensembl-gff3 data/ensembl/Klebsiella_pneumoniae_gca_000764615.ASM76461v1.63.gff3.gz \\
        --gff3 results/pipeline/merged/GCA_000764615.1.klebtrace.gff3 \\
        --out results/ensembl_tracks/GCA_000764615.1.ensembl_names.gff3
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path


def read_fasta_headers(fasta: Path) -> dict[str, tuple[int, str]]:
    """Return {accession: (length, full header)}."""
    out: dict[str, list] = {}
    current = None
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                current = line[1:].split()[0]
                out[current] = [0, line[1:].strip()]
            elif current:
                out[current][0] += len(line.strip())
    return {k: (v[0], v[1]) for k, v in out.items()}


def read_ensembl_regions(gff3: Path) -> dict[str, int]:
    """Return {Ensembl sequence name: length} from the ##sequence-region lines."""
    opener = gzip.open if str(gff3).endswith(".gz") else open
    regions = {}
    with opener(gff3, "rt") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                _, name, _start, end = line.split()
                regions[name] = int(end)
            elif not line.startswith("#"):
                break
    return regions


def mentions(header: str, name: str) -> bool:
    if name.lower() == "chromosome":
        return "plasmid" not in header.lower()
    return re.search(rf"(?<![\w-]){re.escape(name)}(?![\w-])", header) is not None


def match_sequences(fasta: dict[str, tuple[int, str]], ensembl: dict[str, int]) -> dict[str, str]:
    """Return {accession: Ensembl name}, one to one."""
    mapping, used = {}, set()
    for name, length in ensembl.items():
        same_length = [acc for acc, (n, _) in fasta.items() if n == length and acc not in used]
        if len(same_length) > 1:
            same_length = [acc for acc in same_length if mentions(fasta[acc][1], name)]
        if len(same_length) != 1:
            raise ValueError(f"Ensembl sequence {name} ({length} bp) matches {len(same_length)} FASTA sequences")
        mapping[same_length[0]] = name
        used.add(same_length[0])
    missing = sorted(set(fasta) - used)
    if missing:
        raise ValueError(f"FASTA sequences with no Ensembl match: {missing[:5]}")
    return mapping


def rename_gff3(lines: list[str], mapping: dict[str, str]) -> list[str]:
    out = []
    for line in lines:
        if line.startswith("##sequence-region"):
            _, seqid, start, end = line.split()
            out.append(f"##sequence-region {mapping[seqid]} {start} {end}")
        elif line.startswith("#") or not line.strip():
            out.append(line)
        else:
            seqid, rest = line.split("\t", 1)
            out.append(f"{mapping[seqid]}\t{rest}")
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for name in ("fasta", "ensembl-gff3", "gff3", "out"):
        p.add_argument(f"--{name}", type=Path, required=True)
    p.add_argument("--mapping", type=Path, help="also write the accession-to-Ensembl name table here")
    args = p.parse_args(argv)

    mapping = match_sequences(read_fasta_headers(args.fasta), read_ensembl_regions(args.ensembl_gff3))
    lines = rename_gff3(args.gff3.read_text().splitlines(), mapping)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    if args.mapping:
        args.mapping.write_text("accession\tensembl_name\n" + "".join(f"{a}\t{e}\n" for a, e in mapping.items()))
    print(f"{len(mapping)} sequences renamed; written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
