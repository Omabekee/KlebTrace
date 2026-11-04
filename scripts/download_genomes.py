#!/usr/bin/env python3
"""Download the study genomes from ENA and check each one is complete.

Each assembly in samples.tsv is fetched as FASTA and compared with ENA's own
assembly record (sequence count and total length). Every contig named in the
portal's calls must also be present, so portal genes can later be located on
the same sequences. Headers are shortened from 'ENA|X|X.1 description' to
'X.1 description', which matches the portal's contig names and avoids '|' in
sequence IDs.

Usage
-----
    python scripts/download_genomes.py --samples data/samples.tsv \\
        --calls data/portal_calls.tsv --outdir data/genomes \\
        --manifest data/genomes_manifest.tsv
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ENA_API = "https://www.ebi.ac.uk/ena/browser/api"


def fetch(url: str, timeout: int = 300) -> bytes:
    """Return the body of a URL."""
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read()


def parse_assembly_xml(xml: bytes) -> dict[str, str | int]:
    """Return the level, title, sequence count and total length from an ENA assembly record.

    FASTA holds one sequence per scaffold, so scaffold-count is the expected
    number of sequences. count-contig is higher when scaffolds contain gaps.
    """
    root = ET.fromstring(xml)
    attributes = {
        a.findtext("TAG"): a.findtext("VALUE") for a in root.iter("ASSEMBLY_ATTRIBUTE")
    }
    return {
        "assembly_level": root.findtext(".//ASSEMBLY_LEVEL") or "",
        "title": root.findtext(".//TITLE") or "",
        "expected_sequences": int(attributes["scaffold-count"]),
        "expected_bp": int(attributes["total-length"]),
    }


def normalise_header(line: str) -> str:
    """Shorten '>ENA|X|X.1 description' to '>X.1 description'."""
    if line.startswith(">ENA|"):
        first, _, rest = line[1:].partition(" ")
        return f">{first.split('|')[-1]} {rest}".rstrip()
    return line


def normalise_fasta(text: str) -> tuple[str, dict[str, int], int]:
    """Return normalised FASTA text, contig lengths, and the G+C count."""
    lines: list[str] = []
    lengths: dict[str, int] = {}
    gc = 0
    contig = None
    for line in text.splitlines():
        if line.startswith(">"):
            line = normalise_header(line)
            contig = line[1:].split()[0]
            if contig in lengths:
                raise ValueError(f"Duplicate contig name: {contig}")
            lengths[contig] = 0
        elif contig is None:
            raise ValueError("Sequence found before the first FASTA header")
        else:
            lengths[contig] += len(line)
            gc += sum(line.upper().count(base) for base in "GC")
        lines.append(line)
    return "\n".join(lines) + "\n", lengths, gc


def n50(lengths: list[int]) -> int:
    """Return the contig length at which half the assembly is in contigs at least that long."""
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running * 2 >= sum(lengths):
            return length
    return 0


def check_genome(
    lengths: dict[str, int], record: dict[str, str | int], portal_contigs: set[str]
) -> list[str]:
    """Return a list of problems; an empty list means the genome passed."""
    problems = []
    if len(lengths) != record["expected_sequences"]:
        problems.append(f"{len(lengths)} sequences, ENA record says {record['expected_sequences']}")
    if sum(lengths.values()) != record["expected_bp"]:
        problems.append(f"{sum(lengths.values())} bp, ENA record says {record['expected_bp']}")
    missing = sorted(portal_contigs - set(lengths))
    if missing:
        problems.append(f"portal contigs missing from FASTA: {', '.join(missing)}")
    return problems


def download_genome(assembly: str, outdir: Path, portal_contigs: set[str]) -> dict[str, object]:
    """Download one genome, write it to outdir and return its manifest row."""
    record = parse_assembly_xml(fetch(f"{ENA_API}/xml/{assembly}"))
    raw = gzip.decompress(fetch(f"{ENA_API}/fasta/{assembly}?download=true&gzip=true"))
    fasta, lengths, gc = normalise_fasta(raw.decode())
    problems = check_genome(lengths, record, portal_contigs)
    path = outdir / f"{assembly}.fasta"
    path.write_text(fasta)
    total = sum(lengths.values())
    return {
        "assembly": assembly,
        "assembly_level": record["assembly_level"],
        "n_sequences": len(lengths),
        "total_bp": total,
        "n50": n50(list(lengths.values())),
        "gc_percent": round(100 * gc / total, 2) if total else 0,
        "md5": hashlib.md5(fasta.encode()).hexdigest(),
        "status": "PASS" if not problems else "FAIL: " + "; ".join(problems),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--samples", default="data/samples.tsv", type=Path)
    p.add_argument("--calls", default="data/portal_calls.tsv", type=Path)
    p.add_argument("--outdir", default="data/genomes", type=Path)
    p.add_argument("--manifest", default="data/genomes_manifest.tsv", type=Path)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    calls = pd.read_csv(args.calls, sep="\t", dtype=str)
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for assembly in samples["assembly"]:
        contigs = set(calls.loc[calls["assembly_ID"] == assembly, "region"].dropna())
        rows.append(download_genome(assembly, args.outdir, contigs))
        print(rows[-1])

    manifest = pd.DataFrame(rows)
    manifest.to_csv(args.manifest, sep="\t", index=False)
    failed = manifest[manifest["status"] != "PASS"]
    if not failed.empty:
        print(f"{len(failed)} genome(s) failed checks.", file=sys.stderr)
        return 1
    print(f"All {len(manifest)} genomes passed. Manifest written to {args.manifest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
