#!/usr/bin/env python3
"""Screen portal-negative isolates for carbapenemase genes at the DNA level.

Gene-based calls can miss a gene that an assembly has split across contig
breaks, because no single gene model covers it. This script maps every
carbapenemase allele in Kleborate's database directly onto each contig, joins
the pieces per allele, and reports whether each genome carries an intact gene,
a gene split into pieces, a partial match, or none.

A gene counts as split only if every piece runs into a contig end, as expected
when the assembler broke the sequence inside the gene. Pieces that together
cover the gene but sit inside contigs are called "interrupted" (for example, an
insertion inside the gene), since they are not explained by an assembly break.

Usage
-----
    python scripts/screen_carbapenemases.py --amr-db Kleborate_AMRdb.csv \\
        --amr-fasta Kleborate_AMRdb.fasta --minimap2 minimap2 \\
        --outdir data/screen_genomes --out results/portal_negative_screen.tsv
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

import download_genomes as dg

MIN_IDENTITY = 0.95
MIN_HIT_BP = 50
COMPLETE_COVERAGE = 0.95
CONTIG_END_BP = 10  # a piece within this distance of a contig end "runs into" it


def covered_bp(intervals: list[tuple[int, int]]) -> int:
    """Return the number of bases covered by possibly overlapping intervals."""
    total, end = 0, -1
    for start, stop in sorted(intervals):
        start = max(start, end)
        if stop > start:
            total += stop - start
            end = stop
    return total


def classify_allele(hits: list[dict], query_len: int) -> tuple[str, float]:
    """Classify one allele's hits in a genome as intact, split, partial or none."""
    if not hits:
        return "none", 0.0
    coverage = covered_bp([(h["qstart"], h["qend"]) for h in hits]) / query_len
    longest = max(h["qend"] - h["qstart"] for h in hits) / query_len
    if longest >= COMPLETE_COVERAGE:
        return "intact", coverage
    if coverage >= COMPLETE_COVERAGE:
        return ("split" if all(at_contig_end(h) for h in hits) else "interrupted"), coverage
    return "partial", coverage


def at_contig_end(hit: dict) -> bool:
    """True if the aligned piece starts or ends at (or within CONTIG_END_BP of) a contig end."""
    return hit["tstart"] <= CONTIG_END_BP or hit["tlen"] - hit["tend"] <= CONTIG_END_BP


def parse_paf(text: str) -> pd.DataFrame:
    """Parse minimap2 PAF lines, keeping confident hits."""
    rows = []
    for line in text.splitlines():
        f = line.split("\t")
        identity = int(f[9]) / int(f[10])
        hit_bp = int(f[3]) - int(f[2])
        if identity >= MIN_IDENTITY and hit_bp >= MIN_HIT_BP:
            rows.append(
                {
                    "allele": f[0].split("__")[2],
                    "query_len": int(f[1]),
                    "qstart": int(f[2]),
                    "qend": int(f[3]),
                    "contig": f[5],
                    "tlen": int(f[6]),
                    "tstart": int(f[7]),
                    "tend": int(f[8]),
                    "matches": int(f[9]),
                }
            )
    return pd.DataFrame(rows, columns=["allele", "query_len", "qstart", "qend", "contig",
                                       "tlen", "tstart", "tend", "matches"])


def best_call(hits: pd.DataFrame) -> dict[str, object]:
    """Return the best-supported allele for one genome: most matching bases, then fewest pieces."""
    if hits.empty:
        return {"call": "none", "allele": "", "coverage": 0.0, "n_pieces": 0, "contigs": ""}
    rank = {"intact": 0, "split": 1, "interrupted": 2, "partial": 3}
    calls = []
    for allele, group in hits.groupby("allele"):
        status, coverage = classify_allele(group.to_dict("records"), int(group["query_len"].iloc[0]))
        calls.append(
            {
                "call": status,
                "allele": allele,
                "coverage": round(coverage, 3),
                "n_pieces": len(group),
                "contigs": ";".join(sorted(group["contig"].unique())),
                "_matches": int(group["matches"].sum()),
            }
        )
    calls.sort(key=lambda c: (rank[c["call"]], -c["_matches"], c["n_pieces"], c["allele"]))
    best = calls[0]
    best.pop("_matches")
    return best


def carbapenemase_references(db_csv: Path, fasta: Path, out: Path) -> int:
    """Write the database's carbapenemase alleles to out and return how many were kept."""
    # Bla_Carb alone is too broad: at 95% identity, the common ESBL CTX-M-15 matches
    # CTX-M-33, and OXA-10 matches OXA-233, so both would be reported as carbapenemases.
    db = pd.read_csv(db_csv, dtype=str).fillna("")
    alleles = set(
        db.loc[
            (db["bla_class"] == "Bla_Carb")
            & db["bla_description"].str.contains("carbapenem-hydrolyzing")
            & (db["gene"] != "CTX-M"),
            "allele",
        ]
    )
    kept, keep = 0, False
    with open(fasta) as src, open(out, "w") as dst:
        for line in src:
            if line.startswith(">"):
                keep = line[1:].split("__")[2] in alleles
                kept += keep
            if keep:
                dst.write(line)
    return kept


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--isolates", default="data/portal_isolates.tsv", type=Path)
    p.add_argument("--amr-db", required=True, type=Path)
    p.add_argument("--amr-fasta", required=True, type=Path)
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--outdir", default="data/screen_genomes", type=Path)
    p.add_argument("--out", default="results/portal_negative_screen.tsv", type=Path)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    iso = pd.read_csv(args.isolates, sep="\t", dtype=str).fillna("")
    targets = iso[
        (iso["phenotype"] == "resistant")
        & (iso["n_carbapenemase_genes"] == "0")
        & (iso["strict_resistant"] == "True")
        & iso["assembly"].str.match(r"^GC[AF]_")
    ]
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    refs = args.outdir / "carbapenemase_refs.fasta"
    print(f"{carbapenemase_references(args.amr_db, args.amr_fasta, refs)} carbapenemase alleles", flush=True)

    rows = []
    for t in targets.to_dict("records"):
        assembly = t["assembly"]
        row = {k: t[k] for k in ("biosample", "assembly", "country", "collection_year", "region",
                                 "mic_mg_l", "carbapenem_point_mutations")}
        try:
            info = dg.download_genome(assembly, args.outdir, set())
            row.update(n_sequences=info["n_sequences"], n50=info["n50"], download=info["status"])
            paf = subprocess.run(
                [args.minimap2, "-c", "-k", "11", "-w", "3", "-N", "50", "--secondary=yes",
                 str(args.outdir / f"{assembly}.fasta"), str(refs)],
                capture_output=True, text=True, check=True,
            ).stdout
            row.update(best_call(parse_paf(paf)))
        except Exception as error:  # record the failure and continue with the other genomes
            row.update(download=f"ERROR: {error}", call="error")
        rows.append(row)
        print(assembly, row.get("call"), row.get("allele"), row.get("coverage"), row.get("n_pieces"), flush=True)
        pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
