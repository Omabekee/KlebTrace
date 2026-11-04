#!/usr/bin/env python3
"""Score each resistance-call source against the laboratory meropenem result.

Reads the KlebTrace database and, for every genome and source (portal,
AMRFinderPlus, Kleborate, split-gene screen), records whether a carbapenemase
was found and whether the calls explain the laboratory result. Resistance
counts as explained by an acquired carbapenemase, or by changes in both OmpK35
and OmpK36 (the rule used in compare_amr_calls.py); a susceptible isolate is
explained when neither is found.

The split-gene screen only looks for carbapenemases and was run only on the
portal-negative isolates, so it is scored on those genomes alone.

Usage:
    python scripts/assess_accuracy.py --database results/pipeline/database/klebtrace.sqlite \\
        --out results/tool_assessment/accuracy.tsv
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path

import pandas as pd

SOURCES = ["portal", "amrfinderplus", "kleborate", "split_gene_screen"]


def porins(genes: list[str]) -> set[str]:
    """Porin genes with a change, from either naming style (ompK36_D135DGD, ompK36:p.134_135insGD)."""
    return {g.split("_")[0].split(":")[0] for g in genes if g.lower().startswith(("ompk35", "ompk36"))}


def score(con: sqlite3.Connection) -> pd.DataFrame:
    genomes = pd.read_sql(
        "SELECT g.assembly, g.selection_group, p.phenotype, p.mic_mg_l FROM genomes g "
        "JOIN phenotypes p ON p.assembly = g.assembly AND p.antibiotic = 'meropenem'", con)
    calls = pd.read_sql("SELECT assembly, source, gene, is_carbapenemase, element_subtype FROM resistance_calls", con)
    screened = set(calls.loc[calls["source"] == "split_gene_screen", "assembly"]) | set(
        genomes.loc[genomes["selection_group"].str.contains("no_carbapenemase|portal_negative"), "assembly"])
    rows = []
    for g in genomes.itertuples():
        for source in SOURCES:
            if source == "split_gene_screen" and g.assembly not in screened:
                continue
            c = calls[(calls["assembly"] == g.assembly) & (calls["source"] == source)]
            carb = sorted(c.loc[c["is_carbapenemase"] == 1, "gene"])
            changed = porins(c.loc[c["element_subtype"].str.startswith("POINT", na=False), "gene"].tolist())
            predicted = "resistant" if carb or changed == {"ompK35", "ompK36"} else "susceptible"
            rows.append({
                "assembly": g.assembly, "phenotype": g.phenotype, "mic_mg_l": g.mic_mg_l, "source": source,
                "carbapenemases": ";".join(carb) or "none",
                "porins_changed": ";".join(sorted(changed)) or "none",
                "explains_phenotype": "n/a (carbapenemases only)" if source == "split_gene_screen" and not carb
                else ("yes" if predicted == g.phenotype else "no"),
            })
    return pd.DataFrame(rows)


def summarise(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source, group in table.groupby("source", sort=False):
        scored = group[group["explains_phenotype"].isin(["yes", "no"])]
        rows.append({"source": source, "genomes": len(group),
                     "carbapenemase_found": int((group["carbapenemases"] != "none").sum()),
                     "explains_phenotype": f"{int((scored['explains_phenotype'] == 'yes').sum())} of {len(scored)}"})
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--database", type=Path, default=Path("results/pipeline/database/klebtrace.sqlite"))
    p.add_argument("--out", type=Path, default=Path("results/tool_assessment/accuracy.tsv"))
    args = p.parse_args(argv)
    with sqlite3.connect(args.database) as con:
        table = score(con)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False)
    s = summarise(table)
    s.to_csv(args.out.with_name(args.out.stem + "_summary.tsv"), sep="\t", index=False)
    print(table.to_string(index=False))
    print()
    print(s.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
