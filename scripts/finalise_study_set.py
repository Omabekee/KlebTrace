#!/usr/bin/env python3
"""Finalise the study set using the carbapenemase re-check of portal-negative isolates.

query_portal.py selects six genomes from the AMR Portal's calls alone. The
re-check (screen_carbapenemases.py, confirmed with Kleborate) then showed that
both "no carbapenemase" picks carry a whole carbapenemase gene split across
contig breaks. This script:

1. relabels any no-carbapenemase pick whose screen call is intact or split as
   `portal_negative_carbapenemase_on_reannotation`, with its family filled in;
2. adds one `resistant_no_carbapenemase_confirmed` genome: strict resistant,
   no carbapenemase sequence in the screen and none in Kleborate, ranked by
   known region, then fewest sequences (a near-complete assembly is least
   likely to hide a split gene), then most recent and biosample;
3. rewrites data/portal_calls.tsv so it covers all the final genomes
   (skip with --no-portal-calls, for example when offline).

Usage:
    python scripts/finalise_study_set.py \\
        --kleborate results/kleborate_negatives/general__species_output.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

NO_CARBAPENEMASE_GROUPS = (
    "resistant_no_carbapenemase_porin_mutation",
    "resistant_no_carbapenemase_no_porin_mutation",
)
REANNOTATED_GROUP = "portal_negative_carbapenemase_on_reannotation"
CONFIRMED_GROUP = "resistant_no_carbapenemase_confirmed"
WHOLE_GENE_CALLS = ("intact", "split")
FAMILY_PREFIXES = {"NDM": "NDM", "KPC": "KPC", "OXA-": "OXA-48-like"}

KLEBORATE_CARB = "kpsc__amr__Bla_Carb_acquired"
KLEBORATE_PORIN = "kpsc__amr__Omp_mutations"
KLEBORATE_ESBL = "kpsc__amr__Bla_ESBL_acquired"
KLEBORATE_ST = "kpsc__mlst__ST"


def family_of(allele: str) -> str:
    for prefix, family in FAMILY_PREFIXES.items():
        if allele.startswith(prefix):
            return family
    return ""


def relabel_reannotated(samples: pd.DataFrame, screen: pd.DataFrame) -> pd.DataFrame:
    calls = screen.set_index("assembly")[["call", "allele"]]
    samples = samples.copy()
    for i, row in samples.iterrows():
        if row["selection_group"] not in NO_CARBAPENEMASE_GROUPS or row["assembly"] not in calls.index:
            continue
        call, allele = calls.loc[row["assembly"]]
        if call in WHOLE_GENE_CALLS:
            samples.loc[i, "selection_group"] = REANNOTATED_GROUP
            samples.loc[i, "carbapenemase_family"] = family_of(str(allele))
            samples.loc[i, "reannotated_carbapenemase"] = f"bla{allele} ({call})"
    return samples


def pick_confirmed_negative(
    screen: pd.DataFrame, kleborate: pd.DataFrame, exclude: set[str]
) -> pd.DataFrame:
    kl = kleborate.rename(columns={"strain": "assembly"})
    kl = kl[["assembly", KLEBORATE_ST, KLEBORATE_CARB, KLEBORATE_ESBL, KLEBORATE_PORIN]]
    candidates = screen[(screen["call"] == "none") & ~screen["assembly"].isin(exclude)]
    candidates = candidates.merge(kl, on="assembly", how="inner")
    carb = candidates[KLEBORATE_CARB].fillna("-").astype(str)
    candidates = candidates[carb.isin(["-", ""])].copy()
    candidates["unknown_region"] = candidates["region"].fillna("Unknown").str.lower().eq("unknown")
    ranked = candidates.sort_values(
        ["unknown_region", "n_sequences", "collection_year", "biosample"],
        ascending=[True, True, False, True],
        na_position="last",
    )
    return ranked.head(1).drop(columns="unknown_region")


def finalise(
    samples: pd.DataFrame, screen: pd.DataFrame, isolates: pd.DataFrame, kleborate: pd.DataFrame
) -> pd.DataFrame:
    out = relabel_reannotated(samples, screen)
    pick = pick_confirmed_negative(screen, kleborate, set(out["assembly"]))
    if pick.empty:
        print("Warning: no confirmed negative available.", file=sys.stderr)
        return out
    row = isolates[isolates["assembly"].isin(pick["assembly"])].copy()
    row["selection_group"] = CONFIRMED_GROUP
    row["carbapenemase_family"] = ""
    row["reannotated_carbapenemase"] = "none (screen and Kleborate)"
    row["kleborate_porin"] = pick[KLEBORATE_PORIN].iloc[0]
    return pd.concat([out, row[[c for c in row.columns if c in out.columns or c == "kleborate_porin"]]],
                     ignore_index=True)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--samples", type=Path, default=Path("data/samples.tsv"))
    p.add_argument("--isolates", type=Path, default=Path("data/portal_isolates.tsv"))
    p.add_argument("--screen", type=Path, default=Path("results/portal_negative_screen.tsv"))
    p.add_argument("--kleborate", type=Path, required=True,
                   help="Kleborate output table for the screen-negative genomes")
    p.add_argument("--out", type=Path, default=None, help="defaults to --samples (updated in place)")
    p.add_argument("--portal-calls", type=Path, default=Path("data/portal_calls.tsv"))
    p.add_argument("--no-portal-calls", action="store_true", help="do not refresh --portal-calls")
    args = p.parse_args(argv)

    read = lambda path: pd.read_csv(path, sep="\t", dtype={"collection_year": "Int64"})  # noqa: E731
    samples, screen, isolates = read(args.samples), read(args.screen), read(args.isolates)
    kleborate = pd.read_csv(args.kleborate, sep="\t", dtype=str)
    if REANNOTATED_GROUP in set(samples["selection_group"]) or CONFIRMED_GROUP in set(samples["selection_group"]):
        print(f"{args.samples} is already finalised; rerun query_portal.py first.", file=sys.stderr)
        return 1

    final = finalise(samples, screen, isolates, kleborate)
    final.to_csv(args.out or args.samples, sep="\t", index=False)
    if not args.no_portal_calls:
        import query_portal as qp

        geno = qp.RELEASE_URL.format(release=qp.DEFAULT_RELEASE) + "genotype.parquet"
        calls = qp.portal_calls(qp.connect([geno]), geno, final["assembly"].tolist())
        calls.to_csv(args.portal_calls, sep="\t", index=False)
        print(f"Portal calls for {calls['assembly_ID'].nunique()} genomes written to {args.portal_calls}")
    shown = ["selection_group", "carbapenemase_family", "assembly", "country", "collection_year",
             "mic_mg_l", "reannotated_carbapenemase"]
    print(final[[c for c in shown if c in final.columns]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
