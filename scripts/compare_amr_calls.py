#!/usr/bin/env python3
"""Compare KlebTrace AMR calls with the AMR Portal calls and the laboratory result.

For each study genome the script lines up three sources of evidence:

    laboratory   meropenem phenotype (and MIC columns, if present) from samples.tsv
    portal       the AMR Portal's AMRFinderPlus calls, from portal_calls.tsv
    KlebTrace    this project's AMRFinderPlus calls, from amrfinderplus_all.tsv

and writes:

    comparison_carbapenem.tsv   one row per genome: carbapenem mechanisms from
                                each source, whether the two gene sets agree,
                                and whether the mechanisms found explain the
                                laboratory result
    comparison_all_genes.tsv    one row per genome and AMR gene, showing whether
                                the portal, KlebTrace or both reported it

Usage:
    python scripts/compare_amr_calls.py \
        --samples data/samples.tsv \
        --portal data/portal_calls.tsv \
        --klebtrace results/phase2_galaxy/amrfinderplus_all.tsv \
        --outdir results/comparison
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# Column names in each input. Portal names follow the AMR Portal genotype
# table; KlebTrace names follow the AMRFinderPlus report. The portal gene is
# amr_element_symbol: gene_symbol is the assembly's own label and can be
# generic (bla_2 for blaKPC-2), wrong or empty (see data/README.md).
PORTAL = {"assembly": "assembly_ID", "gene": "amr_element_symbol",
          "subclass": "subclass", "subtype": "element_subtype",
          "type": "element_type", "name": "reference_name",
          "contig": "region", "start": "region_start", "end": "region_end"}
KLEBTRACE = {"assembly": "assembly", "gene": "Element symbol",
             "subclass": "Subclass", "subtype": "Subtype",
             "type": "Type", "method": "Method", "name": "Closest reference name",
             "contig": "Contig id", "start": "Start", "end": "Stop"}

# An allele AMRFinderPlus cannot name (for example, a gene split across contig
# ends) is reported with its family symbol and family subclass, such as
# "blaOXA" with BETA-LACTAM. Its closest reference still names the enzyme as
# carbapenem-hydrolyzing, so carbapenemases are matched on either field.
CARBAPENEMASE_REFERENCE = "CARBAPENEM-HYDROLYZING"


def strip_version(accession: str) -> str:
    return re.sub(r"\.\d+$", "", str(accession).strip())


def gene_family(symbol: str) -> str:
    """blaKPC-2 -> blaKPC, blaOXA-48 -> blaOXA-48, ompK36_D135DGD -> ompK36.

    OXA enzymes keep their number, since OXA-48 and OXA-1 are unrelated
    enzymes with different activity.
    """
    symbol = str(symbol)
    if "_" in symbol:
        return symbol.split("_")[0]
    if symbol.startswith("blaOXA"):
        return symbol
    return re.sub(r"-\d+$", "", symbol)


def tidy(df: pd.DataFrame, cols: dict[str, str], source: str) -> pd.DataFrame:
    out = pd.DataFrame({
        "assembly": df[cols["assembly"]].map(strip_version),
        "gene": df[cols["gene"]].astype(str),
        "subclass": df[cols["subclass"]].fillna("").astype(str).str.upper(),
        "subtype": df[cols["subtype"]].fillna("").astype(str).str.upper(),
        "name": df[cols["name"]].fillna("").astype(str).str.upper() if cols["name"] in df else "",
    })
    out["type"] = df[cols["type"]].fillna("").astype(str).str.upper() if "type" in cols and cols["type"] in df else "AMR"
    out["partial"] = (
        df[cols["method"]].fillna("").astype(str).str.startswith("PARTIAL")
        if "method" in cols and cols["method"] in df else False
    )
    out["source"] = source
    # The portal repeats a call once per antibiotic it affects (split_subclass),
    # so keep one row per gene and position; true copies have different positions.
    position = [cols[c] for c in ("contig", "start", "end") if cols[c] in df]
    if position:
        out = out[~pd.concat([out[["assembly", "gene"]], df[position]], axis=1).duplicated()]
    return out


def mechanisms(calls: pd.DataFrame) -> tuple[list[str], list[str]]:
    carb = calls[calls["subclass"].str.contains("CARBAPENEM")
                 | calls["name"].str.contains(CARBAPENEMASE_REFERENCE)]
    genes = carb[~carb["subtype"].str.startswith("POINT")]
    points = carb[carb["subtype"].str.startswith("POINT")]

    def labels(rows: pd.DataFrame) -> list[str]:
        out = []
        for gene, group in rows.groupby("gene", sort=True):
            tag = gene + (f" (x{len(group)})" if len(group) > 1 else "")
            if group["partial"].any():
                tag += " [partial]"
            out.append(tag)
        return out

    return labels(genes), labels(points)


def explains(phenotype: str, genes: list[str], points: list[str]) -> str:
    """Does the genotype account for the laboratory result?

    Resistance is considered explained by an acquired carbapenemase, or by
    loss or alteration of both major porins (OmpK35 and OmpK36).
    """
    families = {gene_family(p.split(" ")[0]) for p in points}
    porin_pair = {"ompK35", "ompK36"} <= families
    predicted = "resistant" if genes or porin_pair else "susceptible"
    phenotype = str(phenotype).lower()
    if phenotype not in ("resistant", "susceptible"):
        return "no laboratory result"
    return "yes" if predicted == phenotype else f"no (genotype predicts {predicted})"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--samples", type=Path, default=Path("data/samples.tsv"))
    p.add_argument("--portal", type=Path, default=Path("data/portal_calls.tsv"))
    p.add_argument("--klebtrace", type=Path, default=Path("results/phase2_galaxy/amrfinderplus_all.tsv"))
    p.add_argument("--outdir", type=Path, default=Path("results/comparison"))
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    samples["key"] = samples["assembly"].map(strip_version)
    portal = tidy(pd.read_csv(args.portal, sep="\t", dtype=str), PORTAL, "portal")
    kleb = tidy(pd.read_csv(args.klebtrace, sep="\t", dtype=str), KLEBTRACE, "klebtrace")

    # Keep the laboratory columns that exist in this version of samples.tsv
    lab_cols = [c for c in samples.columns
                if c in ("phenotype", "selection_group", "country", "region")
                or re.search(r"mic|measurement|unit|standard", c, re.I)]

    rows = []
    for _, s in samples.iterrows():
        key = s["key"]
        pg, pp = mechanisms(portal[portal["assembly"] == key])
        annotated = (kleb["assembly"] == key).any()
        kg, kp = mechanisms(kleb[kleb["assembly"] == key])
        p_fam = {gene_family(g.split(" ")[0]) for g in pg}
        k_fam = {gene_family(g.split(" ")[0]) for g in kg}
        row = {"assembly": s["assembly"]}
        row.update({c: s[c] for c in lab_cols})
        row.update({
            "portal_carbapenemase": ";".join(pg) or "none",
            "klebtrace_carbapenemase": (";".join(kg) or "none") if annotated else "not annotated",
            "carbapenemase_agree": ("yes" if p_fam == k_fam else "no") if annotated else "not annotated",
            "portal_porin_changes": ";".join(pp) or "none",
            "klebtrace_porin_changes": (";".join(kp) or "none") if annotated else "not annotated",
            "klebtrace_explains_phenotype": explains(s.get("phenotype", ""), kg, kp) if annotated else "not annotated",
            "portal_explains_phenotype": explains(s.get("phenotype", ""), pg, pp),
            "n_portal_amr_calls": int((portal["assembly"] == key).sum()),
            "n_klebtrace_amr_calls": int(((kleb["assembly"] == key) & (kleb["type"] == "AMR")).sum()),
        })
        rows.append(row)
    pd.DataFrame(rows).to_csv(args.outdir / "comparison_carbapenem.tsv", sep="\t", index=False)

    # Gene-level presence/absence, AMR genes only
    k_amr = kleb[kleb["type"] == "AMR"]
    both = pd.concat([portal, k_amr])
    both = both[both["assembly"].isin(set(samples["key"]))]
    table = (both.groupby(["assembly", "gene", "source"]).size()
             .unstack("source", fill_value=0).reset_index())
    for col in ("portal", "klebtrace"):
        if col not in table:
            table[col] = 0
    table["found_by"] = table.apply(
        lambda r: "both" if r["portal"] and r["klebtrace"]
        else ("portal only" if r["portal"] else "KlebTrace only"), axis=1)
    table = table.rename(columns={"portal": "portal_rows", "klebtrace": "klebtrace_rows"})
    table.to_csv(args.outdir / "comparison_all_genes.tsv", sep="\t", index=False)

    print(f"Compared {len(rows)} genomes; outputs in {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
