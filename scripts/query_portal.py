#!/usr/bin/env python3
"""Query the EMBL-EBI AMR Portal Parquet release for KlebTrace.

Steps
-----
1. Join laboratory phenotypes (phenotype.parquet) to genome-based resistance
   calls (genotype.parquet) for one species and one antibiotic.
2. Build one row per isolate with an assembly: its meropenem phenotype,
   metadata and any carbapenem-class resistance genes the portal reports.
3. Summarise isolates by phenotype, carbapenemase status and world region.
4. Select a small, deliberately mixed study set:
     - resistant isolates with a carbapenemase gene
     - resistant isolates with no carbapenemase gene, with and without a porin mutation
     - susceptible isolates
   spreading the picks across world regions where possible.
5. Extract every portal genotype call for the selected assemblies, for later
   comparison with KlebTrace's own annotation.

DuckDB reads the Parquet files in place, including over HTTPS, so the full
release does not need to be downloaded.

Usage
-----
    python scripts/query_portal.py --inspect
    python scripts/query_portal.py --outdir data/
    python scripts/query_portal.py --pheno local/phenotype.parquet \\
        --geno local/genotype.parquet --outdir data/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import duckdb
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_RELEASE = "2026-07"  # check the releases folder for the latest
RELEASE_URL = "https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/{release}/"

# Column names follow the AMR Portal developer documentation. Run --inspect
# first and correct any name here if the release uses a different spelling.
COLUMNS = {
    "biosample": "BioSample_ID",
    "assembly": "assembly_ID",
    "species": "species",
    "antibiotic": "antibiotic_name",
    "phenotype": "resistance_phenotype",
    "country": "country",
    "year": "collection_year",
    "region": "geographical_region",
    # amr_element_symbol is the AMRFinderPlus call. gene_symbol is the assembly's
    # own gene label, which can be generic (bla_2), wrong (blaNDM-1 for blaNDM-5) or empty.
    "gene": "amr_element_symbol",
    "subclass": "subclass",
    "subtype": "element_subtype",
    "measurement": "measurement",
    "sign": "measurement_sign",
    "units": "measurement_units",
    "clsi": "Updated_phenotype_CLSI",
    "eucast": "Updated_phenotype_EUCAST",
}

# One study genome is selected per carbapenemase family (regex on the gene name).
CARBAPENEMASE_FAMILIES = {
    "NDM": r"^blaNDM",
    "KPC": r"^blaKPC",
    "OXA-48-like": r"^blaOXA-(?:48|162|181|204|232|244|245|370|484)$",
}

# MIC (mg/L) that a resistant isolate with no carbapenemase must exceed, so that
# group contains only isolates resistant under both CLSI and EUCAST breakpoints.
STRICT_MIC_ABOVE = 8

RESISTANT_VALUES = ("resistant", "r")
SUSCEPTIBLE_VALUES = ("susceptible", "s")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def q(name: str) -> str:
    """Quote a column identifier for SQL."""
    return '"' + name.replace('"', '""') + '"'


def connect(paths: list[str]) -> duckdb.DuckDBPyConnection:
    """Open DuckDB, enabling HTTPS reads when any input is a URL."""
    con = duckdb.connect()
    if any(p.startswith(("http://", "https://")) for p in paths):
        con.execute("INSTALL httpfs; LOAD httpfs;")
    return con


def inspect(con: duckdb.DuckDBPyConnection, path: str) -> pd.DataFrame:
    """Return the column names and types of a Parquet file."""
    return con.execute(
        "DESCRIBE SELECT * FROM read_parquet(?)", [path]
    ).df()


# ---------------------------------------------------------------------------
# Core queries
# ---------------------------------------------------------------------------

def isolate_table(
    con: duckdb.DuckDBPyConnection,
    pheno: str,
    geno: str,
    species: str,
    antibiotic: str,
    subclass: str,
    cols: dict[str, str] = COLUMNS,
) -> pd.DataFrame:
    """One row per isolate with an assembly, its phenotype and carbapenem calls.

    An isolate is called resistant if any test for the antibiotic was
    resistant, and susceptible only if every test was susceptible.
    Intermediate or missing results are excluded.

    Carbapenem-class calls are split in two. Acquired genes (AMRFinderPlus
    element subtype other than POINT), such as blaKPC or blaNDM, are reported
    as carbapenemase genes. Point mutations (subtype POINT), such as changes in
    the porin genes ompK35, ompK36 or ompC, are reported separately, because
    they reduce carbapenem entry into the cell but are not carbapenemases.
    """
    c = {k: q(v) for k, v in cols.items()}
    sql = f"""
    WITH ph AS (
        SELECT
            {c['biosample']}                       AS biosample,
            any_value({c['country']})              AS country,
            any_value({c['year']})                 AS collection_year,
            coalesce(any_value({c['region']}), 'Unknown') AS region,
            max(CASE WHEN lower({c['phenotype']}) IN {RESISTANT_VALUES}
                     THEN 1 ELSE 0 END)            AS resistant,
            string_agg(DISTINCT coalesce({c['sign']}, '') || {c['measurement']}, ';')
                FILTER (WHERE lower({c['units']}) = 'mg/l') AS mic_mg_l,
            -- Every mg/L test must exceed the threshold: '<=16' does not show this, '>8' does.
            coalesce(
                bool_and(
                    CASE WHEN {c['sign']} = '>'
                         THEN try_cast({c['measurement']} AS DOUBLE) >= {STRICT_MIC_ABOVE}
                         WHEN coalesce({c['sign']}, '=') IN ('=', '==', '>=')
                         THEN try_cast({c['measurement']} AS DOUBLE) > {STRICT_MIC_ABOVE}
                         ELSE false END
                ) FILTER (WHERE lower({c['units']}) = 'mg/l'),
                false
            )
            AND bool_and(lower(coalesce({c['clsi']}, '')) = 'resistant'
                         AND lower(coalesce({c['eucast']}, '')) = 'resistant') AS strict_resistant
        FROM read_parquet($pheno)
        WHERE lower({c['species']}) LIKE lower($species) || '%'
          AND lower({c['antibiotic']}) = lower($antibiotic)
          AND lower({c['phenotype']}) IN {RESISTANT_VALUES + SUSCEPTIBLE_VALUES}
        GROUP BY 1
    ),
    asm AS (
        SELECT DISTINCT
            {c['biosample']} AS biosample,
            {c['assembly']}  AS assembly
        FROM read_parquet($geno)
    ),
    carb AS (
        SELECT
            {c['assembly']} AS assembly,
            string_agg(DISTINCT {c['gene']}, ';' ORDER BY {c['gene']})
                FILTER (WHERE upper(coalesce({c['subtype']}, '')) <> 'POINT') AS carbapenemase_genes,
            count(DISTINCT {c['gene']})
                FILTER (WHERE upper(coalesce({c['subtype']}, '')) <> 'POINT') AS n_carbapenemase_genes,
            string_agg(DISTINCT {c['gene']}, ';' ORDER BY {c['gene']})
                FILTER (WHERE upper(coalesce({c['subtype']}, '')) = 'POINT')  AS carbapenem_point_mutations
        FROM read_parquet($geno)
        WHERE upper({c['subclass']}) LIKE '%' || upper($subclass) || '%'
        GROUP BY 1
    )
    SELECT
        ph.biosample,
        asm.assembly,
        ph.country,
        ph.collection_year,
        ph.region,
        CASE WHEN ph.resistant = 1 THEN 'resistant' ELSE 'susceptible' END AS phenotype,
        ph.mic_mg_l,
        ph.strict_resistant,
        coalesce(carb.carbapenemase_genes, '')          AS carbapenemase_genes,
        coalesce(carb.n_carbapenemase_genes, 0)         AS n_carbapenemase_genes,
        coalesce(carb.carbapenem_point_mutations, '')   AS carbapenem_point_mutations
    FROM ph
    JOIN asm USING (biosample)
    LEFT JOIN carb USING (assembly)
    ORDER BY ph.biosample, asm.assembly
    """
    params = {
        "pheno": pheno,
        "geno": geno,
        "species": species,
        "antibiotic": antibiotic,
        "subclass": subclass,
    }
    return con.execute(sql, params).df()


def _strip_version(accession: str) -> str:
    """Return an assembly accession without its version suffix (GCA_1.2 -> GCA_1)."""
    return str(accession).split(".")[0]


def load_ensembl_accessions(path: str, column: str = "assembly_accession") -> set[str]:
    """Read assembly accessions from an Ensembl Bacteria species list.

    The species list is a tab-separated file whose header line may start with
    '#'. Accessions are compared without version suffixes, so GCA_1.1 and
    GCA_1.2 both match GCA_1.
    """
    table = pd.read_csv(path, sep="\t", dtype=str, index_col=False)
    table.columns = [c.lstrip("#").strip() for c in table.columns]
    if column not in table.columns:
        raise ValueError(
            f"Column '{column}' not found in {path}. Available: {', '.join(table.columns)}"
        )
    return {_strip_version(a) for a in table[column].dropna()}


def flag_in_ensembl(isolates: pd.DataFrame, ensembl_accessions: set[str]) -> pd.DataFrame:
    """Add an in_ensembl column marking isolates whose assembly is in Ensembl."""
    return isolates.assign(
        in_ensembl=isolates["assembly"].map(
            lambda a: _strip_version(a) in ensembl_accessions if pd.notna(a) else False
        )
    )


def summarise(isolates: pd.DataFrame) -> pd.DataFrame:
    """Count isolates by region, phenotype and carbapenemase status."""
    df = isolates.assign(
        carbapenemase=lambda d: d["n_carbapenemase_genes"].gt(0).map(
            {True: "present", False: "absent"}
        )
    )
    table = (
        df.groupby(["region", "phenotype", "carbapenemase"])
        .size()
        .rename("n_isolates")
        .reset_index()
        .sort_values(["region", "phenotype", "carbapenemase"])
    )
    return table


def _pick(
    group: pd.DataFrame,
    n: int,
    prefer_regions: list[str] | None = None,
    used_regions: set[str] | None = None,
) -> pd.DataFrame:
    """Pick n isolates, updating used_regions so later picks spread across regions.

    Isolates are ranked by: present in Ensembl (so the genome can be displayed
    against Ensembl's own annotation), known region, region not yet used,
    preferred region, then most recent.
    """
    group = group.copy()
    if "in_ensembl" not in group.columns:
        group["in_ensembl"] = False
    used = used_regions if used_regions is not None else set()
    preferred = {r.lower() for r in (prefer_regions or [])}
    group["not_preferred"] = ~group["region"].str.lower().isin(preferred)
    group["unknown_region"] = group["region"].str.lower().eq("unknown")
    picks = []
    for _ in range(min(n, len(group))):
        group["region_used"] = group["region"].isin(used)
        ranked = group.sort_values(
            ["in_ensembl", "unknown_region", "region_used", "not_preferred",
             "collection_year", "biosample"],
            ascending=[False, True, True, True, False, True],
            na_position="last",
        )
        top = ranked.iloc[[0]]
        picks.append(top)
        used.add(top["region"].iloc[0])
        group = group.drop(top.index)
    helpers = ["not_preferred", "unknown_region", "region_used"]
    return pd.concat(picks).drop(columns=helpers) if picks else group.drop(columns=helpers, errors="ignore")


def select_study_set(
    isolates: pd.DataFrame,
    carbapenemase_families: dict[str, str] | None = None,
    n_susceptible: int = 1,
    prefer_regions: list[str] | None = None,
    require_gca: bool = True,
) -> pd.DataFrame:
    """Select a deliberately mixed study set of isolates.

    Resistant isolates with a carbapenemase contribute one genome per family,
    each with only that one carbapenemase gene so the phenotype can be attributed
    to it. Two further resistant isolates have no carbapenemase, one with a
    carbapenem porin point mutation (from a preferred region when one qualifies)
    and one without; both must meet the strict_resistant rule. The susceptible
    control must have no carbapenemase gene.

    By default, only isolates with a GenBank or RefSeq assembly accession
    (GCA_ or GCF_) are eligible, because these can be downloaded as complete
    assemblies from ENA and NCBI and can be matched to Ensembl.
    """
    families = CARBAPENEMASE_FAMILIES if carbapenemase_families is None else carbapenemase_families
    if require_gca:
        isolates = isolates[isolates["assembly"].astype(str).str.match(r"^GC[AF]_")]
    has_gene = isolates["n_carbapenemase_genes"] > 0
    single_gene = isolates["n_carbapenemase_genes"] == 1
    resistant = isolates["phenotype"] == "resistant"

    groups = [
        ("resistant_with_carbapenemase", family,
         isolates[resistant & single_gene
                  & isolates["carbapenemase_genes"].str.contains(pattern, regex=True)], 1)
        for family, pattern in families.items()
    ]
    no_gene = isolates[resistant & ~has_gene & isolates["strict_resistant"].astype(bool)]
    has_porin = no_gene["carbapenem_point_mutations"].fillna("").ne("")
    porin = no_gene[has_porin]
    in_preferred = porin["region"].str.lower().isin({r.lower() for r in (prefer_regions or [])})
    groups += [
        ("resistant_no_carbapenemase_porin_mutation", "", porin[in_preferred] if in_preferred.any() else porin, 1),
        ("resistant_no_carbapenemase_no_porin_mutation", "", no_gene[~has_porin], 1),
        ("susceptible", "", isolates[~resistant & ~has_gene], n_susceptible),
    ]

    picks = []
    used_regions: set[str] = set()
    for label, family, group, n in groups:
        chosen = _pick(group, n, prefer_regions, used_regions).assign(
            selection_group=label, carbapenemase_family=family
        )
        if len(chosen) < n:
            print(
                f"Warning: only {len(chosen)} of {n} isolates available for '{label}' {family}.",
                file=sys.stderr,
            )
        picks.append(chosen)
    return pd.concat(picks, ignore_index=True)


def portal_calls(
    con: duckdb.DuckDBPyConnection,
    geno: str,
    assemblies: list[str],
    cols: dict[str, str] = COLUMNS,
) -> pd.DataFrame:
    """Return every portal genotype call for the given assemblies."""
    if not assemblies:
        return pd.DataFrame()
    sql = f"""
        SELECT *
        FROM read_parquet($geno)
        WHERE {q(cols['assembly'])} IN (SELECT unnest($assemblies))
        ORDER BY {q(cols['assembly'])}, {q(cols['gene'])}
    """
    return con.execute(sql, {"geno": geno, "assemblies": assemblies}).df()


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--release", default=DEFAULT_RELEASE, help="AMR Portal release folder, e.g. 2025-11")
    p.add_argument("--pheno", help="Path or URL to phenotype.parquet (overrides --release)")
    p.add_argument("--geno", help="Path or URL to genotype.parquet (overrides --release)")
    p.add_argument("--species", default="Klebsiella pneumoniae")
    p.add_argument("--antibiotic", default="meropenem")
    p.add_argument("--subclass", default="CARBAPENEM", help="AMRFinderPlus subclass marking carbapenemase genes")
    p.add_argument(
        "--ensembl-list",
        help="Path or URL to the Ensembl Bacteria species list (species_EnsemblBacteria.txt); "
        "isolates whose assembly is in Ensembl are flagged and preferred",
    )
    p.add_argument("--ensembl-column", default="assembly_accession",
                   help="Column in the species list holding assembly accessions")
    p.add_argument("--prefer-regions", nargs="*", default=[],
                   help="World regions to prioritise when selecting isolates, e.g. Africa")
    p.add_argument("--allow-non-gca", action="store_true",
                   help="Allow isolates whose assembly is not a GCA_/GCF_ accession")
    p.add_argument("--outdir", default="data", type=Path)
    p.add_argument("--inspect", action="store_true", help="Print the Parquet schemas and exit")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    base = RELEASE_URL.format(release=args.release)
    pheno = args.pheno or base + "phenotype.parquet"
    geno = args.geno or base + "genotype.parquet"
    con = connect([pheno, geno])

    if args.inspect:
        for label, path in (("phenotype", pheno), ("genotype", geno)):
            print(f"\n=== {label}: {path}")
            print(inspect(con, path).to_string(index=False))
        return 0

    isolates = isolate_table(con, pheno, geno, args.species, args.antibiotic, args.subclass)
    if isolates.empty:
        print("No isolates matched. Run --inspect and check COLUMNS and filter values.", file=sys.stderr)
        return 1

    if args.ensembl_list:
        accessions = load_ensembl_accessions(args.ensembl_list, args.ensembl_column)
        isolates = flag_in_ensembl(isolates, accessions)
        print(f"Isolates whose assembly is also in Ensembl Bacteria: {int(isolates['in_ensembl'].sum())}")

    summary = summarise(isolates)
    selected = select_study_set(
        isolates, prefer_regions=args.prefer_regions, require_gca=not args.allow_non_gca
    )
    calls = portal_calls(con, geno, selected["assembly"].dropna().tolist())

    args.outdir.mkdir(parents=True, exist_ok=True)
    isolates.to_csv(args.outdir / "portal_isolates.tsv", sep="\t", index=False)
    summary.to_csv(args.outdir / "portal_summary.tsv", sep="\t", index=False)
    selected.to_csv(args.outdir / "samples.tsv", sep="\t", index=False)
    calls.to_csv(args.outdir / "portal_calls.tsv", sep="\t", index=False)

    print(f"Isolates with an assembly and a {args.antibiotic} result: {len(isolates)}")
    print("\nSummary by region, phenotype and carbapenemase status:")
    print(summary.to_string(index=False))
    print("\nSelected study set:")
    shown = ["selection_group", "carbapenemase_family", "biosample", "assembly", "region", "collection_year",
             "phenotype", "mic_mg_l", "carbapenemase_genes", "carbapenem_point_mutations"]
    if "in_ensembl" in selected.columns:
        shown.append("in_ensembl")
    print(selected[shown].to_string(index=False))
    print(f"\nOutputs written to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
