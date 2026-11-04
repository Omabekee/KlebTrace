"""Unit tests for scripts/query_portal.py using small mock Parquet files."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import query_portal as qp  # noqa: E402


@pytest.fixture()
def parquet_files(tmp_path: Path) -> tuple[str, str]:
    """Write mock phenotype and genotype Parquet files."""
    pheno = pd.DataFrame(
        {
            "BioSample_ID": ["S1", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"],
            "species": ["Klebsiella pneumoniae"] * 6
            + ["Klebsiella pneumoniae subsp. pneumoniae", "Escherichia coli", "Klebsiella pneumoniae"],
            "antibiotic_name": ["meropenem", "Meropenem", "meropenem", "meropenem",
                                "meropenem", "meropenem", "meropenem", "meropenem", "meropenem"],
            "resistance_phenotype": ["resistant", "susceptible", "resistant", "resistant",
                                     "susceptible", "intermediate", "resistant", "resistant",
                                     "resistant"],
            "country": ["Nigeria", "Nigeria", "UK", "India", "Brazil", "UK", "Kenya", "UK", "Ghana"],
            "collection_year": [2020, 2020, 2021, 2019, 2022, 2021, 2023, 2020, 2024],
            "geographical_region": ["Africa", "Africa", "Europe", "Asia", "Americas",
                                    "Europe", "Africa", "Europe", "Africa"],
        }
    )
    geno = pd.DataFrame(
        {
            "BioSample_ID": ["S1", "S1", "S2", "S2", "S3", "S4", "S6", "S7", "S8"],
            "assembly_ID": ["GCA_1", "GCA_1", "GCA_2", "GCA_2", "GCA_3", "GCA_4", "GCA_6",
                            "GCA_7", "ERZ8"],
            "amr_element_symbol": ["blaKPC-2", "blaSHV-11", "blaNDM-1", "ompC", "ompK36", "blaSHV-1",
                            "blaOXA-48", "blaKPC-3", "blaVIM-1"],
            "subclass": ["CARBAPENEM", "BETA-LACTAM", "CARBAPENEM", "CARBAPENEM", "CARBAPENEM",
                         "BETA-LACTAM", "CARBAPENEM", "CARBAPENEM", "CARBAPENEM"],
            "element_subtype": ["AMR", "AMR", "AMR", "POINT", "POINT", "AMR", "AMR", "AMR", "AMR"],
        }
    )
    pheno["measurement"] = ["16", "16", "16", "16", "1", "1", "16", "16", "16"]
    pheno["measurement_sign"] = "="
    pheno["measurement_units"] = "mg/L"
    pheno["Updated_phenotype_CLSI"] = pheno["resistance_phenotype"]
    pheno["Updated_phenotype_EUCAST"] = pheno["resistance_phenotype"]
    p = tmp_path / "phenotype.parquet"
    g = tmp_path / "genotype.parquet"
    pheno.to_parquet(p)
    geno.to_parquet(g)
    return str(p), str(g)


def build(parquet_files):
    pheno, geno = parquet_files
    con = qp.connect([pheno, geno])
    iso = qp.isolate_table(con, pheno, geno, "Klebsiella pneumoniae", "meropenem", "CARBAPENEM")
    return con, iso


def test_isolate_table_filters_and_calls_phenotypes(parquet_files):
    _, iso = build(parquet_files)
    iso = iso.set_index("biosample")
    # E. coli (S7) excluded; S5 (intermediate only) excluded; subspecies name kept
    assert set(iso.index) == {"S1", "S2", "S3", "S4", "S6", "S8"}
    # S1 has one resistant and one susceptible test, so it is resistant
    assert iso.loc["S1", "phenotype"] == "resistant"
    assert iso.loc["S4", "phenotype"] == "susceptible"


def test_carbapenemase_genes_exclude_point_mutations(parquet_files):
    _, iso = build(parquet_files)
    iso = iso.set_index("biosample")
    assert iso.loc["S1", "carbapenemase_genes"] == "blaKPC-2"
    assert iso.loc["S6", "carbapenemase_genes"] == "blaOXA-48"
    # S2: ompC point mutation is reported separately, not as a carbapenemase
    assert iso.loc["S2", "carbapenemase_genes"] == "blaNDM-1"
    assert iso.loc["S2", "carbapenem_point_mutations"] == "ompC"
    # S3: porin mutation only, so no carbapenemase
    assert iso.loc["S3", "n_carbapenemase_genes"] == 0
    assert iso.loc["S3", "carbapenem_point_mutations"] == "ompK36"


def test_non_gca_assemblies_excluded_by_default(parquet_files):
    _, iso = build(parquet_files)
    sel = qp.select_study_set(iso)
    assert "ERZ8" not in set(sel["assembly"])
    sel_all = qp.select_study_set(iso, {"VIM": "^blaVIM"}, require_gca=False)
    assert sel_all.loc[sel_all["selection_group"] == "resistant_with_carbapenemase",
                       "assembly"].tolist() == ["ERZ8"]


def test_prefer_regions(parquet_files):
    _, iso = build(parquet_files)
    sel = qp.select_study_set(iso, {"any": "^bla"}, prefer_regions=["Africa"])
    first = sel[sel["selection_group"] == "resistant_with_carbapenemase"]
    assert first["region"].tolist() == ["Africa"]


def test_select_study_set_groups(parquet_files):
    _, iso = build(parquet_files)
    sel = qp.select_study_set(iso)
    counts = sel["selection_group"].value_counts().to_dict()
    assert counts == {
        "resistant_with_carbapenemase": 3,
        "resistant_no_carbapenemase_porin_mutation": 1,
        "susceptible": 1,
    }
    without = sel.loc[sel["selection_group"] == "resistant_no_carbapenemase_porin_mutation", "biosample"]
    assert without.tolist() == ["S3"]


def test_one_genome_per_carbapenemase_family(parquet_files):
    _, iso = build(parquet_files)
    sel = qp.select_study_set(iso)
    with_gene = sel[sel["selection_group"] == "resistant_with_carbapenemase"]
    assert dict(zip(with_gene["carbapenemase_family"], with_gene["biosample"])) == {
        "NDM": "S2", "KPC": "S1", "OXA-48-like": "S6",
    }


def test_susceptible_control_has_no_carbapenemase():
    iso = pd.DataFrame(
        {
            "biosample": ["A", "B", "C"],
            "assembly": ["GCA_A", "GCA_B", "GCA_C"],
            "collection_year": [2024, 2020, 2020],
            "region": ["Africa", "Asia", "Asia"],
            "phenotype": ["susceptible", "susceptible", "resistant"],
            "carbapenemase_genes": ["blaOXA-181", "", ""],
            "n_carbapenemase_genes": [1, 0, 0],
            "strict_resistant": [False, False, False],
            "carbapenem_point_mutations": ["", "", ""],
        }
    )
    sel = qp.select_study_set(iso, {})
    assert sel.loc[sel["selection_group"] == "susceptible", "biosample"].tolist() == ["B"]


def test_no_carbapenemase_subgroups_are_strict_and_prefer_region():
    iso = pd.DataFrame(
        {
            "biosample": ["KE", "BR", "IT", "RU", "XX"],
            "assembly": ["GCA_KE", "GCA_BR", "GCA_IT", "GCA_RU", "GCA_XX"],
            "collection_year": [2015, 2023, 2011, 2017, 2024],
            "region": ["Africa", "Americas", "Europe", "Europe", "Asia"],
            "phenotype": ["resistant"] * 5,
            "carbapenemase_genes": [""] * 5,
            "n_carbapenemase_genes": [0] * 5,
            "strict_resistant": [True, True, True, False, False],
            "carbapenem_point_mutations": ["ompK36_D135DGD", "ompK35_E132K", "", "", "ompK36_G133D"],
        }
    )
    sel = qp.select_study_set(iso, {}, prefer_regions=["Africa"]).set_index("selection_group")
    assert sel.loc["resistant_no_carbapenemase_porin_mutation", "biosample"] == "KE"
    assert sel.loc["resistant_no_carbapenemase_no_porin_mutation", "biosample"] == "IT"


def test_strict_resistant_uses_sign_units_and_both_standards(tmp_path):
    #        id, measurement, sign, units, CLSI, EUCAST, expected
    cases = [
        ("exact", "16", "=", "mg/L", "resistant", "resistant", True),
        ("greater", "8", ">", "mg/L", "resistant", "resistant", True),
        ("at_most", "16", "<=", "mg/L", "resistant", "resistant", False),
        ("at_least_8", "8", ">=", "mg/L", "resistant", "resistant", False),
        ("zone", "25", "=", "mm", "resistant", "resistant", False),
        ("eucast_int", "16", "=", "mg/L", "resistant", "intermediate", False),
        ("no_update", "16", "=", "mg/L", None, None, False),
        ("no_mic", None, None, None, "resistant", "resistant", False),
    ]
    rows = [c[:6] for c in cases] + [("repeat", "16", "=", "mg/L", "resistant", "resistant"),
                                     ("repeat", "1", "=", "mg/L", "resistant", "resistant")]
    pheno = pd.DataFrame(rows, columns=["BioSample_ID", "measurement", "measurement_sign", "measurement_units",
                                        "Updated_phenotype_CLSI", "Updated_phenotype_EUCAST"])
    pheno = pheno.assign(species="Klebsiella pneumoniae", antibiotic_name="meropenem",
                         resistance_phenotype="resistant", country="X", collection_year=2020,
                         geographical_region="Africa")
    ids = pheno["BioSample_ID"].unique()
    geno = pd.DataFrame({"BioSample_ID": ids, "assembly_ID": ["GCA_" + i for i in ids],
                         "amr_element_symbol": "blaSHV-1", "subclass": "BETA-LACTAM", "element_subtype": "AMR"})
    p, g = tmp_path / "p.parquet", tmp_path / "g.parquet"
    pheno.to_parquet(p)
    geno.to_parquet(g)
    con = qp.connect([str(p), str(g)])
    iso = qp.isolate_table(con, str(p), str(g), "Klebsiella pneumoniae", "meropenem", "CARBAPENEM")
    result = iso.set_index("biosample")["strict_resistant"].to_dict()
    expected = {c[0]: c[6] for c in cases} | {"repeat": False}
    assert {k: bool(v) for k, v in result.items()} == expected


def test_unknown_region_ranks_below_a_used_region():
    iso = pd.DataFrame(
        {
            "biosample": ["A", "B"],
            "assembly": ["GCA_A", "GCA_B"],
            "collection_year": [2024, 2020],
            "region": ["Unknown", "Asia"],
        }
    )
    assert qp._pick(iso, 1, used_regions={"Asia"})["biosample"].tolist() == ["B"]


def test_selection_spreads_regions(parquet_files):
    _, iso = build(parquet_files)
    sel = qp.select_study_set(iso)
    with_gene = sel[sel["selection_group"] == "resistant_with_carbapenemase"]
    # Two Africa isolates and one Europe isolate exist; Europe must be included
    assert "Europe" in set(with_gene["region"])


def test_portal_calls_returns_all_genes_for_selected(parquet_files):
    con, iso = build(parquet_files)
    pheno, geno = parquet_files
    calls = qp.portal_calls(con, geno, ["GCA_1"])
    assert sorted(calls["amr_element_symbol"]) == ["blaKPC-2", "blaSHV-11"]


def test_summary_counts(parquet_files):
    _, iso = build(parquet_files)
    summary = qp.summarise(iso)
    assert summary["n_isolates"].sum() == len(iso)


@pytest.fixture()
def ensembl_list(tmp_path: Path) -> str:
    """Mock Ensembl Bacteria species list with a '#'-prefixed header."""
    path = tmp_path / "species_EnsemblBacteria.txt"
    path.write_text(
        "#name\tspecies\tassembly_accession\n"
        "klebsiella_x\tKlebsiella pneumoniae\tGCA_2.1\n"
        "ecoli_y\tEscherichia coli\tGCA_999.1\n"
    )
    return str(path)


def test_ensembl_flag_ignores_versions(parquet_files, ensembl_list):
    _, iso = build(parquet_files)
    flagged = qp.flag_in_ensembl(iso, qp.load_ensembl_accessions(ensembl_list)).set_index("biosample")
    assert bool(flagged.loc["S2", "in_ensembl"]) is True
    assert bool(flagged.loc["S6", "in_ensembl"]) is False


def test_selection_prefers_isolates_in_ensembl(parquet_files, ensembl_list):
    _, iso = build(parquet_files)
    flagged = qp.flag_in_ensembl(iso, qp.load_ensembl_accessions(ensembl_list))
    sel = qp.select_study_set(flagged, {"any": "^bla"})
    first = sel[sel["selection_group"] == "resistant_with_carbapenemase"]
    # S6 is more recent, but S2 is in Ensembl and should be preferred
    assert first["biosample"].tolist() == ["S2"]


def test_missing_ensembl_column_raises(ensembl_list):
    with pytest.raises(ValueError):
        qp.load_ensembl_accessions(ensembl_list, column="not_a_column")


def test_main_writes_outputs(parquet_files, tmp_path):
    pheno, geno = parquet_files
    out = tmp_path / "out"
    assert qp.main(["--pheno", pheno, "--geno", geno, "--outdir", str(out)]) == 0
    for name in ("portal_isolates.tsv", "portal_summary.tsv", "samples.tsv", "portal_calls.tsv"):
        assert (out / name).exists()
