"""Unit tests for scripts/finalise_study_set.py."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import finalise_study_set as fs  # noqa: E402

SAMPLES = pd.DataFrame({
    "assembly": ["GCA_1", "GCA_2", "GCA_3"],
    "selection_group": ["resistant_with_carbapenemase",
                        "resistant_no_carbapenemase_porin_mutation",
                        "resistant_no_carbapenemase_no_porin_mutation"],
    "carbapenemase_family": ["NDM", "", ""],
})
SCREEN = pd.DataFrame({
    "biosample": ["S2", "S3", "S4", "S5", "S6", "S7"],
    "assembly": ["GCA_2", "GCA_3", "GCA_4", "GCA_5", "GCA_6", "GCA_7"],
    "region": ["Africa", "Europe", "Unknown", "Americas", "Asia", "Asia"],
    "collection_year": [2015, 2011, None, 2012, 2015, 2016],
    "n_sequences": [106, 304, 3, 5, 50, 4],
    "call": ["split", "partial", "none", "none", "none", "none"],
    "allele": ["OXA-232", "KPC-170", "", "", "", ""],
})
KLEBORATE = pd.DataFrame({
    "strain": ["GCA_4", "GCA_5", "GCA_6", "GCA_7"],
    fs.KLEBORATE_ST: ["ST1161", "ST17", "ST15", "ST11"],
    fs.KLEBORATE_CARB: ["-", "-", "-", "KPC-2"],
    fs.KLEBORATE_ESBL: ["CTX-M-2", "-", "CTX-M-15", "-"],
    fs.KLEBORATE_PORIN: ["OmpK35-25%", "OmpK36-c.T108del", "-", "-"],
})


def test_whole_gene_relabels_but_partial_does_not():
    out = fs.relabel_reannotated(SAMPLES, SCREEN)
    groups = dict(zip(out["assembly"], out["selection_group"]))
    assert groups["GCA_2"] == fs.REANNOTATED_GROUP
    assert out.loc[out["assembly"] == "GCA_2", "carbapenemase_family"].item() == "OXA-48-like"
    assert out.loc[out["assembly"] == "GCA_2", "reannotated_carbapenemase"].item() == "blaOXA-232 (split)"
    assert groups["GCA_3"] == "resistant_no_carbapenemase_no_porin_mutation"
    assert groups["GCA_1"] == "resistant_with_carbapenemase"


def test_confirmed_negative_prefers_known_region_then_fewest_sequences():
    # GCA_7 has the fewest sequences of the known regions but Kleborate finds KPC-2;
    # GCA_4 has fewer sequences than GCA_5 but an unknown region.
    pick = fs.pick_confirmed_negative(SCREEN, KLEBORATE, exclude={"GCA_2", "GCA_3"})
    assert pick["assembly"].tolist() == ["GCA_5"]


def test_genome_without_kleborate_result_is_not_picked():
    pick = fs.pick_confirmed_negative(SCREEN, KLEBORATE[KLEBORATE["strain"] != "GCA_5"], exclude=set())
    assert pick["assembly"].tolist() == ["GCA_6"]
