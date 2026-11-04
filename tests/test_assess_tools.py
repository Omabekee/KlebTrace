"""Unit tests for scripts/assess_added_value.py and scripts/assess_accuracy.py."""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import assess_accuracy as acc  # noqa: E402
import assess_added_value as av  # noqa: E402
import load_database as ld  # noqa: E402


def test_base_name_drops_prokka_copy_number():
    assert av.base_name("blaNDM-1_2") == "blaNDM-1"
    assert av.base_name("dsbD_3") == "dsbD"
    assert av.base_name("bla") == "bla"


def test_classify_same_other_none_and_not_a_gene():
    genes = pd.DataFrame({"seqid": ["p1"] * 2, "strand": ["+"] * 2, "start": [100, 1000], "end": [900, 1900],
                          "name": ["blaNDM-1", "bla"], "description": ["MBL", "Beta-lactamase CTX-M-1"]})
    feats = pd.DataFrame({
        "source": ["AMRFinderPlus"] * 3 + ["PlasmidFinder"], "type": ["gene"] * 3 + ["origin_of_replication"],
        "seqid": ["p1"] * 4, "start": [100, 1000, 5000, 10], "end": [897, 1897, 5300, 60],
        "strand": ["+", "+", "+", "."], "name": ["blaNDM-1", "blaCTX-M-15", "qnrB1", "IncFII"],
        "subclass": [None] * 4, "element_subtype": ["AMR"] * 3 + [None]})
    status = av.classify(feats, genes)["status"].tolist()
    assert status == ["same_name", "other_name", "no_ensembl_gene", "not_a_gene"]


def test_porins_reads_both_naming_styles():
    assert acc.porins(["ompK36_D135DGD", "ompK35:p.Y76Ter", "gyrA_S83F"]) == {"ompK35", "ompK36"}


def test_accuracy_scores_portal_miss(tmp_path):
    con = ld.create(tmp_path / "k.sqlite", ROOT / "database" / "schema.sql")
    con.execute("INSERT INTO genomes (assembly, biosample, selection_group) "
                "VALUES ('K.1', 'S', 'portal_negative_carbapenemase_on_reannotation')")
    con.execute("INSERT INTO phenotypes VALUES ('K.1', 'meropenem', 'resistant', '>=16', 1)")
    con.executemany("INSERT INTO resistance_calls (assembly, source, gene, element_subtype, is_carbapenemase) "
                    "VALUES ('K.1', ?, ?, ?, ?)",
                    [("portal", "ompK36_D135DGD", "POINT", 0), ("amrfinderplus", "blaOXA", "AMR", 1)])
    con.execute("COMMIT")
    table = acc.score(sqlite3.connect(tmp_path / "k.sqlite")).set_index("source")
    assert table.loc["portal", "explains_phenotype"] == "no"
    assert table.loc["amrfinderplus", "explains_phenotype"] == "yes"
