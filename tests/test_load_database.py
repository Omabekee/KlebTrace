"""Unit tests for scripts/load_database.py with a one-genome example."""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import load_database as ld  # noqa: E402


@pytest.fixture
def inputs(tmp_path):
    pd.DataFrame([{"biosample": "S1", "assembly": "G.1", "country": "Kenya", "collection_year": "2015",
                   "region": "Africa", "phenotype": "resistant", "mic_mg_l": ">=16",
                   "strict_resistant": "True", "in_ensembl": "False",
                   "selection_group": "portal_negative_carbapenemase_on_reannotation",
                   "carbapenemase_family": "OXA-48-like"}]).to_csv(tmp_path / "samples.tsv", sep="\t", index=False)
    pd.DataFrame([{"assembly": "G.1", "n_sequences": 1, "total_bp": 1000, "n50": 1000,
                   "gc_percent": 50.0, "md5": "x"}]).to_csv(tmp_path / "manifest.tsv", sep="\t", index=False)
    (tmp_path / "G.1.fasta").write_text(">c1 contig\n" + "A" * 1000 + "\n")
    (tmp_path / "G.1.klebtrace.gff3").write_text(
        "##gff-version 3\n"
        "c1\tBakta\tCDS\t1\t300\t.\t+\t0\tID=B_1;Name=x;product=a%3Bb\n"
        "c1\tAMRFinderPlus\tgene\t3\t737\t.\t-\t.\tID=amrfinderplus_1;Name=blaOXA\n".replace("737", "700"))
    portal_row = {"assembly_ID": "G.1", "amr_element_symbol": "ompK36_D135DGD", "element_subtype": "POINT",
                  "subclass": "CARBAPENEM", "reference_name": "outer membrane porin OmpK36",
                  "amrfinderplus_method": "POINTP", "region": "c1", "region_start": "10",
                  "region_end": "20", "reference_sequence_coverage": "100.0",
                  "reference_sequence_identity": "98.9"}
    # The same call repeated for a second antibiotic, as in the portal's split_subclass rows
    pd.DataFrame([portal_row, portal_row]).to_csv(tmp_path / "portal.tsv", sep="\t", index=False)
    pd.DataFrame([{"Contig id": "c1", "Start": "3", "Stop": "700", "Element symbol": "blaOXA", "Type": "AMR",
                   "Subtype": "AMR", "Subclass": "BETA-LACTAM", "Method": "BLASTX",
                   "% Coverage of reference": "92.45", "% Identity to reference": "100.00",
                   "Closest reference name": "OXA-48 family carbapenem-hydrolyzing class D beta-lactamase OXA-232"}]
                 ).to_csv(tmp_path / "G.1.tsv", sep="\t", index=False)
    pd.DataFrame([{"strain": "G.1", "kpsc__amr__Bla_Carb_acquired": "OXA-232?",
                   "kpsc__amr__Omp_mutations": "ompK35:p.Y76Ter;ompK36:p.134_135insGD"}]
                 ).to_csv(tmp_path / "kleborate.tsv", sep="\t", index=False)
    return tmp_path


def load(d: Path, **extra) -> Path:
    args = ["--schema", str(ROOT / "database" / "schema.sql"), "--samples", str(d / "samples.tsv"),
            "--manifest", str(d / "manifest.tsv"), "--portal-calls", str(d / "portal.tsv"),
            "--screen", str(d / "no_screen.tsv"), "--fasta", str(d / "G.1.fasta"),
            "--gff3", str(d / "G.1.klebtrace.gff3"), "--amrfinderplus", str(d / "G.1.tsv"),
            "--kleborate", str(d / "kleborate.tsv"), "--out", str(d / "k.sqlite")]
    assert ld.main(args) == 0
    return d / "k.sqlite"


def test_everything_loads(inputs):
    con = sqlite3.connect(load(inputs))
    assert con.execute("SELECT count(*) FROM features").fetchone()[0] == 2
    assert con.execute("SELECT product FROM features WHERE gff_id = 'B_1'").fetchone()[0] == "a;b"
    calls = dict(con.execute("SELECT source, count(*) FROM resistance_calls GROUP BY source").fetchall())
    assert calls == {"portal": 1, "amrfinderplus": 1, "kleborate": 3}   # portal repeat collapsed


def test_unnamed_oxa_counts_as_carbapenemase(inputs):
    con = sqlite3.connect(load(inputs))
    flags = dict(con.execute("SELECT source, is_carbapenemase FROM resistance_calls "
                             "WHERE gene IN ('blaOXA', 'OXA-232?')").fetchall())
    assert flags == {"amrfinderplus": 1, "kleborate": 1}
    assert con.execute("SELECT is_partial FROM resistance_calls WHERE gene = 'OXA-232?'").fetchone()[0] == 1


def test_bad_feature_rolls_back_and_leaves_no_database(inputs):
    gff = inputs / "G.1.klebtrace.gff3"
    gff.write_text(gff.read_text().replace("c1\tBakta", "contig_1\tBakta"))
    with pytest.raises(sqlite3.IntegrityError):
        load(inputs)
    assert not (inputs / "k.sqlite").exists()
