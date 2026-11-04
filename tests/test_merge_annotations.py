"""Unit tests for scripts/merge_annotations.py, using a tiny two-sequence genome."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import merge_annotations as ma  # noqa: E402

BAKTA = """##gff-version 3
##sequence-region CP1.1 1 1000
CP1.1\tBakta\tregion\t1\t1000\t.\t+\t.\tID=CP1.1;Name=CP1.1
CP1.1\tProdigal\tCDS\t101\t400\t.\t+\t0\tID=ABC_00005;Name=metallo-beta-lactamase NDM-1;product=subclass B1 metallo-beta-lactamase NDM-1;gene=blaNDM-1
CP2.1\tBLAST+\toriC\t5\t50\t.\t?\t.\tID=ABC_00010;Name=origin of replication
##FASTA
>CP1.1
ACGT
"""
AMR_COLUMNS = ["Contig id", "Start", "Stop", "Strand", "Element symbol", "Element name", "Type", "Subtype",
               "Class", "Subclass", "Method", "% Coverage of reference", "% Identity to reference",
               "Closest reference accession", "Closest reference name"]


@pytest.fixture
def genome(tmp_path):
    (tmp_path / "g.fasta").write_text(">CP1.1 chromosome\n" + "A" * 1000 + "\n>CP2.1 plasmid\n" + "C" * 200 + "\n")
    (tmp_path / "bakta.gff3").write_text(BAKTA)
    pd.DataFrame([["CP1.1", 101, 397, "+", "blaNDM-1", "metallo-beta-lactamase; NDM=1", "AMR", "AMR",
                   "BETA-LACTAM", "CARBAPENEM", "ALLELEX", "100.00", "100.00", "WP_1.1", "NDM-1"]],
                 columns=AMR_COLUMNS).to_csv(tmp_path / "amr.tsv", sep="\t", index=False)
    pd.DataFrame([{"Database": "enterobacteriales", "Plasmid": "IncFII(K)", "Identity": "100.0",
                   "Query / Template length": "148 / 148", "Contig": "CP2.1 Klebsiella plasmid",
                   "Position in contig": "10..157", "Note": "", "Accession number": "CP000648"}]
                 ).to_csv(tmp_path / "pf.tsv", sep="\t", index=False)
    pd.DataFrame([{"strain": "G.1", "general__species__species": "Klebsiella pneumoniae",
                   "kpsc__mlst__ST": "ST15", "kpsc__amr__Bla_Carb_acquired": "NDM-1"}]
                 ).to_csv(tmp_path / "kl.tsv", sep="\t", index=False)
    return tmp_path


def run(d: Path) -> list[str]:
    return ma.merge("G.1", d / "g.fasta", d / "bakta.gff3", d / "amr.tsv", d / "kl.tsv", d / "pf.tsv")


def features(lines):
    return [line.split("\t") for line in lines if not line.startswith("#")]


def test_merge_harmonises_sources_types_and_order(genome):
    lines = run(genome)
    assert "##sequence-region CP2.1 1 200" in lines
    assert "# kleborate ST=ST15" in lines
    assert not any(line.startswith(">") for line in lines)
    rows = features(lines)
    assert [r[0] for r in rows] == ["CP1.1", "CP1.1", "CP1.1", "CP2.1", "CP2.1"]
    assert {r[1] for r in rows} == {"Bakta", "AMRFinderPlus", "PlasmidFinder"}
    cds = next(r for r in rows if r[2] == "CDS")
    assert cds[8].endswith(";predicted_by=Prodigal")
    assert next(r for r in rows if r[1] == "Bakta" and r[0] == "CP2.1")[2] == "origin_of_replication"


def test_amrfinderplus_attributes_are_encoded(genome):
    amr = next(r for r in features(run(genome)) if r[1] == "AMRFinderPlus")
    assert amr[3:5] == ["101", "397"]
    assert "product=metallo-beta-lactamase%3B NDM%3D1" in amr[8]
    assert "Dbxref=RefSeq:WP_1.1" in amr[8]


def test_feature_on_unknown_contig_is_rejected(genome):
    (genome / "bakta.gff3").write_text(BAKTA.replace("CP2.1\tBLAST+", "contig_7\tBLAST+"))
    with pytest.raises(ValueError, match="not in the FASTA"):
        run(genome)


def test_feature_past_sequence_end_is_rejected(genome):
    (genome / "bakta.gff3").write_text(BAKTA.replace("\t5\t50\t", "\t5\t500\t"))
    with pytest.raises(ValueError, match="beyond"):
        run(genome)
