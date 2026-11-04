"""Unit tests for scripts/to_ensembl_seqids.py."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import to_ensembl_seqids as te  # noqa: E402

PMK1 = {
    "CP1.1": (5000, "CP1.1 Klebsiella pneumoniae strain PMK1, complete genome."),
    "CP2.1": (300, "CP2.1 Klebsiella pneumoniae strain PMK1 plasmid pPMK1-A, complete sequence."),
    "CP3.1": (300, "CP3.1 Klebsiella pneumoniae strain PMK1 plasmid pPMK1-NDM, complete sequence."),
}


def test_unique_lengths_and_ties_broken_by_name():
    # Two plasmids of equal length are told apart by the name in the FASTA header
    mapping = te.match_sequences(PMK1, {"Chromosome": 5000, "pPMK1-NDM": 300, "pPMK1-A": 300})
    assert mapping == {"CP1.1": "Chromosome", "CP3.1": "pPMK1-NDM", "CP2.1": "pPMK1-A"}


def test_unmatched_sequence_is_an_error():
    with pytest.raises(ValueError, match="no Ensembl match"):
        te.match_sequences(PMK1, {"Chromosome": 5000, "pPMK1-A": 300})
    with pytest.raises(ValueError, match="matches 0"):
        te.match_sequences(PMK1, {"Chromosome": 4999, "pPMK1-A": 300, "pPMK1-NDM": 300})


def test_name_must_be_a_whole_word():
    # "contig000001" must not match a header that names contig0000011
    header = "X.1 genome assembly, contig: ERScontig0000011"
    assert not te.mentions(header, "ERScontig000001")
    assert te.mentions("Y.1 genome assembly, contig: ERScontig000001", "ERScontig000001")


def test_rename_gff3_changes_seqids_and_pragmas_only():
    lines = ["##gff-version 3", "##sequence-region CP3.1 1 300", "# kleborate ST=ST15",
             "CP3.1\tAMRFinderPlus\tgene\t1\t100\t.\t+\t.\tID=a1;Name=blaNDM-1"]
    out = te.rename_gff3(lines, {"CP3.1": "pPMK1-NDM"})
    assert out == ["##gff-version 3", "##sequence-region pPMK1-NDM 1 300", "# kleborate ST=ST15",
                   "pPMK1-NDM\tAMRFinderPlus\tgene\t1\t100\t.\t+\t.\tID=a1;Name=blaNDM-1"]
