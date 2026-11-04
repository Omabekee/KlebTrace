"""Unit tests for scripts/download_genomes.py (no network access)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import download_genomes as dg  # noqa: E402

RECORD = {"expected_sequences": 2, "expected_bp": 10}


def test_normalise_fasta_shortens_headers_and_counts():
    text = ">ENA|CP1|CP1.1 Klebsiella plasmid\nACGT\nGG\n>ENA|CP2|CP2.1 chromosome\nAATT\n"
    fasta, lengths, gc = dg.normalise_fasta(text)
    assert fasta.startswith(">CP1.1 Klebsiella plasmid\nACGT\nGG\n>CP2.1 chromosome")
    assert lengths == {"CP1.1": 6, "CP2.1": 4}
    assert gc == 4


def test_duplicate_contig_is_rejected():
    with pytest.raises(ValueError):
        dg.normalise_fasta(">ENA|A|A.1 x\nAC\n>ENA|A|A.1 y\nGT\n")


def test_n50():
    assert dg.n50([50, 30, 20]) == 50
    assert dg.n50([40, 30, 20, 10]) == 30


def test_check_genome_passes_when_complete():
    assert dg.check_genome({"A.1": 6, "B.1": 4}, RECORD, {"A.1"}) == []


@pytest.mark.parametrize(
    "lengths, portal, expected",
    [
        ({"A.1": 10}, set(), "1 sequences"),
        ({"A.1": 6, "B.1": 3}, set(), "9 bp"),
        ({"A.1": 6, "B.1": 4}, {"C.1"}, "C.1"),
    ],
)
def test_check_genome_reports_problems(lengths, portal, expected):
    problems = dg.check_genome(lengths, RECORD, portal)
    assert any(expected in p for p in problems)


def test_parse_assembly_xml():
    xml = b"""<ROOT><ASSEMBLY><TITLE>T</TITLE><ASSEMBLY_LEVEL>contig</ASSEMBLY_LEVEL>
    <ASSEMBLY_ATTRIBUTES>
    <ASSEMBLY_ATTRIBUTE><TAG>scaffold-count</TAG><VALUE>119</VALUE></ASSEMBLY_ATTRIBUTE>
    <ASSEMBLY_ATTRIBUTE><TAG>count-contig</TAG><VALUE>144</VALUE></ASSEMBLY_ATTRIBUTE>
    <ASSEMBLY_ATTRIBUTE><TAG>total-length</TAG><VALUE>6024081</VALUE></ASSEMBLY_ATTRIBUTE>
    </ASSEMBLY_ATTRIBUTES></ASSEMBLY></ROOT>"""
    assert dg.parse_assembly_xml(xml) == {
        "assembly_level": "contig", "title": "T", "expected_sequences": 119, "expected_bp": 6024081,
    }
