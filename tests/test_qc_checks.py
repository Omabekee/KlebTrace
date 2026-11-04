"""Unit tests for scripts/qc_checks.py."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import qc_checks as qc  # noqa: E402

HEADER = ["##gff-version 3"]


def row(seqid="c1", source="Bakta", type_="CDS", start=1, end=300, strand="+", phase="0", attrs="ID=f1"):
    return [seqid, source, type_, str(start), str(end), ".", strand, phase, attrs]


@pytest.mark.parametrize("bad, message", [
    (row(strand=""), "strand ''"),                                  # Phase 0 ISEScan fault
    (row(start=0), "below 1"),
    (row(start=500, end=100), "after end"),
    (row(phase="."), "CDS phase"),
    (row(attrs="ID=f1;product=a;b"), "without '='"),                # unencoded ';'
    (row()[:8], "8 columns"),
])
def test_gff3_valid_catches_each_fault(bad, message):
    status, detail = qc.check_gff3_valid("G", HEADER, [bad])[2:]
    assert status == "FAIL" and message in detail


def test_gff3_valid_passes_good_file():
    assert qc.check_gff3_valid("G", HEADER, [row()])[2] == "PASS"


def test_seqids_must_be_in_fasta_and_within_length():
    gff = qc.gff3_frame([row(seqid="contig_1"), row(attrs="ID=f2", end=2000)])
    status, detail = qc.check_seqids_match_fasta("G", gff, {"c1": 1000})[2:]
    assert status == "FAIL" and "contig_1 not in FASTA" in detail and "ends at 2000" in detail


def test_duplicate_ids_fail():
    gff = qc.gff3_frame([row(), row()])
    assert qc.check_unique_ids("G", gff)[2] == "FAIL"


def test_cds_density_scales_with_genome_size():
    # 5,813 CDS is above a fixed 5,600 limit, but normal for a 6.13 Mb genome (Saint Lucia)
    gff = qc.gff3_frame([row(attrs=f"ID=f{i}") for i in range(5813)])
    assert qc.check_cds_density("G", gff, 6_133_985)[2] == "PASS"
    assert qc.check_cds_density("G", gff, 9_000_000)[2] == "WARN"


@pytest.mark.parametrize("strand, d_start, d_end, expected", [
    ("+", 0, 3, (0, 3)),      # same start, Bakta includes the stop codon
    ("-", -3, 0, (0, 3)),     # minus strand: the stop codon is at the lower coordinate
    ("+", -6, 3, (6, 3)),     # Bakta chose a start codon two codons upstream
])
def test_compare_ends(strand, d_start, d_end, expected):
    assert qc.compare_ends(strand, d_start, d_end) == expected


def amr(start, end, strand="+", name="blaNDM-1", attrs=""):
    return row(source="AMRFinderPlus", type_="gene", start=start, end=end, strand=strand, phase=".",
               attrs=f"ID=a{start};Name={name};element_subtype=AMR;method=ALLELEX{attrs}")


def test_coordinates_agree_accepts_stop_codon_and_new_start():
    gff = qc.gff3_frame([row(start=47947, end=48759), amr(47947, 48756),          # NDM-1 in PMK1
                         row(attrs="ID=f2", start=997, end=2003, strand="-"),
                         amr(1000, 1997, strand="-", name="blaX")])               # other start codon
    status, detail = qc.check_coordinates_agree("G", gff)[2:]
    assert status == "PASS" and "1 same start, 1 different start codon" in detail


def test_coordinates_agree_flags_missing_bakta_gene():
    # Kenya: Bakta has no CDS for the OXA-232 split across contig ends
    gff = qc.gff3_frame([amr(3, 737, strand="-", name="blaOXA")])
    status, detail = qc.check_coordinates_agree("G", gff)[2:]
    assert status == "WARN" and "has no Bakta CDS" in detail
