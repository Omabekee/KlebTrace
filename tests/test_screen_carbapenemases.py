"""Unit tests for scripts/screen_carbapenemases.py (no network access or minimap2)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import screen_carbapenemases as sc  # noqa: E402


def hit(start, stop, tstart=0, tend=None, tlen=5000):
    # By default the piece starts at the beginning of its contig.
    return {"qstart": start, "qend": stop, "tstart": tstart,
            "tend": tstart + (stop - start) if tend is None else tend, "tlen": tlen}


def paf(allele, qlen, start, stop, contig, matches, aligned):
    return "\t".join([f"1__X__{allele}__9", str(qlen), str(start), str(stop), "+", contig,
                      "1000", "0", str(stop - start), str(matches), str(aligned), "60"])


def test_covered_bp_merges_overlaps():
    assert sc.covered_bp([(0, 360), (125, 831), (765, 882)]) == 882
    assert sc.covered_bp([(0, 10), (20, 30)]) == 20


def test_classify_intact_split_partial_none():
    assert sc.classify_allele([hit(0, 882)], 882)[0] == "intact"
    assert sc.classify_allele([hit(0, 360), hit(125, 831), hit(765, 882)], 882)[0] == "split"
    assert sc.classify_allele([hit(0, 737), hit(737, 798)], 798)[0] == "split"
    assert sc.classify_allele([hit(0, 300)], 882)[0] == "partial"
    assert sc.classify_allele([], 882)[0] == "none"


def test_parse_paf_filters_weak_hits():
    text = "\n".join([paf("KPC-2", 882, 0, 882, "c1", 882, 882),
                      paf("OXA-48", 798, 0, 798, "c2", 700, 798),
                      paf("NDM-1", 813, 0, 40, "c3", 40, 40)])
    assert sc.parse_paf(text)["allele"].tolist() == ["KPC-2"]


def test_best_call_prefers_exact_allele_and_reports_pieces():
    hits = sc.parse_paf("\n".join([paf("KPC-3", 882, 0, 700, "c1", 699, 700),
                                   paf("KPC-3", 882, 700, 882, "c2", 180, 182),
                                   paf("KPC-2", 882, 0, 700, "c1", 700, 700),
                                   paf("KPC-2", 882, 700, 882, "c2", 182, 182)]))
    call = sc.best_call(hits)
    assert (call["call"], call["allele"], call["n_pieces"]) == ("split", "KPC-2", 2)
    assert call["contigs"] == "c1;c2"


def test_best_call_with_no_hits():
    assert sc.best_call(sc.parse_paf(""))["call"] == "none"


def test_references_exclude_non_carbapenemases(tmp_path):
    db = tmp_path / "db.csv"
    db.write_text(
        "allele,gene,bla_class,bla_description\n"
        "KPC-2,KPC,Bla_Carb,carbapenem-hydrolyzing class A beta-lactamase\n"
        "CTX-M-33,CTX-M,Bla_Carb,class A carbapenem-hydrolyzing beta-lactamase\n"
        "OXA-233,OXA,Bla_Carb,OXA-10 family class D beta-lactamase\n"
        "TEM-1,TEM,Bla,beta-lactamase\n"
    )
    fasta = tmp_path / "db.fasta"
    fasta.write_text("".join(f">1__X__{a}__1\nACGT\n" for a in ("KPC-2", "CTX-M-33", "OXA-233", "TEM-1")))
    out = tmp_path / "refs.fasta"
    assert sc.carbapenemase_references(db, fasta, out) == 1
    assert out.read_text() == ">1__X__KPC-2__1\nACGT\n"


def test_split_needs_every_piece_at_a_contig_end():
    # Kenya: OXA-232 at the two ends of one 6,141 bp contig.
    kenya = [hit(0, 737, tstart=0, tlen=6141), hit(737, 798, tstart=6080, tlen=6141)]
    assert sc.classify_allele(kenya, 798)[0] == "split"
    # The same pieces in the middle of a contig: an interrupted gene, not an assembly break.
    inside = [hit(0, 737, tstart=1000, tlen=6141), hit(737, 798, tstart=4000, tlen=6141)]
    assert sc.classify_allele(inside, 798)[0] == "interrupted"
