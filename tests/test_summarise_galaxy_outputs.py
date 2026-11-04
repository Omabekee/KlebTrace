"""Unit tests for scripts/summarise_galaxy_outputs.py, using small mock Galaxy outputs."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import summarise_galaxy_outputs as sg  # noqa: E402

AMR_COLUMNS = ["Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol", "Element name",
               "Scope", "Type", "Subtype", "Class", "Subclass", "Method", "% Coverage of reference",
               "Closest reference name"]


def amr_row(symbol, subtype, subclass, method, closest, type_="AMR"):
    return ["NA", "c1", "1", "100", "+", symbol, "", "core", type_, subtype, "BETA-LACTAM", subclass,
            method, "100", closest]


def write_genome(indir: Path, acc: str, rows: list[list[str]]) -> None:
    for sub in ("amrfinderplus_report", "integronfinder_summary", "plasmidfinder", "isescan_summary"):
        (indir / sub).mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=AMR_COLUMNS).to_csv(indir / "amrfinderplus_report" / f"{acc}.tsv",
                                                  sep="\t", index=False)
    (indir / "integronfinder_summary" / f"{acc}.tsv").write_text(
        "# integron_finder 2.0.5\nID_replicon\tCALIN\tcomplete\tIn0\ttopology\tsize\nc1\t0\t0\t0\tlin\t5000\n"
    )
    (indir / "plasmidfinder" / f"{acc}.tsv").write_text("Database\tPlasmid\nenterobacteriales\tIncFII(K)\n")
    (indir / "isescan_summary" / f"{acc}.tsv").write_text(
        "# seqid family nIS %Genome bps4IS dnaLen\nc1 IS3 2 1.0 2000 5000\nc1 total 2 1.0 2000 5000\n"
    )


MANIFEST = pd.DataFrame(
    {"assembly": ["G1"], "n_sequences": [250], "total_bp": [5000], "n50": [5000]}
).set_index("assembly")


def test_generic_carbapenemase_symbol_is_still_counted(tmp_path):
    # A split OXA-232 is reported as "blaOXA" with subclass BETA-LACTAM (as in GCA_019928365.1).
    write_genome(tmp_path, "G1", [
        amr_row("blaOXA", "AMR", "BETA-LACTAM", "BLASTX",
                "OXA-48 family carbapenem-hydrolyzing class D beta-lactamase OXA-232"),
        amr_row("blaOXA-1", "AMR", "CEPHALOSPORIN", "ALLELEX", "oxacillin-hydrolyzing class D beta-lactamase OXA-1"),
    ])
    row, _ = sg.summarise_genome("G1", tmp_path, MANIFEST)
    assert row["carbapenemase_genes"] == "blaOXA [closest OXA-232]"
    assert row["other_esbl_ampc"] == "blaOXA-1"
    assert "carbapenemase allele not named" in row["qc_flags"]


def test_partial_copies_porins_and_manifest_stats(tmp_path):
    write_genome(tmp_path, "G1", [
        amr_row("blaNDM-1", "AMR", "CARBAPENEM", "ALLELEX", "subclass B1 metallo-beta-lactamase NDM-1"),
        amr_row("blaNDM-1", "AMR", "CARBAPENEM", "PARTIAL_CONTIG_ENDX", "subclass B1 metallo-beta-lactamase NDM-1"),
        amr_row("ompK35_Y76Ter", "POINT_DISRUPT", "CARBAPENEM", "POINTX", "outer membrane porin OmpK35"),
    ])
    row, _ = sg.summarise_genome("G1", tmp_path, MANIFEST)
    assert row["carbapenemase_genes"] == "blaNDM-1 (x2) [partial]"
    assert row["carbapenem_porin_changes"] == "ompK35_Y76Ter"
    assert row["n_sequences"] == 250
    assert row["n_is_elements"] == 2
    assert "fragmented assembly (250 sequences)" in row["qc_flags"]
    assert "carbapenemase allele not named" not in row["qc_flags"]
