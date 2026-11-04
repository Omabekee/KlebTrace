"""Unit tests for the Phase 2 helper scripts."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import compare_amr_calls as cmp  # noqa: E402
import make_gff3_tracks as tracks  # noqa: E402

AMR_COLUMNS = ["Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol",
               "Element name", "Scope", "Type", "Subtype", "Class", "Subclass", "Method",
               "Target length", "Reference sequence length", "% Coverage of reference",
               "% Identity to reference", "Alignment length", "Closest reference accession",
               "Closest reference name", "HMM accession", "HMM description"]


def amr_row(contig, start, stop, symbol, subclass="CARBAPENEM", subtype="AMR",
            method="ALLELEX", gtype="AMR", name="test product"):
    row = dict.fromkeys(AMR_COLUMNS, "NA")
    row.update({"Contig id": contig, "Start": str(start), "Stop": str(stop), "Strand": "+",
                "Element symbol": symbol, "Element name": name, "Type": gtype,
                "Subtype": subtype, "Class": "BETA-LACTAM", "Subclass": subclass,
                "Method": method, "% Coverage of reference": "100.00",
                "% Identity to reference": "100.00"})
    return row


# --- make_gff3_tracks -------------------------------------------------------

def test_amrfinder_coordinates_unchanged_and_attributes_encoded():
    df = pd.DataFrame([amr_row("CP008933.1", 47947, 48756, "blaNDM-1",
                               name="metallo-beta-lactamase; NDM=1")])
    gff = tracks.amrfinder_to_gff(df)
    assert gff.loc[0, ["seqid", "start", "end", "type"]].tolist() == ["CP008933.1", 47947, 48756, "gene"]
    # ';' and '=' inside a value must be percent-encoded in GFF3
    assert "product=metallo-beta-lactamase%3B NDM%3D1" in gff.loc[0, "attributes"]


def test_point_mutation_type():
    df = pd.DataFrame([amr_row("c1", 10, 20, "ompK36_D135DGD", subtype="POINT")])
    assert tracks.amrfinder_to_gff(df).loc[0, "type"] == "point_mutation"


def test_plasmidfinder_contig_trimmed_and_positions_split():
    df = pd.DataFrame([{
        "Database": "enterobacteriales", "Plasmid": "IncHI1B(pNDM-MAR)", "Identity": "100.0",
        "Query / Template length": "570 / 570",
        "Contig": "CP008933.1 Klebsiella pneumoniae strain PMK1 plasmid pPMK1-NDM, complete sequence",
        "Position in contig": "184417..184986", "Note": "", "Accession number": "JN420336",
    }])
    gff = tracks.plasmidfinder_to_gff(df)
    assert gff.loc[0, ["seqid", "start", "end", "strand"]].tolist() == ["CP008933.1", 184417, 184986, "."]


# --- compare_amr_calls ------------------------------------------------------

def test_gene_family():
    assert cmp.gene_family("blaKPC-2") == "blaKPC"
    assert cmp.gene_family("blaKPC") == "blaKPC"
    assert cmp.gene_family("blaOXA-48") == "blaOXA-48"
    assert cmp.gene_family("ompK36_D135DGD") == "ompK36"


def test_explains_phenotype():
    assert cmp.explains("resistant", ["blaNDM-1"], []) == "yes"
    assert cmp.explains("resistant", [], ["ompK35_Y76Ter", "ompK36_D135DGD"]) == "yes"
    assert cmp.explains("resistant", [], ["ompK36_D135DGD"]).startswith("no")
    assert cmp.explains("susceptible", [], []) == "yes"


def test_main_compares_portal_and_klebtrace(tmp_path):
    samples = pd.DataFrame({"assembly": ["GCA_1.1", "GCA_2.1"],
                            "phenotype": ["resistant", "susceptible"]})
    portal = pd.DataFrame({
        "assembly_ID": ["GCA_1.1", "GCA_1.1"],
        "gene_symbol": ["bla_2", "bla_1"],
        "amr_element_symbol": ["blaKPC-2", "blaSHV-11"],
        "subclass": ["CARBAPENEM", "BETA-LACTAM"],
        "element_subtype": ["AMR", "AMR"],
        "element_type": ["AMR", "AMR"],
    })
    kleb = pd.DataFrame([amr_row("c1", 1, 700, "blaKPC", method="PARTIAL_CONTIG_ENDX"),
                         amr_row("c1", 900, 1300, "fosA", subclass="FOSFOMYCIN")])
    kleb.insert(0, "assembly", ["GCA_1.1", "GCA_2.1"])
    for name, df in (("s", samples), ("p", portal), ("k", kleb)):
        df.to_csv(tmp_path / f"{name}.tsv", sep="\t", index=False)

    out = tmp_path / "out"
    assert cmp.main(["--samples", str(tmp_path / "s.tsv"), "--portal", str(tmp_path / "p.tsv"),
                     "--klebtrace", str(tmp_path / "k.tsv"), "--outdir", str(out)]) == 0
    res = pd.read_csv(out / "comparison_carbapenem.tsv", sep="\t").set_index("assembly")
    assert res.loc["GCA_1.1", "carbapenemase_agree"] == "yes"
    assert res.loc["GCA_1.1", "klebtrace_carbapenemase"] == "blaKPC [partial]"
    assert res.loc["GCA_2.1", "klebtrace_explains_phenotype"] == "yes"
    genes = pd.read_csv(out / "comparison_all_genes.tsv", sep="\t")
    assert set(genes["found_by"]) == {"portal only", "KlebTrace only"}
    assert res.loc["GCA_2.1", "klebtrace_carbapenemase"] == "none"


def test_portal_uses_element_symbol_and_drops_antibiotic_repeats():
    # The same SHV-12 call repeated for two antibiotics (split_subclass) is one gene;
    # two OXA-48 copies at different positions stay two.
    portal = pd.DataFrame({
        "assembly_ID": ["G.1"] * 4,
        "gene_symbol": ["bla_1", "bla_1", "bla_4", "bla_7"],
        "amr_element_symbol": ["blaSHV-12", "blaSHV-12", "blaOXA-48", "blaOXA-48"],
        "subclass": ["CEPHALOSPORIN", "CEPHALOSPORIN", "CARBAPENEM", "CARBAPENEM"],
        "element_subtype": ["AMR"] * 4, "element_type": ["AMR"] * 4,
        "reference_name": [""] * 4,
        "region": ["c1", "c1", "c2", "c3"], "region_start": ["1", "1", "5", "5"],
        "region_end": ["9", "9", "20", "20"],
    })
    tidy = cmp.tidy(portal, cmp.PORTAL, "portal")
    assert tidy["gene"].tolist() == ["blaSHV-12", "blaOXA-48", "blaOXA-48"]
    assert cmp.mechanisms(tidy)[0] == ["blaOXA-48 (x2)"]


def test_unnamed_carbapenemase_allele_is_counted():
    # GCA_019928365.1: OXA-232 split across contig ends is reported as "blaOXA",
    # subclass BETA-LACTAM, closest reference OXA-232.
    kleb = pd.DataFrame([amr_row("c1", 3, 737, "blaOXA", subclass="BETA-LACTAM", method="BLASTX",
                                 name="OXA-48 family class D beta-lactamase")])
    kleb["Closest reference name"] = "OXA-48 family carbapenem-hydrolyzing class D beta-lactamase OXA-232"
    kleb.insert(0, "assembly", "G.1")
    genes, _ = cmp.mechanisms(cmp.tidy(kleb, cmp.KLEBTRACE, "klebtrace"))
    assert genes == ["blaOXA"]


def test_genome_without_galaxy_output_is_not_annotated(tmp_path):
    samples = pd.DataFrame({"assembly": ["G.1"], "phenotype": ["resistant"]})
    portal = pd.DataFrame(columns=list(cmp.PORTAL.values()))
    kleb = pd.DataFrame([amr_row("c1", 1, 700, "blaKPC-2")])
    kleb.insert(0, "assembly", "OTHER.1")
    for name, df in (("s", samples), ("p", portal), ("k", kleb)):
        df.to_csv(tmp_path / f"{name}.tsv", sep="\t", index=False)
    cmp.main(["--samples", str(tmp_path / "s.tsv"), "--portal", str(tmp_path / "p.tsv"),
              "--klebtrace", str(tmp_path / "k.tsv"), "--outdir", str(tmp_path)])
    res = pd.read_csv(tmp_path / "comparison_carbapenem.tsv", sep="\t")
    assert res.loc[0, "klebtrace_explains_phenotype"] == "not annotated"
