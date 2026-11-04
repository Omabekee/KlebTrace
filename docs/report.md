# KlebTrace: Re-annotating *Klebsiella pneumoniae* Genomes from the AMR Portal to Find Missing Carbapenem Resistance Calls

#### Author: Chioma B Oselu (@Omabekee)

#### GitHub Code: [KlebTrace](https://github.com/Omabekee/KlebTrace)

This project was carried out as part of my EMBL-UNESCO fellowship at EMBL-EBI, where I work on the AMR Portal. It uses the portal's public data to test how well its resistance calls explain laboratory results, and what specialist annotation tools could add to Ensembl.

---

## Table of Contents
1. [Introduction](#1-introduction)
2. [Data](#2-data)
3. [Methods](#3-methods)
4. [Results](#4-results)
5. [Conclusion and Future Directions](#5-conclusion-and-future-directions)
6. [References](#references)

---

## 1. Introduction

*Klebsiella pneumoniae* is a common cause of hospital infections, and strains resistant to carbapenems leave very few treatment options. Carbapenem resistance in this species usually comes from an acquired carbapenemase gene such as *bla*KPC, *bla*NDM or *bla*OXA-48, often carried on a plasmid. It can also come from changes in the outer membrane porins OmpK35 and OmpK36, which reduce how much drug gets into the cell.

The EMBL-EBI AMR Portal links laboratory antibiotic test results to resistance genes found in each isolate's genome by AMRFinderPlus. This makes it possible to ask whether the genes in the portal explain the laboratory result, and what is going on when they do not.

### 1.1 Project Aim

To re-annotate a small, deliberately chosen set of *K. pneumoniae* genomes from the AMR Portal, compare the resistance calls with the portal's calls and the laboratory meropenem result, and judge what specialist tools add to Ensembl's own annotation.

### 1.2 Questions

1. **Annotation agreement.** When the genomes are annotated again with Bakta, AMRFinderPlus and Kleborate, are the portal's resistance genes found?
2. **Genotype and phenotype.** Do the genes found explain the laboratory result, and what explains resistant isolates with no carbapenemase in the portal?
3. **Genomic context.** Do the carbapenemase genes sit on plasmids or on the chromosome?

## 2. Data

### 2.1 AMR Portal

AMR Portal release 2026-07 was read directly from its Parquet files. It holds 4,824 *K. pneumoniae* isolates with a genome assembly and a meropenem test result. An isolate was called resistant if any meropenem test was resistant, and susceptible only if every test was susceptible.

For isolates with no carbapenemase in the portal, a stricter rule was used: the MIC had to be above 8 mg/L and the isolate had to be resistant under both the updated CLSI and EUCAST breakpoints. Of the 109 resistant isolates with no carbapenemase, 47 met this rule, and 37 of these had a GenBank assembly.

### 2.2 Study Genomes

Seven genomes were chosen so that each one tests a different question.

| Group | Assembly | Country | Year | Meropenem MIC (mg/L) | Portal carbapenemase |
|---|---|---|---|---|---|
| Carbapenemase, NDM | GCA_000764615.1 | Nepal | 2011 | resistant, no MIC | *bla*NDM-1 |
| Carbapenemase, KPC | GCA_049383265.1 | Brazil | 2023 | >128 | *bla*KPC-2 |
| Carbapenemase, OXA-48-like | GCA_043905605.1 | Tunisia | 2022 | >32 | *bla*OXA-48 |
| No carbapenemase in the portal | GCA_019928365.1 | Kenya | 2015 | >=16 | none |
| No carbapenemase in the portal | GCA_001182445.1 | Italy | 2011 | >32 | none |
| No carbapenemase, confirmed | GCA_000567225.2 | USA | 2012 | >32 | none |
| Susceptible control | GCA_900774885.1 | Saint Lucia | 2018 | susceptible, no MIC | none |

The assemblies were downloaded from ENA. Two of them, Nepal and Saint Lucia, are also in Ensembl Bacteria release 63, which gave a baseline annotation to compare against.

## 3. Methods

### 3.1 Selecting the Genomes

The AMR Portal Parquet files were queried with DuckDB in Python, reading only the needed columns over HTTPS. Carbapenemase genes were taken from the AMRFinderPlus element symbol, and porin point mutations were kept separate from acquired genes. One genome was picked for each carbapenemase family, preferring genomes that are in Ensembl. The seventh genome was added after the re-screen in section 3.3.

### 3.2 Annotation in Galaxy

The genomes were annotated on usegalaxy.org with a workflow that runs Bakta 1.9.4, AMRFinderPlus 4.2.7 (database 2026-05-15, organism *K. pneumoniae*, nucleotide input), PlasmidFinder 2.1.6, IntegronFinder 2.0.5 and ISEScan 1.7.2.3 on all genomes in one run. Kleborate 3.3.0b1 is not on usegalaxy.org and was run locally.

### 3.3 Re-screen for Split Carbapenemase Genes

The 37 strictly resistant isolates with no carbapenemase in the portal were screened again. The 568 carbapenemase alleles in Kleborate's database were aligned to each assembly with minimap2, keeping matches of at least 50 bp and 95% identity. For each allele, the pieces found on different contigs were added together. A gene was called split when the pieces covered at least 95% of it and every piece ended at a contig edge, and partial when they covered less. For one split gene, public Illumina reads were mapped across the join between the two contig ends to confirm it. The genomes with no match were also run through Kleborate.

### 3.4 Pipeline and Database

The Galaxy workflow was rebuilt as a Nextflow pipeline using the same tool versions, run in Apptainer containers. The Bakta step reuses the Galaxy output, because Bakta's full database (39.4 GB) was too large for the workstation. The pipeline merges the outputs into one GFF3 file per genome, loads genomes, features, laboratory results and resistance calls from every source into an SQLite database, and runs nine quality checks.

### 3.5 Comparison with Ensembl

Ensembl's GFF3 files for Nepal and Saint Lucia were downloaded. Ensembl uses different sequence names from ENA (for example "pPMK1-NDM" instead of CP008933.1), so the KlebTrace files were renamed by matching sequence lengths. Each tool's features were then matched to Ensembl's genes, and the merged annotation for Nepal was uploaded to the Ensembl browser as a custom track.

## 4. Results

### 4.1 The Portal Misses Carbapenemase Genes Split Across Contigs

The portal and KlebTrace agreed on carbapenemase status for 5 of the 7 genomes. In the other two, the portal reports no carbapenemase but the genome carries one.

| Genome | Portal | AMRFinderPlus (DNA) | Kleborate | Split-gene screen |
|---|---|---|---|---|
| Kenya | none | *bla*OXA, closest OXA-232 | OXA-232? | OXA-232, split in 2 pieces |
| Italy | none | *bla*KPC, partial at contig end | none | KPC-2, split in 3 pieces |

In Kenya, OXA-232 is split between the two ends of a 6.1 kb contig. 291 sequencing reads cross the join between the two ends, which shows that the contig is a small circular plasmid with the gene intact. This contig has about 5.5 times the read depth of the large contigs, as expected for a small plasmid with several copies per cell. In Italy, the whole KPC-2 gene is present at 100% identity in three pieces on three contigs, but no reads were available to confirm it.

The likely reason is the way the portal runs AMRFinderPlus. All its calls come from predicted proteins, while the Galaxy run searched the DNA directly. A gene broken across contigs does not give a complete protein, so a protein search can miss it. Bakta did not predict a gene at Kenya's OXA-232 either.

Across the 37 strictly resistant isolates with no carbapenemase in the portal, 2 carry a whole carbapenemase gene split across contigs and 3 more carry part of one: 5 of 37, about 14%.

### 4.2 Genotype and Phenotype

| Source | Explains the laboratory result |
|---|---|
| Portal calls | 4 of 7 |
| AMRFinderPlus on the DNA | 6 of 7 |
| Kleborate | 5 of 7 |

Resistance was counted as explained by a carbapenemase, or by changes in both OmpK35 and OmpK36. The portal also misses some porin changes, such as the OmpK35 stop codon in Kenya. One genome stays unexplained by every source: GCA_000567225.2 has a meropenem MIC above 32 mg/L, no carbapenemase, no extended-spectrum beta-lactamase and a truncated OmpK36 only.

### 4.3 Where the Carbapenemase Genes Sit

In Nepal, the only complete genome in the set, *bla*NDM-1 sits on the plasmid pPMK1-NDM (304,526 bp) in three copies. Each copy is part of an identical block of eight genes (*bla*NDM-1, *ble*, *trpF*, *dsbD*, *cutA*, *groES*, *groL* and an IS3000 transposase), and the block repeats every 8,588 bp. The same plasmid carries CTX-M-15, OXA-1, a class 1 integron, *armA*, and mercury and tellurium resistance genes. In Kenya, OXA-232 is on a small plasmid. The other genomes are in many pieces, so the location of their genes cannot be read directly.

<figure>
  <img src="../galaxy/screenshots/p2_jbrowse_PMK1_NDM_tandem.png" alt="Figure 1: The three NDM-1 copies on plasmid pPMK1-NDM in JBrowse" width="800">
  <figcaption>Figure 1: The three NDM-1 copies on plasmid pPMK1-NDM (45 to 72 kb) in JBrowse, with the Bakta, AMRFinderPlus and PlasmidFinder tracks</figcaption>
</figure>

### 4.4 Other Findings

- **Resistance and virulence together.** The Tunisian OXA-48 genome also carries the aerobactin genes, *rmpA*, *rmpA2* and *peg-344*, which are linked to hypervirulent *K. pneumoniae*. Kleborate gives it the highest virulence score in the set.
- **A duplicated assembly.** In the same genome, contigs 130 to 147 repeat sequence already present earlier in the assembly, so many genes are counted twice, in KlebTrace and in the portal.

### 4.5 Pipeline and Quality Checks

The Nextflow pipeline gave AMRFinderPlus and PlasmidFinder tables identical to the Galaxy outputs, line for line, for all seven genomes. It ran in 13 min 54 s for seven genomes on a 4-CPU workstation, with Bakta excluded. The quality checks gave 50 passes, 13 warnings and no failures. The warnings flag the Kenya and Italy carbapenemases for review, and the duplicated genes in the Tunisian genome.

### 4.6 What the Tools Add to Ensembl

Ensembl's genes for these genomes come from Prokka, and Bakta's genes match them closely: 99.5% have exactly the same coordinates. AMRFinderPlus found no new genes, but it gave the right name to genes that Ensembl names loosely or wrongly.

| AMRFinderPlus | Ensembl |
|---|---|
| *bla*CTX-M-15 | bla, described as CTX-M-1 |
| *bla*SHV-28 | bla, described as SHV-1 |
| *armA* | rmtB |
| *tet(D)* | tetA, class B |
| *qnrB1* | no name |

PlasmidFinder adds plasmid replicons, which Ensembl does not show at all. Kleborate's results have no positions, so they suit a genome summary rather than a gene track.

<figure>
  <img src="../results/tool_assessment/figures/p4_ensembl_PMK1_NDM_klebtrace_track.png" alt="Figure 2: KlebTrace custom track in the Ensembl browser" width="800">
  <figcaption>Figure 2: The KlebTrace track (brown) under Ensembl's genes (red) on pPMK1-NDM in the Ensembl Bacteria release 63 browser</figcaption>
</figure>

The track lines up with every Ensembl gene in the NDM region. In the scorecard, AMRFinderPlus scored 17 of 18 and Kleborate 10 of 18, while Bakta and PlasmidFinder scored 12 of 15 each, because the accuracy criterion does not apply to them. The full scorecard is in [`docs/tool_assessment.md`](tool_assessment.md). A short proposal for how these features could be shown in Ensembl is in [`docs/representation_note.md`](representation_note.md).

## 5. Conclusion and Future Directions

This project re-annotated seven *K. pneumoniae* genomes from the AMR Portal and compared the results with the portal's calls, the laboratory results and Ensembl. The main finding is that carbapenemase genes split across contigs are missing from the portal's calls. Running AMRFinderPlus on the DNA, together with a simple check that adds up gene pieces across contigs, found them in 2 of the 7 study genomes and in 5 of 37 strictly resistant isolates with no carbapenemase in the portal. AMRFinderPlus on the DNA explained the most laboratory results (6 of 7), and it would also give Ensembl more accurate resistance gene names. PlasmidFinder would add plasmid context at little cost.

The study has limits. It used only seven genomes, the split KPC-2 in Italy could not be confirmed with reads, and the explanation for the portal's missing calls comes from its method codes rather than from its pipeline. One resistant genome remains unexplained.

Future work includes running the split-gene check on the whole AMR Portal release, confirming split genes with reads or long-read assemblies, merging the Bakta and AMRFinderPlus features into one feature per gene for display, and looking for other causes of resistance in GCA_000567225.2, such as efflux pumps or AmpC enzymes.

---

## References

1. Carattoli, A., Zankari, E., García-Fernández, A., Voldby Larsen, M., Lund, O., Villa, L., Aarestrup, F.M. & Hasman, H. (2014). In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrobial Agents and Chemotherapy*, 58(7), 3895-3903.
2. Di Tommaso, P., Chatzou, M., Floden, E.W., Barja, P.P., Palumbo, E. & Notredame, C. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4), 316-319.
3. EMBL-EBI AMR Portal. Release 2026-07. https://www.ebi.ac.uk/amr
4. Feldgarden, M., Brover, V., Gonzalez-Escalona, N., Frye, J.G., Haendiges, J., Haft, D.H., ... & Klimke, W. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*, 11, 12728.
5. Lam, M.M.C., Wick, R.R., Watts, S.C., Cerdeira, L.T., Wyres, K.L. & Holt, K.E. (2021). A genomic surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its related species complex. *Nature Communications*, 12, 4188.
6. Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094-3100.
7. Schwengers, O., Jelonek, L., Dieckmann, M.A., Beyvers, S., Blom, J. & Goesmann, A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics*, 7(11), 000685.
8. Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. *Bioinformatics*, 30(14), 2068-2069.
9. The Galaxy Community (2024). The Galaxy platform for accessible, reproducible, and collaborative data analyses: 2024 update. *Nucleic Acids Research*, 52(W1), W83-W94.
