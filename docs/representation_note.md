# How resistance and virulence features could be shown in Ensembl

A short proposal based on what KlebTrace found when its annotation was compared with Ensembl Bacteria's (see [`tool_assessment.md`](tool_assessment.md)).

## The problem

Ensembl's bacterial genes come from Prokka. For resistance genes this gives a gene in the right place, but often with a loose or wrong name: CTX-M-15 is "bla" with the description "Beta-lactamase CTX-M-1", and *armA* is called "rmtB". Nothing tells the user that a gene confers resistance, to which drugs, or how sure the call is. Plasmid replicons and genome-level typing are not shown at all.

## Proposal

1. **Add resistance information to the existing gene, not a second gene.** AMRFinderPlus found no genes that Ensembl lacks: every resistance gene sat on an Ensembl gene. So the resistance call should be attached to that gene as extra attributes, rather than drawn as a separate feature. In the KlebTrace track, each resistance gene appears twice (once from Bakta, once from AMRFinderPlus), which is harder to read.

2. **Carry these attributes for each resistance gene:**

   | Attribute | Example | Why |
   |---|---|---|
   | Allele-level name | blaCTX-M-15 | The exact allele matters for treatment and surveillance |
   | Drug class and subclass | BETA-LACTAM, CEPHALOSPORIN | Tells the user what the gene does |
   | Evidence type | exact allele, partial, partial at contig end | Separates a sure call from a fragment |
   | Coverage and identity | 100% and 100% | Shows how close the match is |
   | Closest reference | RefSeq WP_000239590.1 | Lets the user check the source |
   | Tool and database version | AMRFinderPlus 4.2.7, database 2026-05-15.1 | Makes the call traceable and repeatable |
   | Link to the AMR Portal record | BioSample ID | Connects the gene to the laboratory result |

3. **Show point mutations as variants on the gene.** AMRFinderPlus reports porin and quinolone-target changes (for example *ompK36* D135DGD) with the position of the whole gene. They are better shown as a variant at the changed codon, attached to the gene, than as a gene of their own.

4. **Add plasmid replicons as a new feature type.** PlasmidFinder replicons are short regions (114 to 661 bp in these genomes) that show which plasmid family a sequence belongs to. Shown as their own track, they let a user see which plasmid carries each resistance gene.

5. **Put genome-level typing on the genome page.** Kleborate's species, sequence type, capsule and O type, and virulence and resistance scores have no position, so they belong with the genome's summary information, not in a gene track.

6. **Flag genes that are split or partial.** In two of the seven genomes, a carbapenemase gene was split across contig ends. The AMR Portal's protein-based calls missed both, and in the Kenyan genome Bakta predicted no gene at the split OXA-232 either. A nucleotide search with AMRFinderPlus found the pieces. These should be shown with a clear "partial" or "split across contigs" label, so that a missing gene is not read as a missing resistance mechanism.

7. **Use the same sequence names as Ensembl.** Ensembl names sequences "Chromosome" and by plasmid or contig name, while ENA and most tools use accessions (CP008933.1). A loader must map one to the other, or the features will not appear.

## What this would need

AMRFinderPlus on nucleotide input as the main tool, because it covers all bacteria, has a small database (246 MB) and explained the most laboratory results here; PlasmidFinder alongside it for replicons; a conversion step from their tables to GFF3 (as in `scripts/merge_annotations.py`); and a check that records the tool and database version with every call.
