# KlebTrace

**A reproducible annotation and data-integration pipeline for investigating carbapenem resistance in *Klebsiella pneumoniae***

> Status: pilot in progress (seven genomes), built to scale to the full EMBL-EBI AMR Portal release.

**Report:** [docs/report.md](docs/report.md)

## Overview

*Klebsiella pneumoniae* is a WHO critical-priority pathogen, and carbapenem resistance in this species leaves few treatment options. The [EMBL-EBI AMR Portal](https://www.ebi.ac.uk/amr) links laboratory antimicrobial susceptibility results to genome-based resistance calls for the same isolates.

KlebTrace uses that resource to ask three questions:

1. **Annotation agreement.** When the same genomes are re-annotated independently (Bakta, AMRFinderPlus, Kleborate), are the resistance genes reported by the AMR Portal recovered?
2. **Genotype-phenotype concordance.** Do the detected genes explain each isolate's laboratory meropenem result, and what explains resistant isolates with no known carbapenemase gene?
3. **Genomic context.** Are carbapenemase genes located on plasmids or on the chromosome?

## Workflow

```
AMR Portal (Parquet)  -->  select isolates + extract portal calls   (scripts/query_portal.py)
        |
        v
ENA assemblies (FASTA) -->  annotate: Bakta, AMRFinderPlus, Kleborate, PlasmidFinder
        |                        (prototyped in galaxy/, implemented in nextflow/)
        v
compare calls + phenotypes -->  relational database -->  results/ and docs/
```

## Repository structure

| Folder | Contents |
|---|---|
| `data/` | Selected isolates (`samples.tsv`) and the portal's resistance calls (`portal_calls.tsv`) |
| `galaxy/` | Galaxy workflow prototype, tool versions and screenshots |
| `nextflow/` | Reproducible pipeline (in development) |
| `scripts/` | Python scripts for querying, downloading, parsing, comparing and loading |
| `database/` | Database schema and example queries |
| `tests/` | Unit tests and small test inputs |
| `results/` | Summary tables and figures |
| `docs/` | Project report (`report.md`) and tool assessment (`tool_assessment.md`) |

## Quick start

```bash
conda env create -f environment.yml
conda activate klebtrace

# 1. Inspect the AMR Portal Parquet schema (check column names first)
python scripts/query_portal.py --inspect

# 2. Summarise K. pneumoniae meropenem isolates and select six study genomes
python scripts/query_portal.py --outdir data/ \
    --ensembl-list https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/species_EnsemblBacteria.txt \
    --prefer-regions Africa

# 3. Download the genomes, re-check the portal-negative isolates for split or
#    partial carbapenemase genes, and add one confirmed negative (seventh genome)
python scripts/download_genomes.py
python scripts/screen_carbapenemases.py --amr-db <Kleborate_AMRdb.csv> --amr-fasta <Kleborate_AMRdb.fasta>
kleborate -a data/screen_genomes/*.fasta -o results/kleborate_negatives \
    -m general__species,general__contig_stats,kpsc__mlst,kpsc__amr
python scripts/finalise_study_set.py --kleborate results/kleborate_negatives/general__species_output.txt
python scripts/download_genomes.py

# 4. Summarise the Galaxy annotation outputs (downloaded to results/phase2_galaxy/)
python scripts/summarise_galaxy_outputs.py

# 5. Run tests
pytest tests/
```

## Data sources

- EMBL-EBI AMR Portal release files: https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/
- Genome assemblies: European Nucleotide Archive (ENA)

## Author

Chioma B Oselu

## Licence

Code is released under the MIT Licence. Data remain subject to the terms of their original sources.
