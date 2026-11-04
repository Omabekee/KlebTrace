# KlebTrace Nextflow pipeline

Rebuilds the Phase 2 Galaxy workflow (v0.2) as a reproducible pipeline, using the same tool versions, then merges, loads and checks the results.

```
samples.tsv + FASTA --+--> BAKTA -----------+
                      +--> AMRFINDERPLUS ---+
                      +--> KLEBORATE -------+--> MERGE_ANNOTATIONS --> LOAD_DATABASE --> QC_CHECKS
                      +--> PLASMIDFINDER ---+    (one GFF3/genome)     (SQLite)          (qc_report.tsv)
```

## Status

| Part | State |
|---|---|
| `main.nf`, `nextflow.config`, tool modules | Run end to end on all seven genomes (4 November 2026) |
| Tool containers | Pinned to the Galaxy versions (below); AMRFinderPlus and PlasmidFinder pulled as Apptainer images |
| `database/schema.sql`, `load_database.py` | Database built from the pipeline's own outputs |
| `merge_annotations.py`, `qc_checks.py` | QC report: 50 PASS, 13 WARN, 0 FAIL |
| Reproducibility | AMRFinderPlus and PlasmidFinder tables identical, line for line, to the Galaxy outputs for all seven genomes |
| Run time | 13 min 54 s for seven genomes on a 4-CPU laptop, 0.5 GB peak memory (`results/efficiency/`), Bakta excluded |
| Reference databases | AMRFinderPlus 2026-05-15.1 and PlasmidFinder 81c11f4 downloaded (865 MB with the images); Bakta reused from Galaxy |

## Tools and containers

| Process | Container (quay.io/biocontainers) | Galaxy run |
|---|---|---|
| BAKTA | `bakta:1.9.4--pyhdfd78af_0` | 1.9.4+galaxy3, DB v5.1 full |
| AMRFINDERPLUS | `ncbi-amrfinderplus:4.2.7--hf69ffd2_0` | 4.2.7+galaxy0, DB 2026-05-15.1 |
| KLEBORATE | `kleborate:3.3.0b1--pyhdfd78af_0` (the `laptop` profile uses the local Bioconda install instead) | local Bioconda 3.3.0b1 |
| PLASMIDFINDER | `plasmidfinder:2.1.6--py314hdfd78af_2` | 2.1.6+galaxy1, same container |

AMRFinderPlus runs on the nucleotide sequence (`--nucleotide`), as in Galaxy, because the AMR Portal's protein-based calls miss genes split at contig ends.

## Running

From the repository root:

```bash
conda activate klebtrace
bash scripts/download_databases.sh          # once, about 0.9 GB

# Wiring check: every step runs with placeholder outputs, no tools or databases needed
nextflow run nextflow/main.nf -profile test -stub-run

# Real run on a laptop: Apptainer images for AMRFinderPlus and PlasmidFinder,
# Kleborate from the local conda environment (avoids its 1.5 GB image)
nextflow run nextflow/main.nf -profile laptop
```

Profiles: `laptop`, `singularity`, `docker`, `conda`, and `test` (two genomes: the complete NDM genome and the susceptible control). Run time per step is written to `results/pipeline/pipeline_info/` (timeline, report and trace).

## Reference databases

Downloaded into `databases/` (not in Git) by `scripts/download_databases.sh`, with paths set in `nextflow.config`.

| Database | Version | Size on disk |
|---|---|---|
| AMRFinderPlus | `2026-05-15.1`, indexed with `amrfinder_index` | 246 MB |
| PlasmidFinder | commit `81c11f4` (2023-12-04), as in Galaxy | under 1 MB |
| AMRFinderPlus and PlasmidFinder Apptainer images | as in the table above | 620 MB |
| Bakta | v5.1 full, as in Galaxy | not downloaded (39.4 GB) |

**Bakta.** The full database is too large for the development laptop, so by default the pipeline reuses the Bakta GFF3 files from the Galaxy run (same Bakta version and database). They are not in Git: download the "Annotation and sequences (GFF3)" collection from the Galaxy history and place one file per genome in `results/phase2_galaxy/bakta_gff3/<assembly>.gff3`. To run Bakta itself, for example on a cluster, use `--bakta_gff3_dir false --bakta_db <path>`.
