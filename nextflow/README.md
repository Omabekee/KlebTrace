# KlebTrace Nextflow pipeline

Rebuilds the Phase 2 Galaxy workflow (v0.2) as a reproducible pipeline, using the same tool versions, then merges, loads and checks the results.

```
samples.tsv + FASTA --+--> BAKTA -----------+
                      +--> AMRFINDERPLUS ---+
                      +--> KLEBORATE -------+--> MERGE_ANNOTATIONS --> LOAD_DATABASE --> QC_CHECKS
                      +--> PLASMIDFINDER ---+    (one GFF3/genome)     (SQLite)          (qc_report.tsv)
```

## Status (Phase 3 skeleton)

| Part | State |
|---|---|
| `main.nf`, `nextflow.config`, tool modules | Written; the wiring runs end to end with `-stub-run` |
| Tool containers | Pinned to the Galaxy versions (below) |
| `database/schema.sql` | Written and tested |
| `merge_annotations.py`, `load_database.py`, `qc_checks.py` | Written and tested; run on all seven genomes from the Galaxy outputs (`results/qc_report.tsv`: 50 PASS, 13 WARN, 0 FAIL) |
| Reference databases | Bakta reused from Galaxy (no 40 GB download); AMRFinderPlus and PlasmidFinder databases not yet downloaded |

## Tools and containers

| Process | Container (quay.io/biocontainers) | Galaxy run |
|---|---|---|
| BAKTA | `bakta:1.9.4--pyhdfd78af_0` | 1.9.4+galaxy3, DB v5.1 full |
| AMRFINDERPLUS | `ncbi-amrfinderplus:4.2.7--hf69ffd2_0` | 4.2.7+galaxy0, DB 2026-05-15.1 |
| KLEBORATE | `kleborate:3.3.0b1--pyhdfd78af_0` | local Bioconda 3.3.0b1 |
| PLASMIDFINDER | `plasmidfinder:2.1.6--py314hdfd78af_2` | 2.1.6+galaxy1, same container |

AMRFinderPlus runs on the nucleotide sequence (`--nucleotide`), as in Galaxy, because the AMR Portal's protein-based calls miss genes split at contig ends.

## Running

From the repository root:

```bash
conda activate klebtrace

# Wiring check: every step runs with placeholder outputs, no tools or databases needed
nextflow run nextflow/main.nf -profile test -stub-run

# Real run (WSL: Docker is not connected to WSL here, so use Apptainer/Singularity)
nextflow run nextflow/main.nf -profile singularity
```

Profiles: `singularity`, `docker`, `conda`, and `test` (two genomes: the complete NDM genome and the susceptible control). Run time per step is written to `results/pipeline/pipeline_info/` (timeline, report and trace), which gives the Phase 3 efficiency comparison with Galaxy.

## Reference databases

Stored in `databases/` (not in Git), with paths set in `nextflow.config`.

- **Bakta: not downloaded.** The full v5.1 database is a 39.4 GB download, too large for the development laptop. By default the pipeline reuses the Bakta GFF3 files from the Galaxy run (same Bakta version and database), placed in `results/phase2_galaxy/bakta_gff3/<assembly>.gff3`. To run Bakta itself, for example on a cluster, use `--bakta_gff3_dir false --bakta_db <path>`.
- **AMRFinderPlus database `2026-05-15.1`:** about 120 MB to download.
- **PlasmidFinder database** at commit `81c11f4` (2023-12-04), as in Galaxy: a few MB.
