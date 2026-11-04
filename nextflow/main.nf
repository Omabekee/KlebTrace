#!/usr/bin/env nextflow
// KlebTrace pipeline: annotate every genome in samples.tsv with Bakta,
// AMRFinderPlus, Kleborate and PlasmidFinder, merge the results into one GFF3
// per genome, load everything into SQLite and run the quality checks.
//
//   nextflow run nextflow/main.nf -profile docker      (or singularity, conda)
//   nextflow run nextflow/main.nf -profile test -stub-run   (wiring check, no tools)

nextflow.enable.dsl = 2

include { BAKTA }         from './modules/bakta'
include { AMRFINDERPLUS } from './modules/amrfinderplus'
include { KLEBORATE }     from './modules/kleborate'
include { PLASMIDFINDER } from './modules/plasmidfinder'
include { MERGE_ANNOTATIONS; LOAD_DATABASE; QC_CHECKS } from './modules/klebtrace'

workflow {
    // One (assembly, fasta) pair per row of samples.tsv
    genomes = Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.assembly, file("${params.genomes}/${row.assembly}.fasta", checkIfExists: true)) }

    // Bakta's full database is about 40 GB, too large for this laptop. By
    // default the pipeline reuses the Bakta GFF3 files from the Galaxy run
    // (same Bakta version and database); set --bakta_gff3_dir false and
    // --bakta_db to run Bakta here instead (for example on an HPC cluster).
    if (params.bakta_gff3_dir) {
        bakta_gff3 = genomes.map { assembly, fasta ->
            tuple(assembly, file("${params.bakta_gff3_dir}/${assembly}.gff3", checkIfExists: !workflow.stubRun))
        }
    } else {
        BAKTA(genomes, file(params.bakta_db, checkIfExists: !workflow.stubRun))
        bakta_gff3 = BAKTA.out.gff3
    }
    AMRFINDERPLUS(genomes, file(params.amrfinder_db, checkIfExists: !workflow.stubRun))
    KLEBORATE(genomes)
    PLASMIDFINDER(genomes, file(params.plasmidfinder_db, checkIfExists: !workflow.stubRun))

    merged = MERGE_ANNOTATIONS(
        genomes
            .join(bakta_gff3)
            .join(AMRFINDERPLUS.out.report)
            .join(KLEBORATE.out.report)
            .join(PLASMIDFINDER.out.report)
    )

    db = LOAD_DATABASE(
        file("${projectDir}/../database/schema.sql", checkIfExists: true),
        file(params.samples),
        file(params.manifest, checkIfExists: true),
        file(params.portal_calls, checkIfExists: true),
        file(params.screen, checkIfExists: true),
        genomes.map { it[1] }.collect(),
        merged.gff3.map { it[1] }.collect(),
        AMRFINDERPLUS.out.report.map { it[1] }.collect(),
        KLEBORATE.out.report.map { it[1] }.collect()
    )

    QC_CHECKS(db.db, genomes.map { it[1] }.collect(), merged.gff3.map { it[1] }.collect())
}
