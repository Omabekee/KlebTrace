// KlebTrace's own Python steps. They run in the klebtrace conda environment
// (environment.yml), using the scripts in ../scripts.

process MERGE_ANNOTATIONS {
    tag "$assembly"
    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    tuple val(assembly), path(fasta), path(bakta_gff3), path(amrfinder), path(kleborate), path(plasmidfinder)

    output:
    tuple val(assembly), path("${assembly}.klebtrace.gff3"), emit: gff3

    script:
    """
    merge_annotations.py --assembly $assembly --fasta $fasta --bakta $bakta_gff3 \\
        --amrfinderplus $amrfinder --kleborate $kleborate --plasmidfinder $plasmidfinder \\
        --out ${assembly}.klebtrace.gff3
    """

    stub:
    """
    printf '##gff-version 3\\n' > ${assembly}.klebtrace.gff3
    """
}

process LOAD_DATABASE {
    publishDir "${params.outdir}/database", mode: 'copy'

    input:
    path schema
    path samples
    path manifest
    path portal_calls
    path screen
    path fastas
    path gff3s
    path amrfinder_reports
    path kleborate_reports

    output:
    path "klebtrace.sqlite", emit: db

    script:
    """
    load_database.py --schema $schema --samples $samples --manifest $manifest \\
        --portal-calls $portal_calls --screen $screen --fasta $fastas \\
        --gff3 $gff3s --amrfinderplus $amrfinder_reports --kleborate $kleborate_reports \\
        --out klebtrace.sqlite
    """

    stub:
    """
    touch klebtrace.sqlite
    """
}

process QC_CHECKS {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path db
    path fastas
    path gff3s

    output:
    path "qc_report.tsv", emit: report

    script:
    """
    qc_checks.py --database $db --fasta $fastas --gff3 $gff3s --out qc_report.tsv
    """

    stub:
    """
    printf 'check\\tassembly\\tstatus\\tdetail\\n' > qc_report.tsv
    """
}
