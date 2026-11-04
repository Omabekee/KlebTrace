// Kleborate: species, MLST, resistance and virulence for the K. pneumoniae species complex.
// Not available on usegalaxy.org; same beta version as the local Phase 2 run.
// kpsc__cgmlst is left out because its database is a separate download.
process KLEBORATE {
    tag "$assembly"
    container 'quay.io/biocontainers/kleborate:3.3.0b1--pyhdfd78af_0'
    conda 'bioconda::kleborate=3.3.0b1'
    publishDir "${params.outdir}/kleborate", mode: 'copy'

    input:
    tuple val(assembly), path(fasta)

    output:
    tuple val(assembly), path("${assembly}.kleborate.tsv"), emit: report

    script:
    """
    kleborate -a $fasta -o out -m ${params.kleborate_modules}
    mv out/*_output.txt ${assembly}.kleborate.tsv
    """

    stub:
    """
    touch ${assembly}.kleborate.tsv
    """
}
