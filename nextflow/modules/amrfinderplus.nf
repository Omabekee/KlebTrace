// AMRFinderPlus on the nucleotide sequence (not predicted proteins), so a gene
// broken at a contig end is still reported as partial. Settings match Galaxy.
process AMRFINDERPLUS {
    tag "$assembly"
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0'
    conda 'bioconda::ncbi-amrfinderplus=4.2.7'
    publishDir "${params.outdir}/amrfinderplus", mode: 'copy'

    input:
    tuple val(assembly), path(fasta)
    path db

    output:
    tuple val(assembly), path("${assembly}.amrfinderplus.tsv"), emit: report
    tuple val(assembly), path("${assembly}.mutation_all.tsv"),  emit: mutation_all

    script:
    """
    amrfinder --nucleotide $fasta --database $db \\
        --organism Klebsiella_pneumoniae --plus --report_common \\
        --mutation_all ${assembly}.mutation_all.tsv \\
        --coverage_min 0.5 --translation_table 11 --threads $task.cpus \\
        --output ${assembly}.amrfinderplus.tsv
    """

    stub:
    """
    touch ${assembly}.amrfinderplus.tsv ${assembly}.mutation_all.tsv
    """
}
