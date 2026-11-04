// Bakta: general genome annotation (genes, RNAs, ORFs) with GFF3 output.
// Same version as the Galaxy run (1.9.4); contig headers are kept, as in Galaxy.
process BAKTA {
    tag "$assembly"
    container 'quay.io/biocontainers/bakta:1.9.4--pyhdfd78af_0'
    conda 'bioconda::bakta=1.9.4'
    publishDir "${params.outdir}/bakta", mode: 'copy'

    input:
    tuple val(assembly), path(fasta)
    path db

    output:
    tuple val(assembly), path("${assembly}/${assembly}.gff3"), emit: gff3
    tuple val(assembly), path("${assembly}/${assembly}.tsv"),  emit: tsv
    path "${assembly}/${assembly}.txt", emit: summary

    script:
    """
    bakta --db $db --keep-contig-headers --translation-table 11 \\
        --threads $task.cpus --prefix $assembly --output $assembly $fasta
    """

    stub:
    """
    mkdir -p $assembly
    printf '##gff-version 3\\n' > $assembly/${assembly}.gff3
    touch $assembly/${assembly}.tsv $assembly/${assembly}.txt
    """
}
