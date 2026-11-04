// PlasmidFinder: plasmid replicons. 2.1.6, as in Galaxy (3.0.3 failed there).
process PLASMIDFINDER {
    tag "$assembly"
    container 'quay.io/biocontainers/plasmidfinder:2.1.6--py314hdfd78af_2'
    conda 'bioconda::plasmidfinder=2.1.6'
    publishDir "${params.outdir}/plasmidfinder", mode: 'copy'

    input:
    tuple val(assembly), path(fasta)
    path db

    output:
    tuple val(assembly), path("${assembly}.plasmidfinder.tsv"), emit: report

    script:
    """
    mkdir -p out
    plasmidfinder.py -i $fasta -o out -p $db -l 0.6 -t 0.95 -x
    mv out/results_tab.tsv ${assembly}.plasmidfinder.tsv
    """

    stub:
    """
    touch ${assembly}.plasmidfinder.tsv
    """
}
