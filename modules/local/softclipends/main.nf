process SOFTCLIPENDS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(bam), path(bai)
    val telomere_repeat

    output:
    tuple val(meta), path("*.terminal_softclips.tsv")   , emit: softclips
    tuple val(meta), path("*.terminal_ends_summary.tsv") , emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), val('0.23.3'), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Tool: pde_terminal_softclips.py | Script-SHA: c349e9ba65b9c698
    # softclipends v1.2 - FastaReader fallback for plain-gzip FASTA
    pde_terminal_softclips.py \\
        --bam "${bam}" \\
        --assembly "${assembly}" \\
        --motif "${telomere_repeat}" \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.terminal_softclips.tsv
    touch ${prefix}.terminal_ends_summary.tsv
    """
}
