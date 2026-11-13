process SEGMENTQUANTIZED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(summary_tsv)
    tuple val(meta2), path(quantized_bed)
    tuple val(meta3), path(gaps_bed)
    tuple val(meta4), path(red_bed)
    tuple val(meta5), path(regions_bed)
    tuple val(meta6), path(breakpoints_tsv)
    val chromosomes

    output:
    tuple val(meta), path("*.pde_segments.bed"), emit: bed
    tuple val(meta), path("*.grs_summary.tsv"),  emit: summary
    tuple val("${task.process}"), val('python'), val('3.10'), emit: versions_segmentquantized, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def red_arg  = (red_bed && red_bed != "" && red_bed != "null") ? "--red-bed ${red_bed}" : ""
    def chr_val  = meta.chromosomes ?: chromosomes
    def chr_arg  = (chr_val && chr_val != "" && chr_val != "null") ? "--chromosomes \"${chr_val}\"" : ""

    """
    # Tool: pde_segment_quantized.py | Script-SHA: 645d9ae26b722a89
    python3 ${moduleDir}/resources/usr/bin/pde_segment_quantized.py --version
    python3 ${moduleDir}/resources/usr/bin/pde_segment_quantized.py \\
        --summary-tsv ${summary_tsv} \\
        --quantized-bed ${quantized_bed} \\
        --gaps-bed ${gaps_bed} \\
        ${red_arg} \\
        --regions-bed ${regions_bed} \\
        --breakpoints-tsv ${breakpoints_tsv} \\
        ${chr_arg} \\
        --prefix ${prefix} \\
        --outdir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pde_segments.bed
    touch ${prefix}.grs_summary.tsv
    """
}
