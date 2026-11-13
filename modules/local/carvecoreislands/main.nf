process CARVECOREISLANDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_minimap2_python:c9b2a6b2ba2d6d79'

    input:
    tuple val(meta), path(summary_tsv)
    tuple val(meta2), path(quantized_bed)
    tuple val(meta3), path(gaps_bed)
    tuple val(meta4), path(red_bed)
    tuple val(meta5), path(core_coverage_tsv)

    output:
    tuple val(meta), path("*.pde_segments.bed"), emit: bed
    tuple val(meta), path("*.grs_summary.tsv"),  emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Tool: pde_carve_core_islands.py | Script-SHA: 1e139669e3db09f0
    pde_carve_core_islands.py \\
        --summary-tsv ${summary_tsv} \\
        --quantized-bed ${quantized_bed} \\
        --gaps-bed ${gaps_bed} \\
        --red-bed ${red_bed} \\
        --core-coverage-tsv ${core_coverage_tsv} \\
        --prefix ${prefix} \\
        --outdir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.pde_segments.bed"
    touch "${prefix}.grs_summary.tsv"
    """
}
