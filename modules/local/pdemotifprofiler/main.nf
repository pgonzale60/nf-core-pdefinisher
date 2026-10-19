process PDEMOTIFPROFILER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(fimo_tsv), path(breakpoints_bed), path(pde_segments_bed)
    val flank_tolerance

    output:
    tuple val(meta), path("*.fimo_annotated_hits.tsv")     , emit: annotated_hits
    tuple val(meta), path("*.motif_distribution_summary.tsv"), emit: summary_tsv
    tuple val(meta), path("*.motif_distribution_report.md") , emit: report_md
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def tolerance = flank_tolerance ?: 500

    """
    # Tool: pde_motif_profiler.py | Script-SHA: 1dfde1711505874d
    python3 ${moduleDir}/resources/usr/bin/pde_motif_profiler.py \\
        --fimo-tsv "${fimo_tsv}" \\
        --breakpoints-bed "${breakpoints_bed}" \\
        --pde-segments-bed "${pde_segments_bed}" \\
        --flank-tolerance ${tolerance} \\
        --prefix "${prefix}" \\
        --outdir "." \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.fimo_annotated_hits.tsv"
    touch "${prefix}.motif_distribution_summary.tsv"
    touch "${prefix}.motif_distribution_report.md"
    """
}
