process EXTRACTPDEBORDERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(pde_segments_bed), path(breakpoints_bed), path(assembly)
    val window_small
    val window_large

    output:
    tuple val(meta), path("*.pde_border_breakpoints.bed")  , emit: border_breakpoints
    tuple val(meta), path("*.pde_border_flanks_150bp.fasta"), emit: flanks_150bp
    tuple val(meta), path("*.pde_border_flanks_350bp.fasta"), emit: flanks_350bp
    tuple val(meta), path("*.pde_border_flanks.bed")        , emit: flanks_bed
    tuple val(meta), path("*.break_windows_300bp.fasta")    , emit: break_windows_fasta
    tuple val(meta), path("*.census.json")                  , emit: census_json
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def win_s  = window_small ?: 150
    def win_l  = window_large ?: 350
    def bp_arg = breakpoints_bed ? "--breakpoints-bed ${breakpoints_bed}" : ""

    """
    # Tool: pde_extract_borders.py | Script-SHA: b2447cc0959be3d8
    python3 ${moduleDir}/resources/usr/bin/pde_extract_borders.py \\
        --pde-segments-bed "${pde_segments_bed}" \\
        --assembly "${assembly}" \\
        ${bp_arg} \\
        --window-small ${win_s} \\
        --window-large ${win_l} \\
        --prefix "${prefix}" \\
        --outdir "." \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.pde_border_breakpoints.bed"
    touch "${prefix}.pde_border_flanks_150bp.fasta"
    touch "${prefix}.pde_border_flanks_350bp.fasta"
    touch "${prefix}.pde_border_flanks.bed"
    touch "${prefix}.break_windows_300bp.fasta"
    echo '{"n_breaks": 0}' > "${prefix}.census.json"
    """
}
