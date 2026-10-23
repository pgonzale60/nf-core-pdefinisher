process IDENTIFYPUTATIVEGRS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(assembly)
    tuple val(meta2), path(delfies_bed)
    tuple val(meta3), path(regions_bed)
    tuple val(meta4), path(red_bed)

    output:
    tuple val(meta), path("*.pde_segments.bed")         , emit: bed
    tuple val(meta), path("*.grs_summary.tsv")          , emit: summary
    tuple val(meta), path("*.core_coverage.tsv")        , emit: core_coverage_tsv
    tuple val(meta), path("*.core_coverage.txt")        , emit: core_coverage
    tuple val("${task.process}"), val('python'), val('3.10'), emit: versions_identifyputativegrs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def bp_arg   = delfies_bed ? "--breakpoints-bed ${delfies_bed}" : ""
    def cov_arg  = regions_bed ? "--regions-bed ${regions_bed}" : ""
    def red_arg  = (red_bed && red_bed != "" && red_bed != "null") ? "--red-bed ${red_bed}" : ""

    """
    # Tool: pde_identify_grs.py | Script-SHA: ab690418e8a5ec4b
    python3 ${moduleDir}/resources/usr/bin/pde_identify_grs.py --version
    python3 ${moduleDir}/resources/usr/bin/pde_identify_grs.py \\
        --assembly "${assembly}" \\
        ${bp_arg} \\
        ${cov_arg} \\
        ${red_arg} \\
        --prefix "${prefix}" \\
        --outdir "." \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pde_segments.bed
    touch ${prefix}.grs_summary.tsv
    touch ${prefix}.core_coverage.tsv
    echo "0.0" > ${prefix}.core_coverage.txt
    """
}
