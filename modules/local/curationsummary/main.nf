process CURATIONSUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(scorecards), path(validations)

    output:
    tuple val(meta), path("*.curation_summary.md"), emit: summary_md
    tuple val(meta), path("*.curation_summary.tsv"), emit: summary_tsv
    tuple val(meta), path("*.curation_summary.json"), emit: summary_json
    tuple val("${task.process}"), val('python'), val('3.10.4'), emit: versions_summary, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: (meta.id ?: 'cohort')
    def val_arg = validations ? "--validations ${validations}" : ""
    """
    # Tool: pde_curation_summary.py | Script-SHA: 551cedaa6b731a90
    pde_curation_summary.py --version
    pde_curation_summary.py \\
        --scorecards ${scorecards} \\
        ${val_arg} \\
        --prefix "${prefix}" \\
        --outdir "." \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: (meta.id ?: 'cohort')
    """
    touch "${prefix}.curation_summary.md"
    touch "${prefix}.curation_summary.tsv"
    touch "${prefix}.curation_summary.json"
    """
}
