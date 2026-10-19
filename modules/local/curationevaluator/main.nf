process CURATIONEVALUATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(pre_ref), path(post_ref), path(blueprint_json), path(curated_agp), path(curation_audit), path(gap_manifest), path(walker_manifest), path(rdna_manifest), path(terminal_summary), path(trans_manifest), path(red_bed)

    output:
    tuple val(meta), path("*.curation_scorecard.json"), emit: scorecard_json
    tuple val(meta), path("*.curation_report.md"), emit: report_md
    tuple val("${task.process}"), val('python'), val('3.10.4'), emit: versions_evaluator, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def trans_arg = trans_manifest ? "--trans-manifest \"${trans_manifest}\"" : ""
    def red_arg   = red_bed ? "--red-bed \"${red_bed}\"" : ""
    """
    # Tool: pde_evaluate_curation.py | Script-SHA: 846f6145bcd0e4bc
    pde_evaluate_curation.py --version
    pde_evaluate_curation.py \\
        --pre-ref "${pre_ref}" \\
        --post-ref "${post_ref}" \\
        --blueprint "${blueprint_json}" \\
        --agp "${curated_agp}" \\
        --audit "${curation_audit}" \\
        --gap-manifest "${gap_manifest}" \\
        --walker-manifest "${walker_manifest}" \\
        --rdna-manifest "${rdna_manifest}" \\
        --terminal-summary "${terminal_summary}" \\
        ${trans_arg} \\
        ${red_arg} \\
        --sample "${prefix}" \\
        --outdir "." \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.curation_scorecard.json"
    touch "${prefix}.curation_report.md"
    """
}
