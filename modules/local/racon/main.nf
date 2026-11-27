process RACON_CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/racon:1.5.0--h077b44d_8':
        'quay.io/biocontainers/racon:1.5.0--h077b44d_8' }"

    input:
    tuple val(meta), path(prefixes), path(assigned_paf), path(representatives)

    output:
    tuple val(meta), path("consensus.fasta"), emit: consensus
    tuple val("${task.process}"), val('racon'), eval("racon --version 2>&1 | sed 's/^v//'"), topic: versions, emit: versions_racon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    if [[ -s "${assigned_paf}" && -s "${representatives}" ]]; then
        racon ${args} -t ${task.cpus} "${prefixes}" "${assigned_paf}" "${representatives}" > consensus.fasta
    else
        cp "${representatives}" consensus.fasta
    fi
    """

    stub:
    """
    : > consensus.fasta
    """
}
