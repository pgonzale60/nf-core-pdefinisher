process MODDOTPLOT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/moddotplot:0.9.9--pyhdfd78af_0' :
        'quay.io/biocontainers/moddotplot:0.9.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_moddotplot"), emit: plots
    tuple val("${task.process}"), val('moddotplot'), eval("moddotplot --version 2>&1 | head -n 1 | sed 's/moddotplot //'"), topic: versions, emit: versions_moddotplot

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    moddotplot static \\
        -f "${fasta}" \\
        -o "${prefix}_moddotplot" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_moddotplot"
    touch "${prefix}_moddotplot/dummy.png"
    """
}
