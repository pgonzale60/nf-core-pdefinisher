process TIARA {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiara:1.0.3':
        'quay.io/biocontainers/tiara:1.0.3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tiara.txt"), emit: classifications
    tuple val("${task.process}"), val('tiara'), eval("tiara --version 2>&1 | awk '{print \$NF}'"), topic: versions, emit: versions_tiara

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    def input_fa = fasta.name.endsWith('.gz') ? fasta.name.replaceAll(/\.gz$/, '') : fasta.name
    def gunzip   = fasta.name.endsWith('.gz') ? "gzip -dc \"${fasta}\" > \"${input_fa}\"" : ""
    """
    ${gunzip}

    tiara \\
        -i "${input_fa}" \\
        -o "${prefix}.tiara.txt" \\
        ${args} \\
        -t ${task.cpus}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sequence_id\\tclass_fst_stage\\tclass_snd_stage\\n' > "${prefix}.tiara.txt"
    """
}
