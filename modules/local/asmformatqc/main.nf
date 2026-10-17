process ASMFORMATQC {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/sanger-tol/agp-tpf-utils:1.3.4-c2'

    input:
    tuple val(meta), path(agp)

    output:
    tuple val(meta), path("*.validated.agp"), emit: agp
    tuple val("${task.process}"), val('agp-tpf-utils'), val('1.3.4'), topic: versions, emit: versions_agp_tpf_utils

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    asm-format \\
        --qc-overlaps \\
        --output-file "${prefix}.validated.agp" \\
        ${args} \\
        "${agp}"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.validated.agp"
    """
}
