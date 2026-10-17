process TRANSCHROMOSOMALLINKRESOLVER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_minimap2_python:c9b2a6b2ba2d6d79'

    input:
    tuple val(meta), path(blueprint_json), path(bam), path(bai), path(ref_fasta), path(gaps_bed), path(policy_json)

    output:
    tuple val(meta), path("*_trans_manifest.tsv"), emit: manifest
    tuple val("${task.process}"), val('python'), val('3.10.4'), emit: versions_trans, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def gaps_arg = gaps_bed ? "--gaps-bed \"${gaps_bed}\"" : ""
    def pol_arg  = policy_json ? "--policy \"${policy_json}\"" : ""
    """
    # Tool: pde_transchromosomal_link_resolver.py | Script-SHA: 7306433275f90d6d
    pde_transchromosomal_link_resolver.py \\
        --blueprint "${blueprint_json}" \\
        --bam "${bam}" \\
        --ref "${ref_fasta}" \\
        ${gaps_arg} \\
        ${pol_arg} \\
        --sample "${prefix}" \\
        --manifest "${prefix}_trans_manifest.tsv" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_trans_manifest.tsv"
    """
}
