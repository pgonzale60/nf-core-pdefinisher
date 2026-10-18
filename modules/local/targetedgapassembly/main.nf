process TARGETEDGAPASSEMBLY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/hifiasm_abpoa_minimap2_samtools_mindthegap_python:4e99aa7a05e907bf'

    input:
    tuple val(meta), path(blueprint_json), path(bam), path(bai), path(ref_fasta), path(ref_gzi), path(policy_json)

    output:
    tuple val(meta), path("*_gap_patches"), path("*_patches_manifest.tsv"), emit: manifest
    tuple val("${task.process}"), val('hifiasm'), val('0.24.0'), emit: versions_assembly, topic: versions
    tuple val("${task.process}"), val('abpoa'), val('1.5.7'), emit: versions_consensus, topic: versions
    tuple val("${task.process}"), val('mindthegap'), val('2.3.0'), emit: versions_mindthegap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def policy_opt = policy_json ? "--policy ${policy_json}" : ''
    """
    # Tool: pde_targeted_gap_assembly.py | Script-SHA: 30e0bdbfa772e492
    pde_targeted_gap_assembly.py \\
        --blueprint "${blueprint_json}" \\
        --bam "${bam}" \\
        --ref "${ref_fasta}" \\
        --sample "${prefix}" \\
        --outdir "${prefix}_gap_patches" \\
        --manifest "${prefix}_patches_manifest.tsv" \\
        --threads ${task.cpus} \\
        ${policy_opt} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_gap_patches"
    touch "${prefix}_patches_manifest.tsv"
    """
}
