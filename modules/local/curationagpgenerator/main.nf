process CURATIONAGPGENERATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple(
        val(meta),
        path(blueprint_json),
        path(ref_fasta),
        path(gap_patches),
        path(gap_manifest),
        path(walker_patches),
        path(walker_manifest),
        path(rdna_cassettes),
        path(rdna_manifest),
        path(trans_manifest)
    )

    output:
    tuple val(meta), path("*.curated.agp"), emit: agp
    tuple val(meta), path("*.components.fasta"), emit: components_fasta
    tuple val(meta), path("*.curated.fasta"), emit: curated_fasta
    tuple val(meta), path("*.haplotigs.fasta"), path("*.purged_rdna.fasta"), optional: true, emit: segregated_fastas
    tuple val(meta), path("*.curation_audit.tsv"), emit: audit
    tuple val("${task.process}"), val('python'), val('3.10.4'), emit: versions_agp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def gap_arg     = gap_manifest ? "--gap-manifest \"${gap_manifest}\"" : ""
    def walker_arg  = walker_manifest ? "--walker-manifest \"${walker_manifest}\"" : ""
    def rdna_arg    = rdna_manifest ? "--rdna-manifest \"${rdna_manifest}\"" : ""
    def trans_arg   = trans_manifest ? "--trans-manifest \"${trans_manifest}\"" : ""
    """
    # Tool: pde_generate_agp.py | Script-SHA: 8d7d12b2fec32630
    pde_generate_agp.py \\
        --blueprint "${blueprint_json}" \\
        --ref "${ref_fasta}" \\
        ${gap_arg} \\
        ${walker_arg} \\
        ${rdna_arg} \\
        ${trans_arg} \\
        --sample "${prefix}" \\
        --outdir "." \\
        --synthesize-fasta \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.curated.agp"
    touch "${prefix}.components.fasta"
    touch "${prefix}.curated.fasta"
    touch "${prefix}.haplotigs.fasta"
    touch "${prefix}.purged_rdna.fasta"
    touch "${prefix}.curation_audit.tsv"
    """
}
