process RDNACASSETTERESOLVER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(blueprint_json), path(ref_fasta), path(flye_fasta, stageAs: "rdna_asm/*"), path(flye_info), path(flye_paf), path(telomeric_unitigs, stageAs: "telo_utg/*")

    output:
    tuple val(meta), path("*_rdna_cassettes"), path("*_rdna_manifest.tsv"), emit: manifest
    tuple val("${task.process}"), val('python'), val('3.10.4'), emit: versions_resolver, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def telo_arg = telomeric_unitigs ? "--telomeric-unitigs \"${telomeric_unitigs}\"" : ""
    """
    # Tool: pde_rdna_cassette_resolver.py | Script-SHA: 697077972e7da2e9
    pde_rdna_cassette_resolver.py \\
        --blueprint "${blueprint_json}" \\
        --ref "${ref_fasta}" \\
        --flye-fasta "${flye_fasta}" \\
        --flye-info "${flye_info}" \\
        --flye-paf "${flye_paf}" \\
        ${telo_arg} \\
        --sample "${prefix}" \\
        --outdir "${prefix}_rdna_cassettes" \\
        --manifest "${prefix}_rdna_manifest.tsv" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_rdna_cassettes"
    touch "${prefix}_rdna_manifest.tsv"
    """
}
