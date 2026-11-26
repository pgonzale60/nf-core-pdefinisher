process RDNACASSETTERESOLVER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(blueprint_json), path(ref_fasta), path(flye_fasta, stageAs: "rdna_asm/*"), path(flye_info), path(flye_paf), path(telomeric_unitigs, stageAs: "telo_utg/*"), path(rdna_gfa)

    output:
    tuple val(meta), path("*_rdna_cassettes"), path("*_rdna_manifest.tsv"), emit: manifest
    tuple val("${task.process}"), val('python'), val('3.10.4'), emit: versions_resolver, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def fasta_arg = flye_fasta ? "--flye-fasta \"${flye_fasta}\"" : ""
    def info_arg = flye_info ? "--flye-info \"${flye_info}\"" : ""
    def paf_arg = flye_paf ? "--flye-paf \"${flye_paf}\"" : ""
    def telo_arg = telomeric_unitigs ? "--telomeric-unitigs \"${telomeric_unitigs}\"" : ""
    def gfa_arg  = rdna_gfa ? "--rdna-gfa \"${rdna_gfa}\"" : ""
    """
    # Tool: pde_rdna_cassette_resolver.py | Script-SHA: dc53fafd7c32c60e
    pde_rdna_cassette_resolver.py \\
        --blueprint "${blueprint_json}" \\
        --ref "${ref_fasta}" \\
        ${fasta_arg} \\
        ${info_arg} \\
        ${paf_arg} \\
        ${telo_arg} \\
        ${gfa_arg} \\
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
