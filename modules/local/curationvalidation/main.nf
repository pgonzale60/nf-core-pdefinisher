process CURATIONVALIDATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(audit), path(agp), path(sniffles_vcf), path(sniffles_tbi), path(bam), path(bai)

    output:
    tuple val(meta), path("*.curation_validation.json"), emit: validation_json
    tuple val(meta), path("*.curation_validation.md")  , emit: validation_md
    tuple val("${task.process}"), val('python'), val('3.10.4'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), val('0.23.3'), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def agp_arg = agp ? "--agp \"${agp}\"" : ""
    def bam_arg = bam ? "--bam \"${bam}\"" : ""
    """
    # Tool: pde_validate_curation.py | Script-SHA: 67254f539848ab04
    pde_validate_curation.py --version
    pde_validate_curation.py \\
        --audit "${audit}" \\
        ${agp_arg} \\
        --vcf "${sniffles_vcf}" \\
        ${bam_arg} \\
        --sample "${prefix}" \\
        --outdir "." \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.curation_validation.json"
    touch "${prefix}.curation_validation.md"
    """
}
