process DELFIES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(bam), path(bai)

    output:
    tuple val(meta), path("*.breakpoint_locations.bed") , emit: bed
    tuple val(meta), path("*.breakpoint_sequences.fasta"), emit: fasta, optional: true
    tuple val("${task.process}"), val('delfies'), eval("delfies --version | sed 's/delfies, version //'"), topic: versions, emit: versions_delfies

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_assembly = assembly =~ /\.gz$/ ? assembly.name.take(assembly.name.lastIndexOf('.')) : assembly
    def gunzip_assembly = assembly =~ /\.gz$/ ? "gunzip -c ${assembly} > ${input_assembly}" : ""

    """
    $gunzip_assembly

    delfies \\
        --threads ${task.cpus} \\
        ${args} \\
        ${input_assembly} \\
        ${bam} \\
        ${prefix}_delfies

    mv ${prefix}_delfies/breakpoint_locations.bed ${prefix}.breakpoint_locations.bed
    if [ -f "${prefix}_delfies/breakpoint_sequences.fasta" ]; then
        mv ${prefix}_delfies/breakpoint_sequences.fasta ${prefix}.breakpoint_sequences.fasta
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.breakpoint_locations.bed
    touch ${prefix}.breakpoint_sequences.fasta
    """
}
