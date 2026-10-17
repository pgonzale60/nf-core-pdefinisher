process BGZIPFASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_minimap2_python:c9b2a6b2ba2d6d79'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa.gz"),     emit: fasta
    // Neither index is optional. `bgzip --index` always writes the .gzi, and `samtools
    // faidx` always writes the .fai, because the .gzi it needs is already there. The
    // resolution engines read regions through samtools and have no fallback, so a missing
    // index must fail here, at the process that should have produced it.
    tuple val(meta), path("*.fa.gz.gzi"), emit: gzi
    tuple val(meta), path("*.fa.gz.fai"), emit: fai
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_bgzip
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Assemblies arrive as plain gzip (or uncompressed). The engines that read regions do
    # so through `samtools faidx`, which can only index a bgzipped file, and random access
    # to that file also needs the .gzi that bgzip writes beside it. Normalising once here
    # means no engine has to decompress the assembly into its own work directory.
    #
    # gzi comes from bgzip; fai from samtools faidx, which is possible precisely because
    # the .gzi exists by then.
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c "${fasta}" > "${prefix}.fa"
    else
        cat "${fasta}" > "${prefix}.fa"
    fi
    bgzip --threads ${task.cpus} --index ${args} "${prefix}.fa"
    samtools faidx "${prefix}.fa.gz"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.fa.gz"
    touch "${prefix}.fa.gz.gzi"
    touch "${prefix}.fa.gz.fai"
    """
}
