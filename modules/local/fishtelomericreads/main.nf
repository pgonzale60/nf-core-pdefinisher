process FISHTELOMERICREADS {
    tag "$meta.id"
    label 'process_low'

    container 'community.wave.seqera.io/library/samtools_bbmap_pigz:c8f695696758ff77'

    input:
    tuple val(meta), path(reads), val(telomere_repeat)
    val seed_copies
    val max_k
    val hdist
    val bam_exclude_flags

    output:
    tuple val(meta), path("${prefix}.broad_telomeric_reads.fasta.gz"), emit: reads
    tuple val(meta), path("${prefix}.input_headers.sam"), emit: headers
    tuple val(meta), path("${prefix}.fishing.log"), emit: log
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n 1 | sed 's/samtools //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('bbmap'), eval("bbduk.sh --version 2>&1 | head -n 1 | sed 's/BBMap version //'"), topic: versions, emit: versions_bbmap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def is_bam = reads.name.endsWith('.bam')
    def mem_gb = task.memory ? Math.max(1, (task.memory.toGiga() - 1) as int) : 2
    def seed   = telomere_repeat * (seed_copies as int)
    def k_len  = Math.min(seed.length(), max_k as int)
    def bbduk_cmd = """bbduk.sh -Xmx${mem_gb}g \\
        outm=${prefix}.broad_telomeric_reads.fasta.gz \\
        literal="${seed}" \\
        k=${k_len} \\
        rcomp=t \\
        hdist=${hdist} \\
        ordered=t \\
        threads=${task.cpus} \\
        overwrite=t \\
        ${args}"""
    """
    ${is_bam ? """
    samtools view -H "${reads}" > "${prefix}.input_headers.sam"
    samtools fasta -@ 1 -n -F ${bam_exclude_flags} "${reads}" 2> "${prefix}.fishing.log" \\
        | ${bbduk_cmd} in=stdin.fa 2>> "${prefix}.fishing.log"
    """ : """
    : > "${prefix}.input_headers.sam"
    ${bbduk_cmd} in="${reads}" 2> "${prefix}.fishing.log"
    """}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '' | gzip -c > "${prefix}.broad_telomeric_reads.fasta.gz"
    : > "${prefix}.input_headers.sam"
    : > "${prefix}.fishing.log"
    """
}
