process EXTRACTTELOMERICREADS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_seqkit_gawk:5f1679612f236815'

    input:
    tuple val(meta), path(bam), path(bai)
    val telomere_repeat

    output:
    tuple val(meta), path("*.telomeric_reads.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n 1 | sed 's/samtools //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit //'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fwd_repeat=\$(printf "%0.s${telomere_repeat}" {1..3})
    rev_repeat=\$(echo "\$fwd_repeat" | tr 'ACGTacgt' 'TGCAtgca' | rev)

    samtools fastq \\
        -@ ${task.cpus} \\
        ${args} \\
        ${bam} \\
        | seqkit grep -j ${task.cpus} -s -i -p "\$fwd_repeat" -p "\$rev_repeat" \\
        | gzip -c > ${prefix}.telomeric_reads.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip -c > ${prefix}.telomeric_reads.fastq.gz
    """
}
