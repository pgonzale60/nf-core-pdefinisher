process EXTRACTTELOMERICREADS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_bbmap_pigz:c8f695696758ff77'

    input:
    tuple val(meta), path(bam), path(bai)
    val telomere_repeat

    output:
    tuple val(meta), path("*.telomeric_reads.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n 1 | sed 's/samtools //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('bbmap'), eval("bbduk.sh --version 2>&1 | head -n 1 | sed 's/BBMap version //'"), topic: versions, emit: versions_bbmap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fwd_repeat=\$(printf "%0.s${telomere_repeat}" {1..3})
    k_len=\${#fwd_repeat}
    if [ "\$k_len" -gt 18 ]; then
        k_len=18
    fi

    samtools fastq \\
        -@ ${task.cpus} \\
        ${args} \\
        ${bam} \\
        | bbduk.sh \\
            in=stdin.fq \\
            outm=${prefix}.telomeric_reads.fastq.gz \\
            literal="\$fwd_repeat" \\
            k=\$k_len \\
            rcomp=t \\
            hdist=0 \\
            threads=${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip -c > ${prefix}.telomeric_reads.fastq.gz
    """
}
