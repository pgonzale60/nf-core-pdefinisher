process TELOCLIP_EXTRACT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/teloclip_samtools:2e39bcbf9e597d36'

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)
    val telomere_repeat

    output:
    tuple val(meta), path("*_teloclip_overhangs"), emit: overhangs
    tuple val(meta), path("*.teloclip.bam"), path("*.teloclip.bam.bai"), emit: bam, optional: true
    tuple val("${task.process}"), val('teloclip'), eval("teloclip --version 2>&1 | head -n 1 | sed 's/teloclip //'"), topic: versions, emit: versions_teloclip
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n 1 | sed 's/samtools //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    fai="${fasta}.fai"
    if [ ! -f "\$fai" ]; then
        samtools faidx "${fasta}"
    fi

    teloclip filter \\
        --ref-idx "\$fai" \\
        --min-clip 100 \\
        --max-break 50 \\
        --motifs "${telomere_repeat}" \\
        "${bam}" \\
        | samtools view -bS - \\
        | samtools sort -@ ${task.cpus} -o "${prefix}.teloclip.bam"

    samtools index -@ ${task.cpus} "${prefix}.teloclip.bam"

    mkdir -p "${prefix}_teloclip_overhangs"
    teloclip extract \\
        --ref-idx "\$fai" \\
        --min-clip 100 \\
        --out-dir "${prefix}_teloclip_overhangs" \\
        "${prefix}.teloclip.bam"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_teloclip_overhangs"
    touch "${prefix}.teloclip.bam"
    touch "${prefix}.teloclip.bam.bai"
    """
}
