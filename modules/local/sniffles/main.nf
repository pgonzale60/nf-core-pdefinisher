process SNIFFLES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:2.8.0--pyhdfd78af_0' :
        'quay.io/biocontainers/sniffles:2.8.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(tandem_file)
    val(vcf_output)
    val(snf_output)

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf, optional: true
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.snf")       , emit: snf, optional: true
    tuple val("${task.process}"), val('sniffles'), val('2.8.0'), emit: versions_sniffles, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def input_ref      = fasta =~ /\.gz$/ ? fasta.name.take(fasta.name.lastIndexOf('.')) : fasta
    def gunzip_ref     = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ${input_ref}" : ""
    def reference      = fasta ? "--reference ${input_ref}" : ""
    def tandem_repeats = tandem_file ? "--tandem-repeats ${tandem_file}" : ''
    def vcf            = vcf_output ? "--vcf ${prefix}.vcf.gz" : ''
    def snf            = snf_output ? "--snf ${prefix}.snf" : ''

    """
    $gunzip_ref

    sniffles \\
        --input $input \\
        $reference \\
        -t $task.cpus \\
        $tandem_repeats \\
        $vcf \\
        $snf \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf    = vcf_output ? "echo \"\" | gzip > ${prefix}.vcf.gz; touch ${prefix}.vcf.gz.tbi" : ''
    def snf    = snf_output ? "touch ${prefix}.snf" : ''

    """
    ${vcf}
    ${snf}
    """
}
