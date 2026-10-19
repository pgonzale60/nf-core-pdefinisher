process FIMO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meme:5.5.9--pl5321he99cc7f_1' :
        'quay.io/biocontainers/meme:5.5.9--pl5321he99cc7f_1' }"

    input:
    tuple val(meta), path(motif_meme), path(assembly)
    val pvalue_thresh

    output:
    tuple val(meta), path("*.fimo_genome_hits.tsv"), emit: tsv
    tuple val(meta), path("*_fimo_out")            , emit: outdir
    tuple val("${task.process}"), val('fimo'), eval("fimo -version 2>&1 | head -n 1 | sed 's/.*version //'"), topic: versions, emit: versions_fimo

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def thresh          = pvalue_thresh != null ? pvalue_thresh : 1e-4
    def input_assembly  = assembly =~ /\.gz$/ ? assembly.name.take(assembly.name.lastIndexOf('.')) : assembly
    def gunzip_assembly = assembly =~ /\.gz$/ ? "gunzip -c ${assembly} > ${input_assembly}" : ""

    """
    $gunzip_assembly

    fimo \\
        --thresh ${thresh} \\
        --oc "${prefix}_fimo_out" \\
        ${args} \\
        "${motif_meme}" \\
        "${input_assembly}"

    cp "${prefix}_fimo_out/fimo.tsv" "${prefix}.fimo_genome_hits.tsv"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_fimo_out"
    touch "${prefix}_fimo_out/fimo.tsv"
    touch "${prefix}.fimo_genome_hits.tsv"
    """
}
