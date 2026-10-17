process TELOMERICJUNCTIONPROFILER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_minimap2_python:c9b2a6b2ba2d6d79'

    input:
    tuple val(meta), path(reads), path(ref_fasta), path(delfies_bed), path(policy_json)
    val telomere_repeat

    output:
    tuple val(meta), path("*.telomeric_junctions.json"), emit: json
    tuple val(meta), path("*.telomeric_junctions.tsv") , emit: tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def delfies_arg = delfies_bed ? "--delfies-bed \"${delfies_bed}\"" : ""

    """
    # Tool: pde_telomeric_junction_profiler.py | Script-SHA: c2295a79d9399aa6
    pde_telomeric_junction_profiler.py \\
        --reads "${reads}" \\
        --ref "${ref_fasta}" \\
        ${delfies_arg} \\
        --policy "${policy_json}" \\
        --motif "${telomere_repeat}" \\
        --sample "${prefix}" \\
        --threads ${task.cpus} \\
        --output-json "${prefix}.telomeric_junctions.json" \\
        --output-tsv "${prefix}.telomeric_junctions.tsv" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"sample": "${prefix}", "clusters": []}' > ${prefix}.telomeric_junctions.json
    touch ${prefix}.telomeric_junctions.tsv
    """
}
