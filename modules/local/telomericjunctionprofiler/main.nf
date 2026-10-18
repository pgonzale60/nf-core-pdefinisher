process TELOMERICJUNCTIONPROFILER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_minimap2_vsearch_python:1aa370c14225a820'

    input:
    tuple val(meta), path(reads), path(ref_fasta), path(delfies_bed), path(policy_json)
    val telomere_repeat

    output:
    tuple val(meta), path("*.telomeric_junctions.json"), emit: json
    tuple val(meta), path("*.telomeric_junctions.tsv") , emit: tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('vsearch'), eval("vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //' | sed 's/,.*//'"), topic: versions, emit: versions_vsearch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def delfies_arg = delfies_bed ? "--delfies-bed \"${delfies_bed}\"" : ""

    """
    # Tool: pde_telomeric_junction_profiler.py | Script-SHA: e1485524244791fd
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
