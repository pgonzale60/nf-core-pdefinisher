process GFASTATSGAPS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/gfastats:1.3.6--hdcf5f25_3'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.gaps.bed"), emit: gaps_bed
    tuple val("${task.process}"), val('gfastats'), eval("gfastats --version 2>&1 | sed '1!d;s/.* //'"), topic: versions, emit: versions_gfastats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # gfastats cannot open gzipped paths, so the assembly is streamed via stdin.
    # Output: headerless BED (chrom start end N_gap length strand), sorted. An assembly
    # without N-runs yields an empty file, which is a valid result, not an error.
    zcat -f "${assembly}" \\
        | gfastats /dev/stdin -j ${task.cpus} -b g ${args} \\
        | sort -k1,1 -k2,2n \\
        | awk 'BEGIN { OFS="\\t" } { print \$1, \$2, \$3, "N_gap", \$3 - \$2, "." }' \\
        > "${prefix}.gaps.bed"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.gaps.bed"
    """
}
