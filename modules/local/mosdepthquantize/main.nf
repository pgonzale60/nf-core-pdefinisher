process MOSDEPTHQUANTIZE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a3/a3dc5ea2ce788c24079d24d1721ed28086874152c43b5e7dde3f638dcf64336a/data' :
        'community.wave.seqera.io/library/htslib_mosdepth_gzip:4108dd38be84e40a'}"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(red_bed)
    tuple val(meta3), path(core_coverage_tsv)

    output:
    tuple val(meta), path("*.quantized.red_masked.bed"), emit: quantized_bed
    tuple val("${task.process}"), val('mosdepth'), eval("mosdepth --version | sed 's/mosdepth //g'"), topic: versions, emit: versions_mosdepth
    tuple val("${task.process}"), val('gzip'), eval("gzip -V 2>&1 | sed 's/gzip \\([0-9.]*\\).*/\\1/;q'"), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # T_GRS from the classifier (grs_threshold column, row 1 after header)
    T_GRS=\$(awk -F'\\t' 'NR==2{printf "%.6f", \$3}' ${core_coverage_tsv})
    T_FLOOR=\$(awk -v t="\$T_GRS" 'BEGIN{printf "%d", int(t)}')
    T_CEIL=\$(awk -v t="\$T_GRS" 'BEGIN{printf "%d", (t == int(t)) ? int(t) : int(t) + 1}')
    if [ "\$T_CEIL" -le 0 ]; then
        echo "Non-positive T_GRS (\$T_GRS) for ${prefix}: cannot quantize" >&2
        exit 1
    fi

    # Three classes: 0:T_floor (low), T_floor:T_ceil (bracket), T_ceil:inf (high).
    # Integer T_GRS: bracket is empty, emit 0:T and T:inf only.
    if [ "\$T_FLOOR" -eq "\$T_CEIL" ]; then
        QUANT="0:\${T_FLOOR}:"
    else
        QUANT="0:\${T_FLOOR}:\${T_CEIL}:"
    fi

    # MAPQ 0 kept on purpose: repeat reads must not create artificial drops
    mosdepth \\
        --threads $task.cpus \\
        --no-per-base \\
        --mapq 0 \\
        --quantize "\${QUANT}" \\
        ${prefix} \\
        ${bam}

    # Flag (not remove) quantized intervals overlapping RED repeats; column 5 = 1 means repeat-overlapping
    awk 'BEGIN { FS = OFS = "\\t" }
        FNR == NR { n[\$1]++; s[\$1, n[\$1]] = \$2; e[\$1, n[\$1]] = \$3; next }
        {
            hit = 0
            for (i = 1; i <= n[\$1]; i++) if (\$2 < e[\$1, i] && \$3 > s[\$1, i]) { hit = 1; break }
            print \$1, \$2, \$3, \$4, hit
        }' ${red_bed} <(zcat ${prefix}.quantized.bed.gz) > ${prefix}.quantized.red_masked.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.quantized.red_masked.bed
    """
}
