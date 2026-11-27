process NCRFTELOMERES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncrf:1.01.02--h7b50bb2_6':
        'quay.io/biocontainers/ncrf:1.01.02--h7b50bb2_6' }"

    input:
    tuple val(meta), path(reads), val(telomere_repeat)
    val ncrf_scoring
    val ncrf_min_length
    val ncrf_max_noise

    output:
    tuple val(meta), path("${prefix}.telomeric_repeats.ncrf_summary.tsv.gz"), emit: summary
    tuple val(meta), path("${prefix}.telomeric_repeats.ncrf.txt.gz"), emit: raw
    tuple val("${task.process}"), val('ncrf'), eval("NCRF --version 2>&1 | head -n 1 | awk '{print \$NF}'"), topic: versions, emit: versions_ncrf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    if ( set +o pipefail; gzip -cd "${reads}" | grep -q '^>' ); then
        gzip -cd "${reads}" \\
            | NCRF "${telomere_repeat}" \\
                --scoring="${ncrf_scoring}" \\
                --minlength="${ncrf_min_length}" \\
                --maxnoise="${ncrf_max_noise}" \\
                --stats=events \\
                --positionalevents \\
                ${args} \\
            | gzip -c > "${prefix}.telomeric_repeats.ncrf.txt.gz"

        if ! ( set +o pipefail; gzip -cd "${prefix}.telomeric_repeats.ncrf.txt.gz" | tail -n 5 | grep -q '^# ncrf end-of-file\$' ); then
            echo "ERROR: Missing NCRF end-of-file marker" >&2
            exit 1
        fi

        gzip -cd "${prefix}.telomeric_repeats.ncrf.txt.gz" \\
            | ncrf_summary.py \\
            | gzip -c > "${prefix}.telomeric_repeats.ncrf_summary.tsv.gz"
    else
        printf '# ncrf end-of-file\\n' | gzip -c > "${prefix}.telomeric_repeats.ncrf.txt.gz"
        printf 'line\\tmotif\\tseq\\tstart\\tend\\tstrand\\tseqLen\\tquerybp\\tmRatio\\tm\\tmm\\ti\\td\\n' \\
            | gzip -c > "${prefix}.telomeric_repeats.ncrf_summary.tsv.gz"
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '# ncrf end-of-file\\n' | gzip -c > "${prefix}.telomeric_repeats.ncrf.txt.gz"
    printf 'line\\tmotif\\tseq\\tstart\\tend\\tstrand\\tseqLen\\tquerybp\\tmRatio\\tm\\tmm\\ti\\td\\n' \\
        | gzip -c > "${prefix}.telomeric_repeats.ncrf_summary.tsv.gz"
    """
}
