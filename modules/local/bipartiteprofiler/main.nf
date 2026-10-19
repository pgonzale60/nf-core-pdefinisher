process BIPARTITEPROFILER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meme:5.5.9--pl5321he99cc7f_1' :
        'quay.io/biocontainers/meme:5.5.9--pl5321he99cc7f_1' }"

    input:
    tuple val(meta), path(break_windows_fasta), path(best_core_motif), path(census_json), path(core_metrics_json), path(best_background)

    output:
    tuple val(meta), path("*.best_auxiliary_motif.meme")     , emit: best_auxiliary_motif
    tuple val(meta), path("*.bipartite_spacing_summary.tsv") , emit: bipartite_spacing_tsv
    tuple val(meta), path("*.cbm_discovery_report.md")       , emit: discovery_report_md
    tuple val(meta), path("*.core_masked.fasta")             , emit: masked_fasta
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('meme')  , eval("meme -version"), topic: versions, emit: versions_meme

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Tool: cbm_bipartite_profiler.py | Script-SHA: c400fc5f00391126
    cbm_bipartite_profiler.py \\
        --break-windows "${break_windows_fasta}" \\
        --core-meme "${best_core_motif}" \\
        --census-json "${census_json}" \\
        --core-metrics-json "${core_metrics_json}" \\
        --background-bfile "${best_background}" \\
        --sample-id "${prefix}" \\
        --output-aux-meme "${prefix}.best_auxiliary_motif.meme" \\
        --output-spacing-tsv "${prefix}.bipartite_spacing_summary.tsv" \\
        --output-report "${prefix}.cbm_discovery_report.md" \\
        --output-masked-fasta "${prefix}.core_masked.fasta" \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.best_auxiliary_motif.meme"
    touch "${prefix}.bipartite_spacing_summary.tsv"
    echo "# CBM Discovery Report: ${prefix}" > "${prefix}.cbm_discovery_report.md"
    touch "${prefix}.core_masked.fasta"
    """
}
