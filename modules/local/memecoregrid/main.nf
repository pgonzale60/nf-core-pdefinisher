process MEMECOREGRID {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meme:5.5.9--pl5321he99cc7f_1' :
        'quay.io/biocontainers/meme:5.5.9--pl5321he99cc7f_1' }"

    input:
    tuple val(meta), path(break_windows_fasta)

    output:
    tuple val(meta), path("*.best_core_motif.meme")   , emit: best_core_motif
    tuple val(meta), path("*.core_grid_summary.tsv")  , emit: grid_summary_tsv
    tuple val(meta), path("*.core_motif_metrics.json"), emit: motif_metrics_json
    tuple val(meta), path("*.best_background.bfile")  , emit: best_background
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('meme')  , eval("meme -version"), topic: versions, emit: versions_meme

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Tool: cbm_meme_core_grid.py | Script-SHA: 69c6c7b1fc19b4f1
    cbm_meme_core_grid.py \\
        --break-windows "${break_windows_fasta}" \\
        --output-dir "meme_grid_work" \\
        --output-meme "${prefix}.best_core_motif.meme" \\
        --output-tsv "${prefix}.core_grid_summary.tsv" \\
        --output-metrics "${prefix}.core_motif_metrics.json" \\
        --output-background "${prefix}.best_background.bfile" \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.best_core_motif.meme"
    touch "${prefix}.best_background.bfile"
    echo "phase\tcondition_id\tradius\tmarkov_order\tmod\trevcomp\tw_label\twidth\tnsites\tevalue\tconsensus\ttotal_ic\tmean_ic\tconsensus_entropy\toccupancy\tmedian_offset_to_pos0\tpositional_iqr\tcomposite_score" > "${prefix}.core_grid_summary.tsv"
    echo '{"consensus": "A", "composite_score": 0.0}' > "${prefix}.core_motif_metrics.json"
    """
}
