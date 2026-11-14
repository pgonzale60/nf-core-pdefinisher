process DELFIESBREAKWINDOWS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_minimap2_python:c9b2a6b2ba2d6d79'

    input:
    tuple val(meta), path(s2g_tsv), path(pde_segments_bed), path(assembly)

    output:
    tuple val(meta), path("*.breakwindows.bed"), emit: bed
    tuple val(meta), path("*.breakwindows.fa") , emit: fasta
    tuple val(meta), path("*.census.json")     , emit: census_json
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Tool: pde_delfies_break_windows.py | Script-SHA: e3151cbeef57ff38
    pde_delfies_break_windows.py \\
        --s2g-tsv "${s2g_tsv}" \\
        --pde-segments-bed "${pde_segments_bed}" \\
        --assembly "${assembly}" \\
        --prefix "${prefix}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.breakwindows.bed"
    touch "${prefix}.breakwindows.fa"
    echo '{"n_breaks": 0}' > "${prefix}.census.json"
    """
}
