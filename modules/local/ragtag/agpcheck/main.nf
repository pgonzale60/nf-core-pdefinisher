process RAGTAG_AGPCHECK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/ragtag:2.1.0--6125126fe9d553ab'

    input:
    tuple val(meta), path(agp)

    output:
    tuple val(meta), path(agp), emit: agp
    tuple val("${task.process}"), val('ragtag'), eval("ragtag.py --version | sed 's/RagTag //'"), topic: versions, emit: versions_ragtag

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ragtag.py agpcheck ${args} "${agp}"
    """

    stub:
    """
    echo "ragtag.py agpcheck ok"
    """
}
