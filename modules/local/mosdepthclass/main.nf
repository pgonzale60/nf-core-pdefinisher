process MOSDEPTHCLASS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a3/a3dc5ea2ce788c24079d24d1721ed28086874152c43b5e7dde3f638dcf64336a/data' :
        'community.wave.seqera.io/library/htslib_mosdepth_gzip:4108dd38be84e40a'}"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(segments_bed)

    output:
    tuple val(meta), path("*.class.regions.bed.gz")     , emit: class_regions_bed
    tuple val(meta), path("*.class.regions.bed.gz.csi") , emit: class_regions_csi
    tuple val("${task.process}"), val('mosdepth'), eval("mosdepth --version | sed 's/mosdepth //g'"), topic: versions, emit: versions_mosdepth
    tuple val("${task.process}"), val('gzip'), eval("gzip -V 2>&1 | sed 's/gzip \\([0-9.]*\\).*/\\1/;q'"), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '--mapq 20'
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Regions are the final classification segments; column 4 (putative_core, putative_GRS, uncertain) names each region
    mosdepth \\
        --threads $task.cpus \\
        --no-per-base \\
        --by ${segments_bed} \\
        $args \\
        ${prefix}.class \\
        ${bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.class.regions.bed.gz
    touch ${prefix}.class.regions.bed.gz.csi
    """
}
