process TELOMEREWALKER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/hifiasm_abpoa_minimap2_samtools_gfatools_python:b88b84d24c308348'

    input:
    tuple val(meta), path(blueprint_json), path(telomere_profile), path(bam), path(bai), path(ref_fasta), path(ref_gzi), path(telo_unitigs), path(unplaced_paf)
    val telomere_repeat

    output:
    tuple val(meta), path("*_telomere_extensions"), path("*_walker_manifest.tsv"), emit: manifest
    tuple val("${task.process}"), val('python'), val('3.10.4'), emit: versions_walker, topic: versions
    tuple val("${task.process}"), val('abpoa'), eval("abpoa --version | head -n 1"), topic: versions, emit: versions_abpoa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    // meta.chromosomes is a custom meta key, which nf-core lint discourages. It is
    // kept because the cohort run definition sets it per sample and this is the
    // walker's only consumer of it: dropping it silently widens the walk from the
    // chromosome set to every sequence in the assembly. Replace it with a proper
    // val input when the module interface is next revised.
    def chrom_arg = meta['chromosomes'] ? "--chromosomes \"${meta['chromosomes']}\"" : ""
    // Optional, and passed empty by this pipeline: a PROFILE_TELOMERES table of unitig alignment
    // coordinates is not produced anywhere here, and the reconnaissance's soft-clipped-read summary
    // is not one, despite having been handed over as one. Unitig intervals arrive on the Module 2
    // card instead. Standalone callers that do hold such a table can still supply it, and the engine
    // validates the header before using it.
    def profile_arg = telomere_profile ? "--telomere-profile \"${telomere_profile}\"" : ""
    def unitigs_arg = telo_unitigs ? "--unitigs \"${telo_unitigs}\"" : ""
    def unplaced_arg = unplaced_paf ? "--unplaced-paf \"${unplaced_paf}\"" : ""

    def blueprint_arg = blueprint_json ? "--blueprint \"${blueprint_json}\"" : ""
    """
    # Tool: pde_telomere_walker.py | Script-SHA: feb97e7085726a5c
    pde_telomere_walker.py \\
        --ref "${ref_fasta}" \\
        --bam "${bam}" \\
        ${blueprint_arg} \\
        ${profile_arg} \\
        ${unitigs_arg} \\
        ${unplaced_arg} \\
        ${chrom_arg} \\
        --telo-motif "${telomere_repeat}" \\
        --sample "${prefix}" \\
        --outdir "${prefix}_telomere_extensions" \\
        --manifest "${prefix}_walker_manifest.tsv" \\
        --threads ${task.cpus} \
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_telomere_extensions"
    touch "${prefix}_walker_manifest.tsv"
    """
}
