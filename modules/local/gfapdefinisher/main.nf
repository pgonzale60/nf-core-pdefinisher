process GFAPDEFINISHER {
    tag "$meta.id"
    label 'process_low'

    // The installed graph tool can change independently of the executable path.
    cache false

    input:
    tuple val(meta), path(assembly), path(gfa), path(end_archive, stageAs: 'end_archive'),
        path(raw_gfa), path(tiara), path(barrnap), path(ont_paf), path(ont_reads), path(hifi_reads), path(recipe),
        path(hifi_patch_approvals), path(bam_evidence), path(gap_patches, stageAs: 'hifi_gap_payloads')
    path handoff_script
    path package_script
    path edit_plan_script
    path hifi_patch_script
    val executable

    output:
    tuple val(meta), path("${meta.id}.curated.fasta*"), emit: curated_fasta
    tuple val(meta), path("${meta.id}.curated.agp"), emit: curated_agp
    tuple val(meta), path("${meta.id}.components.fasta"), emit: components_fasta
    tuple val(meta), path('graph_evidence'), emit: evidence

    script:
    def inputs = [raw_gfa: raw_gfa, tiara: tiara, barrnap: barrnap, ont_paf: ont_paf,
                  ont_reads: ont_reads, hifi_reads: hifi_reads, recipe: recipe,
                  hifi_patch_approvals: hifi_patch_approvals, gap_patches: gap_patches]
        .findAll { _key, value -> value }
        .collect { key, value -> "--${key.replace('_', '-')} \"${value}\"" }
        .join(' ')
    def bam_args = bam_evidence ? '--bam-evidence ' + bam_evidence.collect { "\"${it}\"" }.join(' ') : ''
    """
    PYTHONPATH="\$PWD" python3 "${handoff_script}" \\
        --curated "${assembly}" \\
        --gfa "${gfa}" \\
        --end-archive "${end_archive}" \\
        --chromosomes "${meta.chromosomes}" \\
        --telomere "${meta.telomere_repeat}" \\
        --gfapdefinisher "${executable}" \\
        --threads ${task.cpus} \\
        --outdir graph_evidence \\
        ${inputs} ${bam_args}
    if [[ -f graph_evidence/assembly/assembly.fasta ]]; then
        ln graph_evidence/assembly/assembly.fasta "${meta.id}.curated.fasta"
    else
        ln graph_evidence/assembly/assembly.fasta.gz "${meta.id}.curated.fasta.gz"
    fi
    ln graph_evidence/assembly/assembly.agp "${meta.id}.curated.agp"
    ln graph_evidence/assembly/components.fasta "${meta.id}.components.fasta"
    """

    stub:
    def fasta_suffix = recipe || !assembly.name.endsWith('.gz') ? '' : '.gz'
    """
    mkdir graph_evidence
    touch "${meta.id}.curated.fasta${fasta_suffix}" "${meta.id}.curated.agp" "${meta.id}.components.fasta"
    touch graph_evidence/inputs.json
    """
}
