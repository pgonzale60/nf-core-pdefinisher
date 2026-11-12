process CURATIONSYNTHESIS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'quay.io/biocontainers/python:3.10.4' }"

    input:
    tuple(
        val(meta),
        path(ref_fasta),
        path(grs_summary),
        path(mosdepth_regions),
        path(terminal_ends_summary),
        path(terminal_softclips),
        path(delfies_bed),
        path(rdna_regions_bed),
        path(barrnap_gff),
        path(flye_info),
        path(unplaced_paf, stageAs: 'unplaced.paf'),
        path(red_bed),
        path(telo_paf, stageAs: 'telo.paf'),
        path(gap_evidence_json),
        path(telo_junction_json)
    )
    path policy_file, stageAs: "curation_policy.json"
    val telo_motif
    val chromosomes

    output:
    tuple val(meta), path("*.curation_manifest.json"), emit: blueprint_json
    tuple val(meta), path("*.curation_blueprint.md"), emit: blueprint_md
    tuple val("${task.process}"), val('python'), val('3.10.4'), emit: versions_synthesis, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def gap_arg    = gap_evidence_json ? "--gap-evidence-json \"${gap_evidence_json}\"" : ""
    def junc_arg   = telo_junction_json ? "--telo-junctions \"${telo_junction_json}\"" : ""
    def eff_chroms = meta.chromosomes ?: chromosomes
    def chr_arg    = eff_chroms ? "--chromosomes \"${eff_chroms}\"" : ""
    def pol_arg    = policy_file ? "--policy \"${policy_file}\"" : ""
    def eff_telo   = meta.telomere_repeat ?: telo_motif
    """
    # Tool: pde_curation_synthesis.py | Script-SHA: 49d0f0a64a54bb98
    pde_curation_synthesis.py --version
    pde_curation_synthesis.py \\
        --sample "${prefix}" \\
        --ref "${ref_fasta}" \\
        --grs-summary "${grs_summary}" \\
        --mosdepth-regions "${mosdepth_regions}" \\
        --terminal-summary "${terminal_ends_summary}" \\
        --terminal-softclips "${terminal_softclips}" \\
        --delfies-bed "${delfies_bed}" \\
        --rdna-bed "${rdna_regions_bed}" \\
        --barrnap-gff "${barrnap_gff}" \\
        --flye-info "${flye_info}" \\
        --unplaced-paf "unplaced.paf" \\
        --red-bed "${red_bed}" \\
        --telo-paf "telo.paf" \\
        ${gap_arg} \\
        ${junc_arg} \\
        ${chr_arg} \\
        ${pol_arg} \\
        --telo-motif "${eff_telo}" \\
        --output-json "${prefix}.curation_manifest.json" \\
        --output-md "${prefix}.curation_blueprint.md" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.curation_manifest.json"
    touch "${prefix}.curation_blueprint.md"
    """
}
