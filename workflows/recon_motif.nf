/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PDE_RECONNAISSANCE  } from '../subworkflows/local/pde_reconnaissance/main'
include { PDE_MOTIF_DISCOVERY } from '../subworkflows/local/pde_motif_discovery/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RECON + MOTIF WORKFLOW (no assembly surgery, evaluation, or validation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RECON_MOTIF_PIPELINE {

    take:
    ch_samplesheet // channel: [ val(meta), path(assembly), path(bam), path(bai) ]

    main:

    ch_assembly = ch_samplesheet.map { meta, assembly, _bam, _bai -> tuple(meta, assembly) }
    ch_bam      = ch_samplesheet.map { meta, _assembly, bam, bai  -> tuple(meta, bam, bai) }

    def val_telomere_repeat = ch_samplesheet.map { meta, _assembly, _bam, _bai -> meta.telomere_repeat ?: (params.telomere_repeat ?: 'TTAGGC') }.first()
    def val_chromosomes     = ch_samplesheet.map { meta, _assembly, _bam, _bai -> meta.chromosomes ?: (params.chromosomes ?: '') }.first()

    PDE_RECONNAISSANCE (
        ch_assembly,
        ch_bam,
        val_telomere_repeat,
        val_chromosomes
    )

    PDE_MOTIF_DISCOVERY (
        ch_assembly,
        PDE_RECONNAISSANCE.out.pde_segments_bed,
        PDE_RECONNAISSANCE.out.delfies_foci_s2g,
        params.fimo_pvalue_threshold ?: 1e-4
    )

    emit:
    grs_summary_tsv            = PDE_RECONNAISSANCE.out.grs_summary_tsv
    pde_segments_bed           = PDE_RECONNAISSANCE.out.pde_segments_bed
    delfies_bed                = PDE_RECONNAISSANCE.out.delfies_bed
    delfies_foci_s2g           = PDE_RECONNAISSANCE.out.delfies_foci_s2g
    grs_quantized_bed          = PDE_RECONNAISSANCE.out.grs_quantized_bed
    terminal_ends_summary      = PDE_RECONNAISSANCE.out.terminal_ends_summary
    best_core_motif            = PDE_MOTIF_DISCOVERY.out.best_core_motif
    core_motif_metrics_json    = PDE_MOTIF_DISCOVERY.out.core_motif_metrics_json
    grid_summary_tsv           = PDE_MOTIF_DISCOVERY.out.grid_summary_tsv
    best_auxiliary_motif       = PDE_MOTIF_DISCOVERY.out.best_auxiliary_motif
    bipartite_spacing_tsv      = PDE_MOTIF_DISCOVERY.out.bipartite_spacing_tsv
    motif_discovery_report_md  = PDE_MOTIF_DISCOVERY.out.discovery_report_md
    fimo_genome_hits           = PDE_MOTIF_DISCOVERY.out.fimo_genome_hits
    fimo_annotated_hits        = PDE_MOTIF_DISCOVERY.out.fimo_annotated_hits
    motif_distribution_summary = PDE_MOTIF_DISCOVERY.out.motif_distribution_summary
    motif_distribution_report  = PDE_MOTIF_DISCOVERY.out.motif_distribution_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
