#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/pdefinisher
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/pdefinisher
    Website: https://nf-co.re/pdefinisher
    Slack  : https://nfcore.slack.com/channels/pdefinisher
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PDEFINISHER             } from './workflows/pdefinisher'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_pdefinisher_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_pdefinisher_pipeline'
include { RECON_MOTIF_PIPELINE    } from './workflows/recon_motif'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PDEFINISHER {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    PDEFINISHER (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )

    emit:
    curated_fasta              = PDEFINISHER.out.curated_fasta
    curated_agp                = PDEFINISHER.out.curated_agp
    components_fasta           = PDEFINISHER.out.components_fasta
    segregated_fastas          = PDEFINISHER.out.segregated_fastas
    curation_audit             = PDEFINISHER.out.curation_audit
    curated_pde_segments       = PDEFINISHER.out.curated_pde_segments
    curated_repeats            = PDEFINISHER.out.curated_repeats
    blueprint_json             = PDEFINISHER.out.blueprint_json
    blueprint_md               = PDEFINISHER.out.blueprint_md
    scorecard_json             = PDEFINISHER.out.scorecard_json
    report_md                  = PDEFINISHER.out.report_md
    validation_json            = PDEFINISHER.out.validation_json
    validation_md              = PDEFINISHER.out.validation_md
    curated_bam                = PDEFINISHER.out.curated_bam
    curated_bai                = PDEFINISHER.out.curated_bai
    best_core_motif            = PDEFINISHER.out.best_core_motif
    core_motif_metrics_json    = PDEFINISHER.out.core_motif_metrics_json
    grid_summary_tsv           = PDEFINISHER.out.grid_summary_tsv
    best_auxiliary_motif       = PDEFINISHER.out.best_auxiliary_motif
    bipartite_spacing_tsv      = PDEFINISHER.out.bipartite_spacing_tsv
    motif_discovery_report_md  = PDEFINISHER.out.motif_discovery_report_md
    fimo_genome_hits           = PDEFINISHER.out.fimo_genome_hits
    fimo_annotated_hits        = PDEFINISHER.out.fimo_annotated_hits
    motif_distribution_summary = PDEFINISHER.out.motif_distribution_summary
    motif_distribution_report  = PDEFINISHER.out.motif_distribution_report
    summary_md                 = PDEFINISHER.out.summary_md
    summary_tsv                = PDEFINISHER.out.summary_tsv
    summary_json               = PDEFINISHER.out.summary_json
    multiqc_report             = PDEFINISHER.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions                   = PDEFINISHER.out.versions
}

//
// WORKFLOW: Reconnaissance and motif discovery only (run with -entry RECON_MOTIF)
//
workflow RECON_MOTIF {

    main:
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    RECON_MOTIF_PIPELINE (
        PIPELINE_INITIALISATION.out.samplesheet
    )
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PDEFINISHER (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_PDEFINISHER.out.multiqc_report
    )

    publish:
    curated_fasta              = NFCORE_PDEFINISHER.out.curated_fasta
    curated_agp                = NFCORE_PDEFINISHER.out.curated_agp
    components_fasta           = NFCORE_PDEFINISHER.out.components_fasta
    segregated_fastas          = NFCORE_PDEFINISHER.out.segregated_fastas
    curation_audit             = NFCORE_PDEFINISHER.out.curation_audit
    curated_pde_segments       = NFCORE_PDEFINISHER.out.curated_pde_segments
    curated_repeats            = NFCORE_PDEFINISHER.out.curated_repeats
    blueprint_json             = NFCORE_PDEFINISHER.out.blueprint_json
    blueprint_md               = NFCORE_PDEFINISHER.out.blueprint_md
    scorecard_json             = NFCORE_PDEFINISHER.out.scorecard_json
    report_md                  = NFCORE_PDEFINISHER.out.report_md
    validation_json            = NFCORE_PDEFINISHER.out.validation_json
    validation_md              = NFCORE_PDEFINISHER.out.validation_md
    curated_bam                = NFCORE_PDEFINISHER.out.curated_bam
    curated_bai                = NFCORE_PDEFINISHER.out.curated_bai
    best_core_motif            = NFCORE_PDEFINISHER.out.best_core_motif
    core_motif_metrics_json    = NFCORE_PDEFINISHER.out.core_motif_metrics_json
    grid_summary_tsv           = NFCORE_PDEFINISHER.out.grid_summary_tsv
    best_auxiliary_motif       = NFCORE_PDEFINISHER.out.best_auxiliary_motif
    bipartite_spacing_tsv      = NFCORE_PDEFINISHER.out.bipartite_spacing_tsv
    motif_discovery_report_md  = NFCORE_PDEFINISHER.out.motif_discovery_report_md
    fimo_genome_hits           = NFCORE_PDEFINISHER.out.fimo_genome_hits
    fimo_annotated_hits        = NFCORE_PDEFINISHER.out.fimo_annotated_hits
    motif_distribution_summary = NFCORE_PDEFINISHER.out.motif_distribution_summary
    motif_distribution_report  = NFCORE_PDEFINISHER.out.motif_distribution_report
    summary_md                 = NFCORE_PDEFINISHER.out.summary_md
    summary_tsv                = NFCORE_PDEFINISHER.out.summary_tsv
    summary_json               = NFCORE_PDEFINISHER.out.summary_json
    multiqc_report             = NFCORE_PDEFINISHER.out.multiqc_report
    versions                   = NFCORE_PDEFINISHER.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUT PUBLISHING (NEXTFLOW 26+)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


output {
    curated_fasta {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curated_assembly") }
        mode params.publish_dir_mode
    }
    curated_agp {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curated_assembly") }
        mode params.publish_dir_mode
    }
    components_fasta {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curated_assembly") }
        mode params.publish_dir_mode
    }
    segregated_fastas {
        path { Object[] args -> PdePublish.pathFor(params, args[0], "curated_assembly") }
        mode params.publish_dir_mode
    }
    curation_audit {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curated_assembly") }
        mode params.publish_dir_mode
    }
    curated_pde_segments {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curated_assembly") }
        mode params.publish_dir_mode
    }
    curated_repeats {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curated_assembly") }
        mode params.publish_dir_mode
    }
    blueprint_json {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curation_synthesis") }
        mode params.publish_dir_mode
    }
    blueprint_md {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curation_synthesis") }
        mode params.publish_dir_mode
    }
    scorecard_json {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curation_evaluation") }
        mode params.publish_dir_mode
    }
    report_md {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curation_evaluation") }
        mode params.publish_dir_mode
    }
    validation_json {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curation_validation") }
        mode params.publish_dir_mode
    }
    validation_md {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curation_validation") }
        mode params.publish_dir_mode
    }
    curated_bam {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curation_validation") }
        mode params.publish_dir_mode
    }
    curated_bai {
        path { meta, _file -> PdePublish.pathFor(params, meta, "curation_validation") }
        mode params.publish_dir_mode
    }
    best_core_motif {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    core_motif_metrics_json {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    grid_summary_tsv {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    best_auxiliary_motif {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    bipartite_spacing_tsv {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    motif_discovery_report_md {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    fimo_genome_hits {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    fimo_annotated_hits {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    motif_distribution_summary {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    motif_distribution_report {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
        mode params.publish_dir_mode
    }
    // Cohort-level outputs:
    summary_md {
        path { Object[] _args -> "curation_summary" }
        mode params.publish_dir_mode
    }
    summary_tsv {
        path { Object[] _args -> "curation_summary" }
        mode params.publish_dir_mode
    }
    summary_json {
        path { Object[] _args -> "curation_summary" }
        mode params.publish_dir_mode
    }
    // Pipeline info outputs:
    multiqc_report {
        path { "pipeline_info" }
        mode params.publish_dir_mode
    }
    versions {
        path { "pipeline_info" }
        mode params.publish_dir_mode
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
