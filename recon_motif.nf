#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/pdefinisher: RECON + MOTIF entry
    Runs PDE_RECONNAISSANCE and PDE_MOTIF_DISCOVERY only. No assembly surgery, evaluation or validation.
    Usage: nextflow run recon_motif.nf --input samplesheet.csv --outdir <dir>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RECON_MOTIF_PIPELINE    } from './workflows/recon_motif'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_pdefinisher_pipeline'

workflow {

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

    publish:
    grs_summary_tsv            = RECON_MOTIF_PIPELINE.out.grs_summary_tsv
    pde_segments_bed           = RECON_MOTIF_PIPELINE.out.pde_segments_bed
    delfies_bed                = RECON_MOTIF_PIPELINE.out.delfies_bed
    terminal_ends_summary      = RECON_MOTIF_PIPELINE.out.terminal_ends_summary
    best_core_motif            = RECON_MOTIF_PIPELINE.out.best_core_motif
    core_motif_metrics_json    = RECON_MOTIF_PIPELINE.out.core_motif_metrics_json
    grid_summary_tsv           = RECON_MOTIF_PIPELINE.out.grid_summary_tsv
    best_auxiliary_motif       = RECON_MOTIF_PIPELINE.out.best_auxiliary_motif
    bipartite_spacing_tsv      = RECON_MOTIF_PIPELINE.out.bipartite_spacing_tsv
    motif_discovery_report_md  = RECON_MOTIF_PIPELINE.out.motif_discovery_report_md
    fimo_genome_hits           = RECON_MOTIF_PIPELINE.out.fimo_genome_hits
    fimo_annotated_hits        = RECON_MOTIF_PIPELINE.out.fimo_annotated_hits
    motif_distribution_summary = RECON_MOTIF_PIPELINE.out.motif_distribution_summary
    motif_distribution_report  = RECON_MOTIF_PIPELINE.out.motif_distribution_report
}

output {
    grs_summary_tsv {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    pde_segments_bed {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    delfies_bed {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    terminal_ends_summary {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
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
}
