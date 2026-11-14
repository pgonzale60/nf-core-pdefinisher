#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/pdefinisher: RECON + MOTIF entry
    Runs PDE_RECONNAISSANCE_CORE and PDE_MOTIF_DISCOVERY only (mosdepth coverage, delfies S2G, GRS segments,
    CBM motifs). No assembly surgery, rDNA assembly, read fishing or softclip analysis.
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
    red_bed                    = RECON_MOTIF_PIPELINE.out.red_bed
    gaps_bed                   = RECON_MOTIF_PIPELINE.out.gaps_bed
    delfies_bed                = RECON_MOTIF_PIPELINE.out.delfies_bed
    delfies_foci_s2g           = RECON_MOTIF_PIPELINE.out.delfies_foci_s2g
    mosdepth_summary           = RECON_MOTIF_PIPELINE.out.mosdepth_summary
    mosdepth_regions           = RECON_MOTIF_PIPELINE.out.mosdepth_regions
    grs_quantized_bed          = RECON_MOTIF_PIPELINE.out.grs_quantized_bed
    class_regions_bed          = RECON_MOTIF_PIPELINE.out.class_regions_bed
    class_regions_csi          = RECON_MOTIF_PIPELINE.out.class_regions_csi
    break_windows_fasta        = RECON_MOTIF_PIPELINE.out.break_windows_fasta
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
    gaps_bed {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    red_bed {
        path { meta, file -> "${PdePublish.pathFor(params, meta, "pde_reconnaissance")}/${file.name.replaceFirst(/\.bed$/, '.red.bed')}" }
        mode params.publish_dir_mode
    }
    delfies_bed {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    delfies_foci_s2g {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    mosdepth_summary {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    mosdepth_regions {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    grs_quantized_bed {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    class_regions_bed {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    class_regions_csi {
        path { meta, _file -> PdePublish.pathFor(params, meta, "pde_reconnaissance") }
        mode params.publish_dir_mode
    }
    break_windows_fasta {
        path { meta, _file -> PdePublish.pathFor(params, meta, "motif_discovery") }
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
