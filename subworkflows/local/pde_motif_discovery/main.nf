//
// SUBWORKFLOW: PDE_MOTIF_DISCOVERY
// Discovers chromosomal breakage motifs and auxiliary bipartite signals from authentic PDE break sites,
// and profiles genome-wide motif occurrences across assembly compartments.
//

include { EXTRACTPDEBORDERS } from '../../../modules/local/extractpdeborders/main'
include { MEMECOREGRID      } from '../../../modules/local/memecoregrid/main'
include { BIPARTITEPROFILER } from '../../../modules/local/bipartiteprofiler/main'
include { FIMO              } from '../../../modules/local/fimo/main'
include { PDEMOTIFPROFILER  } from '../../../modules/local/pdemotifprofiler/main'

workflow PDE_MOTIF_DISCOVERY {

    take:
    ch_assembly         // channel: [ val(meta), path(assembly) ]
    ch_pde_segments     // channel: [ val(meta), path(pde_segments_bed) ]
    ch_breakpoints_bed  // channel: [ val(meta), path(breakpoints_bed) ]
    ch_grs_summary      // channel: [ val(meta), path(grs_summary) ]
    ch_ends_summary     // channel: [ val(meta), path(terminal_ends_summary) ]
    pvalue_thresh       // val: e.g. 1e-4

    main:
    ch_versions = Channel.empty()

    // 1. Extract authentic GRS/Core border transitions directly from pde_segments.bed
    ch_borders_input = ch_pde_segments
        .join(ch_breakpoints_bed, by: 0, remainder: true)
        .map { meta, segs, bps -> [ meta, segs, bps ?: [] ] }
        .join(ch_assembly, by: 0)
        .join(ch_grs_summary, by: 0, remainder: true)
        .map { meta, segs, bps, asm, grs -> [ meta, segs, bps, asm, grs ?: [] ] }
        .join(ch_ends_summary, by: 0, remainder: true)
        .map { meta, segs, bps, asm, grs, ends -> tuple(meta, segs, bps, asm, grs ?: [], ends ?: []) }

    EXTRACTPDEBORDERS (
        ch_borders_input,
        150,
        350
    )
    ch_versions = ch_versions.mix(EXTRACTPDEBORDERS.out.versions_python)

    // 2. Tier-1 Core Motif Discovery via Markov-backed parameter grid
    MEMECOREGRID (
        EXTRACTPDEBORDERS.out.break_windows_fasta
    )
    ch_versions = ch_versions.mix(MEMECOREGRID.out.versions_meme)

    // 3. Tier-2 Auxiliary Motif Discovery & Bipartite Spacing Profiler
    ch_bipartite_input = EXTRACTPDEBORDERS.out.break_windows_fasta
        .join(MEMECOREGRID.out.best_core_motif, by: 0)
        .join(EXTRACTPDEBORDERS.out.census_json, by: 0)
        .join(MEMECOREGRID.out.motif_metrics_json, by: 0)
        .join(MEMECOREGRID.out.best_background, by: 0)

    BIPARTITEPROFILER (
        ch_bipartite_input
    )
    ch_versions = ch_versions.mix(BIPARTITEPROFILER.out.versions_meme)

    // 4. Whole-Genome FIMO Scanning across raw assembly
    ch_fimo_input = MEMECOREGRID.out.best_core_motif
        .join(ch_assembly, by: 0)

    FIMO (
        ch_fimo_input,
        pvalue_thresh
    )
    ch_versions = ch_versions.mix(FIMO.out.versions_fimo)

    // 5. Profiling of motif occurrences across authentic break sites, GRS, and core
    ch_profiler_input = FIMO.out.tsv
        .join(EXTRACTPDEBORDERS.out.border_breakpoints, by: 0)
        .join(ch_pde_segments, by: 0)
        .map { meta, hits, bps, segs -> tuple(meta, hits, bps, segs) }

    PDEMOTIFPROFILER (
        ch_profiler_input,
        500
    )
    ch_versions = ch_versions.mix(PDEMOTIFPROFILER.out.versions_python)

    emit:
    border_breakpoints         = EXTRACTPDEBORDERS.out.border_breakpoints
    break_windows_fasta        = EXTRACTPDEBORDERS.out.break_windows_fasta
    census_json                = EXTRACTPDEBORDERS.out.census_json
    best_core_motif            = MEMECOREGRID.out.best_core_motif
    grid_summary_tsv           = MEMECOREGRID.out.grid_summary_tsv
    core_motif_metrics_json    = MEMECOREGRID.out.motif_metrics_json
    best_auxiliary_motif       = BIPARTITEPROFILER.out.best_auxiliary_motif
    bipartite_spacing_tsv      = BIPARTITEPROFILER.out.bipartite_spacing_tsv
    discovery_report_md        = BIPARTITEPROFILER.out.discovery_report_md
    fimo_genome_hits           = FIMO.out.tsv
    fimo_annotated_hits        = PDEMOTIFPROFILER.out.annotated_hits
    motif_distribution_summary = PDEMOTIFPROFILER.out.summary_tsv
    motif_distribution_report  = PDEMOTIFPROFILER.out.report_md
    versions                   = ch_versions
}
