//
// SUBWORKFLOW: CURATION_EVALUATION
// Module 5: Evaluates the impact of automated curation and generates the final scorecard
//

include { CURATIONEVALUATOR } from '../../../modules/local/curationevaluator/main'

workflow CURATION_EVALUATION {

    take:
    ch_pre_ref         // channel: [ val(meta), path(pre_ref_fasta) ]
    ch_post_ref        // channel: [ val(meta), path(post_ref_fasta) ]
    ch_blueprint       // channel: [ val(meta), path(blueprint_json) ]
    ch_curated_agp     // channel: [ val(meta), path(curated_agp) ]
    ch_curation_audit  // channel: [ val(meta), path(curation_audit_tsv) ]
    ch_gap_manifest    // channel: [ val(meta), path(gap_manifest_tsv) ]
    ch_walker_manifest // channel: [ val(meta), path(walker_manifest_tsv) ]
    ch_rdna_manifest   // channel: [ val(meta), path(rdna_manifest_tsv) ]
    ch_terminal_summary// channel: [ val(meta), path(terminal_ends_summary_tsv) ]
    ch_trans_manifest  // channel: [ val(meta), path(trans_manifest_tsv) ]
    ch_red_bed         // channel: [ val(meta), path(red_bed) ]

    main:
    ch_versions = Channel.empty()

    ch_eval_input = ch_pre_ref
        .join(ch_post_ref, by: 0)
        .join(ch_blueprint, by: 0)
        .join(ch_curated_agp, by: 0)
        .join(ch_curation_audit, by: 0)
        .join(ch_gap_manifest, by: 0)
        .join(ch_walker_manifest, by: 0)
        .join(ch_rdna_manifest, by: 0)
        .join(ch_terminal_summary, by: 0)
        .join(ch_trans_manifest, by: 0)
        .join(ch_red_bed, by: 0)

    CURATIONEVALUATOR (
        ch_eval_input
    )
    ch_versions = ch_versions.mix(CURATIONEVALUATOR.out.versions_evaluator)

    emit:
    scorecard_json = CURATIONEVALUATOR.out.scorecard_json
    report_md      = CURATIONEVALUATOR.out.report_md
    versions       = ch_versions
}
