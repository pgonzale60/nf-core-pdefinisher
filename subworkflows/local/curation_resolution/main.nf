//
// SUBWORKFLOW: CURATION_RESOLUTION
// Module 3: Orchestrates Curation Resolution Engines (A: Gap Assembly, B: Telomere Walker, C: rDNA Cassette)
//

include { TARGETEDGAPASSEMBLY         } from '../../../modules/local/targetedgapassembly/main'
include { TELOMEREWALKER               } from '../../../modules/local/telomerewalker/main'
include { RDNACASSETTERESOLVER        } from '../../../modules/local/rdnacassetteresolver/main'
include { TRANSCHROMOSOMALLINKRESOLVER } from '../../../modules/local/transchromosomallinkresolver/main'

workflow CURATION_RESOLUTION {

    take:
    ch_blueprint        // channel: [ val(meta), path(blueprint_json) ]
    ch_ref_fasta        // channel: [ val(meta), path(ref_fasta) ]
    ch_bam              // channel: [ val(meta), path(bam), path(bai) ]
    ch_gaps_bed         // channel: [ val(meta), path(gaps_bed) ]
    ch_policy_json      // channel: [ val(meta), path(policy_json) ]
    ch_telo_unitigs     // channel: [ val(meta), path(telo_unitigs_fasta) ]
    ch_unplaced_paf     // channel: [ val(meta), path(unplaced_paf) ]
    ch_flye_assembly    // channel: [ val(meta), path(flye_fasta), path(flye_info), path(flye_paf) ]
    ch_indexed_ref      // channel: [ val(meta), path(ref_bgzipped), path(ref_gzi) ]
    val_telomere_repeat // val: e.g. 'TTAGGC'

    main:
    ch_versions = Channel.empty()

    // -------------------------------------------------------------------------
    // ENGINE A: TARGETED GAP ASSEMBLY
    // -------------------------------------------------------------------------
    ch_gap_input = ch_blueprint
        .join(ch_bam, by: 0)
        .join(ch_indexed_ref, by: 0)
        .join(ch_policy_json, by: 0)
        // join() flattens rather than nesting, so the indexed reference's two paths arrive
        // as two separate arguments. Every element of the joined tuple gets a name here:
        // there is no "is it a list this time" guesswork, and a channel whose width changes
        // fails loudly instead of destructuring into something loose.
        .map { meta, blueprint, bam, bai, ref_gz, ref_gzi, policy ->
            tuple(meta, blueprint, bam, bai, ref_gz, ref_gzi, policy)
        }

    TARGETEDGAPASSEMBLY (
        ch_gap_input
    )
    ch_versions = ch_versions.mix(TARGETEDGAPASSEMBLY.out.versions_assembly)

    // -------------------------------------------------------------------------
    // ENGINE B: TELOMERE WALKER
    // -------------------------------------------------------------------------
    ch_walker_input = ch_blueprint
        .join(ch_bam, by: 0)
        .join(ch_indexed_ref, by: 0)
        .join(ch_telo_unitigs, by: 0)
        .join(ch_unplaced_paf, by: 0)
        // Five joins, one name per resulting element (see the note on the gap input above):
        // meta, blueprint, bam, bai, ref_gz, ref_gzi, unitigs, unplaced_paf.
        //
        // The walker also accepts a PROFILE_TELOMERES table of unitig alignment coordinates, and that
        // slot is passed empty here on purpose. Nothing in this pipeline produces that table: the
        // reconnaissance emits a soft-clipped-read summary, which is a different file with a
        // different shape, and handing it over as a profile is how a cap came to report a whole
        // unitig as its extension. A unitig's interval now arrives on the Module 2 card
        // (evidence.unitig_extension_region), which is the stage that measured it.
        .map { meta, blueprint, bam, bai, ref_gz, ref_gzi, unitigs, unplaced_paf ->
            tuple(meta, blueprint, [], bam, bai, ref_gz, ref_gzi, unitigs, unplaced_paf)
        }

    TELOMEREWALKER (
        ch_walker_input,
        val_telomere_repeat
    )
    ch_versions = ch_versions.mix(TELOMEREWALKER.out.versions_walker)

    // -------------------------------------------------------------------------
    // ENGINE C: RDNA CASSETTE RESOLVER
    // -------------------------------------------------------------------------
    ch_rdna_input = ch_blueprint
        .join(ch_ref_fasta, by: 0)
        .join(ch_flye_assembly, by: 0)
        .join(ch_telo_unitigs, by: 0)
        .map { meta, blueprint, ref, flye_fa, flye_info, flye_paf, telo_utgs ->
            tuple(meta, blueprint, ref, flye_fa, flye_info, flye_paf, telo_utgs)
        }

    RDNACASSETTERESOLVER (
        ch_rdna_input
    )
    ch_versions = ch_versions.mix(RDNACASSETTERESOLVER.out.versions_resolver)

    // -------------------------------------------------------------------------
    // ENGINE D: TRANSCHROMOSOMAL LINK RESOLVER
    // -------------------------------------------------------------------------
    ch_trans_input = ch_blueprint
        .join(ch_bam, by: 0)
        .join(ch_ref_fasta, by: 0)
        .join(ch_gaps_bed, by: 0)
        .join(ch_policy_json, by: 0)
        .map { meta, blueprint, bam, bai, ref, gaps_bed, policy ->
            tuple(meta, blueprint, bam, bai, ref, gaps_bed, policy)
        }

    TRANSCHROMOSOMALLINKRESOLVER (
        ch_trans_input
    )
    ch_versions = ch_versions.mix(TRANSCHROMOSOMALLINKRESOLVER.out.versions_trans)

    emit:
    gap_manifest    = TARGETEDGAPASSEMBLY.out.manifest.map { meta, dir, tsv -> [ meta, tsv ] }
    gap_patches     = TARGETEDGAPASSEMBLY.out.manifest.map { meta, dir, tsv -> [ meta, dir ] }
    walker_manifest = TELOMEREWALKER.out.manifest.map { meta, dir, tsv -> [ meta, tsv ] }
    walker_patches  = TELOMEREWALKER.out.manifest.map { meta, dir, tsv -> [ meta, dir ] }
    rdna_manifest   = RDNACASSETTERESOLVER.out.manifest.map { meta, dir, tsv -> [ meta, tsv ] }
    rdna_cassettes  = RDNACASSETTERESOLVER.out.manifest.map { meta, dir, tsv -> [ meta, dir ] }
    trans_manifest  = TRANSCHROMOSOMALLINKRESOLVER.out.manifest
    versions        = ch_versions
}
