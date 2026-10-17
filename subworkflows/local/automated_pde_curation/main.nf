//
// SUBWORKFLOW: AUTOMATED_PDE_CURATION
// Master End-to-End Orchestrator for Sanger Automated PDE Curation
//
// Chains:
// Module 1: Reconnaissance
// Module 2: Synthesis (Blueprint)
// Module 3: Resolution Engines (A: Gap, B: Telomere, C: rDNA)
// Module 4: Assembly Surgery & AGP Validation
// Module 5: Post-Surgery Evaluation & Impact Scorecard
//

include { PDE_RECONNAISSANCE   } from '../pde_reconnaissance/main'
include { BGZIPFASTA          } from '../../../modules/local/bgzipfasta/main'
include { CURATIONSYNTHESIS   } from '../../../modules/local/curationsynthesis/main'
include { CURATION_RESOLUTION  } from '../curation_resolution/main'
include { ASSEMBLY_SURGERY     } from '../assembly_surgery/main'
include { CURATION_EVALUATION  } from '../curation_evaluation/main'
include { CURATION_VALIDATION  } from '../curation_validation/main'

workflow AUTOMATED_PDE_CURATION {

    take:
    ch_assembly         // channel: [ val(meta), path(assembly_fasta) ]
    ch_bam              // channel: [ val(meta), path(bam), path(bai) ]
    val_telomere_repeat // val: e.g. 'TTAGGC'
    val_chromosomes     // val: e.g. 'I,II,III,IV,V,X'

    main:
    ch_versions = Channel.empty()

    // =========================================================================
    // MODULE 1: PDE RECONNAISSANCE
    // =========================================================================
    PDE_RECONNAISSANCE (
        ch_assembly,
        ch_bam,
        val_telomere_repeat,
        val_chromosomes
    )

    // =========================================================================
    // MODULE 2: CURATION SYNTHESIS (BLUEPRINT)
    // =========================================================================
    ch_s1  = ch_assembly.map { meta, fasta -> [ meta.id, meta, fasta ] }
    ch_s2  = PDE_RECONNAISSANCE.out.grs_summary_tsv.map { meta, f -> [ meta.id, f ] }
    ch_s3  = PDE_RECONNAISSANCE.out.mosdepth_regions.map { meta, f -> [ meta.id, f ] }
    ch_s4  = PDE_RECONNAISSANCE.out.terminal_ends_summary.map { meta, f -> [ meta.id, f ] }
    ch_s5  = PDE_RECONNAISSANCE.out.terminal_softclips.map { meta, f -> [ meta.id, f ] }
    ch_s6  = PDE_RECONNAISSANCE.out.delfies_bed.map { meta, f -> [ meta.id, f ] }
    ch_s7  = PDE_RECONNAISSANCE.out.rdna_regions_bed.map { meta, f -> [ meta.id, f ] }
    ch_s8  = PDE_RECONNAISSANCE.out.barrnap_gff.map { meta, f -> [ meta.id, f ] }
    ch_s9  = PDE_RECONNAISSANCE.out.rdna_assembly_info.map { meta, f -> [ meta.id, f ] }
    ch_s10 = PDE_RECONNAISSANCE.out.unplaced_paf.map { meta, f -> [ meta.id, f ] }
    ch_s11 = PDE_RECONNAISSANCE.out.red_bed.map { meta, f -> [ meta.id, f ] }
    ch_s12 = PDE_RECONNAISSANCE.out.telo_paf.map { meta, f -> [ meta.id, f ] }
    ch_s13 = PDE_RECONNAISSANCE.out.gap_evidence_json.map { meta, f -> [ meta.id, f ] }
    ch_s14 = PDE_RECONNAISSANCE.out.telo_junction_json.map { meta, f -> [ meta.id, f ] }

    ch_synthesis_input = ch_s1
        .join(ch_s2, by: 0)
        .join(ch_s3, by: 0)
        .join(ch_s4, by: 0)
        .join(ch_s5, by: 0)
        .join(ch_s6, by: 0)
        .join(ch_s7, by: 0)
        .join(ch_s8, by: 0)
        .join(ch_s9, by: 0)
        .join(ch_s10, by: 0)
        .join(ch_s11, by: 0)
        .join(ch_s12, by: 0)
        .join(ch_s13, by: 0)
        .join(ch_s14, by: 0)
        .map { id, meta, fasta, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14 ->
            [ meta, fasta, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14 ]
        }

    CURATIONSYNTHESIS (
        ch_synthesis_input,
        val_telomere_repeat
    )
    ch_versions = ch_versions.mix(CURATIONSYNTHESIS.out.versions_synthesis)

    //
    // MODULE: BGZIPFASTA — normalise the assembly for the resolution engines.
    //
    // Assemblies arrive as plain gzip, which `samtools faidx` cannot index in any
    // directory ("Cannot index files compressed with gzip, please use bgzip") and which
    // the engines' text-mode fallback could not read either, so each engine used to
    // decompress the assembly into its own work directory just to resolve a region. Do it
    // once here instead. `remainder: true` keeps the sample in the channel even if bgzip
    // wrote no .gzi, so a missing index degrades to the engine's fallback rather than
    // silently dropping the sample from the pipeline.
    //
    BGZIPFASTA (
        ch_assembly
    )
    // The .gzi is not optional and there is no fallback: it is what makes the bgzipped
    // reference randomly accessible to the engines below. `remainder: true` keeps the
    // sample in the channel so a missing index is *reported* rather than silently
    // dropping the sample from the pipeline.
    //
    ch_indexed_assembly = BGZIPFASTA.out.fasta
        .join(BGZIPFASTA.out.gzi, by: 0, remainder: true)
        .map { meta, fasta, gzi ->
            if (!gzi) {
                error("BGZIPFASTA produced no .gzi for ${meta.id} (${fasta}). The resolution " +
                      "engines require an indexed reference and there is no fallback.")
            }
            tuple(meta, fasta, gzi)
        }

    // =========================================================================
    // MODULE 3: CURATION RESOLUTION ENGINES (A, B, C)
    // =========================================================================
    // The reconnaissance's soft-clipped-read summary stays out of Module 3. It is consumed by
    // CURATIONSYNTHESIS and CURATION_EVALUATION, but the walker's optional profile slot expects a
    // PROFILE_TELOMERES table -- a different file, which this pipeline does not produce -- and handing
    // the summary over as one is how a unitig cap came to report an entire unitig as its extension.
    // A unitig's interval now travels on the Module 2 card as evidence.unitig_extension_region.

    ch_telo_unitigs  = PDE_RECONNAISSANCE.out.telo_unitigs

    ch_flye_evidence = PDE_RECONNAISSANCE.out.rdna_assembly_fasta
        .join(PDE_RECONNAISSANCE.out.rdna_assembly_info, by: 0)
        .join(PDE_RECONNAISSANCE.out.rdna_assembly_paf, by: 0)
        .map { meta, fa, info, paf ->
            tuple(meta, fa, info, paf)
        }

    ch_gaps_bed      = PDE_RECONNAISSANCE.out.gaps_bed
    ch_policy_json   = ch_assembly.map { meta, fasta -> [ meta, [] ] }

    CURATION_RESOLUTION (
        CURATIONSYNTHESIS.out.blueprint_json,
        ch_assembly,
        ch_bam,
        ch_gaps_bed,
        ch_policy_json,
        ch_telo_unitigs,
        PDE_RECONNAISSANCE.out.unplaced_paf,
        ch_flye_evidence,
        ch_indexed_assembly,
        val_telomere_repeat
    )
    ch_versions = ch_versions.mix(CURATION_RESOLUTION.out.versions)

    // =========================================================================
    // MODULE 4: ASSEMBLY SURGERY & AGP VALIDATION
    // =========================================================================
    ASSEMBLY_SURGERY (
        CURATIONSYNTHESIS.out.blueprint_json,
        ch_assembly,
        CURATION_RESOLUTION.out.gap_manifest,
        CURATION_RESOLUTION.out.gap_patches,
        CURATION_RESOLUTION.out.walker_manifest,
        CURATION_RESOLUTION.out.walker_patches,
        CURATION_RESOLUTION.out.rdna_manifest,
        CURATION_RESOLUTION.out.rdna_cassettes,
        CURATION_RESOLUTION.out.trans_manifest
    )
    ch_versions = ch_versions.mix(ASSEMBLY_SURGERY.out.versions)

    // =========================================================================
    // MODULE 5: POST-SURGERY EVALUATION & SCORECARD
    // =========================================================================
    CURATION_EVALUATION (
        ch_assembly,
        ASSEMBLY_SURGERY.out.curated_fasta,
        CURATIONSYNTHESIS.out.blueprint_json,
        ASSEMBLY_SURGERY.out.curated_agp,
        ASSEMBLY_SURGERY.out.audit,
        CURATION_RESOLUTION.out.gap_manifest,
        CURATION_RESOLUTION.out.walker_manifest,
        CURATION_RESOLUTION.out.rdna_manifest,
        PDE_RECONNAISSANCE.out.terminal_ends_summary,
        CURATION_RESOLUTION.out.trans_manifest,
        PDE_RECONNAISSANCE.out.red_bed
    )
    ch_versions = ch_versions.mix(CURATION_EVALUATION.out.versions)

    // =========================================================================
    // MODULE 6: POST-CURATION MAPPING & SNIFFLES VALIDATION
    // =========================================================================
    ch_val_reads = ch_bam.map { meta, bam, bai -> tuple(meta, bam) }
    CURATION_VALIDATION (
        ch_val_reads,
        ASSEMBLY_SURGERY.out.curated_fasta,
        ASSEMBLY_SURGERY.out.curated_agp,
        ASSEMBLY_SURGERY.out.audit
    )
    ch_versions = ch_versions.mix(CURATION_VALIDATION.out.versions)

    emit:
    blueprint_json    = CURATIONSYNTHESIS.out.blueprint_json
    blueprint_md      = CURATIONSYNTHESIS.out.blueprint_md
    curated_agp       = ASSEMBLY_SURGERY.out.curated_agp
    curated_fasta     = ASSEMBLY_SURGERY.out.curated_fasta
    components_fasta  = ASSEMBLY_SURGERY.out.components_fasta
    segregated_fastas = ASSEMBLY_SURGERY.out.segregated_fastas
    curation_audit    = ASSEMBLY_SURGERY.out.audit
    scorecard_json    = CURATION_EVALUATION.out.scorecard_json
    report_md         = CURATION_EVALUATION.out.report_md
    validation_json   = CURATION_VALIDATION.out.validation_json
    validation_md     = CURATION_VALIDATION.out.validation_md
    curated_bam       = CURATION_VALIDATION.out.curated_bam
    curated_bai       = CURATION_VALIDATION.out.curated_bai
    versions          = ch_versions
}
