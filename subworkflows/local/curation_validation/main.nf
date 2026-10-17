//
// SUBWORKFLOW: CURATION_VALIDATION
// Re-maps HiFi reads to the post-surgery curated assembly, runs Sniffles SV calling,
// and evaluates concordance across surgical edits and discovers flagrant errors.
//

include { MINIMAP2_ALIGN     } from '../../../modules/local/minimap2/align/main'
include { SNIFFLES           } from '../../../modules/local/sniffles/main'
include { CURATIONVALIDATION } from '../../../modules/local/curationvalidation/main'

workflow CURATION_VALIDATION {

    take:
    ch_reads        // channel: [ val(meta), path(reads) ] (or [val(meta), path(bam), path(bai)])
    ch_curated_fa   // channel: [ val(meta), path(curated_fasta) ]
    ch_curated_agp  // channel: [ val(meta), path(curated_agp) ]
    ch_curation_audit // channel: [ val(meta), path(curation_audit_tsv) ]

    main:
    ch_versions = Channel.empty()

    // 1. Map HiFi reads to curated assembly
    ch_map_input = ch_reads
        .join(ch_curated_fa, by: 0)
        .multiMap { meta, reads, fasta ->
            reads: tuple(meta, reads)
            fasta: tuple(meta, fasta)
        }

    MINIMAP2_ALIGN (
        ch_map_input.reads,
        ch_map_input.fasta,
        true,   // bam_format
        'bai',  // bam_index_extension
        false,  // cigar_paf_format
        false,  // cigar_bam
        false   // bed_bool
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions_minimap2)

    // 2. Call Structural Variants with Sniffles
    ch_sniffles_input = MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index, by: 0)
        .join(ch_curated_fa, by: 0)
        .multiMap { meta, bam, bai, fasta ->
            bam:   tuple(meta, bam, bai)
            fasta: tuple(meta, fasta)
        }

    SNIFFLES (
        ch_sniffles_input.bam,
        ch_sniffles_input.fasta,
        [[:], []],       // tandem_repeats BED
        true,            // vcf_output
        false            // snf_output
    )
    ch_versions = ch_versions.mix(SNIFFLES.out.versions_sniffles)

    // 3. Evaluate Locus Concordance and Flagrant Assembly Errors
    ch_val_input = ch_curation_audit
        .join(ch_curated_agp, by: 0)
        .join(SNIFFLES.out.vcf, by: 0)
        .join(SNIFFLES.out.tbi, by: 0)
        .join(MINIMAP2_ALIGN.out.bam, by: 0)
        .join(MINIMAP2_ALIGN.out.index, by: 0)
        .map { meta, audit, agp, vcf, tbi, bam, bai ->
            tuple(meta, audit, agp, vcf, tbi, bam, bai)
        }

    CURATIONVALIDATION (
        ch_val_input
    )
    ch_versions = ch_versions.mix(CURATIONVALIDATION.out.versions_python)

    emit:
    curated_bam     = MINIMAP2_ALIGN.out.bam
    curated_bai     = MINIMAP2_ALIGN.out.index
    sniffles_vcf    = SNIFFLES.out.vcf
    sniffles_tbi    = SNIFFLES.out.tbi
    validation_json = CURATIONVALIDATION.out.validation_json
    validation_md   = CURATIONVALIDATION.out.validation_md
    versions        = ch_versions
}
