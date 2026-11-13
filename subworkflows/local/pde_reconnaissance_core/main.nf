//
// SUBWORKFLOW: PDE_RECONNAISSANCE_CORE
// Reconnaissance subset required by RECON_MOTIF: GRS classification, delfies breakpoints and mosdepth coverage (MOSDEPTH, MOSDEPTHQUANTIZE, MOSDEPTHCLASS)
//

include { RED                  } from '../../../modules/local/red/main'
include { DELFIES              } from '../../../modules/local/delfies/main'
include { MOSDEPTH             } from '../../../modules/nf-core/mosdepth/main'
include { IDENTIFYPUTATIVEGRS  } from '../../../modules/local/identifyputativegrs/main'
include { MOSDEPTHQUANTIZE     } from '../../../modules/local/mosdepthquantize/main'
include { MOSDEPTHCLASS        } from '../../../modules/local/mosdepthclass/main'
include { GFASTATSGAPS                   } from '../../../modules/local/gfastatsgaps/main'
include { SEGMENTQUANTIZED               } from '../../../modules/local/segmentquantized/main'

workflow PDE_RECONNAISSANCE_CORE {

    take:
    ch_assembly          // channel: [ val(meta), path(assembly) ]
    ch_bam               // channel: [ val(meta), path(bam), path(bai) ]
    val_chromosomes      // string or list: e.g. 'I,II,III,IV,V,X'

    main:

    RED (
        ch_assembly
    )

    ch_delfies_input = ch_assembly
        .combine(ch_bam, by: 0)
        .map { meta, assembly, bam, bai ->
            tuple(meta, assembly, bam, bai)
        }

    DELFIES (
        ch_delfies_input
    )

    ch_mosdepth_input = ch_bam
        .join(ch_assembly, by: 0)
        .multiMap { meta, bam, bai, fasta ->
            bam:   tuple(meta, bam, bai, [])
            fasta: tuple(meta, fasta)
        }

    MOSDEPTH (
        ch_mosdepth_input.bam,
        ch_mosdepth_input.fasta,
        []
    )

    ch_identify_grs_input = ch_assembly
        .join(DELFIES.out.foci_s2g, by: 0)
        .join(MOSDEPTH.out.regions_bed, by: 0)
        .join(RED.out.bed, by: 0)
        .multiMap { meta, assembly, foci, mosdepth, red_bed ->
            assembly: tuple(meta, assembly)
            foci:     tuple(meta, foci)
            mosdepth: tuple(meta, mosdepth)
            red:      tuple(meta, red_bed)
        }

    IDENTIFYPUTATIVEGRS (
        ch_identify_grs_input.assembly,
        ch_identify_grs_input.foci,
        ch_identify_grs_input.mosdepth,
        ch_identify_grs_input.red,
        val_chromosomes
    )

    ch_quantize_input = ch_bam
        .join(RED.out.bed, by: 0)
        .join(IDENTIFYPUTATIVEGRS.out.core_coverage_tsv, by: 0)
        .multiMap { meta, bam, bai, red_bed, cov_tsv ->
            bam: tuple(meta, bam, bai)
            red: tuple(meta, red_bed)
            cov: tuple(meta, cov_tsv)
        }

    MOSDEPTHQUANTIZE (
        ch_quantize_input.bam,
        ch_quantize_input.red,
        ch_quantize_input.cov
    )

    // MODULE: GFASTATSGAPS — assembly N-gap intervals, input to SEGMENTQUANTIZED
    GFASTATSGAPS (
        ch_assembly
    )

    // MODULE: SEGMENTQUANTIZED — final three-class call from quantized depth, anchored by DELFIES breakpoints
    ch_segment_input = IDENTIFYPUTATIVEGRS.out.summary
        .join(MOSDEPTHQUANTIZE.out.quantized_bed, by: 0)
        .join(GFASTATSGAPS.out.gaps_bed, by: 0)
        .join(RED.out.bed, by: 0)
        .join(MOSDEPTH.out.regions_bed, by: 0)
        .join(IDENTIFYPUTATIVEGRS.out.breakpoints_tsv, by: 0)
        .multiMap { meta, summary, quantized, gaps, red_bed, regions, breakpoints ->
            summary:     tuple(meta, summary)
            quantized:   tuple(meta, quantized)
            gaps:        tuple(meta, gaps)
            red:         tuple(meta, red_bed)
            regions:     tuple(meta, regions)
            breakpoints: tuple(meta, breakpoints)
        }

    SEGMENTQUANTIZED (
        ch_segment_input.summary,
        ch_segment_input.quantized,
        ch_segment_input.gaps,
        ch_segment_input.red,
        ch_segment_input.regions,
        ch_segment_input.breakpoints,
        val_chromosomes
    )

    // MODULE: MOSDEPTHCLASS — per-region depth over the final classification segments (core / GRS / uncertain)
    ch_class_input = ch_bam
        .join(SEGMENTQUANTIZED.out.bed, by: 0)
        .multiMap { meta, bam, bai, segments_bed ->
            bam:      tuple(meta, bam, bai)
            segments: tuple(meta, segments_bed)
        }

    MOSDEPTHCLASS (
        ch_class_input.bam,
        ch_class_input.segments
    )

    emit:
    red_bed                    = RED.out.bed
    delfies_bed                = DELFIES.out.bed
    delfies_foci_s2g           = DELFIES.out.foci_s2g
    mosdepth_summary           = MOSDEPTH.out.summary_txt
    mosdepth_regions           = MOSDEPTH.out.regions_bed
    grs_quantized_bed          = MOSDEPTHQUANTIZE.out.quantized_bed
    class_regions_bed          = MOSDEPTHCLASS.out.class_regions_bed
    class_regions_csi          = MOSDEPTHCLASS.out.class_regions_csi
    grs_summary_tsv            = SEGMENTQUANTIZED.out.summary
    pde_segments_bed           = SEGMENTQUANTIZED.out.bed
}
