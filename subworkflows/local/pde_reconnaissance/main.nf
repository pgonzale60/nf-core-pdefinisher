//
// SUBWORKFLOW: PDE_RECONNAISSANCE
// Automated pre-curation reconnaissance for Programmed DNA Elimination (PDE) assemblies
//

include { DELFIES                } from '../../../modules/local/delfies/main'
include { MOSDEPTH               } from '../../../modules/nf-core/mosdepth/main'
include { BARRNAP                } from '../../../modules/local/barrnap/main'
include { EXTRACTTELOMERICREADS  } from '../../../modules/local/extracttelomericreads/main'
include { EXTRACTRDNAREADS       } from '../../../modules/local/extractrdnareads/main'
include { SOFTCLIPENDS           } from '../../../modules/local/softclipends/main'
include { TELOMERICJUNCTIONPROFILER } from '../../../modules/local/telomericjunctionprofiler/main'
include { MINIMAP2_ALIGN as ALIGN_UNPLACED } from '../../../modules/local/minimap2/align/main'

include { SNIFFLES                        } from '../../../modules/local/sniffles/main'
include { FILTERCHROMOSOMES               } from '../../../modules/local/filterchromosomes/main'
include { RED                             } from '../../../modules/local/red/main'
include { IDENTIFYPUTATIVEGRS             } from '../../../modules/local/identifyputativegrs/main'
include { MOSDEPTHQUANTIZE               } from '../../../modules/local/mosdepthquantize/main'
include { MOSDEPTHCLASS                  } from '../../../modules/local/mosdepthclass/main'
include { EXTRACTGAPEVIDENCE            } from '../../../modules/local/extractgapevidence/main'
include { TELOCLIP_EXTRACT              } from '../../../modules/local/teloclip/extract/main'

workflow PDE_RECONNAISSANCE {

    take:
    ch_assembly         // channel: [ val(meta), path(assembly) ]
    ch_bam              // channel: [ val(meta), path(bam), path(bai) ]
    val_telomere_repeat // string: Telomeric repeat motif string (e.g. 'TTAGGC')
    val_chromosomes     // string or list: e.g. 'I,II,III,IV,V,X'

    main:

    //
    // MODULE: FILTERCHROMOSOMES — Separate canonical chromosomes from non-chromosomal sequences
    //
    FILTERCHROMOSOMES (
        ch_assembly,
        val_chromosomes
    )

    //
    // MODULE: RED — Repeat detection and masking across the entire assembly
    //
    RED (
        ch_assembly
    )

    //
    // MODULE: MINIMAP2_ALIGN (ALIGN_UNPLACED) — Align non-chromosomal contigs against chromosomes
    //
    ch_align_unplaced_input = FILTERCHROMOSOMES.out.non_chromosomal
        .combine(FILTERCHROMOSOMES.out.chromosomes, by: 0)
        .multiMap { meta, non_chr, chr ->
            reads: tuple(meta, non_chr)
            reference: tuple(meta, chr)
        }

    ALIGN_UNPLACED (
        ch_align_unplaced_input.reads,
        ch_align_unplaced_input.reference,
        false,  // bam_format
        '',     // bam_index_extension
        true,   // cigar_paf_format
        false,  // cigar_bam
        false   // bed_bool
    )

    //
    // MODULE: DELFIES — Somatic diminution breakpoint detection with neo-telomere addition
    //
    ch_delfies_input = ch_assembly
        .combine(ch_bam, by: 0)
        .map { meta, assembly, bam, bai ->
            tuple(meta, assembly, bam, bai)
        }

    DELFIES (
        ch_delfies_input
    )

    //
    // MODULE: MOSDEPTH — 500 bp sliding window coverage and chromosome summary
    //
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

    //
    // MODULE: BARRNAP — Predict ribosomal RNA genes across the assembly (eukaryota)
    //
    ch_barrnap_input = ch_assembly
        .map { meta, fasta -> tuple(meta, fasta, 'euk') }

    BARRNAP (
        ch_barrnap_input
    )

    //
    // MODULE: EXTRACTTELOMERICREADS — Dynamically extract telomeric reads from BAM
    //
    EXTRACTTELOMERICREADS (
        ch_bam,
        val_telomere_repeat
    )

    //
    // MODULE: SOFTCLIPENDS — Annotate softclipped alignments and overhangs at sequence ends
    //
    ch_softclipends_input = ch_assembly
        .combine(ch_bam, by: 0)
        .map { meta, assembly, bam, bai ->
            tuple(meta, assembly, bam, bai)
        }

    SOFTCLIPENDS (
        ch_softclipends_input,
        val_telomere_repeat
    )

    //
    // MODULE: TELOCLIP_EXTRACT — Dedicated teloclip filtering and soft-clipped overhang extraction
    //
    ch_teloclip_input = ch_bam
        .join(ch_assembly, by: 0)
        .map { meta, bam, bai, fasta ->
            tuple(meta, bam, bai, fasta)
        }

    TELOCLIP_EXTRACT (
        ch_teloclip_input,
        val_telomere_repeat
    )

    //
    // MODULE: TELOMERICJUNCTIONPROFILER — PDE-aware subtelomeric junction deconvolution
    //
    def policy_candidates = [
        params.curation_policy,
        "${projectDir}/assets/curation_policy.json",
        "${workflow.projectDir}/assets/curation_policy.json"
    ]
    def resolved_policy = policy_candidates.find { it && file(it).exists() }
    if (!resolved_policy) {
        error("Could not find assets/curation_policy.json in any expected location (projectDir: ${projectDir}, workflow.projectDir: ${workflow.projectDir}, or params.curation_policy).")
    }
    ch_policy_file = file(resolved_policy, checkIfExists: true)
    ch_profiler_input = EXTRACTTELOMERICREADS.out.reads
        .join(ch_assembly, by: 0)
        .join(DELFIES.out.bed, by: 0)
        .map { meta, reads, fasta, delfies ->
            tuple(meta, reads, fasta, delfies, ch_policy_file)
        }

    TELOMERICJUNCTIONPROFILER (
        ch_profiler_input,
        val_telomere_repeat
    )

    //
    // MODULE: SNIFFLES — Structural variant calling and somatic mosaicism detection
    //
    ch_sniffles_input = ch_bam
        .join(ch_assembly, by: 0)
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

    //
    // MODULE: EXTRACTGAPEVIDENCE — Granular physical read observations around assembly gaps
    //
    ch_gap_ev_input = ch_assembly
        .combine(ch_bam, by: 0)
        .join(SNIFFLES.out.vcf, by: 0)
        .map { meta, assembly, bam, bai, vcf ->
            tuple(meta, assembly, bam, bai, vcf)
        }

    EXTRACTGAPEVIDENCE (
        ch_gap_ev_input
    )

    //
    // MODULE: EXTRACTRDNAREADS — Extract HiFi reads covering 18S/28S rDNA arrays and flanking anchors
    //
    ch_extract_rdna_input = BARRNAP.out.gff
        .join(ch_bam, by: 0)
        .multiMap { meta, gff, bam, bai ->
            gff: tuple(meta, gff)
            bam: tuple(meta, bam, bai)
        }

    EXTRACTRDNAREADS (
        ch_extract_rdna_input.gff,
        ch_extract_rdna_input.bam,
        50000 // flank_bp: 50 kb flanking anchor into chromosome core
    )

    //
    // MODULE: IDENTIFYPUTATIVEGRS — Classify regions into putative_core vs putative_GRS
    //
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

    //
    // MODULE: MOSDEPTHQUANTIZE — Quantize coverage at the inferred T_GRS (MAPQ 0), flag RED repeat overlaps
    //
    ch_quantize_input = ch_bam
        .join(RED.out.bed, by: 0)
        .join(IDENTIFYPUTATIVEGRS.out.core_coverage_tsv, by: 0)
        .multiMap { meta, bam, bai, red_bed, cov_tsv ->
            bam:      tuple(meta, bam, bai)
            red:      tuple(meta, red_bed)
            cov:      tuple(meta, cov_tsv)
        }

    MOSDEPTHQUANTIZE (
        ch_quantize_input.bam,
        ch_quantize_input.red,
        ch_quantize_input.cov
    )

    // MODULE: MOSDEPTHCLASS — per-region depth over the final classification segments (core / GRS / uncertain)
    ch_class_input = ch_bam
        .join(IDENTIFYPUTATIVEGRS.out.bed, by: 0)
        .multiMap { meta, bam, bai, segments_bed ->
            bam:      tuple(meta, bam, bai)
            segments: tuple(meta, segments_bed)
        }

    MOSDEPTHCLASS (
        ch_class_input.bam,
        ch_class_input.segments
    )

    ch_meta_with_coverage = IDENTIFYPUTATIVEGRS.out.core_coverage
        .map { meta, cov_file ->
            def txt = cov_file.text.trim()
            if (!txt) {
                error "Empty core_coverage file for '${meta.id}': cannot establish the somatic core coverage baseline"
            }
            meta + [coverage: txt.toFloat()]
        }

    emit:
    chromosomes_fasta          = FILTERCHROMOSOMES.out.chromosomes
    non_chromosomal_fasta      = FILTERCHROMOSOMES.out.non_chromosomal
    red_bed                    = RED.out.bed
    unplaced_paf               = ALIGN_UNPLACED.out.paf
    delfies_bed                = DELFIES.out.bed
    delfies_fasta              = DELFIES.out.fasta
    delfies_foci_s2g           = DELFIES.out.foci_s2g
    grs_quantized_bed          = MOSDEPTHQUANTIZE.out.quantized_bed
    class_regions_bed          = MOSDEPTHCLASS.out.class_regions_bed
    mosdepth_summary           = MOSDEPTH.out.summary_txt
    mosdepth_regions           = MOSDEPTH.out.regions_bed
    barrnap_gff                = BARRNAP.out.gff
    rdna_reads                 = EXTRACTRDNAREADS.out.reads
    rdna_regions_bed           = EXTRACTRDNAREADS.out.bed

    core_coverage              = ch_meta_with_coverage
    core_coverage_tsv          = IDENTIFYPUTATIVEGRS.out.core_coverage_tsv
    pde_segments_bed           = IDENTIFYPUTATIVEGRS.out.bed
    grs_summary_tsv            = IDENTIFYPUTATIVEGRS.out.summary
    telomeric_reads            = EXTRACTTELOMERICREADS.out.reads
    terminal_softclips         = SOFTCLIPENDS.out.softclips
    terminal_ends_summary      = SOFTCLIPENDS.out.summary
    teloclip_overhangs         = TELOCLIP_EXTRACT.out.overhangs
    telo_junction_json         = TELOMERICJUNCTIONPROFILER.out.json
    telo_junction_tsv          = TELOMERICJUNCTIONPROFILER.out.tsv
    sniffles_vcf               = SNIFFLES.out.vcf
    sniffles_tbi               = SNIFFLES.out.tbi
    gaps_bed                   = EXTRACTGAPEVIDENCE.out.gaps_bed
    gap_evidence_json          = EXTRACTGAPEVIDENCE.out.gap_evidence_json
}
