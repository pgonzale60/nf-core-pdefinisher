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
include { HIFIASM                         } from '../../../modules/local/hifiasm/main'
include { HIFIASM as HIFIASM_RDNA         } from '../../../modules/local/hifiasm/main'
include { MINIMAP2_ALIGN                  } from '../../../modules/local/minimap2/align/main'
include { MINIMAP2_ALIGN as ALIGN_UNPLACED } from '../../../modules/local/minimap2/align/main'
include { MINIMAP2_ALIGN as ALIGN_RDNA_TO_ASSEMBLY } from '../../../modules/local/minimap2/align/main'

include { SNIFFLES                        } from '../../../modules/local/sniffles/main'
include { FILTERCHROMOSOMES               } from '../../../modules/local/filterchromosomes/main'
include { RED                             } from '../../../modules/local/red/main'
include { IDENTIFYPUTATIVEGRS             } from '../../../modules/local/identifyputativegrs/main'
include { EXTRACTGAPEVIDENCE            } from '../../../modules/local/extractgapevidence/main'

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
    // MODULE: TELOMERICJUNCTIONPROFILER — PDE-aware subtelomeric junction deconvolution
    //
    def policy_candidates = [
        params.curation_policy,
        "${projectDir}/assets/curation_policy.json",
        "${workflow.projectDir}/assets/curation_policy.json",
        "/nfs/users/nfs_p/pg17/tests/nf-core-modules/assets/curation_policy.json"
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
    // MODULE: HIFIASM — De novo assembly of telomeric unitigs from extracted reads
    //
    ch_hifiasm_input = EXTRACTTELOMERICREADS.out.reads
        .map { meta, reads ->
            def read_list = (reads instanceof List) ? reads : [reads]
            tuple(meta, read_list, [])
        }

    HIFIASM (
        ch_hifiasm_input,
        [[:], []],
        [[:], [], []],
        [[:], []],
        [[:], []]
    )

    //
    // MODULE: MINIMAP2_ALIGN — Map assembled telomeric unitigs against baseline assembly
    //
    ch_minimap2_reads = HIFIASM.out.assembly_fasta
        .filter { _meta, fasta -> fasta != null && fasta.size() > 0 }

    ch_minimap2_input = ch_minimap2_reads
        .combine(ch_assembly, by: 0)
        .multiMap { meta, unitigs, assembly ->
            reads: tuple(meta, unitigs)
            reference: tuple(meta, assembly)
        }

    MINIMAP2_ALIGN (
        ch_minimap2_input.reads,
        ch_minimap2_input.reference,
        false,  // bam_format
        '',     // bam_index_extension
        true,   // cigar_paf_format
        false,  // cigar_bam
        false   // bed_bool
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
        .join(DELFIES.out.bed, by: 0)
        .join(MOSDEPTH.out.regions_bed, by: 0)
        .join(RED.out.bed, by: 0)
        .multiMap { meta, assembly, delfies, mosdepth, red_bed ->
            assembly: tuple(meta, assembly)
            delfies:  tuple(meta, delfies)
            mosdepth: tuple(meta, mosdepth)
            red:      tuple(meta, red_bed)
        }

    IDENTIFYPUTATIVEGRS (
        ch_identify_grs_input.assembly,
        ch_identify_grs_input.delfies,
        ch_identify_grs_input.mosdepth,
        ch_identify_grs_input.red
    )

    ch_meta_with_coverage = IDENTIFYPUTATIVEGRS.out.core_coverage
        .map { meta, cov_file ->
            def txt = cov_file.text.trim()
            if (!txt) {
                error "Empty core_coverage file for '${meta.id}': cannot establish the somatic core coverage baseline"
            }
            meta + [coverage: txt.toFloat()]
        }

    //
    // MODULE: HIFIASM (HIFIASM_RDNA) — High-speed string-graph assembly of the rDNA repeat array and anchors
    //
    ch_hifiasm_rdna_input = EXTRACTRDNAREADS.out.reads
        .map { meta, reads ->
            def read_list = (reads instanceof List) ? reads : [reads]
            tuple(meta, read_list, [])
        }

    HIFIASM_RDNA (
        ch_hifiasm_rdna_input,
        [[:], []],
        [[:], [], []],
        [[:], []],
        [[:], []]
    )

    ch_rdna_assembly_reads = HIFIASM_RDNA.out.assembly_fasta
        .map { meta, fasta ->
            def target_fa = fasta
            if (fasta instanceof List) {
                target_fa = fasta.find { it.name.endsWith('.p_utg.fa') } ?:
                            fasta.find { it.name.endsWith('.p_ctg.fa') } ?:
                            fasta[0]
            }
            tuple(meta, target_fa)
        }
        .filter { _meta, fasta -> fasta != null }

    //
    // MODULE: MINIMAP2_ALIGN (ALIGN_RDNA_TO_ASSEMBLY) — Map the rDNA assembly back onto
    // the baseline assembly. A terminal rDNA locus is resolved as
    // anchor -> rDNA model -> distal cap, and the single-copy core anchor is identified from
    // this mapping: a contig aligning to the locus chromosome at high MAPQ is the core flank.
    // Without it the cassette resolver can only fall back to guessing an anchor.
    //
    ch_rdna_align_input = ch_rdna_assembly_reads
        .map { meta, f -> [ meta.subMap(['id', 'chromosomes', 'telomere_repeat']), f ] }
        .join(ch_assembly, by: 0)
        .multiMap { meta, rdna_fa, asm ->
            reads:     tuple(meta, rdna_fa)
            reference: tuple(meta, asm)
        }

    ALIGN_RDNA_TO_ASSEMBLY (
        ch_rdna_align_input.reads,
        ch_rdna_align_input.reference,
        false,  // bam_format
        '',     // bam_index_extension
        true,   // cigar_paf_format
        false,  // cigar_bam
        false   // bed_bool
    )

    emit:
    chromosomes_fasta          = FILTERCHROMOSOMES.out.chromosomes
    non_chromosomal_fasta      = FILTERCHROMOSOMES.out.non_chromosomal
    red_bed                    = RED.out.bed
    unplaced_paf               = ALIGN_UNPLACED.out.paf
    delfies_bed                = DELFIES.out.bed
    delfies_fasta              = DELFIES.out.fasta
    mosdepth_summary           = MOSDEPTH.out.summary_txt
    mosdepth_regions           = MOSDEPTH.out.regions_bed
    barrnap_gff                = BARRNAP.out.gff
    rdna_reads                 = EXTRACTRDNAREADS.out.reads
    rdna_regions_bed           = EXTRACTRDNAREADS.out.bed
    rdna_assembly_fasta        = ch_rdna_assembly_reads.map { meta, f -> [ meta.subMap(['id', 'chromosomes', 'telomere_repeat']), f ] }
    rdna_assembly_gfa          = HIFIASM_RDNA.out.assembly_graphs.map { meta, f -> [ meta.subMap(['id', 'chromosomes', 'telomere_repeat']), f ] }
    rdna_assembly_info         = HIFIASM_RDNA.out.log.map { meta, f -> [ meta.subMap(['id', 'chromosomes', 'telomere_repeat']), f ] }
    rdna_assembly_paf          = ALIGN_RDNA_TO_ASSEMBLY.out.paf.map { meta, f -> [ meta.subMap(['id', 'chromosomes', 'telomere_repeat']), f ] }

    core_coverage              = ch_meta_with_coverage
    core_coverage_tsv          = IDENTIFYPUTATIVEGRS.out.core_coverage_tsv
    pde_segments_bed           = IDENTIFYPUTATIVEGRS.out.bed
    grs_summary_tsv            = IDENTIFYPUTATIVEGRS.out.summary
    telomeric_reads            = EXTRACTTELOMERICREADS.out.reads
    terminal_softclips         = SOFTCLIPENDS.out.softclips
    terminal_ends_summary      = SOFTCLIPENDS.out.summary
    telo_unitigs               = HIFIASM.out.assembly_fasta
    telo_paf                   = MINIMAP2_ALIGN.out.paf
    telo_junction_json         = TELOMERICJUNCTIONPROFILER.out.json
    telo_junction_tsv          = TELOMERICJUNCTIONPROFILER.out.tsv
    sniffles_vcf               = SNIFFLES.out.vcf
    sniffles_tbi               = SNIFFLES.out.tbi
    gaps_bed                   = EXTRACTGAPEVIDENCE.out.gaps_bed
    gap_evidence_json          = EXTRACTGAPEVIDENCE.out.gap_evidence_json
}
