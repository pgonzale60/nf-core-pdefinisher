//
// Subworkflow: PDE_CURATED_VALIDATION
// Automated decontamination, mitochondrial segregation, and residue validation
// for curated chromosome-level assemblies.
//

include { EXTRACT_UNMAPPED_HIFI } from '../../../modules/local/pde_curated_validation/main'
include { MASH_DECONTAMINATE     } from '../../../modules/local/pde_curated_validation/main'
include { MITOHIFI_SEGREGATE     } from '../../../modules/local/pde_curated_validation/main'
include { HIFIASM_RESIDUE        } from '../../../modules/local/pde_curated_validation/main'
include { CHARACTERIZE_RESIDUE   } from '../../../modules/local/pde_curated_validation/main'
include { BUILD_VALIDATION_PACKAGE } from '../../../modules/local/pde_curated_validation/main'

workflow PDE_CURATED_VALIDATION {
    take:
    ch_bam                // channel: [ val(meta), path(bam), path(bai) ]
    val_curated_fasta     // value / path: /path/to/curated.fa
    ch_mash_sketch        // channel: path(sketch_msh)
    ch_mito_fasta         // channel: path(mito_fasta)
    ch_mito_gbk           // channel: path(mito_gbk)
    ch_blast_nt           // channel: path(blast_nt)
    ch_taxonomy           // channel: path(taxonomy_dir)

    main:
    ch_versions = Channel.empty()

    // Step 1: Extract unmapped reads
    EXTRACT_UNMAPPED_HIFI(ch_bam)
    ch_versions = ch_versions.mix(EXTRACT_UNMAPPED_HIFI.out.versions)

    // Step 2: Mash decontamination
    MASH_DECONTAMINATE(
        EXTRACT_UNMAPPED_HIFI.out.fastq,
        ch_mash_sketch
    )
    ch_versions = ch_versions.mix(MASH_DECONTAMINATE.out.versions)

    // Step 3: MitoHiFi assembly & segregation
    MITOHIFI_SEGREGATE(
        MASH_DECONTAMINATE.out.decontaminated_fastq,
        ch_mito_fasta,
        ch_mito_gbk
    )
    ch_versions = ch_versions.mix(MITOHIFI_SEGREGATE.out.versions)

    // Step 4: De novo assemble residue reads
    HIFIASM_RESIDUE(
        MITOHIFI_SEGREGATE.out.residue_fastq
    )
    ch_versions = ch_versions.mix(HIFIASM_RESIDUE.out.versions)

    // Step 5: Characterize residue contigs
    ch_residue_combined = HIFIASM_RESIDUE.out.residue_contigs
        .join(MITOHIFI_SEGREGATE.out.residue_fastq)

    CHARACTERIZE_RESIDUE(
        ch_residue_combined,
        ch_blast_nt,
        ch_taxonomy
    )
    ch_versions = ch_versions.mix(CHARACTERIZE_RESIDUE.out.versions)

    // Step 6: Build Validation Artifact
    ch_package_inputs = EXTRACT_UNMAPPED_HIFI.out.fastq
        .join(MASH_DECONTAMINATE.out.contaminant_fastq)
        .join(MITOHIFI_SEGREGATE.out.mito_fastq)
        .join(MITOHIFI_SEGREGATE.out.residue_fastq)
        .join(HIFIASM_RESIDUE.out.residue_contigs)
        .join(CHARACTERIZE_RESIDUE.out.summary)

    BUILD_VALIDATION_PACKAGE(
        ch_package_inputs,
        val_curated_fasta
    )

    emit:
    validation_package = BUILD_VALIDATION_PACKAGE.out.validation_package
    unmapped_fastq     = EXTRACT_UNMAPPED_HIFI.out.fastq
    decontaminated_fq = MASH_DECONTAMINATE.out.decontaminated_fastq
    contaminant_fq    = MASH_DECONTAMINATE.out.contaminant_fastq
    mito_fasta        = MITOHIFI_SEGREGATE.out.mito_fasta
    mito_fastq        = MITOHIFI_SEGREGATE.out.mito_fastq
    residue_fastq     = MITOHIFI_SEGREGATE.out.residue_fastq
    residue_contigs   = HIFIASM_RESIDUE.out.residue_contigs
    blob_summary      = CHARACTERIZE_RESIDUE.out.summary
    versions          = ch_versions
}
