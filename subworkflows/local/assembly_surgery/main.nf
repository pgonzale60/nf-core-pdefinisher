//
// SUBWORKFLOW: ASSEMBLY_SURGERY
// Module 4: Generates curated AGP v2.1, synthesizes curated FASTA, and performs official Sanger TOL AGP QC
//

include { CURATIONAGPGENERATOR } from '../../../modules/local/curationagpgenerator/main'
include { ASMFORMATQC          } from '../../../modules/local/asmformatqc/main'
include { RAGTAG_AGPCHECK      } from '../../../modules/local/ragtag/agpcheck/main'

workflow ASSEMBLY_SURGERY {

    take:
    ch_blueprint       // channel: [ val(meta), path(blueprint_json) ]
    ch_ref_fasta       // channel: [ val(meta), path(ref_fasta) ]
    ch_gap_manifest    // channel: [ val(meta), path(gap_manifest) ]
    ch_gap_patches     // channel: [ val(meta), path(gap_patches) ]
    ch_walker_manifest // channel: [ val(meta), path(walker_manifest) ]
    ch_walker_patches  // channel: [ val(meta), path(walker_patches) ]
    ch_rdna_manifest   // channel: [ val(meta), path(rdna_manifest) ]
    ch_rdna_cassettes  // channel: [ val(meta), path(rdna_cassettes) ]
    ch_trans_manifest  // channel: [ val(meta), path(trans_manifest) ]

    main:
    ch_versions = Channel.empty()

    ch_generator_input = ch_blueprint
        .join(ch_ref_fasta, by: 0)
        .join(ch_gap_patches, by: 0)
        .join(ch_gap_manifest, by: 0)
        .join(ch_walker_patches, by: 0)
        .join(ch_walker_manifest, by: 0)
        .join(ch_rdna_cassettes, by: 0)
        .join(ch_rdna_manifest, by: 0)
        .join(ch_trans_manifest, by: 0)

    CURATIONAGPGENERATOR (
        ch_generator_input
    )
    ch_versions = ch_versions.mix(CURATIONAGPGENERATOR.out.versions_agp)

    // Validate AGP syntax via RagTag agpcheck
    RAGTAG_AGPCHECK (
        CURATIONAGPGENERATOR.out.agp
    )
    ch_versions = ch_versions.mix(RAGTAG_AGPCHECK.out.versions_ragtag)

    // Validate and QC AGP through official Sanger TOL agp-tpf-utils
    ASMFORMATQC (
        CURATIONAGPGENERATOR.out.agp
    )
    ch_versions = ch_versions.mix(ASMFORMATQC.out.versions_agp_tpf_utils)

    emit:
    curated_agp       = ASMFORMATQC.out.agp
    components_fasta  = CURATIONAGPGENERATOR.out.components_fasta
    curated_fasta     = CURATIONAGPGENERATOR.out.curated_fasta
    segregated_fastas = CURATIONAGPGENERATOR.out.segregated_fastas
    audit             = CURATIONAGPGENERATOR.out.audit
    versions          = ch_versions
}
