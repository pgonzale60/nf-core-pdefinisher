//
// SUBWORKFLOW: GRAPH_PREREQUISITES
// Computes missing gfapdefinisher prerequisites from primary GFA and HiFi reads:
//   1. GFA unitig FASTA extraction (when Tiara or Barrnap is not pre-supplied)
//   2. Tiara taxonomic/contaminant classification on primary GFA unitigs
//   3. Barrnap eukaryotic rDNA prediction on primary GFA unitigs
//   4. Telomere end-family archive (BBDuk fishing -> NCRF -> minimap2 clustering -> Racon consensus -> summary)
// Pre-supplied samplesheet artifacts (tiara, barrnap, end_archive) are reused as-is.
//

include { GFATOFASTA                  } from '../../../modules/local/gfatofasta/main'
include { TIARA                       } from '../../../modules/local/tiara/main'
include { BARRNAP as BARRNAP_GFA      } from '../../../modules/local/barrnap/main'
include { FISHTELOMERICREADS          } from '../../../modules/local/fishtelomericreads/main'
include { NCRFTELOMERES               } from '../../../modules/local/ncrftelomeres/main'
include { RACON_CONSENSUS             } from '../../../modules/local/racon/main'
include { TELOMEREENDCLUSTERS_CLUSTER   } from '../../../modules/local/telomereendclusters/main'
include { TELOMEREENDCLUSTERS_SUMMARIZE } from '../../../modules/local/telomereendclusters/main'

workflow GRAPH_PREREQUISITES {

    take:
    ch_graph_samples     // channel: [ val(delivery_meta), path(assembly), path(gfa), val(inputs_map), path(bam_or_empty) ]
    end_clusters_config  // path: telomeric_end_clusters.json
    cluster_script       // path: scripts/telomere_end_clusters.py

    main:

    // -------------------------------------------------------------------------
    // 1. Export primary GFA unitigs to FASTA if Tiara or Barrnap is needed
    // -------------------------------------------------------------------------
    ch_need_unitigs = ch_graph_samples
        .filter { _meta, _asm, _gfa, inputs, _bam -> !inputs.tiara || !inputs.barrnap }
        .map { meta, _asm, gfa, _inputs, _bam -> tuple(meta, gfa) }

    GFATOFASTA(ch_need_unitigs)
    ch_unitigs_by_id = GFATOFASTA.out.fasta.map { meta, fasta -> tuple(meta.id, fasta) }

    // -------------------------------------------------------------------------
    // 2. Tiara on primary GFA unitigs (or reuse precomputed tiara)
    // -------------------------------------------------------------------------
    ch_tiara_provided = ch_graph_samples
        .filter { _meta, _asm, _gfa, inputs, _bam -> inputs.tiara }
        .map { meta, _asm, _gfa, inputs, _bam -> tuple(meta.id, file(inputs.tiara, checkIfExists: true)) }

    ch_need_tiara = ch_graph_samples
        .filter { _meta, _asm, _gfa, inputs, _bam -> !inputs.tiara }
        .map { meta, _asm, _gfa, _inputs, _bam -> tuple(meta.id, meta) }
        .join(ch_unitigs_by_id, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, fasta -> tuple(meta, fasta) }

    TIARA(ch_need_tiara)
    ch_tiara = ch_tiara_provided.mix(TIARA.out.classifications.map { meta, f -> tuple(meta.id, f) })

    // -------------------------------------------------------------------------
    // 3. Barrnap on primary GFA unitigs (or reuse precomputed barrnap)
    // -------------------------------------------------------------------------
    ch_barrnap_provided = ch_graph_samples
        .filter { _meta, _asm, _gfa, inputs, _bam -> inputs.barrnap }
        .map { meta, _asm, _gfa, inputs, _bam -> tuple(meta.id, file(inputs.barrnap, checkIfExists: true)) }

    ch_need_barrnap = ch_graph_samples
        .filter { _meta, _asm, _gfa, inputs, _bam -> !inputs.barrnap }
        .map { meta, _asm, _gfa, _inputs, _bam -> tuple(meta.id, meta) }
        .join(ch_unitigs_by_id, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, fasta -> tuple(meta, fasta, 'euk') }

    BARRNAP_GFA(ch_need_barrnap)
    ch_barrnap = ch_barrnap_provided.mix(BARRNAP_GFA.out.gff.map { meta, f -> tuple(meta.id, f) })

    // Single source of truth: read fishing & NCRF settings from telomeric_end_clusters.json
    def end_cfg = new groovy.json.JsonSlurper().parse(end_clusters_config)

    // -------------------------------------------------------------------------
    // 4. Telomere end-family archive (or reuse precomputed end_archive)
    // -------------------------------------------------------------------------
    ch_archive_provided = ch_graph_samples
        .filter { _meta, _asm, _gfa, inputs, _bam -> inputs.end_archive }
        .map { meta, _asm, _gfa, inputs, _bam -> tuple(meta.id, file(inputs.end_archive, checkIfExists: true)) }

    ch_need_archive = ch_graph_samples
        .filter { _meta, _asm, _gfa, inputs, _bam -> !inputs.end_archive }

    ch_fish_input = ch_need_archive.map { meta, _asm, _gfa, inputs, bam ->
        def read_source = inputs.hifi_reads ? file(inputs.hifi_reads, checkIfExists: true) : bam
        tuple(meta, read_source, meta.telomere_repeat)
    }

    FISHTELOMERICREADS(
        ch_fish_input,
        end_cfg.bbduk_seed_copies,
        end_cfg.bbduk_max_k,
        end_cfg.bbduk_hdist,
        end_cfg.bam_exclude_flags
    )

    ch_ncrf_input = FISHTELOMERICREADS.out.reads.map { meta, reads ->
        tuple(meta, reads, meta.telomere_repeat)
    }
    NCRFTELOMERES(
        ch_ncrf_input,
        end_cfg.ncrf_scoring,
        end_cfg.ncrf_min_length,
        end_cfg.ncrf_max_noise
    )

    ch_cluster_input = FISHTELOMERICREADS.out.reads
        .map { meta, reads -> tuple(meta.id, meta, reads) }
        .join(NCRFTELOMERES.out.summary.map { meta, summary -> tuple(meta.id, summary) }, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .join(FISHTELOMERICREADS.out.headers.map { meta, headers -> tuple(meta.id, headers) }, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .join(FISHTELOMERICREADS.out.log.map { meta, log -> tuple(meta.id, log) }, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, reads, summary, headers, log -> tuple(meta, reads, summary, headers, log) }

    TELOMEREENDCLUSTERS_CLUSTER(
        ch_cluster_input,
        end_clusters_config,
        cluster_script
    )

    RACON_CONSENSUS(TELOMEREENDCLUSTERS_CLUSTER.out.racon_inputs)

    ch_asm_by_id = ch_need_archive.map { meta, asm, _gfa, _inputs, _bam -> tuple(meta.id, asm) }
    ch_summarize_input = FISHTELOMERICREADS.out.reads
        .map { meta, reads -> tuple(meta.id, meta, reads) }
        .join(NCRFTELOMERES.out.summary.map { meta, summary -> tuple(meta.id, summary) }, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_asm_by_id, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .join(RACON_CONSENSUS.out.consensus.map { meta, consensus -> tuple(meta.id, consensus) }, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .join(
            TELOMEREENDCLUSTERS_CLUSTER.out.cluster_state.map { meta, inward, end_rec, fate, nonterm, exposure, memb, sens, reps, assigned ->
                tuple(meta.id, inward, end_rec, fate, nonterm, exposure, memb, sens, reps, assigned)
            },
            by: 0, failOnMismatch: true, failOnDuplicate: true
        )
        .map { _id, meta, reads, summary, asm, consensus, inward, end_rec, fate, nonterm, exposure, memb, sens, reps, assigned ->
            tuple(meta, reads, summary, asm, consensus, inward, end_rec, fate, nonterm, exposure, memb, sens, reps, assigned)
        }

    TELOMEREENDCLUSTERS_SUMMARIZE(
        ch_summarize_input,
        end_clusters_config,
        cluster_script
    )

    ch_end_archive = ch_archive_provided.mix(
        TELOMEREENDCLUSTERS_SUMMARIZE.out.end_archive.map { meta, dir -> tuple(meta.id, dir) }
    )

    ch_resolved = ch_end_archive
        .join(ch_tiara, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_barrnap, by: 0, failOnMismatch: true, failOnDuplicate: true)

    emit:
    resolved    = ch_resolved // [ val(sample_id), path(end_archive), path(tiara), path(barrnap) ]
    tiara       = TIARA.out.classifications
    barrnap     = BARRNAP_GFA.out.gff
    end_archive = TELOMEREENDCLUSTERS_SUMMARIZE.out.end_archive
}
