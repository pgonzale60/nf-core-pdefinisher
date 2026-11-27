#!/usr/bin/env nextflow

// Reference-first graph finishing pilot. The worker environment provides the graph tools.
params.curated = null
params.gfa = null
params.raw_gfa = null
params.end_archive = null
params.ont_reads = null
params.ont_paf = null
params.hifi_reads = null
params.tiara = null
params.barrnap = null
params.recipe = null
params.gfapdefinisher = 'gfapdefinisher'

process GRAPH_FINISH {
    label 'process_low'
    // An executable on PATH can change independently of the staged files.
    cache false

    input:
    path curated
    path gfa
    path end_archive, stageAs: 'end_archive'
    path raw_gfa
    path ont_reads
    path ont_paf
    path hifi_reads
    path tiara
    path barrnap
    path recipe
    path pilot_script
    path package_script
    val chromosomes
    val telomere_repeat
    val gfapdefinisher

    output:
    path 'graph_finish', emit: results

    script:
    def optional = [raw_gfa: raw_gfa, ont_reads: ont_reads, ont_paf: ont_paf, hifi_reads: hifi_reads,
                    tiara: tiara, barrnap: barrnap, recipe: recipe]
        .findAll { _key, value -> value }
        .collect { key, value -> "--${key.replace('_', '-')} \"${value}\"" }
        .join(' ')
    """
    python3 "${pilot_script}" \\
        --curated "${curated}" \\
        --gfa "${gfa}" \\
        --end-archive "${end_archive}" \\
        --chromosomes "${chromosomes}" \\
        --telomere "${telomere_repeat}" \\
        --gfapdefinisher "${gfapdefinisher}" \\
        --threads ${task.cpus} \\
        --outdir graph_finish \\
        ${optional}
    """

    stub:
    """
    mkdir graph_finish
    touch graph_finish/inputs.json
    """
}

workflow {
    main:
    if (!params.curated || !params.gfa || !params.end_archive || !params.chromosomes || !params.outdir) {
        error 'Supply --curated, --gfa, --end_archive, --chromosomes and --outdir'
    }
    GRAPH_FINISH(
        file(params.curated, checkIfExists: true),
        file(params.gfa, checkIfExists: true),
        file(params.end_archive, checkIfExists: true),
        params.raw_gfa ? file(params.raw_gfa, checkIfExists: true) : [],
        params.ont_reads ? file(params.ont_reads, checkIfExists: true) : [],
        params.ont_paf ? file(params.ont_paf, checkIfExists: true) : [],
        params.hifi_reads ? file(params.hifi_reads, checkIfExists: true) : [],
        params.tiara ? file(params.tiara, checkIfExists: true) : [],
        params.barrnap ? file(params.barrnap, checkIfExists: true) : [],
        params.recipe ? file(params.recipe, checkIfExists: true) : [],
        file("${projectDir}/scripts/graph_finish.py", checkIfExists: true),
        file("${projectDir}/scripts/package_reviewed_assembly.py", checkIfExists: true),
        params.chromosomes,
        params.telomere_repeat,
        params.gfapdefinisher,
    )

    publish:
    graph_finish = GRAPH_FINISH.out.results
}

output {
    graph_finish {
        path '.'
        mode params.publish_dir_mode
    }
}
