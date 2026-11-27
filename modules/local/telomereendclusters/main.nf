process TELOMEREENDCLUSTERS_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    container 'community.wave.seqera.io/library/samtools_minimap2_vsearch_python:1aa370c14225a820'

    input:
    tuple val(meta), path(reads), path(ncrf_summary), path(headers), path(fishing_log)
    path config
    path cluster_script

    output:
    tuple val(meta), path("prefixes.fasta"), path("assigned.paf"), path("representatives.fasta"), emit: racon_inputs
    tuple val(meta), path("inward_full.fasta.gz"), path("end_records.tsv.gz"), path("read_fate.tsv.gz"),
        path("nonterminal_hits.tsv.gz"), path("movie_exposure.tsv"), path("membership.tsv.gz"),
        path("identity_sensitivity.tsv"), path("representatives.fasta"), path("assigned.paf"), emit: cluster_state
    tuple val("${task.process}"), val('minimap2'), eval("minimap2 --version"), topic: versions, emit: versions_minimap2

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    read -r preset secondary overlap_ratio < <(python3 -c 'import json,sys;c=json.load(open(sys.argv[1]));print(c["minimap2_preset"],c["secondary_alignments"],c["overlap_min_ratio"])' "${config}")

    python3 "${cluster_script}" extract \\
        --config "${config}" \\
        --out . \\
        --reads "${reads}" \\
        --ncrf "${ncrf_summary}" \\
        --headers "${headers}" \\
        --fishing-log "${fishing_log}"

    if [[ -s prefixes.fasta ]]; then
        minimap2 -x "\$preset" -c --cs=short --for-only --secondary=yes -N "\$secondary" -p "\$overlap_ratio" -t ${task.cpus} \\
            prefixes.fasta prefixes.fasta | gzip -c > overlaps.paf.gz
    else
        : | gzip -c > overlaps.paf.gz
    fi

    python3 "${cluster_script}" cluster \\
        --config "${config}" \\
        --out . \\
        --overlaps overlaps.paf.gz

    if [[ -s representatives.fasta && -s prefixes.fasta ]]; then
        minimap2 -x "\$preset" -c --cs=short --for-only --secondary=yes -N "\$secondary" -p "\$overlap_ratio" -t ${task.cpus} \\
            representatives.fasta prefixes.fasta | gzip -c > members.paf.gz
    else
        : | gzip -c > members.paf.gz
    fi

    python3 "${cluster_script}" assigned \\
        --config "${config}" \\
        --out . \\
        --alignments members.paf.gz
    """

    stub:
    """
    : > prefixes.fasta
    : > assigned.paf
    : > representatives.fasta
    : | gzip -c > inward_full.fasta.gz
    printf 'end_id\\tread_id\\tmovie\\tread_end\\tread_length\\ttract_start\\ttract_end\\ttract_span\\ttract_parts\\tend_gap\\tstrict_gap\\tflank_length\\tprefix_length\\tstatus\\n' | gzip -c > end_records.tsv.gz
    printf 'read_id\\tmovie\\tread_length\\tncrf_hits\\tterminal_ends\\tstatus\\n' | gzip -c > read_fate.tsv.gz
    printf 'read_id\\tstart\\tend\\tstrand\\tread_length\\tmatch_ratio\\tstatus\\n' | gzip -c > nonterminal_hits.tsv.gz
    printf 'movie\\tinput_reads\\n' > movie_exposure.tsv
    printf 'end_id\\tcluster_id\\tseed_end_id\\tprefix_length\\toverlap_degree\\n' | gzip -c > membership.tsv.gz
    printf 'identity\\tsequence_families\\n' > identity_sensitivity.tsv
    """
}

process TELOMEREENDCLUSTERS_SUMMARIZE {
    tag "$meta.id"
    label 'process_low'

    container 'community.wave.seqera.io/library/samtools_minimap2_vsearch_python:1aa370c14225a820'

    input:
    tuple val(meta), path(reads), path(ncrf_summary), path(assembly), path(consensus),
        path(inward_full), path(end_records), path(read_fate), path(nonterminal_hits),
        path(movie_exposure), path(membership), path(identity_sensitivity),
        path(representatives), path(assigned_paf)
    path config
    path cluster_script

    output:
    tuple val(meta), path("end_archive"), emit: end_archive
    tuple val("${task.process}"), val('minimap2'), eval("minimap2 --version"), topic: versions, emit: versions_minimap2

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    read -r preset secondary genome_ratio < <(python3 -c 'import json,sys;c=json.load(open(sys.argv[1]));print(c["minimap2_preset"],c["secondary_alignments"],c["genome_min_ratio"])' "${config}")

    mkdir -p end_archive
    cp -L "${consensus}" end_archive/consensus.fasta
    cp -L "${representatives}" end_archive/representatives.fasta
    cp -L "${end_records}" end_archive/end_records.tsv.gz
    cp -L "${read_fate}" end_archive/read_fate.tsv.gz
    cp -L "${nonterminal_hits}" end_archive/nonterminal_hits.tsv.gz
    cp -L "${movie_exposure}" end_archive/movie_exposure.tsv
    cp -L "${membership}" end_archive/membership.tsv.gz
    cp -L "${identity_sensitivity}" end_archive/identity_sensitivity.tsv
    cp -L "${assigned_paf}" end_archive/assigned.paf

    if ( set +o pipefail; gzip -cd "${inward_full}" | grep -q '^>' ); then
        minimap2 -x "\$preset" -c --cs=short --secondary=yes -N "\$secondary" -p "\$genome_ratio" -t ${task.cpus} \\
            "${assembly}" "${inward_full}" | gzip -c > genome.paf.gz
    else
        : | gzip -c > genome.paf.gz
    fi

    python3 "${cluster_script}" summarize \\
        --config "${config}" \\
        --out end_archive \\
        --assigned end_archive/assigned.paf \\
        --genome genome.paf.gz \\
        --reference "${assembly}" \\
        --reads "${reads}" \\
        --ncrf "${ncrf_summary}"

    rm -f end_archive/assigned.paf
    """

    stub:
    """
    mkdir -p end_archive
    : > end_archive/consensus.fasta
    printf 'cluster_id\\tend_observations\\tunique_reads\\n' > end_archive/families.tsv
    printf 'cluster_id\\tpattern_id\\treads\\talleles\\n' | gzip -c > end_archive/haplotype_patterns.tsv.gz
    printf 'end_id\\tcluster_id\\tseed_end_id\\tprefix_length\\toverlap_degree\\n' | gzip -c > end_archive/membership.tsv.gz
    printf 'end_id\\tread_id\\tmovie\\tread_end\\tread_length\\ttract_start\\ttract_end\\ttract_span\\ttract_parts\\tend_gap\\tstrict_gap\\tflank_length\\tprefix_length\\tstatus\\n' | gzip -c > end_archive/end_records.tsv.gz
    printf 'cluster_id\\treference\\tstrand\\tlocus_group\\tmedian_coordinate\\tnative_sides\\tall_candidate_reads\\tunique_mapping_reads\\tambiguity_status\\n' > end_archive/locus_support.tsv
    cp -L "${config}" end_archive/MANIFEST.yaml
    """
}
