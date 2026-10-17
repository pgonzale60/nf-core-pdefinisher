process FILTERCHROMOSOMES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_seqkit_gawk:5f1679612f236815'

    input:
    tuple val(meta), path(assembly)
    val chr_names // e.g. "I,II,III,IV,V,X" or list of chromosome names

    output:
    tuple val(meta), path("*.chromosomes.fasta"), emit: chromosomes
    tuple val(meta), path("*.non_chromosomal.fasta"), emit: non_chromosomal
    tuple val("${task.process}"), val('seqkit'), val('2.10.0'), emit: versions_filterchromosomes, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chr_str = meta.chromosomes ? meta.chromosomes : ((chr_names instanceof List) ? chr_names.join(",") : "${chr_names}")

    """
    # Create pattern file for exact anchored header matching
    echo "${chr_str}" | tr ',' '\\n' | grep . | awk '{print "^" \$0 "\$"}' > chr_patterns.txt

    # Extract canonical chromosomes
    seqkit grep -r -f chr_patterns.txt ${assembly} > ${prefix}.chromosomes.fasta

    # Extract non-chromosomal sequences (invert match)
    seqkit grep -v -r -f chr_patterns.txt ${assembly} > ${prefix}.non_chromosomal.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.chromosomes.fasta
    touch ${prefix}.non_chromosomal.fasta
    """
}
