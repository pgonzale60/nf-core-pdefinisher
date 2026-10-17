process EXTRACTGAPEVIDENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delfies:0.11.0--pyhdfd78af_0' :
        'quay.io/biocontainers/delfies:0.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(bam), path(bai), path(sniffles_vcf)

    output:
    tuple val(meta), path("*.gaps.bed"),          emit: gaps_bed
    tuple val(meta), path("*.gap_evidence.json"), emit: gap_evidence_json
    tuple val("${task.process}"), val('python'), val('3.10.4'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), val('0.23.3'), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Tool: pde_extract_gap_evidence.py | Script-SHA: 0876e6489ea84794
    python3 - << 'PYEOF'
import re, gzip

open_fn = gzip.open if "${assembly}".endswith('.gz') else open
with open_fn("${assembly}", 'rt') as f, open("${prefix}.gaps.bed", 'w') as out:
    c_name = None
    seq_chunks = []
    def dump_gaps(name, seq):
        for m in re.finditer(r'[Nn]+', seq):
            out.write(f"{name}\\t{m.start()}\\t{m.end()}\\n")
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if c_name:
                dump_gaps(c_name, "".join(seq_chunks))
            c_name = line[1:].split()[0]
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if c_name:
        dump_gaps(c_name, "".join(seq_chunks))
PYEOF

    pde_extract_gap_evidence.py \\
        --gaps "${prefix}.gaps.bed" \\
        --sniffles "${sniffles_vcf}" \\
        --bam "${bam}" \\
        --sample "${prefix}" \\
        ${args} > "${prefix}.gap_evidence.json"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.gaps.bed"
    touch "${prefix}.gap_evidence.json"
    """
}
