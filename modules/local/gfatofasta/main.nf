process GFATOFASTA {
    tag "$meta.id"
    label 'process_single'

    container 'community.wave.seqera.io/library/samtools_minimap2_vsearch_python:1aa370c14225a820'

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("${prefix}.unitigs.fasta"), emit: fasta
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${gfa}" "${prefix}.unitigs.fasta" <<'PY'
import gzip
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
opener = gzip.open if src.suffix == ".gz" else open
count = 0
with opener(src, "rt") as handle, dst.open("w") as out:
    for line in handle:
        if line.startswith("S\\t"):
            parts = line.rstrip("\\n").split("\\t")
            name, seq = parts[1], parts[2]
            if seq == "*":
                raise SystemExit(f"Graph segment {name} in {src} has no sequence ('*')")
            out.write(f">{name}\\n{seq}\\n")
            count += 1
if count == 0:
    raise SystemExit(f"No S-line sequences found in {src}")
PY
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>utg000001l\\nACGT\\n' > "${prefix}.unitigs.fasta"
    """
}
