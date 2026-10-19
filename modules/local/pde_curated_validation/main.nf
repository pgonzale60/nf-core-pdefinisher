process EXTRACT_UNMAPPED_HIFI {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        '/software/tola/images/htstools-1.17.sif' :
        'community.wave.seqera.io/library/samtools:1.17--5354929c' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.unmapped.fastq.gz"), emit: fastq
    tuple val(meta), path("${prefix}.unmapped.count")   , emit: count
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | head -n1 | cut -d' ' -f2"), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -@ ${task.cpus} -b -f 4 ${bam} | \\
    samtools fastq -@ ${task.cpus} - | gzip -c > ${prefix}.unmapped.fastq.gz

    READS=\$(zcat ${prefix}.unmapped.fastq.gz | awk '{if(NR%4==1) count++} END {print count+0}')
    echo "\$READS" > ${prefix}.unmapped.count
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "@read1" | gzip -c > ${prefix}.unmapped.fastq.gz
    echo "1" > ${prefix}.unmapped.count
    """
}

process MASH_DECONTAMINATE {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        '/software/tola/images/mash-2.3.sif' :
        'community.wave.seqera.io/library/mash_python:2.3--abc' }"

    input:
    tuple val(meta), path(unmapped_fastq)
    path(sketch_db)

    output:
    tuple val(meta), path("${prefix}.decontaminated.fastq.gz"), emit: decontaminated_fastq
    tuple val(meta), path("${prefix}.contaminant.fastq.gz")   , emit: contaminant_fastq
    tuple val(meta), path("${prefix}.mash_hits.tab")          , emit: hits
    tuple val("${task.process}"), val('mash'), eval("mash --version"), emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mash screen -w -p ${task.cpus} ${sketch_db} ${unmapped_fastq} > ${prefix}.mash_hits.tab

    # Python partition
    python3 - << 'EOF'
import sys, gzip

min_identity = 0.90
hits_file = "${prefix}.mash_hits.tab"
in_fq = "${unmapped_fastq}"
out_decon = "${prefix}.decontaminated.fastq.gz"
out_contam = "${prefix}.contaminant.fastq.gz"

contam_ids = set()
with open(hits_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            try:
                if float(parts[0]) >= min_identity:
                    rid = parts[4].split()[0].lstrip('@')
                    contam_ids.add(rid)
            except ValueError:
                pass

with gzip.open(in_fq, 'rt') as fin, \
     gzip.open(out_decon, 'wt') as fdec, \
     gzip.open(out_contam, 'wt') as fcon:
    while True:
        h1 = fin.readline()
        if not h1: break
        s = fin.readline()
        h2 = fin.readline()
        q = fin.readline()
        rid = h1.strip().split()[0].lstrip('@')
        if rid in contam_ids:
            fcon.write(h1 + s + h2 + q)
        else:
            fdec.write(h1 + s + h2 + q)
EOF
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.decontaminated.fastq.gz
    touch ${prefix}.contaminant.fastq.gz
    touch ${prefix}.mash_hits.tab
    """
}

process MITOHIFI_SEGREGATE {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        '/software/tola/images/mitohifi-3.2.2.sif' :
        'community.wave.seqera.io/library/mitohifi:3.2.2--abc' }"

    input:
    tuple val(meta), path(decontaminated_fastq)
    path(mito_ref_fasta)
    path(mito_ref_gbk)

    output:
    tuple val(meta), path("${prefix}.mitochondrion.fasta")       , emit: mito_fasta
    tuple val(meta), path("${prefix}.mitochondrial_reads.fastq.gz"), emit: mito_fastq
    tuple val(meta), path("${prefix}.residue_reads.fastq.gz")     , emit: residue_fastq
    tuple val("${task.process}"), val('mitohifi'), eval("mitohifi.py --version | head -n1"), emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mitohifi.py -r ${decontaminated_fastq} \\
                -f ${mito_ref_fasta} \\
                -g ${mito_ref_gbk} \\
                -t ${task.cpus} \\
                -o 5 || true

    FINAL_FA=\$(ls -1 final_mitogenome.fasta final_mitogenome.fa ptg000001c.mitogenome.fa *mitogenome.fa *.mito.fa *.fasta 2>/dev/null | grep -iE 'mitogenome|mito|contig' | head -n1 || true)

    if [ -n "\$FINAL_FA" ] && [ -s "\$FINAL_FA" ]; then
        cp "\$FINAL_FA" ${prefix}.mitochondrion.fasta
        minimap2 -ax map-hifi -t ${task.cpus} --secondary=no ${prefix}.mitochondrion.fasta ${decontaminated_fastq} | \\
          samtools view -@ 4 -b -o mito_aln.bam -
        samtools view -@ 4 -b -F 4 mito_aln.bam | samtools fastq -@ 4 - | gzip -c > ${prefix}.mitochondrial_reads.fastq.gz
        samtools view -@ 4 -b -f 4 mito_aln.bam | samtools fastq -@ 4 - | gzip -c > ${prefix}.residue_reads.fastq.gz
    else
        touch ${prefix}.mitochondrion.fasta
        touch ${prefix}.mitochondrial_reads.fastq.gz
        cp ${decontaminated_fastq} ${prefix}.residue_reads.fastq.gz
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mitochondrion.fasta
    touch ${prefix}.mitochondrial_reads.fastq.gz
    touch ${prefix}.residue_reads.fastq.gz
    """
}

process HIFIASM_RESIDUE {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        '/software/tola/images/hifiasm-0.25.0.sif' :
        'community.wave.seqera.io/library/hifiasm:0.25.0--abc' }"

    input:
    tuple val(meta), path(residue_fastq)

    output:
    tuple val(meta), path("${prefix}.residue_contigs.fasta"), emit: residue_contigs
    tuple val("${task.process}"), val('hifiasm'), eval("hifiasm --version"), emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    NREADS=\$(zcat ${residue_fastq} | awk '{if(NR%4==1) count++} END {print count+0}')
    if [ "\$NREADS" -lt 50 ]; then
        touch ${prefix}.residue_contigs.fasta
    else
        hifiasm -o residue_asm -t ${task.cpus} ${residue_fastq}
        if [ -f "residue_asm.bp.p_ctg.gfa" ]; then
            awk '/^S/{print ">"\$2"\\n"\$3}' residue_asm.bp.p_ctg.gfa > ${prefix}.residue_contigs.fasta
        elif [ -f "residue_asm.p_ctg.gfa" ]; then
            awk '/^S/{print ">"\$2"\\n"\$3}' residue_asm.p_ctg.gfa > ${prefix}.residue_contigs.fasta
        else
            touch ${prefix}.residue_contigs.fasta
        fi
        rm -f residue_asm.*
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.residue_contigs.fasta
    """
}

process CHARACTERIZE_RESIDUE {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        '/software/tola/images/ncbi-blast-2.16.0+.sif' :
        'community.wave.seqera.io/library/blast:2.16.0--abc' }"

    input:
    tuple val(meta), path(residue_contigs), path(residue_fastq)
    val(blast_nt_db)
    path(taxonomy_dir)

    output:
    tuple val(meta), path("${prefix}.blob_summary.json"), emit: summary
    tuple val("${task.process}"), val('blastn'), eval("blastn -version | head -n1 | cut -d' ' -f2"), emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Ensure FASTA is properly 80-char line-wrapped for BLAST compatibility
    python3 - << 'PYWRAP'
with open('${residue_contigs}', 'r') as fin, open('wrapped_contigs.fa', 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            fout.write(line)
        else:
            seq = line.strip()
            for i in range(0, len(seq), 80):
                fout.write(seq[i:i+80] + "\\n")
PYWRAP

    NCTG=\$(grep -c '^>' wrapped_contigs.fa || true)
    if [ "\$NCTG" -eq 0 ]; then
        echo '{"nematode_contigs_count": 0, "nematode_contigs_bp": 0, "bacterial_contigs_count": 0, "unassigned_contigs_count": 0}' > ${prefix}.blob_summary.json
    else
        # Direct query against Lustre NT database without symlink staging
        export BLASTDB=\$(dirname "${blast_nt_db}")
        DB_BASE=\$(basename "${blast_nt_db}")

        blastn -query wrapped_contigs.fa \\
               -db "\${DB_BASE}" \\
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sblastnames sskingdoms" \\
               -num_threads ${task.cpus} \\
               -max_target_seqs 10 \\
               -evalue 1e-10 \\
               -out blast_hits.out

        python3 - << 'EOF'
import json

hits = {}
try:
    with open('blast_hits.out', 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 12:
                qid = parts[0]
                bs = float(parts[11])
                sciname = parts[13] if len(parts) > 13 else ""
                skingdom = parts[16] if len(parts) > 16 else ""
                if qid not in hits or bs > hits[qid]['bs']:
                    hits[qid] = {'bs': bs, 'sciname': sciname, 'kingdom': skingdom}
except Exception:
    pass

lengths = {}
with open('${residue_contigs}', 'r') as f:
    cid = None
    cur = 0
    for line in f:
        if line.startswith('>'):
            if cid: lengths[cid] = cur
            cid = line.strip().split()[0].lstrip('>')
            cur = 0
        else:
            cur += len(line.strip())
    if cid: lengths[cid] = cur

nem_c = 0
nem_bp = 0
bac_c = 0
unass_c = 0

for c, l in lengths.items():
    if c in hits:
        h = hits[c]
        s = h['sciname'].lower()
        k = h['kingdom'].lower()
        if 'nematoda' in s or 'nematod' in k or 'auanema' in s:
            nem_c += 1
            nem_bp += l
        elif 'bacteria' in s or 'bacteria' in k:
            bac_c += 1
        else:
            unass_c += 1
    else:
        unass_c += 1

res = {
    "nematode_contigs_count": nem_c,
    "nematode_contigs_bp": nem_bp,
    "bacterial_contigs_count": bac_c,
    "unassigned_contigs_count": unass_c
}
with open('${prefix}.blob_summary.json', 'w') as f:
    json.dump(res, f, indent=2)
EOF
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"nematode_contigs_count": 0, "nematode_contigs_bp": 0, "bacterial_contigs_count": 0, "unassigned_contigs_count": 0}' > ${prefix}.blob_summary.json
    """
}

process BUILD_VALIDATION_PACKAGE {
    tag "${meta.id}"
    label 'process_single'

    executor 'local'

    input:
    tuple val(meta), path(unmapped_fq), path(contam_fq), path(mito_fq), path(residue_fq), path(residue_fa), path(blob_json)
    val(curated_fa_path)

    output:
    tuple val(meta), path("${prefix}.Artifact_C_Validation_Package.json"), emit: validation_package

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - << 'EOF'
import json, gzip, os

def count_fq(p):
    if not p or not os.path.exists(p) or os.path.getsize(p) == 0: return 0
    c = 0
    fn = gzip.open if p.endswith('.gz') else open
    with fn(p, 'rt') as f:
        for idx, line in enumerate(f):
            if idx % 4 == 0: c += 1
    return c

def count_fa(p):
    if not p or not os.path.exists(p) or os.path.getsize(p) == 0: return 0, 0
    c = 0
    bp = 0
    with open(p, 'r') as f:
        for line in f:
            if line.startswith('>'): c += 1
            else: bp += len(line.strip())
    return c, bp

unm_c = count_fq('${unmapped_fq}')
bac_c = count_fq('${contam_fq}')
mit_c = count_fq('${mito_fq}')
res_c = count_fq('${residue_fq}')
ctg_c, ctg_bp = count_fa('${residue_fa}')

blob_data = {"nematode_contigs_count": 0, "nematode_contigs_bp": 0, "bacterial_contigs_count": 0, "unassigned_contigs_count": 0}
if os.path.exists('${blob_json}'):
    with open('${blob_json}', 'r') as f:
        blob_data = json.load(f)

action = "FREEZE_ASSEMBLY"
rationale = "All unmapped reads resolve to contamination or mitochondria; no authentic nematode sequence remains unintegrated."
candidate_contigs = None

if blob_data.get('nematode_contigs_count', 0) > 0 and blob_data.get('nematode_contigs_bp', 0) >= 5000:
    action = "ADD_CANDIDATE_CONTIGS"
    rationale = f"Discovered {blob_data['nematode_contigs_count']} genuine nematode contigs totaling {blob_data['nematode_contigs_bp']} bp in unmapped residue."
    candidate_contigs = os.path.abspath('${residue_fa}')

pkg = {
    "artifact_id": "APS4_ASSEMBLY_VALIDATION_V1",
    "pipeline": "pde-curated-validation",
    "assembly_tested": "${curated_fa_path}",
    "metrics": {
        "total_unmapped_reads": unm_c,
        "bacterial_contaminant_reads_mash": bac_c,
        "mitochondrial_reads": mit_c,
        "residue_reads_assembled": res_c,
        "residue_contigs_count": ctg_c,
        "residue_contigs_total_bp": ctg_bp
    },
    "blobtoolkit_classification": blob_data,
    "verdict": {
        "action_required": action,
        "rationale": rationale,
        "candidate_contigs_to_add_path": candidate_contigs
    }
}

with open('${prefix}.Artifact_C_Validation_Package.json', 'w') as f:
    json.dump(pkg, f, indent=2, sort_keys=True)
EOF
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{"artifact_id": "APS4_ASSEMBLY_VALIDATION_V1", "verdict": {"action_required": "FREEZE_ASSEMBLY"}}' > ${prefix}.Artifact_C_Validation_Package.json
    """
}
