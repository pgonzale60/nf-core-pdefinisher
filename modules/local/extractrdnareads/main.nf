process EXTRACTRDNAREADS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/samtools_seqkit_gawk:5f1679612f236815'

    input:
    tuple val(meta), path(gff)
    tuple val(meta2), path(bam), path(bai)
    val flank_bp

    output:
    tuple val(meta), path("*.rdna_reads.fastq.gz"), emit: reads
    tuple val(meta), path("*.rdna_regions.bed")   , emit: bed
    tuple val("${task.process}"), val('samtools'), val('1.20'), emit: versions_extractrdnareads, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def flank  = flank_bp ?: 50000

    """
    # 1. Extract 18S and 28S rRNA feature coordinates from Barrnap GFF
    awk -F'\\t' '\$3 == "rRNA" && (\$9 ~ /18S/ || \$9 ~ /28S/) {
        print \$1 FS \$4-1 FS \$5
    }' ${gff} | sort -k1,1 -k2,2n > raw_features.bed

    if [ -s raw_features.bed ]; then
        # Merge overlapping/adjacent rRNA annotations within 10 kb
        awk 'BEGIN{OFS="\\t"}{
            if (NR==1){c=\$1; s=\$2; e=\$3}
            else {
                if (\$1==c && \$2 <= e+10000) { if (\$3>e) e=\$3 }
                else { print c, s, e; c=\$1; s=\$2; e=\$3 }
            }
        }END{if (NR>0) print c, s, e}' raw_features.bed > merged_raw.bed

        # Filter out singleton barrnap hits whose merged core span is below the policy
        # threshold (rdna_min_core_bp). A single stray 28S/18S unit (~1-2 kb) padded
        # by flank_bp would otherwise masquerade as a terminal rDNA array and trigger
        # the tripartite cassette pathway in curationsynthesis. Genuine multi-copy
        # tandem arrays span tens of kilobases even before flanking and pass this gate.
        MIN_CORE_BP=\$(python3 -c "
import json
try:
    p = json.load(open('${projectDir}/assets/curation_policy.json'))
    print(p.get('repeat_and_rdna', {}).get('rdna_min_core_bp', 10000))
except Exception:
    print(10000)
" 2>/dev/null || echo 10000)
        awk -v min_core="\${MIN_CORE_BP}" 'BEGIN{OFS="\\t"}(\$3-\$2) >= min_core' merged_raw.bed > merged_arrays.bed

        # Expand with flank_bp to capture chromosome anchoring regions
        awk -v flank=${flank} 'BEGIN{OFS="\\t"}{
            start = (\$2 - flank > 0) ? \$2 - flank : 0;
            end = \$3 + flank;
            print \$1, start, end;
        }' merged_arrays.bed > ${prefix}.rdna_regions.bed
    else
        touch ${prefix}.rdna_regions.bed
    fi

    # 2. Extract reads mapped to the rDNA arrays and flanking anchors
    if [ -s ${prefix}.rdna_regions.bed ]; then
        samtools view \\
            -@ ${task.cpus} \\
            -b \\
            -F 2308 \\
            -L ${prefix}.rdna_regions.bed \\
            ${bam} \\
            | samtools fastq -@ ${task.cpus} - \\
            | gzip -c > ${prefix}.rdna_reads.fastq.gz
    else
        echo "" | gzip -c > ${prefix}.rdna_reads.fastq.gz
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rdna_regions.bed
    echo "" | gzip -c > ${prefix}.rdna_reads.fastq.gz
    """
}
