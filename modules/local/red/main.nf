process RED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/red:2018.09.10--h9948957_3' :
        'quay.io/biocontainers/red:2018.09.10--h9948957_3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('red'), val('2018.09.10'), emit: versions_red, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    if [ ! -s "${fasta}" ]; then
        touch ${prefix}.bed
    else
        mkdir -p genome_dir out_rpt
        # Red requires input files with .fa extension in the directory
        if [[ "${fasta}" == *.gz ]]; then
            gzip -dc "${fasta}" > genome_dir/${prefix}.fa
        else
            cp -L "${fasta}" genome_dir/${prefix}.fa
        fi

        Red \\
            -gnm genome_dir \\
            -rpt out_rpt \\
            -frm 2 \\
            -cor ${task.cpus} \\
            ${args}

        # Red writes repeat coordinates in tabular format to out_rpt with .bed extension
        # Ensure column 1 contains strictly the sequence ID (strip header description/whitespace)
        if [ -s out_rpt/${prefix}.bed ]; then
            awk 'BEGIN {FS="\\t"; OFS="\\t"} {split(\$1, a, " "); \$1=a[1]; print}' out_rpt/${prefix}.bed > ${prefix}.bed
        elif [ -s out_rpt/${prefix}.rpt ]; then
            awk 'BEGIN {FS="\\t"; OFS="\\t"} {split(\$1, a, " "); \$1=a[1]; print}' out_rpt/${prefix}.rpt > ${prefix}.bed
        else
            touch ${prefix}.bed
        fi
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
