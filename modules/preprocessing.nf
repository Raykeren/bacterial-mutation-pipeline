process FASTP {
    tag "Quality trimming ${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    publishDir "${params.outdir}/${meta.id}/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("*_trimmed_1.fastq.gz"), path("*_trimmed_2.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastp \
        -i ${reads1} \
        -I ${reads2} \
        -o ${prefix}_trimmed_1.fastq.gz \
        -O ${prefix}_trimmed_2.fastq.gz \
        --json ${prefix}.fastp.json \
        --html ${prefix}.fastp.html \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}