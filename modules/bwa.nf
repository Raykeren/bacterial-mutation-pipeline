#!/usr/bin/env nextflow

process BWA_MEM {
    tag "BWA alignment on ${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/bwa:0.7.17--h7132678_9'

    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy'

    input:
    tuple val(meta), path(reads1), path(reads2)
    path genome

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bwa index ${genome}
    bwa mem ${args} ${genome} ${reads1} ${reads2} | \
        samtools sort -@ ${task.cpus} -o ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -i 'Version' | sed 's/Version: //')
    END_VERSIONS
    """
}