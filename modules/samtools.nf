process SAMTOOLS {
    tag "Processing BAM files for ${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.18--h50ea8bc_1'

    publishDir "${params.outdir}/${meta.id}/bam", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam_sorted
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${bam}
    samtools index ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | grep ^samtools | sed 's/samtools //')
    END_VERSIONS
    """
}