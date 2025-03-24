process VARIANT_CALL {
    tag "Variant calling for ${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/bcftools:1.18--h8b25389_0'

    publishDir "${params.outdir}/${meta.id}/variants", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    bcftools mpileup -f ${fasta} ${args} ${bam} | \
        bcftools call -mv -Oz -o ${prefix}.vcf.gz
    bcftools index ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}