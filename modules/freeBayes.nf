process FREEBAYES {
    tag "Variant calling with FreeBayes for ${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/freebayes:1.3.6--h1870644_3'

    publishDir "${params.outdir}/${meta.id}/variants/freebayes", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    # Run FreeBayes
    freebayes \\
        -f ${fasta} \\
        ${args} \\
        ${bam} \\
        > ${prefix}.vcf

    # Compress and index VCF
    bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(freebayes --version | grep -o 'version .*' | cut -f2 -d' ')
        bgzip: \$(bgzip --version | head -n1 | cut -f2 -d' ')
        tabix: \$(tabix --version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: 1.3.6
        bgzip: 1.16
        tabix: 1.16
    END_VERSIONS
    """
}