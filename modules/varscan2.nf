process VARSCAN {
    tag "Variant calling with VarScan2 for ${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/varscan:2.4.4--hdfd78af_1'

    publishDir "${params.outdir}/${meta.id}/variants/varscan2", mode: 'copy'

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
    # Generate mpileup file
    samtools mpileup \
        -f ${fasta} \
        ${bam} \
        > ${prefix}.mpileup

    # Run VarScan2 for SNP calling
    varscan mpileup2snp \
        ${prefix}.mpileup \
        --min-coverage 8 \
        --min-var-freq 0.01 \
        --min-avg-qual 20 \
        --p-value 0.05 \
        ${args} \
        > ${prefix}.vcf

    # Compress and index VCF
    bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2: \$(varscan 2>&1 | grep -i 'VarScan v' | sed 's/VarScan v//g')
        samtools: \$(samtools --version | grep ^samtools | sed 's/samtools //')
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
        varscan2: 2.4.4
        samtools: 1.18
        bgzip: 1.16
        tabix: 1.16
    END_VERSIONS
    """
}