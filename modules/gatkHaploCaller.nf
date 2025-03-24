process GATK_HAPLOTYPE_CALLER {
    tag "HaplotypeCaller on ${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    publishDir "${params.outdir}/${meta.id}/variants/gatk", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(dict)  // Sequence dictionary required by GATK

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def avail_mem = 6
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O ${prefix}.vcf.gz \
        --sample-ploidy 1 \
        --dont-use-soft-clipped-bases true \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}