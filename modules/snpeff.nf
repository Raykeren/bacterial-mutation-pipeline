process ANNOTATE {
    tag "Annotating variants for ${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/snpeff:5.1--hdfd78af_2'

    publishDir "${params.outdir}/${meta.id}/annotated", mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    path(db)

    output:
    tuple val(meta), path("*.annotated.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    snpEff ann ${args} ${db} ${vcf} > ${prefix}.annotated.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(snpEff -version 2>&1 | grep -i 'version' | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}