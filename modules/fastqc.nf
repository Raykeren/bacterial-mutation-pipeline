process FASTQC {
    tag "FastQC on ${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    
    publishDir "${params.outdir}/${meta.id}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("*_fastqc.{zip,html}"), emit: fastqc_reports
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    fastqc \\
        --threads $task.cpus \\
        ${args} \\
        ${reads1} ${reads2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed -e "s/FastQC v//g")
    END_VERSIONS
    """

    stub:
    """
    touch ${reads1.baseName}_fastqc.html
    touch ${reads1.baseName}_fastqc.zip
    touch ${reads2.baseName}_fastqc.html
    touch ${reads2.baseName}_fastqc.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.12.1
    END_VERSIONS
    """
}