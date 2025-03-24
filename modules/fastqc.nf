process FASTQC {
    tag "FastQC on ${sample}"
    label 'process_low'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    
    publishDir "${params.outdir}/${sample}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads1), path(reads2)  // Changed sample to meta for consistency

    output:
    tuple val(meta), path("*_fastqc.{zip,html}"), emit: fastqc_reports
    path "versions.yml", emit: versions           // Added versions output

    script:
    """
    fastqc ${reads1} ${reads2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed -e "s/FastQC v//g")
    END_VERSIONS
    """
}