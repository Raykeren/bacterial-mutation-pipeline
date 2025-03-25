#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/preprocessing'
include { BWA_MEM } from './modules/bwa'
include { SAMTOOLS } from './modules/samtools'
include { VARIANT_CALL as BCFTOOLS_CALL } from './modules/bcftools'
include { FREEBAYES } from './modules/freeBayes'
include { VARSCAN } from './modules/varscan2'
include { GATK_HAPLOTYPE_CALLER } from './modules/gatkHaploCaller'
include { ANNOTATE } from './modules/snpeff'

// Pipeline parameter defaults
workflow {
    // Pipeline parameter defaults
    params {
        def input = "samples.csv"
        def genome = "reference.fasta"
        def outdir = "results"
        def variant_callers = ['bcftools', 'freebayes', 'varscan2', 'gatk']
    }
    
    // Log pipeline info
    log.info """
             BACTERIAL MUTATION PIPELINE
             ==========================
             input           : ${params.input}
             genome         : ${params.genome}
             outdir         : ${params.outdir}
             variant_callers: ${params.variant_callers.join(', ')}
             """
             
    // Create input channel from CSV file
    def input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.sample]
            tuple(meta, file(row.fastq_1), file(row.fastq_2))
        }
    
    // Create reference genome channel
    def genome_ch = Channel.fromPath(params.genome, checkIfExists: true)

    // Quality control
    FASTQC(input_ch)
    
    // Preprocessing
    FASTP(input_ch)
    
    // Alignment
    BWA_MEM(FASTP.out.reads, genome_ch)
    SAMTOOLS(BWA_MEM.out.bam)
    
    // Create output directory
    params.outdir.mkdirs()

    // Variant calling
    def vcf_ch = Channel.empty()
    
    if ('bcftools' in params.variant_callers) {
        BCFTOOLS_CALL(SAMTOOLS.out.bam_sorted, genome_ch)
        vcf_ch = vcf_ch.mix(BCFTOOLS_CALL.out.vcf)
    }
    
    if ('freebayes' in params.variant_callers) {
        FREEBAYES(SAMTOOLS.out.bam_sorted, genome_ch, params.outdir)
        vcf_ch = vcf_ch.mix(FREEBAYES.out.vcf)
    }
    
    if ('varscan2' in params.variant_callers) {
        VARSCAN(SAMTOOLS.out.bam_sorted, genome_ch, params.outdir)
        vcf_ch = vcf_ch.mix(VARSCAN.out.vcf)
    }
    
    if ('gatk' in params.variant_callers) {
        GATK_HAPLOTYPE_CALLER(SAMTOOLS.out.bam_sorted, genome_ch, params.outdir, params.variant_callers)
        vcf_ch = vcf_ch.mix(GATK_HAPLOTYPE_CALLER.out.vcf)
    }

    // Annotation
    ANNOTATE(vcf_ch, genome_ch)
}