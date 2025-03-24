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
params {
    input = "samples.csv"
    genome = "reference.fasta"
    outdir = "results"
    variant_callers = ['bcftools', 'freebayes', 'varscan2', 'gatk']  // Allow user to specify which callers to use
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

// Create channels
workflow {
    // Input channel from CSV
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.sample]
            tuple(meta, file(row.fastq_1), file(row.fastq_2))
        }
        .set { input_samples }

    // Reference genome channel
    genome_ch = Channel.fromPath(params.genome)
    
    // Quality control and preprocessing
    FASTQC(input_samples)
    FASTP(input_samples)
    
    // Alignment
    BWA_MEM(FASTP.out.reads, genome_ch)
    SAMTOOLS(BWA_MEM.out.bam)

    // Variant calling with multiple callers
    if ('bcftools' in params.variant_callers) {
        BCFTOOLS_CALL(SAMTOOLS.out.bam_sorted, genome_ch)
    }
    if ('freebayes' in params.variant_callers) {
        FREEBAYES(SAMTOOLS.out.bam_sorted, genome_ch)
    }
    if ('varscan2' in params.variant_callers) {
        VARSCAN(SAMTOOLS.out.bam_sorted, genome_ch)
    }
    if ('gatk' in params.variant_callers) {
        // Create dictionary and index for GATK
        GATK_HAPLOTYPE_CALLER(SAMTOOLS.out.bam_sorted, genome_ch)
    }

    // Collect and merge VCF outputs from all callers
    vcf_ch = Channel.empty()
    if ('bcftools' in params.variant_callers) {
        vcf_ch = vcf_ch.mix(BCFTOOLS_CALL.out.vcf)
    }
    if ('freebayes' in params.variant_callers) {
        vcf_ch = vcf_ch.mix(FREEBAYES.out.vcf)
    }
    if ('varscan2' in params.variant_callers) {
        vcf_ch = vcf_ch.mix(VARSCAN.out.vcf)
    }
    if ('gatk' in params.variant_callers) {
        vcf_ch = vcf_ch.mix(GATK_HAPLOTYPE_CALLER.out.vcf)
    }

    // Annotation
    ANNOTATE(vcf_ch, genome_ch)
}