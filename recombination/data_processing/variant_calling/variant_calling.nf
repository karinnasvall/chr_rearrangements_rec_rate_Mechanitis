#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reference = '/data/tol/teams/meier/lustre/users/kn9/recombination/rawdata_to_vcf/03_variant_call_lysimnia_rptmask/ref/ilMecLysi212_1.renamed.masked.fa.gz'
params.bam_list = '/data/tol/teams/meier/lustre/users/kn9/recombination/pop_gen/ref_lysimnia/variant_calling/list_bam.txt'
params.chrom_list = '/data/tol/teams/meier/lustre/users/kn9/recombination/rawdata_to_vcf/03_variant_call_lysimnia_rptmask/ref/list_SUPER_chr.txt'
params.outdir = '/lustre/scratch125/tol/teams/meier/users/kn9/recombination/pop_gen/ref_lysimnia/variant_calling/output_vcf'
params.min_mq = 20
params.min_bq = 20

// Process 1: mpileup by chromosome
process BCFTOOLS_MPILEUP {
    tag "$chrom"
    
    publishDir "${params.outdir}/mpileup", mode: 'copy'
    
    module 'bcftools/1.20--h8b25389_0'

    input:
    val chrom
    path reference
    path bam_list
    
    output:
    tuple val(chrom), path("${chrom}.mpileup.bcf"), emit: mpileup_bcf
    
    script:
    """
    bcftools mpileup \\
        --threads 5 \\
        -f ${reference} \\
        -b ${bam_list} \\
        -r ${chrom} \\
        --min-MQ ${params.min_mq} \\
        --min-BQ ${params.min_bq} \\
        -O b \\
        -a FORMAT/DP,FORMAT/AD \\
        --skip-indels \\
        -o ${chrom}.mpileup.bcf
    """
}

// Process 2: variant calling by chromosome
process BCFTOOLS_CALL {
    tag "$chrom"
    
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    
    module 'bcftools/1.20--h8b25389_0'

    input:
    tuple val(chrom), path(mpileup_bcf)
    
    output:
    tuple val(chrom), path("${chrom}.call.MQ_BQ20.vcf.gz"), emit: vcf
    
    script:
    """
    bcftools call \\
        --threads 5 \\
        --annotate GQ,GP \\
        -mO z \\
        -o ${chrom}.call.MQ_BQ20.vcf.gz \\
        ${mpileup_bcf}
    """
}


// Workflow
workflow {
    // Create channel from chromosome list file
    chromosomes_ch = Channel
        .fromPath(params.chrom_list)
        .splitText()
        .map { it.trim() }

    // Stage reference and BAM list
    reference_file = file(params.reference)
    bam_list_file = file(params.bam_list)

    // Run mpileup
    mpileup_results = BCFTOOLS_MPILEUP(
        chromosomes_ch,
        reference_file,
        bam_list_file
    )

    // Run variant calling
    call_results = BCFTOOLS_CALL(mpileup_results.mpileup_bcf)
    
}

workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
}