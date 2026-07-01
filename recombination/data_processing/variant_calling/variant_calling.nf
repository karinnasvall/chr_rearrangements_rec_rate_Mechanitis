#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reference = '/data/tol/teams/meier/lustre/users/kn9/recombination/rawdata_to_vcf/03_variant_call_polymnia_rptmask/ref/ilMecPoly1_1.renamed.masked.fa.gz'
params.bam_list = '/data/tol/teams/meier/lustre/users/kn9/recombination/pedigree/01_rawdata_to_vcf/02_mapping_polymnia/list_dedup_bam.txt'
params.chrom_list = '/data/tol/teams/meier/lustre/users/kn9/recombination/rawdata_to_vcf/03_variant_call_polymnia_rptmask/ref/list_SUPER_chr.txt'
params.outdir = '/lustre/scratch125/tol/teams/meier/users/kn9/recombination/pedigree/01_rawdata_to_vcf/03_variant_call_polymnia/output_vcf'
params.min_mq = 20
params.min_bq = 20
params.plink2_path = '/lustre/scratch125/tol/teams/meier/users/kn9/recombination/rawdata_to_vcf/03_variant_call_polymnia_rptmask/nf_dir/plink.yaml'

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

// Process 3: Filter variants - mark filters
process BCFTOOLS_FILTER_MARK {
    tag "$chrom"

    module 'bcftools/1.20--h8b25389_0'

    input:
    tuple val(chrom), path(vcf)

    output:
    tuple val(chrom), path("${chrom}.markfilter.vcf.gz"), emit: markfilter_vcf

    script:
    """
    bcftools filter \\
        --threads 5 \\
        --include "TYPE!=\\"indel\\" && STRLEN(REF) == 1 && QUAL >= 10 && MQ >= 20 && INFO/DP < 5000 && (FMT/GQ >= 10 & FMT/DP >= 6)" \\
        --SnpGap 5 \\
        --set-GTs . \\
        -O z \\
        -o ${chrom}.markfilter.vcf.gz \\
        ${vcf}
    """
}

// Process 4: Apply PASS filters
process BCFTOOLS_FILTER_PASS {
    tag "$chrom"

    module 'bcftools/1.20--h8b25389_0'

    input:
    tuple val(chrom), path(markfilter_vcf)

    output:
    tuple val(chrom), path("${chrom}.passfilter.vcf.gz"), emit: passfilter_vcf

    script:
    """
    bcftools view \\
        --threads 5 \\
        --apply-filters "PASS" \\
        -O z \\
        -o ${chrom}.passfilter.vcf.gz \\
        ${markfilter_vcf}
    """
}

// Process 5: Normalize variants
process BCFTOOLS_NORMALIZE {
    tag "$chrom"

    module 'bcftools/1.20--h8b25389_0'
    module 'htslib-1.19/perl-5.38.0'

    input:
    tuple val(chrom), path(passfilter_vcf)
    val reference

    output:
    tuple val(chrom), path("${chrom}.normfilter.vcf.gz"), path("${chrom}.normfilter.vcf.gz.tbi"), emit: normfilter_vcf

    script:
    """
    bcftools view \\
        --threads 5 \\
	    -M2 -m1 \\
        -O z \\
        -o ${chrom}.normfilter.vcf.gz \\
        ${passfilter_vcf}

    tabix -f ${chrom}.normfilter.vcf.gz
    """
}

// Process 6: PLINK2 simplification
process PLINK2_SIMPLIFY {
    tag "$chrom"
    conda "${params.plink2_path}"

    input:
    tuple val(chrom), path(normfilter_vcf), path(normfilter_tbi)

    output:
    tuple val(chrom), path("${chrom}.simplify.vcf"), emit: simplify_vcf

    script:
    """
    plink2 \\
        --double-id \\
        --threads 10 \\
        --vcf ${normfilter_vcf} \\
        --geno 0.9 \\
        --allow-extra-chr \\
        --recode vcf-iid \\
        --out ${chrom}.simplify
    """
}

// Process 7: Compress and index final VCF
process COMPRESS_INDEX_FINAL {
    tag "$chrom"
    publishDir "${params.outdir}/filtered", mode: 'copy'

    module 'samtools/1.20--h50ea8bc_0'
    module 'bcftools/1.20--h8b25389_0'
    module 'htslib-1.19/perl-5.38.0'

    input:
    tuple val(chrom), path(simplify_vcf)

    output:
    tuple val(chrom), path("${chrom}.simplify.vcf.gz"), path("${chrom}.simplify.vcf.gz.tbi"), emit: final_vcf

    script:
    """
    bgzip -@ 5 -f ${simplify_vcf}
    tabix -f -p vcf ${chrom}.simplify.vcf.gz
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
    
    // Run filtering pipeline
    markfilter_results = BCFTOOLS_FILTER_MARK(call_results.vcf)
    passfilter_results = BCFTOOLS_FILTER_PASS(markfilter_results.markfilter_vcf)
    normfilter_results = BCFTOOLS_NORMALIZE(passfilter_results.passfilter_vcf, reference_file)
    simplify_results = PLINK2_SIMPLIFY(normfilter_results.normfilter_vcf)
    COMPRESS_INDEX_FINAL(simplify_results.simplify_vcf)
}
