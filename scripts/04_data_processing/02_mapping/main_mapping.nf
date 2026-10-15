#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * WGS mapping and deduplication pipeline
 * Process 1: BWA-MEM2 alignment and sorting
 * Process 2: Collate, fixmate, sort and index
 * Process 3: Mark duplicates
 */


// Parameters
params.reference = "/data/tol/teams/meier/lustre/users/kn9/recombination/rawdata_to_vcf/00_input/ref_genomes/ilMecPoly1.hap1.1.primary.fa.gz"
params.input_dir = "/data/tol/teams/meier/lustre/users/kn9/recombination/pedigree/01_rawdata_to_vcf/01_trim_qc/fastp_polymnia"
params.mapped_dir = "/data/tol/teams/meier/lustre/users/kn9/recombination/pedigree/01_rawdata_to_vcf/02_mapping_polymnia/bam"
params.dedup_dir = "/data/tol/teams/meier/lustre/users/kn9/recombination/pedigree/01_rawdata_to_vcf/02_mapping_polymnia/bam_dedup"


// Process 1: BWA-MEM2 alignment and sorting
process bwaMapping {
    publishDir "${params.mapped_dir}", mode: 'copy'
    
    module 'bwa-mem2/2.2.1--hd03093a_2'
    module 'samtools-1.19/python-3.12.0'
    
    tag "${sample_id}"
    
    input:
    val reference  // Changed from 'path' to 'val'
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: sorted_bam
    
    script:
    def r1 = reads[0]
    def r2 = reads[1]
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:Illumina"
    """
    bwa-mem2 mem -t 25 -M -R "${rg}" ${reference} ${r1} ${r2} | samtools sort -@ 10 -m 4G -o ${sample_id}.sorted.bam
    """
}

// Process 2: Collate, fixmate, sort and index
process collateSortFixmate {
    publishDir "${params.dedup_dir}", mode: 'copy', pattern: "${sample_id}.sorted.sorted.bam*"
    
    module 'samtools-1.19/python-3.12.0'
    
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.sorted.bam"), path("${sample_id}.sorted.sorted.bam.bai"), emit: fixed_bam
    
    script:
    """
    # Create tmp directory for this sample
    mkdir -p tmp_${sample_id}
    
    # Collate → fixmate → sort → index
    samtools collate --threads 5 -T tmp_${sample_id}/D_${sample_id}.collate -O -u ${bam} \
    | samtools fixmate --threads 5 -m -u - - \
    | samtools sort --threads 5 -O bam -o ${sample_id}.sorted.sorted.bam -
    
    samtools index --threads 5 ${sample_id}.sorted.sorted.bam
    
    # Clean up tmp directory
    rm -rf tmp_${sample_id}
    """
}

// Process 3: Mark duplicates
process markDuplicates {
    publishDir "${params.dedup_dir}", mode: 'copy'
    
    module 'samtools-1.19/python-3.12.0'
    
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.dedup.bam"), path("${sample_id}.sorted.dedup.bam.bai"), emit: dedup_bam
    
    script:
    """
    # Mark duplicates → index
    samtools markdup --threads 5 -r -s ${bam} ${sample_id}.sorted.dedup.bam
    
    samtools index --threads 5 ${sample_id}.sorted.dedup.bam
    """
}

// Workflow
workflow {
    // Prepare reference channel - pass as value (path string), not as staged file
    reference_ch = Channel.value(params.reference)
    
    // Prepare read pairs channel for files like: sample.fastp.r1.fq.gz and sample.fastp.r2.fq.gz
    reads_ch = Channel
        .fromFilePairs("${params.input_dir}/*_R{1,2}_fastp_trimmed.fq.gz", checkIfExists: true)
        .map { sample_id, files -> 
            // sample_id will be the part before .fastp.r1 or .fastp.r2
            tuple(sample_id, files)
        }

    // Run BWA-MEM2 mapping
    mapped_bams = bwaMapping(reference_ch, reads_ch)
    
    // Run collate, fixmate, sort and index
    fixed_bams = collateSortFixmate(mapped_bams.sorted_bam)
    
    // Run mark duplicates
    dedup_bams = markDuplicates(fixed_bams.fixed_bam)
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
