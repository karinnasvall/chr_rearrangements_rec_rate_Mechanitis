nextflow.enable.dsl=2

// ✔️ Input/output directory configuration
params.input_dir = "/data/tol/teams/meier/lustre/users/kn9/recombination/pedigree/01_rawdata_to_vcf/00_input/fastq_raw"
params.output_dir = "/data/tol/teams/meier/lustre/users/kn9/recombination/pedigree/01_rawdata_to_vcf/01_trim_qc"

workflow {

    // ✔️ Detect *_R1.fq.gz files and extract sample IDs
    fastp_ch = Channel
    .fromFilePairs("${params.input_dir}/*_R{1,2}.fq.gz")
    .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }

    // 🔁 Run the pipeline steps in order
    qc_reports_raw   = fastqc_raw(fastp_ch)		   // Step 0: fastqc raw reads
    trimmed_fastqs   = process_fastp(fastp_ch)            // Step 1: fastp
    qc_reports       = fastqc(trimmed_fastqs)               // Step 2: fastqc
}


// 📊 FASTQC: Generate per-sample quality reports
///////////////////////////////////////////////////////////////////////////////////////////
process fastqc_raw {

    tag "$sample_id"
    publishDir "${params.output_dir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(r1_gz), path(r2_gz)

    output:
    path("*_fastqc.zip")
    path("*_fastqc.html")

    script:
    """
    module load openjdk-17.0.8.1_1

    /software/team347/ev4/FastQC/fastqc -t 8 $r1_gz $r2_gz
    """
}

///////////////////////////////////////////////////////////////////////////////////////////
// 🔬 FASTP: Quality trimming
///////////////////////////////////////////////////////////////////////////////////////////
process process_fastp {

    tag "$sample_id"
    publishDir "${params.output_dir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1_fastp_trimmed.fq.gz"), path("${sample_id}_R2_fastp_trimmed.fq.gz")

    script:
    """
    module load fastp/0.23.4--hadf994f_3

    fastp \\
        --in1 ${params.input_dir}/${sample_id}_R1.fq.gz \\
        --in2 ${params.input_dir}/${sample_id}_R2.fq.gz \\
        --out1 ${sample_id}_R1_fastp_trimmed.fq.gz \\
        --out2 ${sample_id}_R2_fastp_trimmed.fq.gz \\
        --trim_poly_g --trim_poly_x \\
        --length_required 50 \\
        --thread 16 \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json
    """
}


///////////////////////////////////////////////////////////////////////////////////////////
// 📊 FASTQC: Generate per-sample quality reports
///////////////////////////////////////////////////////////////////////////////////////////
process fastqc {

    tag "$sample_id"
    publishDir "${params.output_dir}/fastqc_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(r1_gz), path(r2_gz)

    output:
    path("*_fastqc.zip")
    path("*_fastqc.html")

    script:
    """
    module load openjdk-17.0.8.1_1

    /software/team347/ev4/FastQC/fastqc -t 8 $r1_gz $r2_gz
    """
}
