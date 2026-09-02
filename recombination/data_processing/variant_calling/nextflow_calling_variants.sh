#!/bin/bash
#BSUB -J variant_calling
#BSUB -o variant_calling-%I-%J.output
#BSUB -e variant_calling-%I-%J.error
#BSUB -n 1
#BSUB -M 2000
#BSUB -R "select[mem>2000] rusage[mem=2000] span[hosts=1]"
#BSUB -q oversubscribed


#get a list of bam files
#ls /data/tol/teams/meier/lustre/users/kn9/recombination/pop_gen/ref_lysimnia/mapping/bam_dedup/*.sorted.dedup.bam > list_bam.txt
# add the previoulsy mapped polymnia
#ls /data/tol/teams/meier/lustre/users/kn9/recombination/rawdata_to_vcf/02_bam_lysimnia/bam_dedup/*.sorted.dedup.bam | grep -f <(cut -f1 -d";" ../../ref_polymnia/list_lysimnia_incl.csv ) >> list_bam.txt


ml nextflow/25.10.0-10289 
module load singularityce-4.1.0/python-3.11.6 

nextflow run variant_calling.nf 
#nextflow run variant_calling.nf -resume
