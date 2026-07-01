#!/bin/bash
#BSUB -J mapping_DNA
#BSUB -o ant_DNA_mapping-%I-%J.output
#BSUB -e ant_DNA_mapping-%I-%J.error
#BSUB -n 1
#BSUB -M 2000
#BSUB -R "select[mem>2000] rusage[mem=2000] span[hosts=1]"
#BSUB -q oversubscribed

ml nextflow/25.10.0-10289 
module load singularityce-4.1.0/python-3.11.6 

nextflow run main_mapping.nf 
#nextflow run main_mapping.nf -resume
