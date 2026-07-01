#!/bin/bash
#BSUB -J calling
#BSUB -o ant_calling-%I-%J.output
#BSUB -e ant_calling-%I-%J.error
#BSUB -n 1
#BSUB -M 5000
#BSUB -R "select[mem>5000] rusage[mem=5000] span[hosts=1]"
#BSUB -q long

ml nextflow/25.10.0-10289 
module load singularityce-4.1.0/python-3.11.6 

nextflow run variant_calling.nf -resume
#nextflow run main.nf -resume
