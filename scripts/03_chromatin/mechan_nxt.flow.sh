#!/bin/bash
#BSUB -J Hic_pro[7]
#BSUB -o /data/tol/users/sn19/lustre/HiC_project/log/script_hic-pro-%I-%J.output
#BSUB -e /data/tol/users/sn19/lustre/HiC_project/log/script_hic-pro-%I-%J.error
#BSUB -n 5
#BSUB -M 20000
#BSUB -R "select[mem>20000] rusage[mem=20000] span[hosts=1]"
#BSUB -q week

module load nextflow/24.10.4-5934 
module load singularityce-4.1.0/python-3.11.6 

# Paths
OUTDIR="/data/tol/users/sn19/lustre/HiC_project/results/high_cr/Mechanitis/"

# Species list file
SPECIES_LIST="/data/tol/users/sn19/lustre/HiC_project/data/high_cr/Mechanitis/sample_sheet/species_list.csv"

# Read the line based on array index
LINE=$(sed -n "${LSB_JOBINDEX}p" $SPECIES_LIST)

# Parse the CSV line into variables
IFS=',' read -r CODE NAME CSV_PATH FASTA_PATH <<< "$LINE"

# Check if variables were parsed correctly
if [ -z "$CODE" ] || [ -z "$CSV_PATH" ] || [ -z "$FASTA_PATH" ]; then
    echo "Error: Failed to parse line $LSB_JOBINDEX from species list"
    exit 1
fi

echo "Processing species: $CODE ($NAME)"
echo "Input CSV: $CSV_PATH"
echo "Reference FASTA: $FASTA_PATH"

#  Run Nextflow with parameters
# the path to the folder where you have downloaded nf-core/hic pipeline https://nf-co.re/hic/2.1.0/

nextflow run /lustre/scratch125/tol/teams/meier/projects/heliconius/08_chromatin_structure/hic/ \
  --input "${CSV_PATH}" \
  --fasta "${FASTA_PATH}" \
  --digestion 'arimaV2' \
  --hicpro_maps \
  --bin_size '1000,2000,5000,10000,20000,40000,50000,60000,80000,100000,150000,200000' \
  --res_compartments '40000,100000' \
  --res_tads '2000,5000' \
  --split_fastq \
  --save_raw_maps \
  --outdir "${OUTDIR}/${CODE}" \
  -w '/data/tol/users/sn19/lustre/HiC_project/tmp' \
  -profile sanger,singularity \
  -name 'HiC_fast_Mechan_9'
    
    echo "Finished $species"
    

echo " All runs completed successfully!"
