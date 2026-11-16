#!/bin/bash
#BSUB -J cram2.fastq_mechanitis[2]
#BSUB -o /data/tol/users/sn19/lustre/HiC_project/log/cram2.fastq_Mechanitis-%I-%J.output
#BSUB -e /data/tol/users/sn19/lustre/HiC_project/log/cram2.fastq_Mechanitis-%I-%J.error
#BSUB -n 5
#BSUB -M 1000
#BSUB -R "select[mem>1000] rusage[mem=1000] span[hosts=1]"
#BSUB -q normal

module load samtools/1.20--h50ea8bc_0

GENOMES="/lustre/scratch125/tol/teams/meier/data/assemblies_Sanger/insects"
OUT_DIR="/data/tol/users/sn19/lustre/HiC_project/data/high_cr/Mechanitis/fastq"
LIST="/data/tol/users/sn19/lustre/HiC_project/data/high_cr/Mechanitis/tolid_list.txt"

module load speciesops/0.6.4
#  TOLID 
read -r tolid species <<< "$(sed -n "${LSB_JOBINDEX}p" "$LIST")"
if [ -z "$tolid" ]; then
    echo "ERROR: No TOLID found for index ${LSB_JOBINDEX}"
    exit 1
fi
echo "Processing $tolid"

# Speciesops to find cram files
GENOMES=$(speciesops getdir --tolid ${tolid} | awk '
/Species directory path:/ { if (NF > 3) print $NF;
    else { getline; print $1 }
    exit
}')

#  Hi-C folder
HIC_DIR="${GENOMES}/genomic_data/${tolid}/hic-arima2"
# CRAM file
CRAM_FILE=$(ls ${HIC_DIR}/*.cram 2>/dev/null | head -n 1)
if [ ! -f "$CRAM_FILE" ]; then
    echo "ERROR: No CRAM found in $HIC_DIR"
    exit 1
fi
# Output names
R1_OUT="${OUT_DIR}/${tolid}_R1.fastq.gz"
R2_OUT="${OUT_DIR}/${tolid}_R2.fastq.gz"

# Convert CRAM to FASTQ
samtools fastq -@ 10 \
  -1 "$R1_OUT" \
  -2 "$R2_OUT" \
  -0 /dev/null \
  -s /dev/null \
  -n "$CRAM_FILE"

