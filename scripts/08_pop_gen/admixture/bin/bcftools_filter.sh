#!/bin/bash
set -euo pipefail

VCF=$(sed -n -e ${LSB_JOBINDEX}p list_vcfs.txt)

# filters
MIN=3
MAX=50
n_ind=110
MIN_AF=0.03
#3/110

MIN_FILTER=330
MAX_FILTER=5000

VCF_PATH=../ref_polymnia/variant_calling/output_vcf/vcfs/
OUT=input_vcf/${VCF}_filtered_DP${MIN_FILTER}_${MAX_FILTER}.vcf.gz

ml bcftools/1.20--h8b25389_0
bcftools view -S samples.txt $VCF_PATH/${VCF}.vcf.gz \
    | bcftools view -m2 -M2 -v snps \
        -i "QUAL>=30 && INFO/DP>=${MIN_FILTER} && INFO/DP<=${MAX_FILTER}"  \
    | bcftools +setGT -Oz -o $OUT -- -t q -n . -i "FMT/DP<${MIN}"
        

bcftools index $OUT
