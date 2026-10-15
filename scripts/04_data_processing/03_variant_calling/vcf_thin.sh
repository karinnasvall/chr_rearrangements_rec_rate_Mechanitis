#!/bin/bash
# script to thin vcf files and concatenate them



#load modules
module load bcftools/1.20--h8b25389_0  

#list with vcf files including path
LIST_VCF=$1
OUT_DIR=$2


mkdir $OUT_DIR
cd $OUT_DIR


for file in $(cat $LIST_VCF)
do
OUT_FILE=$(basename ${file%.*})
echo "Processing file: ${OUT_FILE}"

bcftools view -v snps ${file} | bcftools +prune -w 10000bp -n 1 -N rand -Oz  \
| bcftools view -O z --output ${OUT_FILE%.*}_thin_snps.vcf.gz
bcftools index -t ${OUT_FILE%.*}_thin_snps.vcf.gz
echo ">>> Thinning complete for ${OUT_FILE}"    
done

ls *_thin_snps.vcf.gz | sort -V  > ../list_vcf_sorted.txt
#the vcf files have to be in correct order
bcftools concat -f ../list_vcf_sorted.txt -Oz --output vcf_thin_concat_snps.vcf.gz
echo ">>> VCF thinning and concatenation complete!"

