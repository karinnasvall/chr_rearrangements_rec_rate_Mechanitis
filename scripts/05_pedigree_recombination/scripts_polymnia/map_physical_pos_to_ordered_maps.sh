#!/bin/sh

# script for mapping ordered files to physical position

OUTPUT=$1
INPUT_SNPS=$2

for file in $(ls $OUTPUT/*/ordered.*)
do
echo "Mapping file ${file} to physical position"

awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' $INPUT_SNPS $file > ${file}_mapped

done
#because of first line of snps.txt, we use NR-1 instead of NRll 05  
