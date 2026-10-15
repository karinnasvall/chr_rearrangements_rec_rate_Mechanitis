#!/bin/bash
set -euo pipefail
# summarize bcftools stats outputs across multiple chromosomes, files named chr1.something...more.stats

IN_DIR=$1
OUT_DIR=$2

#check if input directory exists
if [ ! -d "$IN_DIR" ]; then
  echo "Input directory $IN_DIR does not exist. Exiting."
  exit 1
fi
# create output directory
mkdir -p $OUT_DIR

# Output files
echo -e "chr\tn_sites\tn_var" > $OUT_DIR/sites_per_chr.tsv
echo -e "chr\tti_tv" > $OUT_DIR/titv_per_chr.tsv
echo -e "chr\tmutation\tcount" > $OUT_DIR/mutation_spectrum.tsv
echo -e "sample\tchr\tnRefHom\tnNonRefHom\tnHets\tnIndels\tavg_depth\tnMissing\theterozygosity" > $OUT_DIR/per_sample_summary.tsv

for f in ${IN_DIR}/*.stats; do
  chr=$(basename "$f" | sed 's/\..*$//')
  echo "Processing $chr from $f"

  # number of sites
  n_sites=$(grep -P "^SN" "$f" | awk -F'\t' '$3=="number of records:"{print $4}')
  nvar=$(grep -P "^SN" "$f" | awk -F'\t' '$3=="number of SNPs:"{print $4}')

  # Ti/Tv
  titv=$(grep -P "^TSTV" "$f" | awk -F'\t' '{print $5}')

  # mutation spectrum
  grep -P "^ST" "$f" | awk -F'\t' -v chr="$chr" '{print chr "\t" $3 "\t" $4}'  >> $OUT_DIR/mutation_spectrum.tsv

  echo -e "${chr}\t${n_sites}\t${nvar}" >> $OUT_DIR/sites_per_chr.tsv
  echo -e "${chr}\t${titv}" >> $OUT_DIR/titv_per_chr.tsv
  
  # per-sample counts with heterozygosity calculation
  grep -P "^PSC" "$f" | awk -F'\t' -v chr="$chr" '{
    nRefHom=$4; nNonRefHom=$5; nHets=$6; nIndels=$9; avg_depth=$10; nMissing=$14
    total_called = nRefHom + nNonRefHom + nHets
    het = (total_called > 0) ? nHets / total_called : 0
    printf "%s\t%s\t%d\t%d\t%d\t%d\t%.1f\t%d\t%.6f\n", $3, chr, nRefHom, nNonRefHom, nHets, nIndels, avg_depth, nMissing, het
  }' >> $OUT_DIR/per_sample_summary.tsv

done

# Merge summary into one table
join -t $'\t' -1 1 -2 1 <(sort -V $OUT_DIR/sites_per_chr.tsv) <(sort -V $OUT_DIR/titv_per_chr.tsv) > $OUT_DIR/summary_per_chr.tsv
#move header first
grep "chr" $OUT_DIR/summary_per_chr.tsv | cat - <(grep -v "chr" $OUT_DIR/summary_per_chr.tsv) > $OUT_DIR/summary_per_chr_temp.tsv; mv $OUT_DIR/summary_per_chr_temp.tsv $OUT_DIR/summary_per_chr.tsv
# clean up intermediate files
rm $OUT_DIR/sites_per_chr.tsv $OUT_DIR/titv_per_chr.tsv
echo "Done!"
echo "Summary: $OUT_DIR/summary_per_chr.tsv"
echo "Mutation spectrum: $OUT_DIR/mutation_spectrum.tsv"  
echo "Per-sample summary: $OUT_DIR/per_sample_summary.tsv"