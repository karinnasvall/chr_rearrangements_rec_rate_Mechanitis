#!/bin/bash

# Script that evaluates the log likelihood difference between different repetitions of LepMap3 OrderMarkers 

# Retrieve the likelihood from each lg and print the one is best for fitting the step function in the next step (06_fitstepfunction.bsub)


OUT_DIR=$1
OUT_TABLE=$OUT_DIR/likelihoods_table.txt
OUT_BEST=$OUT_DIR/best_likelihood_files.txt


# Collect all: rep  LG  likelihood  path
for f in $OUT_DIR/rep*/ordered.*; do
  rep=$(basename $(dirname $f))
  lg=$(basename $f | sed 's/ordered\.\([0-9]*\).*/\1/')
  lik=$(grep "likelihood" "$f" | awk '{print $NF}')
  [ -n "$lik" ] && printf "%s\t%s\t%s\t%s\n" "$rep" "$lg" "$lik" "$f"
done | sort -k2V > $OUT_TABLE


# Output 2: path to best-likelihood ordered file per LG (for FitStepFunction)
awk -F'\t' '
  { if (!($2 in best) || $3+0 > best[$2]+0) { best[$2]=$3; path[$2]=$4 } }
  END { for (lg in path) print path[lg] }
' $OUT_TABLE | sort -V > "$OUT_BEST"

echo "Table written to:      $OUT_TABLE"
echo "Best rep paths written: $OUT_BEST"
