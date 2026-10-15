# summarize RepeatMasker output
#!/bin/sh


# Usage: ./summarize_repeat_per_chr.sh <input_file> <output_file>
# input shoudl be merged RepeatMasker output file, if overlapping repeats have been merged
# and are in comma separated in the 4th column. They will be summarized as "Overlapping" in the output.
INPUT=$1
OUTPUT=$2


awk -v OFS='\t' '{
  rep = $4
  if (rep ~ /,/) rep = "Overlapping"
  print $1, $2, $3, rep, $3 - $2
}' $INPUT | sort -k1,1 -k4,4 |bedtools groupby -g 1,4 -c 5 -o sum,count > $OUTPUT



