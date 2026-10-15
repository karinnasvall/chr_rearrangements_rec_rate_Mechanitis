#!/usr/bin/env bash
set -euo pipefail

############################################
# Feature Density Pipeline
# Works for genes, repeats, or any GFF feature
############################################

usage() {
    echo ""
    echo "Usage:"
    echo "  $0 -g genome.fa -a annotation.gff3 -f feature_name -w window_size -o output_dir"
    echo ""
    echo "Example:"
    echo "  $0 -g genome.fa -a helixer.gff3 -f gene -w 100000 -o results/"
    echo "  $0 -g genome.fa -a repeats.gff3 -f repeat_region -w 50000 -o results/"
    echo ""
    exit 1
}

############################################
# Parse arguments
############################################

while getopts "g:a:f:w:o:" opt; do
  case ${opt} in
    g ) GENOME=$OPTARG ;;
    a ) GFF=$OPTARG ;;
    f ) FEATURE=$OPTARG ;;
    w ) WINDOW=$OPTARG ;;
    o ) OUT_DIR=$OPTARG ;;
    * ) usage ;;
  esac
done

if [[ -z "${GENOME:-}" || -z "${GFF:-}" || -z "${FEATURE:-}" || -z "${WINDOW:-}" || -z "${OUT_DIR:-}" ]]; then
    usage
fi

############################################
# Setup
############################################

FASTA=$(basename "$GENOME")
PREFIX="${FASTA%%.*}"

mkdir -p "$OUT_DIR"

echo "Genome:      $GENOME"
echo "Annotation:  $GFF"
echo "Feature:     $FEATURE"
echo "Window size: $WINDOW"
echo "Output dir:  $OUT_DIR"
echo ""

############################################
# Check required files
############################################

if [[ ! -f "$GENOME" ]]; then
    echo "ERROR: Genome file not found."
    exit 1
fi

if [[ ! -f "$GFF" ]]; then
    echo "ERROR: GFF file not found."
    exit 1
fi

############################################
# Ensure genome index (.fai) exists
############################################

FAI1="${GENOME}.fai"
FAI2="${GENOME%.*}.fai"

if [[ -f "$FAI1" ]]; then
    FAI="$FAI1"

elif [[ -f "$FAI2" ]]; then
    FAI="$FAI2"

else
    echo "Genome index (.fai) not found."
    echo "Attempting to create index with samtools faidx..."

    if command -v samtools &> /dev/null; then
        samtools faidx "$GENOME"

        # Check again after indexing
        if [[ -f "$FAI1" ]]; then
            FAI="$FAI1"
        elif [[ -f "$FAI2" ]]; then
            FAI="$FAI2"
        else
            echo "ERROR: samtools ran but .fai file was not created."
            exit 1
        fi

        echo "Successfully created index: $FAI"

    else
        echo "ERROR: samtools not found in PATH. Cannot generate .fai index."
        exit 1
    fi
fi




############################################
# 1 Create windows
############################################

WINDOW_BED="$OUT_DIR/${PREFIX}_${FEATURE}_${WINDOW}bp_windows.bed"

echo "Creating genomic windows..."
bedtools makewindows -g "$FAI" -w "$WINDOW" > "$WINDOW_BED"

############################################
# 2 Extract feature from GFF
############################################

FEATURE_GFF="$OUT_DIR/${PREFIX}_${FEATURE}.gff3"

echo "Extracting feature: $FEATURE"
awk -v f="$FEATURE" '$3 == f' "$GFF" > "$FEATURE_GFF"

############################################
# 3 Count features per window
############################################

COUNT_FILE="$OUT_DIR/${PREFIX}_${FEATURE}_density_count.tsv"

echo "Counting features per window..."
bedtools intersect -a "$WINDOW_BED" -b "$FEATURE_GFF" -c \
| bedtools sort \
> "$COUNT_FILE"

############################################
# 4 Merge overlapping features
############################################

MERGED_BED="$OUT_DIR/${PREFIX}_${FEATURE}_merged.bed"

echo "Merging overlapping features..."

bedtools sort -i "$FEATURE_GFF" \
| awk 'BEGIN{FS=OFS="\t"}{print $1,$4,$5,$9}' \
| bedtools merge -c 4 -o distinct \
| awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4}' \
> "$MERGED_BED"

############################################
# 5 Calculate occupancy (bp per window)
############################################

BP_FILE="$OUT_DIR/${PREFIX}_${FEATURE}_density_bp.tsv"

echo "Calculating occupancy per window..."

bedtools intersect -a "$WINDOW_BED" -b "$MERGED_BED" -wao \
| bedtools groupby -g 1,2,3 -c 8 -o sum \
| bedtools sort \
> "$BP_FILE"

############################################
# 6 Merge count + occupancy
############################################

FINAL_TSV="$OUT_DIR/${PREFIX}_${FEATURE}_${WINDOW}bp_density.tsv"

echo "Merging occupancy and count..."

awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {key=$1 FS $2 FS $3; data[key]=$4; next}
     {
        key=$1 FS $2 FS $3
        if (key in data) {
            print $0, data[key]
        } else {
            print $0, 0
        }
     }' "$BP_FILE" "$COUNT_FILE" \
| bedtools sort \
> "$FINAL_TSV"

############################################
# 7 Chromosome summary
############################################

echo "Generating chromosome summary..."

awk 'BEGIN{FS=OFS="\t"}{print $1, $2, $3, $3-$2}' "$MERGED_BED" \
| bedtools groupby -g 1 -c 4 -o sum \
> "$OUT_DIR/${PREFIX}_${FEATURE}_occupancy_chr.tsv"

bedtools groupby -i "$FEATURE_GFF" -g 1 -c 3 -o count \
> "$OUT_DIR/${PREFIX}_${FEATURE}_count_chr.tsv"

join -1 1 -2 1 \
  <(sort -k1,1 "$OUT_DIR/${PREFIX}_${FEATURE}_occupancy_chr.tsv") \
  <(sort -k1,1 "$OUT_DIR/${PREFIX}_${FEATURE}_count_chr.tsv") \
| join -1 1 -2 1 - \
  <(cut -f1,2 "${FAI}" | sort -k1,1) \
  -a1 \
> "$OUT_DIR/${PREFIX}_${FEATURE}_density_chr.tsv"

############################################
# 8 Cleanup
############################################

echo "Cleaning intermediate files..."
rm "$BP_FILE" "$COUNT_FILE" "$WINDOW_BED" "$FEATURE_GFF" "$OUT_DIR/${PREFIX}_${FEATURE}_occupancy_chr.tsv" "$OUT_DIR/${PREFIX}_${FEATURE}_count_chr.tsv"

echo ""
echo "Main output:"
echo "  $FINAL_TSV"
echo "  ${PREFIX}_${FEATURE}_density_chr.tsv"
echo ""
