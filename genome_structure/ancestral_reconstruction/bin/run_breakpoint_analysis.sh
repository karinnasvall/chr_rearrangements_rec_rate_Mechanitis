#!/bin/bash
set -euo pipefail

# Usage: bash run_breakpoint_analysis.sh <ANC_NODE> <Q1> <Q2>
# Example: bash run_breakpoint_analysis.sh N5 ilMecPoly1 ilMecLysi212
# Check node names in species tree output from Orthofinder
# python3 ../bin/plot_tree.py ../results/03_orthofinder/species_tree_node_labels.txt tree_agora.pdf

ANC=${1:?Usage: $0 <ANC_NODE> <Q1> <Q2>}
Q1=${2:?Usage: $0 <ANC_NODE> <Q1> <Q2>}
Q2=${3:?Usage: $0 <ANC_NODE> <Q1> <Q2>}

OUTDIR="${ANC}_${Q1}_${Q2}"
AGORA_DIR="../results/05_agora/output/agora_output/ancGenomes/basic-workflow"
BED_DIR="../../anc_reconstruction_pipeline_260516/analysis_dore_phyl/bed"
BLOCKS_DIR="../../breakpoint_extraction_260321/results/03_conserved_blocks"
BIN_DIR="../bin"

mkdir -p "${OUTDIR}"

echo "[1/3] Detecting breakpoints from AGORA output..."
python3 "${BIN_DIR}/detect_breakpoints_agora.py" \
    -i "${AGORA_DIR}/ancGenome.${ANC}.list.bz2" \
    --bed-dir "${BED_DIR}/" \
    --q1 "${Q1}" \
    --q2 "${Q2}" \
    -o "${OUTDIR}/gene_alignment.tsv" \
    2> "${OUTDIR}/gene_alignment.log"

echo "[2/3] Plotting breakpoints..."
python3 "${BIN_DIR}/plot_breakpoints_agora.py" \
    -i "${OUTDIR}/gene_alignment.tsv" \
    -o "${OUTDIR}/chr_painting_breakpoints.pdf"

echo "[3/3] Refining breakpoints with cactus blocks..."
python3 "${BIN_DIR}/refine_breakpoints.py" \
    --input-bed-q1 "${OUTDIR}/gene_alignment.${Q1}_breakpoints.bed" \
    --input-bed-q2 "${OUTDIR}/gene_alignment.${Q2}_breakpoints.bed" \
    --blocks-q1 "${BLOCKS_DIR}/conserved_blocks_by_query1.tsv" \
    --blocks-q2 "${BLOCKS_DIR}/conserved_blocks_by_query2.tsv" \
    --bed-q1 "${OUTDIR}/refined_${Q1}.bed" \
    --bed-q2 "${OUTDIR}/refined_${Q2}.bed" \
    2> "${OUTDIR}/refined_breakpoints.log"

echo "Done. Results in ${OUTDIR}/"
