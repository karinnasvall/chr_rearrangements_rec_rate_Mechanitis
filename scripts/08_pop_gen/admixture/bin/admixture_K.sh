#!/bin/bash
set -euo pipefail
ADMIX_PATH=/software/team347/kn9/admixture_linux-1.4.0/

K=${LSB_JOBINDEX}
SEED=42

$ADMIX_PATH/admixture --cv -j4 -s $SEED output/pruned.bed ${K} > output_adm/log_K${K}.out
