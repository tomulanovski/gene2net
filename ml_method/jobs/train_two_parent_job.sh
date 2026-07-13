#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1

# Two-parent graft phaser training: WGD detection + unordered two-parent placement.
# Warm-start from the current rooted one-parter; the reshaped partner head (1->2
# slots) reinitializes via filter_compatible_state_dict.
#
# Prereq: relabel the rooted training data with two-parent labels first:
#   for c in ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
#     python scripts/relabel_from_metadata.py --two-parent --training-subdir training_rooted/$c
#   done
#
# Submit with (BASE_DIR, DATA_DIR = the 6 rooted config dirs):
#   sbatch --export=ALL,BASE_DIR=$PWD,\
# DATA_DIR="data/mul_trees_2k/training_rooted/ils_low ... (all 6)",\
# CONFIG=$PWD/configs/reconstruct_two_parent.yaml,\
# OUTPUT_DIR=$PWD/output/reconstruct_two_parent,\
# INIT_FROM=$PWD/output/reconstruct_cladelabels_rooted/best_model.pt,\
# CLADE_LABELS=1 \
#     jobs/train_two_parent_job.sh

set -euo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate final_project

export PYTHONPATH="${PYTHONPATH:-}:${BASE_DIR}"
export PYTHONUNBUFFERED=1

echo "============================================================================"
echo "Two-parent graft phaser: detection + unordered two-parent placement"
echo "Data:   ${DATA_DIR}"
echo "Config: ${CONFIG:-${BASE_DIR}/configs/reconstruct_two_parent.yaml}"
echo "Output: ${OUTPUT_DIR:-${BASE_DIR}/output/reconstruct_two_parent}"
echo "Warm-start: ${INIT_FROM:-<none>}"
echo "Device: $(python -c 'import torch; print("cuda" if torch.cuda.is_available() else "cpu")')"
echo "Date: $(date)"
echo "============================================================================"

python "${BASE_DIR}/scripts/train_reconstruct.py" \
    --data-dir ${DATA_DIR} \
    --config "${CONFIG:-${BASE_DIR}/configs/reconstruct_two_parent.yaml}" \
    --output-dir "${OUTPUT_DIR:-${BASE_DIR}/output/reconstruct_two_parent}" \
    ${INIT_FROM:+--init-from "$INIT_FROM"} \
    ${CLADE_LABELS:+--clade-labels}
