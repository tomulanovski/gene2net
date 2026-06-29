#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1

# Reconstruction training: WGD detection + partner-edge prediction.

set -euo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate final_project

export PYTHONPATH="${PYTHONPATH:-}:${BASE_DIR}"
export PYTHONUNBUFFERED=1

echo "============================================================================"
echo "Reconstruction: WGD detection + partner-edge prediction"
echo "============================================================================"
echo "Data: ${DATA_DIR}"
echo "Config: ${CONFIG:-${BASE_DIR}/configs/reconstruct.yaml}"
echo "Output: ${OUTPUT_DIR:-${BASE_DIR}/output/reconstruct}"
echo "Device: $(python -c 'import torch; print("cuda" if torch.cuda.is_available() else "cpu")')"
echo "Date: $(date)"
echo "============================================================================"

echo "Clade labels: ${CLADE_LABELS:-<no>}"
echo "============================================================================"

python "${BASE_DIR}/scripts/train_reconstruct.py" \
    --data-dir ${DATA_DIR} \
    --config "${CONFIG:-${BASE_DIR}/configs/reconstruct.yaml}" \
    --output-dir "${OUTPUT_DIR:-${BASE_DIR}/output/reconstruct}" \
    ${INIT_FROM:+--init-from "$INIT_FROM"} \
    ${CLADE_LABELS:+--clade-labels}
