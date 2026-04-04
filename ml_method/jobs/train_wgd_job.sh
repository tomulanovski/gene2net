#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1

# WGD Detector training: GIN gene tree encoder + GAT species tree.
# Needs more memory and time than Phase 1 due to gene tree processing.

set -euo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate final_project

export PYTHONPATH="${PYTHONPATH:-}:${BASE_DIR}"
export PYTHONUNBUFFERED=1

echo "============================================================================"
echo "WGD Detector: GIN Gene Tree Encoder + GAT Species Tree"
echo "============================================================================"
echo "Data: ${DATA_DIR}"
echo "Config: ${CONFIG:-${BASE_DIR}/configs/wgd_detector.yaml}"
echo "Output: ${OUTPUT_DIR:-${BASE_DIR}/output/wgd_detector}"
echo "Device: $(python -c 'import torch; print("cuda" if torch.cuda.is_available() else "cpu")')"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Date: $(date)"
echo "============================================================================"

python "${BASE_DIR}/scripts/train_wgd.py" \
    --data-dir "${DATA_DIR}" \
    --config "${CONFIG:-${BASE_DIR}/configs/wgd_detector.yaml}" \
    --output-dir "${OUTPUT_DIR:-${BASE_DIR}/output/wgd_detector}"
