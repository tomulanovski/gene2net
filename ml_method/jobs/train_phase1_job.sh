#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1

# Phase 1 training: features-only GNN for binary WGD detection.
# Much faster than full model — no gene tree encoding needed.

set -euo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate final_project

export PYTHONPATH="${PYTHONPATH:-}:${BASE_DIR}"
export PYTHONUNBUFFERED=1

echo "============================================================================"
echo "Phase 1: Features-only GNN Training"
echo "============================================================================"
echo "Data: ${DATA_DIR}"
echo "Config: ${BASE_DIR}/configs/phase1.yaml"
echo "Output: ${BASE_DIR}/output/phase1"
echo "Device: $(python -c 'import torch; print("cuda" if torch.cuda.is_available() else "cpu")')"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Date: $(date)"
echo "============================================================================"

python "${BASE_DIR}/scripts/train_phase1.py" \
    --data-dir "${DATA_DIR}" \
    --config "${BASE_DIR}/configs/phase1.yaml" \
    --output-dir "${BASE_DIR}/output/phase1"
