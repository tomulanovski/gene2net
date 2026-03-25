#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1

# Gene2Net-GNN training job. Called by submit_train.sh.

set -euo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate final_project

export PYTHONPATH="${PYTHONPATH:-}:${BASE_DIR}"
export PYTHONUNBUFFERED=1

echo "============================================================================"
echo "Gene2Net-GNN Training"
echo "============================================================================"
echo "Data: ${DATA_DIR}"
echo "Config: ${BASE_DIR}/configs/default.yaml"
echo "Output: ${BASE_DIR}/output"
echo "Device: $(python -c 'import torch; print("cuda" if torch.cuda.is_available() else "cpu")')"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Date: $(date)"
echo "============================================================================"

python "${BASE_DIR}/scripts/train.py" \
    --data-dir "${DATA_DIR}" \
    --config "${BASE_DIR}/configs/default.yaml" \
    --output-dir "${BASE_DIR}/output" \
    --device auto \
    --val-split 0.2 \
    --seed 42
