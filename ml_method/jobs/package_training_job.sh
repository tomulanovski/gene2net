#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Package one training example. Called by submit_package_training.sh.

set -euo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
export PYTHONPATH="${PYTHONPATH:-}:${BASE_DIR}"

python "${BASE_DIR}/scripts/package_training_data.py" \
    --index "$SLURM_ARRAY_TASK_ID" \
    --mul-trees-dir "$MUL_TREES_DIR" \
    --config "$CONFIG_NAME"
