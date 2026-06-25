#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Package one training example. Called by submit_package_training.sh.

set -euo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate final_project

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
export PYTHONPATH="${PYTHONPATH:-}:${BASE_DIR}"

N_PER_JOB="${N_PER_JOB:-1}"
TREE_INDEX=$(( SLURM_ARRAY_TASK_ID * N_PER_JOB + ${INDEX_OFFSET:-0} ))

# Optional: re-root the species tree (rooted-tree experiment) and/or a custom
# output dir so rooted data lands separately from the unrooted data.
EXTRA_ARGS=()
[ "${ROOT_SPECIES_TREE:-0}" = "1" ] && EXTRA_ARGS+=(--root-species-tree)
[ -n "${OUTPUT_DIR:-}" ] && EXTRA_ARGS+=(--output-dir "$OUTPUT_DIR")

python "${BASE_DIR}/scripts/package_training_data.py" \
    --index "$TREE_INDEX" \
    --n-batch "$N_PER_JOB" \
    --mul-trees-dir "$MUL_TREES_DIR" \
    --config "$CONFIG_NAME" \
    "${EXTRA_ARGS[@]}"
