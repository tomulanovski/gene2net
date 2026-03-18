#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Generate one MUL-tree per array task (uses SimPhy for species tree).
# Called by submit_generate_mul_trees.sh — do not run directly.

set -euo pipefail

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net

export PYTHONPATH="${PYTHONPATH:-}:/groups/itay_mayrose/tomulanovski/gene2net/ml_method"

python "$SCRIPT" \
    --index "$SLURM_ARRAY_TASK_ID" \
    --output-dir "$OUTPUT_DIR" \
    --tree-height "$TREE_HEIGHT"
