#!/bin/bash
#SBATCH --job-name=g2n_perm
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/perm_%j.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/perm_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Permutation importance for edge AND node features on the validation set.
# Ranks features by their effect on detection F1 (no retrain). Used to justify
# which features matter and, importantly, to test the node-feature design choice.
#
# Submit with:
#   sbatch --export=ALL,MODEL_DIR=output/phase1_feat9,\
# MUL_TREES_DIR=data/mul_trees_2k,FEATURES=both jobs/perm_importance_job.sh

set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"

BASE=/groups/itay_mayrose/tomulanovski/gene2net/ml_method
cd "$BASE"

MODEL_DIR="${MODEL_DIR:-output/phase1_feat9_full}"
MUL_TREES_DIR="${MUL_TREES_DIR:-data/mul_trees_2k}"
FEATURES="${FEATURES:-both}"
THRESHOLD="${THRESHOLD:-0.88}"
MAX_SAMPLES="${MAX_SAMPLES:-1000}"

# All training-config dirs (the val split is taken inside the script, seeded).
if [ -z "${DATA_DIRS:-}" ]; then
    DATA_DIRS=$(ls -d "$MUL_TREES_DIR"/training/*/ 2>/dev/null | tr '\n' ' ')
fi
echo "Data dirs: $DATA_DIRS"

conda activate final_project
python scripts/permutation_importance.py \
    --data-dir $DATA_DIRS \
    --model-dir "$MODEL_DIR" \
    --threshold "$THRESHOLD" \
    --max-samples "$MAX_SAMPLES" \
    --features "$FEATURES"

echo "Done."
