#!/bin/bash
#SBATCH --job-name=partner_import
#SBATCH --partition=itaym-pool
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --output=partner_import_%j.out

# Re-run ONLY the partner-accuracy permutation importance on the shipped model. The previous
# run OOM-killed at 32G because it builds the E x E pairwise partner features; bumped to 128G.
# Detection importance already succeeded, so it is not repeated here.
# If this still OOMs, resubmit with --max-samples 500 (halves the held sample set).

set -o pipefail
CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate final_project
cd /groups/itay_mayrose/tomulanovski/gene2net/ml_method
export PYTHONPATH="$PWD"

DATA_DIR=""
for c in ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
  DATA_DIR="$DATA_DIR data/mul_trees_2k/training_rooted/$c"
done

echo "================ partner-accuracy permutation importance ($(date +%H:%M:%S)) ================"
python _archive/experiments/scripts/partner_permutation_importance.py \
    --data-dir $DATA_DIR --model-dir output/reconstruct_final --config configs/reconstruct_final.yaml \
    || echo "  !! partner permutation importance FAILED"

echo "================ DONE ($(date +%H:%M:%S)) ================"
