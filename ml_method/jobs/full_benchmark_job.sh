#!/bin/bash
#SBATCH --job-name=g2n_fullbench
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/fullbench_%j.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/fullbench_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Full benchmark of the model on the 21 networks across all 12 configs:
# reconstruct (final_project env) then score (gene2net env). One clean table.

set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
BASE=/groups/itay_mayrose/tomulanovski/gene2net/ml_method
cd "$BASE"

MODEL="${MODEL_DIR:-output/reconstruct_allo}"
THRESHOLD="${THRESHOLD:-0.9}"
STRATEGY="${STRATEGY:-raw}"

CONFIGS="conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
conf_dup_loss_low_10M conf_dup_loss_low_10M_ne1M conf_dup_loss_low_10M_ne2M \
conf_dup_loss_medium_10M conf_dup_loss_medium_10M_ne1M conf_dup_loss_medium_10M_ne2M \
conf_dup_loss_high_10M conf_dup_loss_high_10M_ne1M conf_dup_loss_high_10M_ne2M"

echo "=== RECONSTRUCT (final_project) | model=$MODEL threshold=$THRESHOLD strategy=$STRATEGY ==="
conda activate final_project
for C in $CONFIGS; do
    echo "#### BENCHMARK $C ####"
    python scripts/benchmark_networks.py --model-dir "$MODEL" --config "$C" \
        --threshold "$THRESHOLD" --strategies "$STRATEGY" \
        --out-base "$MODEL/full/$C" || echo "  benchmark failed for $C"
done

echo "=== SCORE (gene2net) ==="
conda activate gene2net
for C in $CONFIGS; do
    echo "#### SCORE $C ####"
    python scripts/score_reconstructions.py \
        --recon-dir "$MODEL/full/$C/$STRATEGY" --timeout 90 || echo "  score failed for $C"
done

echo "=== DONE ==="
