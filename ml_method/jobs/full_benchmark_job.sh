#!/bin/bash
#SBATCH --job-name=g2n_fullbench
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/fullbench_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/fullbench_%A_%a.err
#SBATCH --array=0-11
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Benchmark one config per array task (parallel). reconstruct (final_project)
# then score (gene2net). Uses a capped strategy so networks are well-formed and
# actually scoreable (raw over-counts -> malformed -> graph_edit_distance hangs).

set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
BASE=/groups/itay_mayrose/tomulanovski/gene2net/ml_method
cd "$BASE"

MODEL="${MODEL_DIR:-output/reconstruct_allo}"
THRESHOLD="${THRESHOLD:-0.9}"
STRATEGY="${STRATEGY:-cap}"

CONFIGS=(conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
conf_dup_loss_low_10M conf_dup_loss_low_10M_ne1M conf_dup_loss_low_10M_ne2M \
conf_dup_loss_medium_10M conf_dup_loss_medium_10M_ne1M conf_dup_loss_medium_10M_ne2M \
conf_dup_loss_high_10M conf_dup_loss_high_10M_ne1M conf_dup_loss_high_10M_ne2M)

C=${CONFIGS[$SLURM_ARRAY_TASK_ID]}
echo "=== config=$C | model=$MODEL threshold=$THRESHOLD strategy=$STRATEGY ==="

echo ">>> RECONSTRUCT ($C)"
conda activate final_project
python scripts/benchmark_networks.py --model-dir "$MODEL" --config "$C" \
    --threshold "$THRESHOLD" --strategies "$STRATEGY" --out-base "$MODEL/full/$C"

echo ">>> SCORE ($C)"
conda activate gene2net
python scripts/score_reconstructions.py \
    --recon-dir "$MODEL/full/$C/$STRATEGY" --timeout 120

echo "=== DONE $C ==="
