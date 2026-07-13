#!/bin/bash
#SBATCH --job-name=g2n_2parent
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/2parent_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/2parent_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner
#SBATCH --array=0-5

# Phase 0 two-parent construction oracle. Per config, run three builds and score:
#   1. true parents on the TRUE backbone  -> faithfulness (expect edit ~0)
#   2. true parents on the ASTRAL backbone -> construction CEILING
#   3. co-cluster top-2 parents on ASTRAL  -> what the fix buys with today's signal
# Reads the three numbers vs Polyphest to pick Approach A / B / reframe.
#
# Submit with:
#   sbatch jobs/two_parent_oracle_job.sh

set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net

BASE=/groups/itay_mayrose/tomulanovski/gene2net/ml_method
cd "$BASE"

MUL_TREES_DIR="${MUL_TREES_DIR:-data/mul_trees_2k}"
N="${N:-200}"

CONFIGS=(ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M)
CONFIG=${CONFIGS[$SLURM_ARRAY_TASK_ID]}

run () {  # $1=backbone $2=parents-mode $3=root(none|hybrid) $4=out-suffix
    OUT="output/two_parent_oracle/${CONFIG}_$4"
    echo ">>> [$CONFIG] backbone=$1 parents=$2 root=$3 -> $OUT"
    python scripts/two_parent_oracle.py --mul-trees-dir "$MUL_TREES_DIR" \
        --config "$CONFIG" --backbone "$1" --parents "$2" --root-backbone "$3" \
        --n "$N" --out-dir "$OUT"
    python scripts/score_reconstructions.py --recon-dir "$OUT" --timeout 120
    echo "=================== SUMMARY [$CONFIG $4] ==================="
    python scripts/summarize_oracle.py --scores "$OUT/reconstruction_scores.csv" \
        --mapping "$OUT/mapping_scores.csv"
}

echo "============================================================"
echo "Two-parent oracle | config=$CONFIG n=$N"
echo "============================================================"

run true   true    none   true_bb                 # faithfulness: expect edit ~0
run astral true    none   astral_ceiling          # ceiling, UNROOTED astral (has rooting confound)
run astral true    hybrid astral_ceiling_rooted   # ceiling on the ROOTED backbone (pipeline-fair)
run astral coclust hybrid astral_coclust_rooted   # realistic + rooted = closest to a real pipeline

echo "=================== FAITHFULNESS LOCALIZATION [$CONFIG] ==================="
python scripts/localize_faithfulness.py \
    --scores "output/two_parent_oracle/${CONFIG}_true_bb/reconstruction_scores.csv" \
    --mul-trees-dir "$MUL_TREES_DIR"

echo "Done. Read: astral_ceiling_rooted vs astral_coclust_rooted vs current(0.68) vs Polyphest."
