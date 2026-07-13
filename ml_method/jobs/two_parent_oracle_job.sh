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

run () {  # $1=backbone $2=parents $3=root(none|hybrid) $4=build(graft|detach) $5=out-suffix
    OUT="output/two_parent_oracle/${CONFIG}_$5"
    echo ">>> [$CONFIG] backbone=$1 parents=$2 root=$3 build=$4 -> $OUT"
    python scripts/two_parent_oracle.py --mul-trees-dir "$MUL_TREES_DIR" \
        --config "$CONFIG" --backbone "$1" --parents "$2" --root-backbone "$3" \
        --build "$4" --n "$N" --out-dir "$OUT"
    python scripts/score_reconstructions.py --recon-dir "$OUT" --timeout 120
    echo "=================== SUMMARY [$CONFIG $5] ==================="
    python scripts/summarize_oracle.py --scores "$OUT/reconstruction_scores.csv" \
        --mapping "$OUT/mapping_scores.csv"
}

echo "============================================================"
echo "Two-parent oracle | config=$CONFIG n=$N"
echo "============================================================"

# graft = the nested-safe build (default); detach = old build, kept to show the fix.
run true   true    none   detach true_bb_detach          # old floor (~0.159)
run true   true    none   graft  true_bb_graft           # NEW floor: expect << 0.159
run astral true    hybrid graft  astral_ceiling_rooted   # ceiling, rooted + nested-safe build
run astral coclust hybrid graft  astral_coclust_rooted   # realistic: rooted + coclust + nested-safe

echo "=================== FAITHFULNESS LOCALIZATION [$CONFIG] (graft) ==================="
python scripts/localize_faithfulness.py \
    --scores "output/two_parent_oracle/${CONFIG}_true_bb_graft/reconstruction_scores.csv" \
    --mul-trees-dir "$MUL_TREES_DIR"

echo "Done. Read: true_bb_graft vs true_bb_detach(0.159); astral_ceiling_rooted vs Polyphest."
