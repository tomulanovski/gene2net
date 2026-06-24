#!/bin/bash
#SBATCH --job-name=g2n_reroot
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/reroot_%j.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/reroot_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Rooting test: reconstruct predicted events on ASTRAL RE-ROOTED to the true root,
# then score. Compare to the existing astral (~0.76) and true (~0.35) numbers. If
# this lands near 0.35, the ASTRAL backbone gap is just rooting (cheap fix).
#
# Submit: sbatch --export=ALL,CONFIG=ils_low,N=200 jobs/backbone_rerooted_job.sh

set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
BASE=/groups/itay_mayrose/tomulanovski/gene2net/ml_method
cd "$BASE"

MODEL_DIR="${MODEL_DIR:-output/reconstruct_allo}"
MUL_TREES_DIR="${MUL_TREES_DIR:-data/mul_trees_2k}"
CONFIG="${CONFIG:-ils_low}"
THRESHOLD="${THRESHOLD:-0.9}"
N="${N:-200}"
OUT="$MODEL_DIR/backbone_exp/${CONFIG}_rerooted"

echo ">>> Reconstruct on RE-ROOTED ASTRAL (final_project env)"
conda activate final_project
python scripts/reconstruct_mul_tree.py \
    --model-dir "$MODEL_DIR" --mul-trees-dir "$MUL_TREES_DIR" \
    --config "$CONFIG" --start 0 --n "$N" --thresholds "$THRESHOLD" \
    --backbone astral_rerooted --workers "${SLURM_CPUS_PER_TASK:-8}" \
    --out-base "$OUT"

echo ">>> Score (gene2net env)"
conda activate gene2net
python scripts/score_reconstructions.py --recon-dir "$OUT/t${THRESHOLD}" --timeout 120

echo "Done. Compare mean edit_distance_multree to astral (~0.76) and true (~0.35)."
