#!/bin/bash
#SBATCH --job-name=g2n_bbone
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/bbone_%j.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/bbone_%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# True-backbone experiment: run the model's PREDICTED events on (a) the ASTRAL
# backbone (the real pipeline) and (b) the true diploid backbone
# (species_tree_NNNN.nex), on the SAME samples + threshold, then score both.
# This isolates the ASTRAL-backbone cost given our real predictions, closing the
# 2x2: {true,astral} backbone x {true,predicted} events.
#
# Submit with:
#   sbatch --export=ALL,CONFIG=ils_low,THRESHOLD=0.9,N=200,\
# MUL_TREES_DIR=data/mul_trees_2k,MODEL_DIR=output/reconstruct_allo \
#     jobs/backbone_oracle_job.sh

set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"

BASE=/groups/itay_mayrose/tomulanovski/gene2net/ml_method
cd "$BASE"

MODEL_DIR="${MODEL_DIR:-output/reconstruct_allo}"
MUL_TREES_DIR="${MUL_TREES_DIR:-data/mul_trees_2k}"
CONFIG="${CONFIG:-ils_low}"
THRESHOLD="${THRESHOLD:-0.9}"
START="${START:-0}"
N="${N:-200}"
REPLICATE="${REPLICATE:-1}"
NWORKERS="${SLURM_CPUS_PER_TASK:-8}"

OUT_ASTRAL="$MODEL_DIR/backbone_exp/${CONFIG}_astral"
OUT_TRUE="$MODEL_DIR/backbone_exp/${CONFIG}_true"

echo "============================================================"
echo "Backbone experiment | config=$CONFIG threshold=$THRESHOLD n=$N"
echo "model=$MODEL_DIR  mul_trees=$MUL_TREES_DIR"
echo "============================================================"

echo ">>> Reconstruct on ASTRAL + TRUE backbones (final_project env)"
conda activate final_project
for BB in astral true; do
    if [ "$BB" = "astral" ]; then OUT="$OUT_ASTRAL"; else OUT="$OUT_TRUE"; fi
    echo "--- backbone=$BB -> $OUT ---"
    python scripts/reconstruct_mul_tree.py \
        --model-dir "$MODEL_DIR" --mul-trees-dir "$MUL_TREES_DIR" \
        --config "$CONFIG" --replicate "$REPLICATE" \
        --start "$START" --n "$N" --thresholds "$THRESHOLD" \
        --backbone "$BB" --workers "$NWORKERS" --out-base "$OUT"
done

echo ">>> Score both (gene2net env)"
conda activate gene2net
for BB in astral true; do
    if [ "$BB" = "astral" ]; then OUT="$OUT_ASTRAL"; else OUT="$OUT_TRUE"; fi
    echo "=================== SCORE backbone=$BB ==================="
    python scripts/score_reconstructions.py \
        --recon-dir "$OUT/t${THRESHOLD}" --timeout 120
done

echo "Done. Compare mean edit_distance: ASTRAL (real pipeline) vs TRUE (ceiling)."
