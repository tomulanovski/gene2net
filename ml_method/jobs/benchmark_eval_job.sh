#!/bin/bash
#SBATCH --job-name=gene2net_bench
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/bench_%j.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/bench_%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Benchmark the model on the 21 networks across thresholds (reconstruct + score),
# all on one compute node. Two conda envs: final_project for inference, gene2net
# for scoring.
#
# Submit with:
#   sbatch --export=ALL,CONFIG=conf_ils_low_10M,THRESHOLDS="0.9 0.95 0.97 0.99",\
# MODEL_DIR=output/reconstruct_allo jobs/benchmark_eval_job.sh

set -eo pipefail
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"

BASE=/groups/itay_mayrose/tomulanovski/gene2net/ml_method
cd "$BASE"

MODEL_DIR="${MODEL_DIR:-output/reconstruct_allo}"
CONFIG="${CONFIG:-conf_ils_low_10M}"
THRESHOLDS="${THRESHOLDS:-0.9 0.95 0.97 0.99}"
REPLICATE="${REPLICATE:-1}"
NWORKERS="${SLURM_CPUS_PER_TASK:-8}"

echo "============================================================"
echo "Benchmark eval | config=$CONFIG | thresholds=$THRESHOLDS"
echo "model=$MODEL_DIR | replicate=$REPLICATE"
echo "============================================================"

echo ">>> Reconstruction (final_project env)"
conda activate final_project
for T in $THRESHOLDS; do
    echo "--- reconstruct threshold $T ---"
    python scripts/benchmark_networks.py \
        --model-dir "$MODEL_DIR" --config "$CONFIG" --replicate "$REPLICATE" \
        --threshold "$T" --out-dir "$MODEL_DIR/benchmark/${CONFIG}_t${T}"
done

echo ">>> Scoring (gene2net env)"
conda activate gene2net
for T in $THRESHOLDS; do
    echo "=================== SCORE threshold $T ==================="
    python scripts/score_reconstructions.py \
        --recon-dir "$MODEL_DIR/benchmark/${CONFIG}_t${T}" --workers "$NWORKERS"
done

echo "Done."
