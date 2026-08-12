#!/bin/bash
# ==============================================================================
# Finalize the locked model. Trains, all at configs/reconstruct_final.yaml with
# away labels (the HPO winner: gat, hidden 256, depth 4, dropout 0.1, lr 5e-4, wd 1e-5):
#   reconstruct_final  - full data. The shipped model AND the learning-curve 100% point.
#   lc_final_{F}       - learning curve at 10/25/50/75% of the training pool.
# 5 GPU jobs. After reconstruct_final finishes, run the re-benchmark against
# output/reconstruct_final (separate step).
#
# Usage:  ./jobs/submit_finalize.sh
# ==============================================================================
set -euo pipefail

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
DATA_DIR=""
for c in ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
    DATA_DIR="$DATA_DIR ${BASE_DIR}/data/mul_trees_2k/training_rooted/$c"
done

PARTITION="gpu-rotemhsh-pool"; ACCOUNT="--account=itaym-users_v2"; QOS="owner"
LOG_DIR="${BASE_DIR}/logs"; mkdir -p "$LOG_DIR"
CFG="${BASE_DIR}/configs/reconstruct_final.yaml"

submit () {   # $1=tag   $2=TRAIN_FRACTION (empty for full data)
    local tag="$1" frac="$2"
    echo "Submitting: ${tag}${frac:+  (fraction ${frac})}"
    sbatch \
        --job-name="g2n_${tag}" \
        --output="${LOG_DIR}/${tag}_%j.out" --error="${LOG_DIR}/${tag}_%j.err" \
        --time=12:00:00 --mem=96G --cpus-per-task=4 \
        --partition=${PARTITION} ${ACCOUNT} --qos=${QOS} --gres=gpu:1 \
        --export=ALL,DATA_DIR="${DATA_DIR}",BASE_DIR="${BASE_DIR}",CONFIG="${CFG}",OUTPUT_DIR="${BASE_DIR}/output/${tag}",INIT_FROM="",AWAY_LABELS=1${frac:+,TRAIN_FRACTION=$frac} \
        "${BASE_DIR}/jobs/train_reconstruct_job.sh"
}

submit "reconstruct_final" ""          # full data: shipped model + learning-curve 100% point
for F in 0.10 0.25 0.50 0.75; do
    submit "lc_final_${F}" "$F"
done

echo "Submitted 5 finalize jobs (final model + 4 learning-curve fractions)."
