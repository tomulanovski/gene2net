#!/bin/bash
# ==============================================================================
# Submit the partner-feature ablation batch:
#   baseline-from-scratch (reconstruct.yaml) + B (dupcond) + D (neff) + B+D.
# Reuses jobs/train_reconstruct_job.sh. Trains from scratch (NO --init-from) so
# every variant starts identically and differs only in the features. Same 6
# rooted config dirs + --clade-labels as the reconstruct_cladelabels_rooted
# baseline (see train_two_parent_job.sh), and the fixed seed=42 val split makes
# the val partner/allo numbers directly comparable across runs.
#
# USAGE (from anywhere):  ./jobs/submit_ablation_batch.sh
# ==============================================================================
set -euo pipefail

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"

DATA_DIR="\
${BASE_DIR}/data/mul_trees_2k/training_rooted/ils_low \
${BASE_DIR}/data/mul_trees_2k/training_rooted/ils_medium \
${BASE_DIR}/data/mul_trees_2k/training_rooted/ils_high \
${BASE_DIR}/data/mul_trees_2k/training_rooted/dup_loss_low_ne1M \
${BASE_DIR}/data/mul_trees_2k/training_rooted/dup_loss_medium_ne1M \
${BASE_DIR}/data/mul_trees_2k/training_rooted/dup_loss_high_ne1M"

# config basename (in configs/, no .yaml)  ->  output dir name (in output/)
# One-parter (our actual method): baseline + B + D + B+D.
declare -A RUNS=(
    [reconstruct_oneparter]=reconstruct_op_baseline
    [reconstruct_dupcond]=reconstruct_op_dupcond
    [reconstruct_neff]=reconstruct_op_neff
    [reconstruct_dupcond_neff]=reconstruct_op_dupcond_neff
)

PARTITION="gpu-rotemhsh-pool"
ACCOUNT="--account=itaym-users_v2"
QOS="owner"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

for cfg in "${!RUNS[@]}"; do
    out="${RUNS[$cfg]}"
    echo "Submitting: config=${cfg}.yaml  ->  output/${out}"
    sbatch \
        --job-name="g2n_abl_${out}" \
        --output="${LOG_DIR}/abl_${out}_%j.out" \
        --error="${LOG_DIR}/abl_${out}_%j.err" \
        --time=1-00:00:00 \
        --mem=128G \
        --cpus-per-task=4 \
        --partition=${PARTITION} \
        ${ACCOUNT} \
        --qos=${QOS} \
        --gres=gpu:1 \
        --export=ALL,DATA_DIR="${DATA_DIR}",BASE_DIR="${BASE_DIR}",CONFIG="${BASE_DIR}/configs/${cfg}.yaml",OUTPUT_DIR="${BASE_DIR}/output/${out}",INIT_FROM="",CLADE_LABELS="1" \
        "${BASE_DIR}/jobs/train_reconstruct_job.sh"
done

echo "Submitted ${#RUNS[@]} jobs. Watch: squeue -u \$USER ; tail -f ${LOG_DIR}/abl_*_*.out"
