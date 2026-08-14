#!/bin/bash
# ==============================================================================
# Launch all ablation training runs for the method chapter, one sbatch each:
#   learning curve : train on 250 / 500 / 1000 / 2000 networks (data saturation)
#   depth          : n_gat_layers 1 / 2 / 4   (base 3 = reconstruct_op_away_clade)
#   width          : hidden_dim 32 / 128       (base 64 = reconstruct_op_away_clade)
#   learning rate  : 0.0003                     (base 0.001)
# All use the final setup: one-parter + away labels (labels_away.pkl, single+clade).
# The base (depth 3, width 64, full data) is the already-trained reconstruct_op_away_clade,
# so it is NOT re-run here. Prereq: configs generated (python scripts/gen_ablation_configs.py).
#
# USAGE:  ./jobs/submit_ablations.sh
# ==============================================================================
set -euo pipefail

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
DATA_DIR=""
for c in ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
    DATA_DIR="$DATA_DIR ${BASE_DIR}/data/mul_trees_2k/training_rooted/$c"
done

PARTITION="gpu-rotemhsh-pool"; ACCOUNT="--account=itaym-users_v2"; QOS="owner"
LOG_DIR="${BASE_DIR}/logs"; mkdir -p "$LOG_DIR"

submit () {  # $1=tag  $2=config basename (in configs/)  $3=extra --export (comma-sep, no leading comma)
    local tag="$1" cfg="$2" extra="$3"
    echo "Submitting: ${tag}  (config ${cfg}.yaml)"
    sbatch \
        --job-name="g2n_${tag}" \
        --output="${LOG_DIR}/${tag}_%j.out" \
        --error="${LOG_DIR}/${tag}_%j.err" \
        --time=1-00:00:00 --mem=128G --cpus-per-task=4 \
        --partition=${PARTITION} ${ACCOUNT} --qos=${QOS} --gres=gpu:1 \
        --export=ALL,DATA_DIR="${DATA_DIR}",BASE_DIR="${BASE_DIR}",CONFIG="${BASE_DIR}/configs/${cfg}.yaml",OUTPUT_DIR="${BASE_DIR}/output/${tag}",INIT_FROM="",${extra} \
        "${BASE_DIR}/jobs/train_reconstruct_job.sh"
}

# learning curve (data saturation)
for N in 250 500 1000 2000; do
    submit "reconstruct_lc_${N}" reconstruct_oneparter "AWAY_LABELS=1,MAX_TRAIN=${N}"
done
# depth ablation
for d in 1 2 4; do
    submit "reconstruct_abl_depth${d}" "abl_depth${d}" "AWAY_LABELS=1"
done
# width ablation
for h in 32 128; do
    submit "reconstruct_abl_hidden${h}" "abl_hidden${h}" "AWAY_LABELS=1"
done
# learning-rate ablation
submit "reconstruct_abl_lr3e-4" "abl_lr3e-4" "AWAY_LABELS=1"

echo "Submitted 10 ablation jobs. Watch: squeue -u \$USER"
