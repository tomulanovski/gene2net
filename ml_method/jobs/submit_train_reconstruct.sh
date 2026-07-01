#!/bin/bash
# ==============================================================================
# Submit reconstruction training: WGD detection + partner-edge prediction.
#
# USAGE:
#   ./submit_train_reconstruct.sh [DATA_DIR] [--all-configs] [--general]
#                                 [--config CONFIG] [--output-dir DIR] [--init-from CKPT]
#
# EXAMPLES:
#   ./submit_train_reconstruct.sh --all-configs --output-dir output/reconstruct_all9
#   ./submit_train_reconstruct.sh --all-configs --init-from output/phase1_feat9_full/best_f1_model.pt
# ==============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
DEFAULT_DATA="${BASE_DIR}/data/mul_trees_2k/training/ils_low"

ALL_CONFIGS_DIRS="\
${BASE_DIR}/data/mul_trees_2k/training/ils_low \
${BASE_DIR}/data/mul_trees_2k/training/ils_medium \
${BASE_DIR}/data/mul_trees_2k/training/ils_high \
${BASE_DIR}/data/mul_trees_2k/training/dup_loss_low \
${BASE_DIR}/data/mul_trees_2k/training/dup_loss_medium \
${BASE_DIR}/data/mul_trees_2k/training/dup_loss_high \
${BASE_DIR}/data/mul_trees_2k/training/dup_loss_low_ne1M \
${BASE_DIR}/data/mul_trees_2k/training/dup_loss_medium_ne1M \
${BASE_DIR}/data/mul_trees_2k/training/dup_loss_high_ne1M"

DATA_DIR=""
USE_GENERAL=false
USE_ALL=false
CONFIG="${BASE_DIR}/configs/reconstruct.yaml"
OUTPUT_DIR="${BASE_DIR}/output/reconstruct"
INIT_FROM=""
CLADE_LABELS=""
ROOTED=""

while [ $# -gt 0 ]; do
    case "$1" in
        --general) USE_GENERAL=true; shift;;
        --all-configs) USE_ALL=true; shift;;
        --config) CONFIG="${BASE_DIR}/$2"; shift 2;;
        --output-dir) OUTPUT_DIR="${BASE_DIR}/$2"; shift 2;;
        --init-from) INIT_FROM="${BASE_DIR}/$2"; shift 2;;
        --clade-labels) CLADE_LABELS=1; shift;;
        --rooted) ROOTED=1; shift;;
        *) [ -z "$DATA_DIR" ] && DATA_DIR="$1"; shift;;
    esac
done

if [ "$USE_ALL" = true ]; then
    DATA_DIR="$ALL_CONFIGS_DIRS"
else
    DATA_DIR="${DATA_DIR:-$DEFAULT_DATA}"
fi

# --rooted: point at the rooted re-packaged data (training_rooted/) instead of training/.
if [ -n "$ROOTED" ]; then
    DATA_DIR="${DATA_DIR//\/training\//\/training_rooted\/}"
fi

if [ "$USE_GENERAL" = true ]; then
    PARTITION="gpu-general-pool"; ACCOUNT=""; QOS="public"
else
    PARTITION="gpu-rotemhsh-pool"; ACCOUNT="--account=itaym-users_v2"; QOS="owner"
fi

LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "Reconstruction Training (detection + partner)"
echo "Data: ${DATA_DIR}"
echo "Config: ${CONFIG}"
echo "Output: ${OUTPUT_DIR}"
echo "Init from: ${INIT_FROM:-<none>}"
echo "Clade labels: ${CLADE_LABELS:-<no>}"
echo "Partition: ${PARTITION}"
echo "================================================"

sbatch \
    --job-name="gene2net_reconstruct" \
    --output="${LOG_DIR}/reconstruct_%j.out" \
    --error="${LOG_DIR}/reconstruct_%j.err" \
    --time=1-00:00:00 \
    --mem=128G \
    --cpus-per-task=4 \
    --partition=${PARTITION} \
    ${ACCOUNT} \
    --qos=${QOS} \
    --gres=gpu:1 \
    --export=ALL,DATA_DIR="${DATA_DIR}",BASE_DIR="${BASE_DIR}",CONFIG="${CONFIG}",OUTPUT_DIR="${OUTPUT_DIR}",INIT_FROM="${INIT_FROM}",CLADE_LABELS="${CLADE_LABELS}" \
    "${BASE_DIR}/jobs/train_reconstruct_job.sh"
