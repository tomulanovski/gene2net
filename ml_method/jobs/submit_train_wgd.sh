#!/bin/bash
# ==============================================================================
# Submit WGD Detector training: GIN gene tree encoder + GAT species tree
#
# USAGE:
#   ./submit_train_wgd.sh [DATA_DIR] [--all-configs] [--general] [--config CONFIG] [--output-dir DIR]
#
# EXAMPLES:
#   ./submit_train_wgd.sh                              # defaults (ils_low only)
#   ./submit_train_wgd.sh --all-configs                # all 9 configs (~17.9K samples)
#   ./submit_train_wgd.sh --all-configs --output-dir output/wgd_all9
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

# Parse arguments
DATA_DIR=""
USE_GENERAL=false
USE_ALL=false
CONFIG="${BASE_DIR}/configs/wgd_detector.yaml"
OUTPUT_DIR="${BASE_DIR}/output/wgd_detector"

while [ $# -gt 0 ]; do
    case "$1" in
        --general) USE_GENERAL=true; shift;;
        --all-configs) USE_ALL=true; shift;;
        --config) CONFIG="${BASE_DIR}/$2"; shift 2;;
        --output-dir) OUTPUT_DIR="${BASE_DIR}/$2"; shift 2;;
        *) [ -z "$DATA_DIR" ] && DATA_DIR="$1"; shift;;
    esac
done

if [ "$USE_ALL" = true ]; then
    DATA_DIR="$ALL_CONFIGS_DIRS"
else
    DATA_DIR="${DATA_DIR:-$DEFAULT_DATA}"
fi

# Set partition
if [ "$USE_GENERAL" = true ]; then
    PARTITION="gpu-general-pool"
    ACCOUNT=""
    QOS="public"
    GPU_FLAG="gpu:1"
    MEM="48G"
    TIME="2-00:00:00"
else
    PARTITION="gpu-rotemhsh-pool"
    ACCOUNT="--account=itaym-users_v2"
    QOS="owner"
    GPU_FLAG="gpu:1"
    MEM="48G"
    TIME="2-00:00:00"
fi

LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "WGD Detector Training"
echo "================================================"
echo "Data: ${DATA_DIR}"
echo "Config: ${CONFIG}"
echo "Output: ${OUTPUT_DIR}"
echo "Partition: ${PARTITION}"
echo "================================================"

sbatch \
    --job-name="gene2net_wgd" \
    --output="${LOG_DIR}/wgd_%j.out" \
    --error="${LOG_DIR}/wgd_%j.err" \
    --time=${TIME} \
    --mem=${MEM} \
    --cpus-per-task=4 \
    --partition=${PARTITION} \
    ${ACCOUNT} \
    --qos=${QOS} \
    --gres=${GPU_FLAG} \
    --export=ALL,DATA_DIR="${DATA_DIR}",BASE_DIR="${BASE_DIR}",CONFIG="${CONFIG}",OUTPUT_DIR="${OUTPUT_DIR}" \
    "${BASE_DIR}/jobs/train_wgd_job.sh"
