#!/bin/bash
# ==============================================================================
# Submit Phase 1 training: features-only GNN for binary WGD detection
#
# USAGE:
#   ./submit_train_phase1.sh [DATA_DIR] [--general] [--config CONFIG] [--output-dir DIR]
#
# EXAMPLES:
#   ./submit_train_phase1.sh                                              # defaults
#   ./submit_train_phase1.sh --config configs/phase1_cw6.yaml --output-dir output/phase1_cw6
#   ./submit_train_phase1.sh --general
# ==============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
DEFAULT_DATA="${BASE_DIR}/data/mul_trees_2k/training/ils_low"

# Parse arguments
DATA_DIR=""
USE_GENERAL=false
CONFIG="${BASE_DIR}/configs/phase1.yaml"
OUTPUT_DIR="${BASE_DIR}/output/phase1"

while [ $# -gt 0 ]; do
    case "$1" in
        --general) USE_GENERAL=true; shift;;
        --config) CONFIG="${BASE_DIR}/$2"; shift 2;;
        --output-dir) OUTPUT_DIR="${BASE_DIR}/$2"; shift 2;;
        *) [ -z "$DATA_DIR" ] && DATA_DIR="$1"; shift;;
    esac
done

DATA_DIR="${DATA_DIR:-$DEFAULT_DATA}"

# Set partition
if [ "$USE_GENERAL" = true ]; then
    PARTITION="gpu-general-pool"
    ACCOUNT=""
    QOS="public"
    GPU_FLAG="gpu:1"
    MEM="32G"
    TIME="1-00:00:00"
else
    PARTITION="gpu-rotemhsh-pool"
    ACCOUNT="--account=itaym-users_v2"
    QOS="owner"
    GPU_FLAG="gpu:1"
    MEM="32G"
    TIME="1-00:00:00"
fi

LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "Phase 1: Features-only GNN Training"
echo "================================================"
echo "Data: ${DATA_DIR}"
echo "Config: ${CONFIG}"
echo "Output: ${OUTPUT_DIR}"
echo "Partition: ${PARTITION}"
echo "================================================"

sbatch \
    --job-name="gene2net_phase1" \
    --output="${LOG_DIR}/phase1_%j.out" \
    --error="${LOG_DIR}/phase1_%j.err" \
    --time=${TIME} \
    --mem=${MEM} \
    --cpus-per-task=4 \
    --partition=${PARTITION} \
    ${ACCOUNT} \
    --qos=${QOS} \
    --gres=${GPU_FLAG} \
    --export=ALL,DATA_DIR="${DATA_DIR}",BASE_DIR="${BASE_DIR}",CONFIG="${CONFIG}",OUTPUT_DIR="${OUTPUT_DIR}" \
    "${BASE_DIR}/jobs/train_phase1_job.sh"
