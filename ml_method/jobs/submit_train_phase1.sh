#!/bin/bash
# ==============================================================================
# Submit Phase 1 training: features-only GNN for binary WGD detection
#
# USAGE:
#   ./submit_train_phase1.sh [DATA_DIR] [--general]
#
# EXAMPLES:
#   ./submit_train_phase1.sh                                    # default data, rotem GPU
#   ./submit_train_phase1.sh /path/to/training/ils_low          # custom data
#   ./submit_train_phase1.sh --general                          # general GPU pool
# ==============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
DEFAULT_DATA="${BASE_DIR}/data/mul_trees_2k/training/ils_low"

# Parse arguments
DATA_DIR=""
USE_GENERAL=false

for arg in "$@"; do
    if [ "$arg" == "--general" ]; then
        USE_GENERAL=true
    elif [ -z "$DATA_DIR" ]; then
        DATA_DIR="$arg"
    fi
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
echo "Partition: ${PARTITION}"
echo "QoS: ${QOS}"
echo "Config: ${BASE_DIR}/configs/phase1.yaml"
echo "Output: ${BASE_DIR}/output/phase1"
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
    --export=ALL,DATA_DIR="${DATA_DIR}",BASE_DIR="${BASE_DIR}" \
    "${BASE_DIR}/jobs/train_phase1_job.sh"
