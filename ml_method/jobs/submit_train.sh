#!/bin/bash
# ==============================================================================
# Submit Gene2Net-GNN training job
#
# USAGE:
#   ./submit_train.sh [DATA_DIR] [--general]
#
# ARGUMENTS:
#   1. DATA_DIR   : Training data directory (default: .../data/mul_trees/training/ils_low)
#   --general     : Use gpu-general-pool instead of gpu-rotemhsh-pool
#
# EXAMPLES:
#   ./submit_train.sh                                                    # default data, rotem GPU
#   ./submit_train.sh /path/to/training/ils_low                         # custom data, rotem GPU
#   ./submit_train.sh /path/to/training/ils_low --general               # custom data, general GPU
#   ./submit_train.sh --general                                         # default data, general GPU
# ==============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
DEFAULT_DATA="${BASE_DIR}/data/mul_trees/training/ils_low"

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
    MEM="64G"
    TIME="08:00:00"
else
    PARTITION="gpu-rotemhsh-pool"
    ACCOUNT="--account=itaym-users_v2"
    QOS="owner"
    GPU_FLAG="gpu:1"
    MEM="48G"
    TIME="08:00:00"
fi

LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "Gene2Net-GNN Training"
echo "================================================"
echo "Data: ${DATA_DIR}"
echo "Partition: ${PARTITION}"
echo "QoS: ${QOS}"
echo "GPU: ${GPU_FLAG}"
echo "Config: ${BASE_DIR}/configs/default.yaml"
echo "Output: ${BASE_DIR}/output"
echo "Logs: ${LOG_DIR}/train_*.out"
echo "================================================"

sbatch \
    --job-name="gene2net_train" \
    --output="${LOG_DIR}/train_%j.out" \
    --error="${LOG_DIR}/train_%j.err" \
    --time=${TIME} \
    --mem=${MEM} \
    --cpus-per-task=4 \
    --partition=${PARTITION} \
    ${ACCOUNT} \
    --qos=${QOS} \
    --gres=${GPU_FLAG} \
    --export=ALL,DATA_DIR="${DATA_DIR}",BASE_DIR="${BASE_DIR}" \
    "${BASE_DIR}/jobs/train_job.sh"
