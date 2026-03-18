#!/bin/bash
# ==============================================================================
# Package SimPhy + ASTRAL output into Gene2Net training samples
#
# USAGE:
#   ./submit_package_training.sh <MUL_TREES_DIR> <CONFIG_NAME>
#
# EXAMPLES:
#   ./submit_package_training.sh /path/to/mul_trees ils_low
#   ./submit_package_training.sh /path/to/mul_trees ils_high
# ==============================================================================

MUL_TREES_DIR=$1
CONFIG_NAME=${2:-ils_low}

if [ -z "$MUL_TREES_DIR" ]; then
    echo "Error: MUL_TREES_DIR is required."
    exit 1
fi

N_TREES=$(ls -1 "${MUL_TREES_DIR}"/mul_tree_*.nex 2>/dev/null | wc -l)
if [ "$N_TREES" -eq 0 ]; then
    echo "Error: No mul_tree_*.nex files found"
    exit 1
fi

MAX_IDX=$((N_TREES - 1))

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "Packaging training data"
echo "================================================"
echo "MUL-trees: ${N_TREES}"
echo "Config: ${CONFIG_NAME}"
echo "Output: ${MUL_TREES_DIR}/training/${CONFIG_NAME}/"
echo "================================================"

sbatch \
    --array="0-${MAX_IDX}" \
    --job-name="package_${CONFIG_NAME}" \
    --output="${LOG_DIR}/package_${CONFIG_NAME}_%a.out" \
    --error="${LOG_DIR}/package_${CONFIG_NAME}_%a.err" \
    --time=00:30:00 \
    --mem=8G \
    --cpus-per-task=1 \
    --partition=itaym-pool \
    --account=itaym-users_v2 \
    --qos=owner \
    --export=ALL,MUL_TREES_DIR="${MUL_TREES_DIR}",CONFIG_NAME="${CONFIG_NAME}" \
    "${BASE_DIR}/jobs/package_training_job.sh"
