#!/bin/bash
# ==============================================================================
# Submit ASTRAL jobs for ML training MUL-trees (after SimPhy completes)
#
# USAGE:
#   ./submit_astral_training.sh <MUL_TREES_DIR> <CONFIG_NAME> [NUM_REPLICATES]
#
# EXAMPLES:
#   ./submit_astral_training.sh /path/to/mul_trees ils_low
#   ./submit_astral_training.sh /path/to/mul_trees ils_high 1
# ==============================================================================

MUL_TREES_DIR=$1
CONFIG_NAME=${2:-ils_low}
NUM_REPLICATES=${3:-1}

if [ -z "$MUL_TREES_DIR" ]; then
    echo "Error: MUL_TREES_DIR is required."
    echo "Usage: ./submit_astral_training.sh <MUL_TREES_DIR> <CONFIG_NAME> [NUM_REPLICATES]"
    exit 1
fi

# Count MUL-trees
N_TREES=$(ls -1 "${MUL_TREES_DIR}"/mul_tree_*.nex 2>/dev/null | wc -l)
if [ "$N_TREES" -eq 0 ]; then
    echo "Error: No mul_tree_*.nex files found in ${MUL_TREES_DIR}"
    exit 1
fi

MAX_IDX=$((N_TREES - 1))

# Check SimPhy output exists for at least the first tree
if [ ! -d "${MUL_TREES_DIR}/simphy/${CONFIG_NAME}/0000" ]; then
    echo "Error: SimPhy output not found at ${MUL_TREES_DIR}/simphy/${CONFIG_NAME}/"
    echo "Did you run submit_simphy_training.sh first?"
    exit 1
fi

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "ASTRAL for ML training"
echo "================================================"
echo "MUL-trees: ${N_TREES} (indices 0-${MAX_IDX})"
echo "Config: ${CONFIG_NAME}"
echo "Replicates: ${NUM_REPLICATES}"
echo "================================================"

sbatch \
    --array="0-${MAX_IDX}" \
    --job-name="astral_train_${CONFIG_NAME}" \
    --output="${LOG_DIR}/astral_train_${CONFIG_NAME}_%a.out" \
    --error="${LOG_DIR}/astral_train_${CONFIG_NAME}_%a.err" \
    --time=02:00:00 \
    --mem=16G \
    --cpus-per-task=8 \
    --partition=itaym-pool \
    --account=itaym-users_v2 \
    --qos=owner \
    --export=ALL,MUL_TREES_DIR="${MUL_TREES_DIR}",CONFIG_NAME="${CONFIG_NAME}",NUM_REPLICATES="${NUM_REPLICATES}" \
    "${BASE_DIR}/jobs/run_astral_training.sh"
