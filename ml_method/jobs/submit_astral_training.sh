#!/bin/bash
# ==============================================================================
# Submit ASTRAL jobs for ML training MUL-trees (after SimPhy completes)
#
# USAGE:
#   ./submit_astral_training.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [NUM_REPLICATES]
#
# EXAMPLES:
#   ./submit_astral_training.sh 1000 /path/to/mul_trees 0            # indices 0-999
#   ./submit_astral_training.sh 1000 /path/to/mul_trees 1000         # indices 1000-1999
#   ./submit_astral_training.sh 200  /path/to/mul_trees 1000         # indices 1000-1199
# ==============================================================================

N_TREES=${1:?Usage: submit_astral_training.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [NUM_REPLICATES]}
MUL_TREES_DIR=${2:?Usage: submit_astral_training.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [NUM_REPLICATES]}
START_INDEX=${3:-0}
CONFIG_NAME=${4:-ils_low}
NUM_REPLICATES=${5:-1}

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "ASTRAL for ML training"
echo "================================================"
echo "MUL-trees dir: ${MUL_TREES_DIR}"
echo "Batch: ${N_TREES} trees, indices ${START_INDEX}-$((START_INDEX + N_TREES - 1))"
echo "Config: ${CONFIG_NAME}"
echo "Replicates: ${NUM_REPLICATES}"
echo "================================================"

sbatch \
    --array="0-$((N_TREES - 1))" \
    --job-name="astral_train_${CONFIG_NAME}" \
    --output="${LOG_DIR}/astral_train_${CONFIG_NAME}_%a.out" \
    --error="${LOG_DIR}/astral_train_${CONFIG_NAME}_%a.err" \
    --time=02:00:00 \
    --mem=16G \
    --cpus-per-task=8 \
    --partition=itaym-pool \
    --account=itaym-users_v2 \
    --qos=owner \
    --export=ALL,MUL_TREES_DIR="${MUL_TREES_DIR}",CONFIG_NAME="${CONFIG_NAME}",NUM_REPLICATES="${NUM_REPLICATES}",INDEX_OFFSET="${START_INDEX}" \
    "${BASE_DIR}/jobs/run_astral_training.sh"
