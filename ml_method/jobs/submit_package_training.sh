#!/bin/bash
# ==============================================================================
# Package SimPhy + ASTRAL output into Gene2Net training samples
#
# USAGE:
#   ./submit_package_training.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME]
#
# EXAMPLES:
#   ./submit_package_training.sh 1000 /path/to/mul_trees 0            # indices 0-999
#   ./submit_package_training.sh 1000 /path/to/mul_trees 1000         # indices 1000-1999
#   ./submit_package_training.sh 200  /path/to/mul_trees 1000         # indices 1000-1199
# ==============================================================================

N_TREES=${1:?Usage: submit_package_training.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [N_PER_JOB]}
MUL_TREES_DIR=${2:?Usage: submit_package_training.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [N_PER_JOB]}
START_INDEX=${3:-0}
CONFIG_NAME=${4:-ils_low}
N_PER_JOB=${5:-20}
N_TASKS=$(( (N_TREES + N_PER_JOB - 1) / N_PER_JOB ))

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "Packaging training data"
echo "================================================"
echo "MUL-trees dir: ${MUL_TREES_DIR}"
echo "Batch: ${N_TREES} trees, indices ${START_INDEX}-$((START_INDEX + N_TREES - 1))"
echo "Config: ${CONFIG_NAME}"
echo "Bundling: ${N_PER_JOB} samples per job → ${N_TASKS} array tasks"
echo "Output: ${MUL_TREES_DIR}/training/${CONFIG_NAME}/"
echo "================================================"

sbatch \
    --array="0-$((N_TASKS - 1))" \
    --job-name="package_${CONFIG_NAME}" \
    --output="${LOG_DIR}/package_${CONFIG_NAME}_%a.out" \
    --error="${LOG_DIR}/package_${CONFIG_NAME}_%a.err" \
    --partition=itaym-pool \
    --account=itaym-users_v2 \
    --qos=owner \
    --export=ALL,MUL_TREES_DIR="${MUL_TREES_DIR}",CONFIG_NAME="${CONFIG_NAME}",INDEX_OFFSET="${START_INDEX}",N_PER_JOB="${N_PER_JOB}" \
    "${BASE_DIR}/jobs/package_training_job.sh"
