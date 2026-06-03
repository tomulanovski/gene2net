#!/bin/bash
# ==============================================================================
# Append the 5 WGD-detection edge features to ALREADY-packaged samples,
# avoiding a full repackage. Edits sample.pkl files in place (4 -> 9 edge dims).
#
# USAGE:
#   ./submit_augment_features.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [N_PER_JOB]
#
# EXAMPLES:
#   ./submit_augment_features.sh 2000 /path/to/mul_trees_2k 0 ils_low
#   ./submit_augment_features.sh 2000 /path/to/mul_trees_2k 0 ils_medium 20
# ==============================================================================

N_TREES=${1:?Usage: submit_augment_features.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [N_PER_JOB]}
MUL_TREES_DIR=${2:?Usage: submit_augment_features.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [N_PER_JOB]}
START_INDEX=${3:-0}
CONFIG_NAME=${4:-ils_low}
N_PER_JOB=${5:-20}
N_TASKS=$(( (N_TREES + N_PER_JOB - 1) / N_PER_JOB ))

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "Augmenting edge features (4 -> 9) in place"
echo "================================================"
echo "MUL-trees dir: ${MUL_TREES_DIR}"
echo "Batch: ${N_TREES} samples, indices ${START_INDEX}-$((START_INDEX + N_TREES - 1))"
echo "Config: ${CONFIG_NAME}"
echo "Bundling: ${N_PER_JOB} samples per job → ${N_TASKS} array tasks"
echo "Target: ${MUL_TREES_DIR}/training/${CONFIG_NAME}/"
echo "================================================"

sbatch \
    --array="0-$((N_TASKS - 1))" \
    --job-name="augment_${CONFIG_NAME}" \
    --output="${LOG_DIR}/augment_${CONFIG_NAME}_%a.out" \
    --error="${LOG_DIR}/augment_${CONFIG_NAME}_%a.err" \
    --partition=itaym-pool \
    --account=itaym-users_v2 \
    --qos=owner \
    --export=ALL,MUL_TREES_DIR="${MUL_TREES_DIR}",CONFIG_NAME="${CONFIG_NAME}",INDEX_OFFSET="${START_INDEX}",N_PER_JOB="${N_PER_JOB}" \
    "${BASE_DIR}/jobs/augment_features_job.sh"
