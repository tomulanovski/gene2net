#!/bin/bash
# ==============================================================================
# Generate random MUL-trees for Gene2Net-GNN training data
#
# USAGE:
#   ./submit_generate_mul_trees.sh N_TREES OUTPUT_DIR START_INDEX
#
# ARGUMENTS:
#   1. N_TREES     : Number of MUL-trees to generate
#   2. OUTPUT_DIR  : Output directory
#   3. START_INDEX : First tree index (default: 0)
#
# EXAMPLES:
#   ./submit_generate_mul_trees.sh 1000 /path/to/mul_trees 0       # indices 0-999
#   ./submit_generate_mul_trees.sh 1000 /path/to/mul_trees 1000    # indices 1000-1999
#   ./submit_generate_mul_trees.sh 500  /path/to/mul_trees 2000    # indices 2000-2499
#
# Each tree index produces a unique, deterministic tree.
# Files that already exist at a given index will be overwritten (same content).
#
# Output per index:
#     mul_tree_NNNN.nex         # MUL-tree in NEXUS
#     species_tree_NNNN.nex     # Original species tree (backbone)
#     metadata_NNNN.json        # Parameters, events, polyploid info
# ==============================================================================

N_TREES=${1:?Usage: submit_generate_mul_trees.sh N_TREES OUTPUT_DIR [START_INDEX]}
OUTPUT_DIR=${2:?Usage: submit_generate_mul_trees.sh N_TREES OUTPUT_DIR [START_INDEX]}
START_INDEX=${3:-0}

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"
SCRIPT="${BASE_DIR}/scripts/generate_mul_trees.py"

mkdir -p "$LOG_DIR" "$OUTPUT_DIR"

echo "================================================"
echo "Generating ${N_TREES} random MUL-trees"
echo "Output: ${OUTPUT_DIR}"
echo "Indices: ${START_INDEX}-$((START_INDEX + N_TREES - 1))"
echo "Logs: ${LOG_DIR}/gen_mul_trees_*.out"
echo "================================================"

sbatch \
    --array="0-$((N_TREES - 1))" \
    --job-name="gen_mul_trees" \
    --output="${LOG_DIR}/gen_mul_trees_%a.out" \
    --error="${LOG_DIR}/gen_mul_trees_%a.err" \
    --time=00:15:00 \
    --mem=4G \
    --cpus-per-task=1 \
    --partition=itaym-pool \
    --account=itaym-users_v2 \
    --qos=owner \
    --export=ALL,OUTPUT_DIR="${OUTPUT_DIR}",SCRIPT="${SCRIPT}",INDEX_OFFSET="${START_INDEX}" \
    "${BASE_DIR}/jobs/generate_mul_trees_job.sh"
