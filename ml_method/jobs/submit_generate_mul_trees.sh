#!/bin/bash
# ==============================================================================
# Generate random MUL-trees for Gene2Net-GNN training data
#
# USAGE:
#   ./submit_generate_mul_trees.sh [N_TREES] [OUTPUT_DIR] [TREE_HEIGHT]
#
# ARGUMENTS:
#   1. N_TREES    : Number of MUL-trees to generate (default: 500)
#   2. OUTPUT_DIR : Output directory (default: /groups/.../ml_method/data/mul_trees)
#   3. TREE_HEIGHT: Root-to-tip in generations (default: 10000000 = 10M)
#
# EXAMPLES:
#   ./submit_generate_mul_trees.sh                          # 500 trees, default path
#   ./submit_generate_mul_trees.sh 100                      # 100 trees (quick test)
#   ./submit_generate_mul_trees.sh 500 /path/to/output      # Custom output dir
#
# Each array task generates one MUL-tree (species tree via SimPhy birth-death).
# Output:
#   output_dir/
#     mul_tree_NNNN.nex         # MUL-tree in NEXUS (input for SimPhy gene trees via -SR)
#     species_tree_NNNN.nex     # Original species tree (backbone)
#     metadata_NNNN.json        # Parameters, events, polyploid info
# ==============================================================================

N_TREES=${1:-500}
OUTPUT_DIR=${2:-/groups/itay_mayrose/tomulanovski/gene2net/ml_method/data/mul_trees}
TREE_HEIGHT=${3:-10000000}

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"
SCRIPT="${BASE_DIR}/scripts/generate_mul_trees.py"

mkdir -p "$LOG_DIR"

echo "================================================"
echo "Generating ${N_TREES} random MUL-trees"
echo "Output: ${OUTPUT_DIR}"
echo "Tree height: ${TREE_HEIGHT} generations"
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
    --export=ALL,OUTPUT_DIR="${OUTPUT_DIR}",TREE_HEIGHT="${TREE_HEIGHT}",SCRIPT="${SCRIPT}" \
    "${BASE_DIR}/jobs/generate_mul_trees_job.sh"
