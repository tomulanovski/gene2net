#!/bin/bash
# ==============================================================================
# Generate random MUL-trees for Gene2Net-GNN training data
#
# USAGE:
#   ./submit_generate_mul_trees.sh [N_TREES] [OUTPUT_DIR]
#
# ARGUMENTS:
#   1. N_TREES    : Number of MUL-trees to generate (default: 500)
#   2. OUTPUT_DIR : Output directory (default: /groups/.../ml_method/data/mul_trees)
#
# The script auto-detects existing trees in OUTPUT_DIR and continues from the
# next available index. Just run it repeatedly to add more trees.
#
# EXAMPLES:
#   ./submit_generate_mul_trees.sh 1000 /path/to/mul_trees   # first run: indices 0-999
#   ./submit_generate_mul_trees.sh 1000 /path/to/mul_trees   # second run: indices 1000-1999
#   ./submit_generate_mul_trees.sh 1000 /path/to/mul_trees   # third run: indices 2000-2999
#
# Each array task generates one MUL-tree (species tree via SimPhy birth-death).
# Tree index = SLURM_ARRAY_TASK_ID + auto_offset → different index = different tree.
#
# Output:
#   output_dir/
#     mul_tree_NNNN.nex         # MUL-tree in NEXUS (input for SimPhy gene trees via -SR)
#     species_tree_NNNN.nex     # Original species tree (backbone)
#     metadata_NNNN.json        # Parameters, events, polyploid info
# ==============================================================================

N_TREES=${1:-500}
OUTPUT_DIR=${2:-/groups/itay_mayrose/tomulanovski/gene2net/ml_method/data/mul_trees}

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"
SCRIPT="${BASE_DIR}/scripts/generate_mul_trees.py"

mkdir -p "$LOG_DIR" "$OUTPUT_DIR"

# Auto-detect offset: find highest existing index in the directory
# Looks at mul_tree_NNNN.nex filenames and picks max index + 1
INDEX_OFFSET=0
if ls "$OUTPUT_DIR"/mul_tree_*.nex &>/dev/null; then
    HIGHEST=$(ls "$OUTPUT_DIR"/mul_tree_*.nex | sed 's/.*mul_tree_\([0-9]*\)\.nex/\1/' | sort -n | tail -1)
    # Strip leading zeros to avoid bash octal interpretation
    HIGHEST=$(echo "$HIGHEST" | sed 's/^0*//')
    HIGHEST=${HIGHEST:-0}
    INDEX_OFFSET=$((HIGHEST + 1))
fi

echo "================================================"
echo "Generating ${N_TREES} random MUL-trees"
echo "Output: ${OUTPUT_DIR}"
echo "Existing trees detected: ${INDEX_OFFSET}"
echo "New indices: ${INDEX_OFFSET}-$((INDEX_OFFSET + N_TREES - 1))"
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
    --export=ALL,OUTPUT_DIR="${OUTPUT_DIR}",SCRIPT="${SCRIPT}",INDEX_OFFSET="${INDEX_OFFSET}" \
    "${BASE_DIR}/jobs/generate_mul_trees_job.sh"
