#!/bin/bash
# ==============================================================================
# Submit SimPhy gene tree simulations for ML training MUL-trees
#
# USAGE:
#   ./submit_simphy_training.sh <MUL_TREES_DIR> <CONFIG_NAME> [NE] [NUM_GENE_TREES] [NUM_REPLICATES]
#
# ARGUMENTS:
#   1. MUL_TREES_DIR   : Directory containing mul_tree_NNNN.nex files
#   2. CONFIG_NAME     : Configuration name (e.g., "ils_low", "ils_high")
#   3. NE              : Effective population size (default: 200000)
#   4. NUM_GENE_TREES  : Number of gene trees per replicate (default: 500)
#   5. NUM_REPLICATES  : Number of replicates (default: 1)
#
# EXAMPLES:
#   ./submit_simphy_training.sh /path/to/mul_trees ils_low
#   ./submit_simphy_training.sh /path/to/mul_trees ils_low 200000 500 1
#   ./submit_simphy_training.sh /path/to/mul_trees ils_high 2000000 500 1
#
# OUTPUT STRUCTURE:
#   MUL_TREES_DIR/
#     mul_tree_0000.nex              # Input (already exists)
#     simphy/
#       ils_low/
#         0000/replicate_1/1/        # SimPhy output (gene trees)
#         0001/replicate_1/1/
#         ...
# ==============================================================================

MUL_TREES_DIR=$1
CONFIG_NAME=${2:-ils_low}
NE_VAL=${3:-200000}
NUM_GENE_TREES=${4:-500}
NUM_REPLICATES=${5:-1}

if [ -z "$MUL_TREES_DIR" ]; then
    echo "Error: MUL_TREES_DIR is required."
    echo "Usage: ./submit_simphy_training.sh <MUL_TREES_DIR> <CONFIG_NAME> [NE] [NUM_GENE_TREES] [NUM_REPLICATES]"
    exit 1
fi

# Count how many mul_tree_NNNN.nex files exist
N_TREES=$(ls -1 "${MUL_TREES_DIR}"/mul_tree_*.nex 2>/dev/null | wc -l)

if [ "$N_TREES" -eq 0 ]; then
    echo "Error: No mul_tree_*.nex files found in ${MUL_TREES_DIR}"
    exit 1
fi

MAX_IDX=$((N_TREES - 1))

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"

mkdir -p "$LOG_DIR"

echo "================================================"
echo "SimPhy gene trees for ML training"
echo "================================================"
echo "MUL-trees dir: ${MUL_TREES_DIR}"
echo "MUL-trees found: ${N_TREES} (indices 0-${MAX_IDX})"
echo "Config: ${CONFIG_NAME}"
echo "Ne: ${NE_VAL}"
echo "Gene trees per replicate: ${NUM_GENE_TREES}"
echo "Replicates: ${NUM_REPLICATES}"
echo "Logs: ${LOG_DIR}/simphy_train_${CONFIG_NAME}_*.out"
echo "================================================"

sbatch \
    --array="0-${MAX_IDX}" \
    --job-name="simphy_train_${CONFIG_NAME}" \
    --output="${LOG_DIR}/simphy_train_${CONFIG_NAME}_%a.out" \
    --error="${LOG_DIR}/simphy_train_${CONFIG_NAME}_%a.err" \
    --time=04:00:00 \
    --mem=16G \
    --cpus-per-task=1 \
    --partition=itaym-pool \
    --account=itaym-users_v2 \
    --qos=owner \
    --export=ALL,MUL_TREES_DIR="${MUL_TREES_DIR}",CONFIG_NAME="${CONFIG_NAME}",NE_VAL="${NE_VAL}",NUM_GENE_TREES="${NUM_GENE_TREES}",NUM_REPLICATES="${NUM_REPLICATES}" \
    "${BASE_DIR}/jobs/simphy_training_job.sh"
