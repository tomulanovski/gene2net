#!/bin/bash
# ==============================================================================
# Submit SimPhy gene tree simulations for ML training MUL-trees
#
# USAGE:
#   ./submit_simphy_training.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [NE] [NUM_GENE_TREES] [NUM_REPLICATES] [DUP_RATE] [LOSS_RATE]
#
# ARGUMENTS:
#   1. N_TREES         : Number of trees to process in this batch
#   2. MUL_TREES_DIR   : Directory containing mul_tree_NNNN.nex files
#   3. START_INDEX     : First tree index (e.g., 0, 1000)
#   4. CONFIG_NAME     : Configuration name (default: "ils_low")
#   5. NE              : Effective population size (default: 200000)
#   6. NUM_GENE_TREES  : Number of gene trees per replicate (default: 500)
#   7. NUM_REPLICATES  : Number of replicates (default: 1)
#   8. DUP_RATE        : Gene duplication rate (default: 0)
#   9. LOSS_RATE       : Gene loss rate (default: 0)
#
# EXAMPLES:
#   ./submit_simphy_training.sh 2000 /path/to/mul_trees 0                                    # low ILS, no dup/loss
#   ./submit_simphy_training.sh 2000 /path/to/mul_trees 0 ils_medium 1000000                  # medium ILS
#   ./submit_simphy_training.sh 2000 /path/to/mul_trees 0 dup_loss_low 200000 500 1 1e-9 1e-9 # low dup/loss
# ==============================================================================

N_TREES=${1:?Usage: submit_simphy_training.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [NE] [NUM_GENE_TREES] [NUM_REPLICATES] [DUP_RATE] [LOSS_RATE]}
MUL_TREES_DIR=${2:?Usage: submit_simphy_training.sh N_TREES MUL_TREES_DIR START_INDEX [CONFIG_NAME] [NE] [NUM_GENE_TREES] [NUM_REPLICATES] [DUP_RATE] [LOSS_RATE]}
START_INDEX=${3:-0}
CONFIG_NAME=${4:-ils_low}
NE_VAL=${5:-200000}
NUM_GENE_TREES=${6:-500}
NUM_REPLICATES=${7:-1}
DUP_RATE_VAL=${8:-0}
LOSS_RATE_VAL=${9:-0}

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
LOG_DIR="${BASE_DIR}/logs"

mkdir -p "$LOG_DIR"

echo "================================================"
echo "SimPhy gene trees for ML training"
echo "================================================"
echo "MUL-trees dir: ${MUL_TREES_DIR}"
echo "Batch: ${N_TREES} trees, indices ${START_INDEX}-$((START_INDEX + N_TREES - 1))"
echo "Config: ${CONFIG_NAME} | Ne: ${NE_VAL}"
echo "Dup rate: ${DUP_RATE_VAL} | Loss rate: ${LOSS_RATE_VAL}"
echo "Gene trees per replicate: ${NUM_GENE_TREES}"
echo "Replicates: ${NUM_REPLICATES}"
echo "Logs: ${LOG_DIR}/simphy_train_${CONFIG_NAME}_*.out"
echo "================================================"

sbatch \
    --array="0-$((N_TREES - 1))" \
    --job-name="simphy_train_${CONFIG_NAME}" \
    --output="${LOG_DIR}/simphy_train_${CONFIG_NAME}_%a.out" \
    --error="${LOG_DIR}/simphy_train_${CONFIG_NAME}_%a.err" \
    --time=04:00:00 \
    --mem=16G \
    --cpus-per-task=1 \
    --partition=itaym-pool \
    --account=itaym-users_v2 \
    --qos=owner \
    --export=ALL,MUL_TREES_DIR="${MUL_TREES_DIR}",CONFIG_NAME="${CONFIG_NAME}",NE_VAL="${NE_VAL}",NUM_GENE_TREES="${NUM_GENE_TREES}",NUM_REPLICATES="${NUM_REPLICATES}",INDEX_OFFSET="${START_INDEX}",DUP_RATE="${DUP_RATE_VAL}",LOSS_RATE_VAL="${LOSS_RATE_VAL}" \
    "${BASE_DIR}/jobs/simphy_training_job.sh"
