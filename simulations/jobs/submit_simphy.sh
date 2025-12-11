#!/bin/bash
# ==============================================================================
# SCRIPT: submit_simphy.sh
# PURPOSE: Submit SimPhy simulations with dynamic logging and parameters.
#
# USAGE: 
#   ./submit_simphy.sh <CONFIG_NAME> [TREE_HEIGHT] [DUP_RATE] [LOSS_RATE] [NE] [ARRAY_RANGE]
#
# ARGUMENTS:
#   1. CONFIG_NAME : Name of the run (e.g., "conf_dup_loss_low").
#   2. TREE_HEIGHT : "10M" or "50M" (Default: 10M).
#   3. DUP_RATE    : Duplication rate per MY (Default: 0).
#   4. LOSS_RATE   : Loss rate per MY (Default: 0).
#   5. NE          : (Optional) Population size (Default: 200000).
#   6. ARRAY_RANGE : (Optional) Array indices to run. Default: "1-21".
#
# EXAMPLES:
#
#   --- ILS only (no dup/loss, varying Ne) ---
#   ./submit_simphy.sh conf_ils_low 10M 0 0                      # Ne=200,000 (default)
#   ./submit_simphy.sh conf_ils_medium 10M 0 0 1000000           # Ne=1,000,000
#   ./submit_simphy.sh conf_ils_high 10M 0 0 2000000             # Ne=2,000,000
#
#   --- Dup/Loss only (default Ne) ---
#   ./submit_simphy.sh conf_dup_loss_low 10M 0.001 0.001         # Low dup/loss rates
#   ./submit_simphy.sh conf_dup_loss_medium 10M 0.005 0.005      # Medium dup/loss rates
#   ./submit_simphy.sh conf_dup_loss_high 10M 0.01 0.01          # High dup/loss rates
#
#   --- ILS + Dup/Loss combined ---
#   ./submit_simphy.sh conf_ils_low_dup_loss_low 10M 0.001 0.001 200000
#   ./submit_simphy.sh conf_ils_medium_dup_loss_low 10M 0.001 0.001 1000000
#   ./submit_simphy.sh conf_ils_high_dup_loss_low 10M 0.001 0.001 2000000
#   ./submit_simphy.sh conf_ils_high_dup_loss_high 10M 0.01 0.01 2000000
#
#   --- Running specific networks ---
#   ./submit_simphy.sh conf_ils_low 10M 0 0 200000 1,13,14       # Only networks 1, 13, 14
#   ./submit_simphy.sh conf_ils_low 10M 0 0 200000 1-5           # Networks 1 through 5
#
#   --- Different tree heights ---
#   ./submit_simphy.sh conf_ils_low_50M 50M 0 0                  # 50M tree height
# ==============================================================================

# 1. READ ARGUMENTS
CONFIG_NAME=$1
TREE_HEIGHT=${2:-10M}       # Default: 10M
DUP_RATE=${3:-0}            # Default: 0
LOSS_RATE=${4:-0}           # Default: 0
NE_VAL=${5:-200000}         # Default: 200,000
ARRAY_RANGE=${6:-1-21}      # Default: 1-21 (All networks)

if [ -z "$CONFIG_NAME" ]; then
    echo "Error: Configuration name is required."
    echo "Usage: ./submit_simphy.sh <CONFIG_NAME> [TREE_HEIGHT] [DUP_RATE] [LOSS_RATE] [NE] [ARRAY]"
    exit 1
fi

echo "Using Ne: ${NE_VAL}"

# 2. EXPORT VARIABLES
export SIMPHY_CONFIG="$CONFIG_NAME"
export SIMPHY_TREE_HEIGHT="$TREE_HEIGHT"
export SIMPHY_NE="$NE_VAL"
export SIMPHY_DUP="$DUP_RATE"
export SIMPHY_LOSS="$LOSS_RATE"

# 3. PATHS
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations"
LOG_DIR="${BASE_DIR}/logs"
SCRIPT_PATH="${BASE_DIR}/jobs/simphy_reusable.sh" 

# 4. SUBMIT
echo "------------------------------------------------"
echo "Submitting: ${CONFIG_NAME}"
echo "Array Range: ${ARRAY_RANGE}"
echo "Tree: ${TREE_HEIGHT} | Ne: ${NE_VAL}"
echo "Dup: ${DUP_RATE} | Loss: ${LOSS_RATE}"
echo "Logs: ${LOG_DIR}/simphy_${CONFIG_NAME}_%a.out"
echo "------------------------------------------------"

sbatch \
    --array="$ARRAY_RANGE" \
    --job-name="simphy_${CONFIG_NAME}" \
    --output="${LOG_DIR}/simphy_${CONFIG_NAME}_%a.out" \
    --error="${LOG_DIR}/simphy_${CONFIG_NAME}_%a.err" \
    "$SCRIPT_PATH"