#!/bin/bash
# ==============================================================================
# SCRIPT: submit_simphy.sh
# PURPOSE: Submit SimPhy simulations with dynamic logging and parameters.
#
# USAGE: 
#   ./submit_simphy.sh <CONFIG_NAME> [TREE_HEIGHT] [DUP_RATE] [LOSS_RATE] [MANUAL_NE]
#
# ARGUMENTS:
#   1. CONFIG_NAME : Name of the run (e.g., "conf_ils_low_10M").
#                    Used for log filenames and auto-detecting Ne.
#   2. TREE_HEIGHT : "10M" or "50M" (Default: 10M).
#   3. DUP_RATE    : Duplication rate per MY (Default: 0).
#   4. LOSS_RATE   : Loss rate per MY (Default: 0).
#   5. MANUAL_NE   : (Optional) Override population size. 
#                    If empty, defaults to: High=2M, Med=1M, Low=200k.
#
# EXAMPLES:
#   1. Standard Low ILS (10M tree, no dup/loss, Ne=10,000):
#      ./submit_simphy.sh conf_ils_low_10M
#
#   2. Medium ILS (10M tree, no dup/loss, Ne=200,000):
#      ./submit_simphy.sh conf_ils_medium_10M
#
#   3. High ILS with 50M Tree (Ne=2,000,000):
#      ./submit_simphy.sh conf_ils_high_50M 50M
#
#   4. Low ILS with Duplication & Loss (Rate 0.001):
#      ./submit_simphy.sh conf_dup_loss_low 10M 0.001 0.001
#
#   5. Custom Ne (Force Ne to 500,000):
#      ./submit_simphy.sh conf_custom_test 10M 0 0 500000
# ==============================================================================

# 1. READ ARGUMENTS
CONFIG_NAME=$1
TREE_HEIGHT=${2:-10M}   # Default: 10M
DUP_RATE=${3:-0}        # Default: 0
LOSS_RATE=${4:-0}       # Default: 0
MANUAL_NE=${5}          # Optional: Manual Ne

if [ -z "$CONFIG_NAME" ]; then
    echo "Error: Configuration name is required."
    echo "Usage: ./submit_simphy.sh <CONFIG_NAME> [TREE_HEIGHT] [DUP_RATE] [LOSS_RATE] [NE]"
    exit 1
fi

# 2. DETERMINE POPULATION SIZE (NE)
if [ -n "$MANUAL_NE" ]; then
    # Use the manually provided Ne
    NE_VAL="$MANUAL_NE"
    echo "Using manual Ne: ${NE_VAL}"
else
    # Auto-detect based on name
    if [[ "$CONFIG_NAME" == *"ils_high"* ]]; then
        NE_VAL=2000000
        echo "Auto-detected 'high' ILS -> Ne = 2,000,000"
    elif [[ "$CONFIG_NAME" == *"ils_medium"* ]]; then
        NE_VAL=1000000 
        echo "Auto-detected 'medium' ILS -> Ne = 1,000,000"
    else 
        NE_VAL=200000
        echo "Auto-detected 'low' (or default) ILS -> Ne = 200,000"
    fi
fi

# 3. EXPORT VARIABLES
export SIMPHY_CONFIG="$CONFIG_NAME"
export SIMPHY_TREE_HEIGHT="$TREE_HEIGHT"
export SIMPHY_NE="$NE_VAL"
export SIMPHY_DUP="$DUP_RATE"
export SIMPHY_LOSS="$LOSS_RATE"

# 4. PATHS
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations"
LOG_DIR="${BASE_DIR}/logs"
SCRIPT_PATH="${BASE_DIR}/jobs/simphy_reusable.sh" 

# 5. SUBMIT
echo "------------------------------------------------"
echo "Submitting: ${CONFIG_NAME}"
echo "Tree: ${TREE_HEIGHT} | Ne: ${NE_VAL}"
echo "Dup: ${DUP_RATE} | Loss: ${LOSS_RATE}"
echo "Logs: ${LOG_DIR}/simphy_${CONFIG_NAME}_%a.out"
echo "------------------------------------------------"

sbatch \
    --job-name="simphy_${CONFIG_NAME}" \
    --output="${LOG_DIR}/simphy_${CONFIG_NAME}_%a.out" \
    --error="${LOG_DIR}/simphy_${CONFIG_NAME}_%a.err" \
    "$SCRIPT_PATH"