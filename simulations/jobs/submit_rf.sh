#!/bin/bash
# ==============================================================================
# SCRIPT: submit_rf.sh
# PURPOSE: Submit RF Analysis & Aggregation for a specific configuration.
#
# USAGE: 
#   ./submit_rf.sh <CONFIG_NAME>
#
# EXAMPLES:
#   ./submit_rf.sh conf_ils_low_10M
#   ./submit_rf.sh conf_dup_loss_ils_high
# ==============================================================================

# 1. READ ARGUMENTS
CONFIG_NAME=$1

if [ -z "$CONFIG_NAME" ]; then
    echo "Error: Configuration name is required."
    echo "Usage: ./submit_rf.sh <CONFIG_NAME>"
    exit 1
fi

# 2. EXPORT VARIABLES
# This variable is read by the calculation script
export RF_CONFIG="$CONFIG_NAME"

# 3. PATHS
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations"
LOG_DIR="${BASE_DIR}/logs"
SCRIPT_PATH="${BASE_DIR}/jobs/calculate_rf_reusable.sh" 

# 4. SUBMIT
echo "------------------------------------------------"
echo "Submitting RF Analysis for: ${CONFIG_NAME}"
echo "Log File: ${LOG_DIR}/rf_${CONFIG_NAME}_%j.out"
echo "------------------------------------------------"

# Note: We use %j (Job ID) because this is a single job, not an array.
sbatch \
    --job-name="rf_${CONFIG_NAME}" \
    --output="${LOG_DIR}/rf_${CONFIG_NAME}_%j.out" \
    --error="${LOG_DIR}/rf_${CONFIG_NAME}_%j.err" \
    "$SCRIPT_PATH"