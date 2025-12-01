#!/bin/bash
# ==============================================================================
# SCRIPT: submit_dup_loss_extraction.sh
# PURPOSE: Submit Duplication/Loss Extraction & Aggregation.
#
# USAGE: 
#   ./submit_dup_loss_extraction.sh <CONFIG_NAME>
#
# EXAMPLES:
#   ./submit_dup_loss_extraction.sh conf_ils_low_10M
#   ./submit_dup_loss_extraction.sh conf_dup_loss_low_50M
# ==============================================================================

# 1. READ ARGUMENTS
CONFIG_NAME=$1

if [ -z "$CONFIG_NAME" ]; then
    echo "Error: Configuration name is required."
    echo "Usage: ./submit_dup_loss_extraction.sh <CONFIG_NAME>"
    exit 1
fi

# 2. EXPORT VARIABLES
export EXTRACT_CONFIG="$CONFIG_NAME"

# 3. PATHS
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations"
LOG_DIR="${BASE_DIR}/logs"
SCRIPT_PATH="${BASE_DIR}/jobs/extract_dup_loss_reusable.sh" 

# 4. SUBMIT
echo "------------------------------------------------"
echo "Submitting Extraction for: ${CONFIG_NAME}"
echo "Log File: ${LOG_DIR}/extract_${CONFIG_NAME}_%j.out"
echo "------------------------------------------------"

sbatch \
    --job-name="extract_${CONFIG_NAME}" \
    --output="${LOG_DIR}/extract_${CONFIG_NAME}_%j.out" \
    --error="${LOG_DIR}/extract_${CONFIG_NAME}_%j.err" \
    "$SCRIPT_PATH"