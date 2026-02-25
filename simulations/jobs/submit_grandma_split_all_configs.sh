#!/bin/bash

# ============================================================================
# SUBMIT_GRANDMA_SPLIT_ALL_CONFIGS.SH
# ============================================================================
# Submits run_grandma_split.sh for all 9 configurations.
# Each submission is a SLURM array of 105 jobs (21 networks × 5 replicates).
#
# Usage:
#   ./submit_grandma_split_all_configs.sh [--dry-run]
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="${SCRIPT_DIR}/run_grandma_split.sh"

CONFIGS=(
    "conf_ils_low_10M"
    "conf_ils_medium_10M"
    "conf_ils_high_10M"
    "conf_dup_loss_low_10M"
    "conf_dup_loss_medium_10M"
    "conf_dup_loss_high_10M"
    "conf_dup_loss_low_10M_ne1M"
    "conf_dup_loss_medium_10M_ne1M"
    "conf_dup_loss_high_10M_ne1M"
)

DRY_RUN=false
if [ "${1:-}" = "--dry-run" ]; then
    DRY_RUN=true
    echo "MODE: DRY RUN (no jobs will be submitted)"
    echo ""
fi

if [ ! -f "$RUN_SCRIPT" ]; then
    echo "ERROR: Run script not found: $RUN_SCRIPT"
    exit 1
fi

echo "Submitting grandma_split for ${#CONFIGS[@]} configurations..."
echo ""

for config in "${CONFIGS[@]}"; do
    echo -n "  ${config} ... "
    if [ "$DRY_RUN" = true ]; then
        echo "[dry run] sbatch --export=CONFIG=${config} ${RUN_SCRIPT}"
    else
        job_id=$(sbatch --export=CONFIG="${config}" "$RUN_SCRIPT" | grep -oP '\d+')
        echo "submitted (job ID: ${job_id})"
    fi
done

echo ""
echo "Done. Monitor with: squeue -u \$USER"
