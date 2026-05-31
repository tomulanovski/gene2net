#!/bin/bash

# ============================================================================
# SUBMIT_GRANDMA_SPLIT_PRIOR_REAL_ALL_CONFIGS.SH
# ============================================================================
# Submits run_grandma_split_prior_real.sh (GRANDMA_SPLIT with the REAL
# ground-truth ploidy prior) for all 9 dup_loss configurations.
# Each submission is a SLURM array of 105 jobs (21 networks x 5 replicates).
#
# The real ploidy prior is derived in-job from each network's MUL-tree (and
# renamed via the per-replicate taxa_map.txt). Results go to a separate
# grandma_split_prior_real/ dir so the inferred-prior results are left intact.
#
# Usage:
#   ./submit_grandma_split_prior_real_all_configs.sh [--dry-run]
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="${SCRIPT_DIR}/run_grandma_split_prior_real.sh"
LOG_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"

CONFIGS=(
    "conf_dup_loss_low_10M"
    "conf_dup_loss_medium_10M"
    "conf_dup_loss_high_10M"
    "conf_dup_loss_low_10M_ne1M"
    "conf_dup_loss_medium_10M_ne1M"
    "conf_dup_loss_high_10M_ne1M"
    "conf_dup_loss_low_10M_ne2M"
    "conf_dup_loss_medium_10M_ne2M"
    "conf_dup_loss_high_10M_ne2M"
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

echo "Submitting grandma_split_prior_real for ${#CONFIGS[@]} configurations..."
echo ""

for config in "${CONFIGS[@]}"; do
    echo -n "  ${config} ... "
    CMD="sbatch --array=1-105 \
        --export=CONFIG=${config} \
        --output=${LOG_DIR}/run_grandma_split_prior_real_${config}_%A_%a.out \
        --error=${LOG_DIR}/run_grandma_split_prior_real_${config}_%A_%a.err \
        ${RUN_SCRIPT}"
    if [ "$DRY_RUN" = true ]; then
        echo "[dry run] ${CMD}"
    else
        job_id=$(eval "$CMD" | grep -oP '\d+$')
        echo "submitted (job ID: ${job_id})"
    fi
done

echo ""
echo "Done. Monitor with: squeue -u \$USER"
