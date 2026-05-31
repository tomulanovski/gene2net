#!/bin/bash

# ============================================================================
# SUBMIT_POLYPHEST_REAL_PIPELINE.SH
# ============================================================================
# Submits run_polyphest_real.sh (Polyphest with the REAL ground-truth ploidy
# prior) for all 9 dup_loss configurations x percentiles {50, 70, 90}.
# Each submission is a SLURM array of 105 jobs (21 networks x 5 replicates).
#
# The real multiset is generated in-job from each network's MUL-tree, and
# results are written to results/<config>/polyphest_real_p<PCT>/ so the
# existing inferred-prior results (polyphest_p<PCT>/) are left intact.
#
# Usage:
#   ./submit_polyphest_real_pipeline.sh [--dry-run]
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="${SCRIPT_DIR}/run_polyphest_real.sh"
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

PERCENTILES=(50 70 90)

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

echo "Submitting polyphest_real for ${#CONFIGS[@]} configs x ${#PERCENTILES[@]} percentiles..."
echo ""

for config in "${CONFIGS[@]}"; do
    for pct in "${PERCENTILES[@]}"; do
        echo -n "  ${config} p${pct} ... "
        CMD="sbatch --array=1-105 \
            --export=CONFIG=${config},PERCENTILE=${pct} \
            --output=${LOG_DIR}/run_polyphest_real_${config}_p${pct}_%A_%a.out \
            --error=${LOG_DIR}/run_polyphest_real_${config}_p${pct}_%A_%a.err \
            ${RUN_SCRIPT}"
        if [ "$DRY_RUN" = true ]; then
            echo "[dry run] ${CMD}"
        else
            job_id=$(eval "$CMD" | grep -oP '\d+$')
            echo "submitted (job ID: ${job_id})"
        fi
    done
done

echo ""
echo "Done. Monitor with: squeue -u \$USER"
