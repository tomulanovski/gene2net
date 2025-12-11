#!/bin/bash

# ============================================================================
# PADRE PIPELINE HELPER
# ============================================================================
# This script provides easy submission of the PADRE pipeline.
#
# Usage:
#   ./submit_padre_pipeline.sh <config_name> [options]
#
# Examples:
#   # Run both prep and PADRE for a configuration
#   ./submit_padre_pipeline.sh conf_ils_low_10M
#
#   # Run only prep step
#   ./submit_padre_pipeline.sh conf_ils_low_10M --prep-only
#
#   # Run only PADRE step (assumes prep already done)
#   ./submit_padre_pipeline.sh conf_ils_low_10M --run-only
#
#   # With custom number of replicates
#   ./submit_padre_pipeline.sh conf_ils_low_10M --replicates 3
#
#   # Dry run (show commands without executing)
#   ./submit_padre_pipeline.sh conf_ils_low_10M --dry-run
# ============================================================================

set -euo pipefail

# ============================================================================
# DEFAULTS
# ============================================================================

NUM_REPLICATES=5
NUM_NETWORKS=21
PREP_ONLY=false
RUN_ONLY=false
DRY_RUN=false
JAVA_MEM="4g"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_SCRIPT="${SCRIPT_DIR}/prep_padre.sh"
RUN_SCRIPT="${SCRIPT_DIR}/run_padre.sh"

# ============================================================================
# USAGE
# ============================================================================

usage() {
    cat << EOF
Usage: $0 <config_name> [options]

Arguments:
  config_name       Configuration name (e.g., conf_ils_low_10M)

Options:
  --replicates N    Number of replicates (default: 5)
  --prep-only       Only run the preparation step
  --run-only        Only run PADRE (assumes prep is done)
  --java-mem SIZE   Java heap size (default: 4g)
  --dry-run         Show commands without executing
  -h, --help        Show this help message

Examples:
  $0 conf_ils_low_10M
  $0 conf_ils_medium_10M --replicates 3
  $0 conf_ils_high_10M --run-only
  $0 conf_ils_low_10M --dry-run

Configuration names from your simulations:
  - conf_ils_low_10M
  - conf_ils_medium_10M
  - conf_ils_high_10M
EOF
    exit 1
}

# ============================================================================
# PARSE ARGUMENTS
# ============================================================================

if [ $# -lt 1 ]; then
    usage
fi

CONFIG="$1"
shift

while [ $# -gt 0 ]; do
    case "$1" in
        --replicates)
            NUM_REPLICATES="$2"
            shift 2
            ;;
        --prep-only)
            PREP_ONLY=true
            shift
            ;;
        --run-only)
            RUN_ONLY=true
            shift
            ;;
        --java-mem)
            JAVA_MEM="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            usage
            ;;
    esac
done

# Validate
if [ "$PREP_ONLY" = true ] && [ "$RUN_ONLY" = true ]; then
    echo "ERROR: Cannot specify both --prep-only and --run-only"
    exit 1
fi

# Calculate array size for run step
ARRAY_SIZE=$((NUM_NETWORKS * NUM_REPLICATES))

# ============================================================================
# DISPLAY CONFIGURATION
# ============================================================================

echo "============================================================================"
echo "PADRE PIPELINE"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Replicates: ${NUM_REPLICATES}"
echo "Networks: ${NUM_NETWORKS}"
echo "Total PADRE jobs: ${ARRAY_SIZE}"
echo ""
echo "PADRE parameters:"
echo "  Java heap size: ${JAVA_MEM}"
echo ""
echo "Steps to run:"
if [ "$RUN_ONLY" = true ]; then
    echo "  [ ] Prep (skipped - --run-only)"
    echo "  [x] Run PADRE"
elif [ "$PREP_ONLY" = true ]; then
    echo "  [x] Prep"
    echo "  [ ] Run PADRE (skipped - --prep-only)"
else
    echo "  [x] Prep"
    echo "  [x] Run PADRE (with dependency)"
fi
echo ""
if [ "$DRY_RUN" = true ]; then
    echo "MODE: DRY RUN (commands will be shown but not executed)"
fi
echo "============================================================================"
echo ""

# ============================================================================
# SUBMIT JOBS
# ============================================================================

run_cmd() {
    if [ "$DRY_RUN" = true ]; then
        echo "[DRY RUN] Would execute: $*"
    else
        echo "Executing: $*"
        eval "$@"
    fi
}

PREP_JOB_ID=""
RUN_JOB_ID=""

# Log directory
LOG_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"

# Submit prep job
if [ "$RUN_ONLY" = false ]; then
    echo "Submitting prep job..."
    
    PREP_CMD="sbatch --export=CONFIG=${CONFIG},NUM_REPLICATES=${NUM_REPLICATES} \
        --output=${LOG_DIR}/prep_padre_${CONFIG}_%A_%a.out \
        --error=${LOG_DIR}/prep_padre_${CONFIG}_%A_%a.err \
        ${PREP_SCRIPT}"
    
    if [ "$DRY_RUN" = true ]; then
        run_cmd "$PREP_CMD"
        PREP_JOB_ID="XXXXX"
    else
        PREP_OUTPUT=$(run_cmd "$PREP_CMD")
        PREP_JOB_ID=$(echo "$PREP_OUTPUT" | grep -oP '\d+$')
        echo "  Prep job submitted: ${PREP_JOB_ID}"
    fi
    echo ""
fi

# Submit run job
if [ "$PREP_ONLY" = false ]; then
    echo "Submitting PADRE run job..."
    
    EXPORT_VARS="CONFIG=${CONFIG},NUM_REPLICATES=${NUM_REPLICATES},JAVA_MEM=${JAVA_MEM}"
    
    if [ -n "$PREP_JOB_ID" ] && [ "$PREP_JOB_ID" != "XXXXX" ]; then
        # Add dependency on prep job
        RUN_CMD="sbatch --array=1-${ARRAY_SIZE} --dependency=afterok:${PREP_JOB_ID} \
            --export=${EXPORT_VARS} \
            --output=${LOG_DIR}/run_padre_${CONFIG}_%A_%a.out \
            --error=${LOG_DIR}/run_padre_${CONFIG}_%A_%a.err \
            ${RUN_SCRIPT}"
    else
        RUN_CMD="sbatch --array=1-${ARRAY_SIZE} \
            --export=${EXPORT_VARS} \
            --output=${LOG_DIR}/run_padre_${CONFIG}_%A_%a.out \
            --error=${LOG_DIR}/run_padre_${CONFIG}_%A_%a.err \
            ${RUN_SCRIPT}"
    fi
    
    if [ "$DRY_RUN" = true ]; then
        run_cmd "$RUN_CMD"
        RUN_JOB_ID="YYYYY"
    else
        RUN_OUTPUT=$(run_cmd "$RUN_CMD")
        RUN_JOB_ID=$(echo "$RUN_OUTPUT" | grep -oP '\d+$')
        echo "  PADRE run job submitted: ${RUN_JOB_ID}"
    fi
    echo ""
fi

# ============================================================================
# SUMMARY
# ============================================================================

echo "============================================================================"
echo "SUBMISSION SUMMARY"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo ""

if [ -n "$PREP_JOB_ID" ]; then
    echo "Prep job ID: ${PREP_JOB_ID}"
    echo "  - 21 array tasks (one per network)"
    echo "  - Will process ${NUM_REPLICATES} replicates each"
fi

if [ -n "$RUN_JOB_ID" ]; then
    echo "Run job ID: ${RUN_JOB_ID}"
    echo "  - ${ARRAY_SIZE} array tasks (${NUM_NETWORKS} networks ª ${NUM_REPLICATES} replicates)"
    if [ -n "$PREP_JOB_ID" ] && [ "$PREP_JOB_ID" != "XXXXX" ]; then
        echo "  - Depends on prep job ${PREP_JOB_ID}"
    fi
fi

echo ""
echo "Monitor with:"
if [ -n "$PREP_JOB_ID" ]; then
    echo "  squeue -j ${PREP_JOB_ID}"
fi
if [ -n "$RUN_JOB_ID" ]; then
    echo "  squeue -j ${RUN_JOB_ID}"
fi
echo ""
echo "Check logs in:"
echo "  ${LOG_DIR}/"
echo "============================================================================"