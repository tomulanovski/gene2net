#!/bin/bash

# ============================================================================
# GRANDMA_SPLIT PIPELINE HELPER
# ============================================================================
# This script provides easy submission of the GRANDMA_SPLIT pipeline.
#
# Note: GRANDMA_SPLIT uses the same input as GRAMPA (gene trees and species tree)
# So it requires GRAMPA prep and ASTRAL to be run first.
#
# Pipeline dependencies:
#   1. prep_grampa.sh    - Process gene trees (same as GRAMPA)
#   2. run_astral.sh     - Run ASTRAL to create species tree (same as GRAMPA)
#   3. run_grandma_split.sh - Run GRANDMA_SPLIT
#
# Usage:
#   ./submit_grandma_split_pipeline.sh <config_name> [options]
#
# Examples:
#   # Run GRANDMA_SPLIT (assumes GRAMPA prep and ASTRAL already done)
#   ./submit_grandma_split_pipeline.sh conf_ils_low_10M
#
#   # Run only prep step (delegates to GRAMPA prep)
#   ./submit_grandma_split_pipeline.sh conf_ils_low_10M --prep-only
#
#   # Run only GRANDMA_SPLIT (assumes prep and ASTRAL done)
#   ./submit_grandma_split_pipeline.sh conf_ils_low_10M --run-only
#
#   # Dry run
#   ./submit_grandma_split_pipeline.sh conf_ils_low_10M --dry-run
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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_SCRIPT="${SCRIPT_DIR}/prep_grampa.sh"
ASTRAL_SCRIPT="${SCRIPT_DIR}/run_astral.sh"
RUN_SCRIPT="${SCRIPT_DIR}/run_grandma_split.sh"

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
  --prep-only       Only run the preparation step (GRAMPA prep + ASTRAL)
  --run-only        Only run GRANDMA_SPLIT (assumes prep and ASTRAL done)
  --dry-run         Show commands without executing
  -h, --help        Show this help message

Pipeline steps:
  1. prep_grampa.sh       - Process gene trees (21 jobs, one per network)
  2. run_astral.sh        - Run ASTRAL (21 jobs, one per network)
  3. run_grandma_split.sh - Run GRANDMA_SPLIT (105 jobs, one per network×replicate)

Note: GRANDMA_SPLIT uses the same prep and ASTRAL steps as GRAMPA.
      Make sure these are completed before running GRANDMA_SPLIT.

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

# Validate conflicting options
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
echo "GRANDMA_SPLIT PIPELINE"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Replicates: ${NUM_REPLICATES}"
echo "Networks: ${NUM_NETWORKS}"
echo "Total GRANDMA_SPLIT jobs: ${ARRAY_SIZE}"
echo ""
echo "Steps to run:"

if [ "$PREP_ONLY" = true ]; then
    echo "  [x] Prep (process gene trees - same as GRAMPA)"
    echo "  [x] ASTRAL"
    echo "  [ ] GRANDMA_SPLIT (skipped - --prep-only)"
elif [ "$RUN_ONLY" = true ]; then
    echo "  [ ] Prep (skipped - --run-only)"
    echo "  [ ] ASTRAL (skipped - --run-only)"
    echo "  [x] GRANDMA_SPLIT"
else
    echo "  [x] Prep (process gene trees - same as GRAMPA)"
    echo "  [x] ASTRAL (with dependency)"
    echo "  [x] GRANDMA_SPLIT (with dependency)"
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
ASTRAL_JOB_ID=""
RUN_JOB_ID=""

# Log directory
LOG_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"

# Determine which steps to run
DO_PREP=true
DO_ASTRAL=true
DO_RUN=true

if [ "$PREP_ONLY" = true ]; then
    DO_RUN=false
elif [ "$RUN_ONLY" = true ]; then
    DO_PREP=false
    DO_ASTRAL=false
fi

# Submit prep job (using GRAMPA prep script)
if [ "$DO_PREP" = true ]; then
    echo "Submitting prep job (using GRAMPA prep)..."
    
    PREP_CMD="sbatch --export=CONFIG=${CONFIG},NUM_REPLICATES=${NUM_REPLICATES} \
        --output=${LOG_DIR}/prep_grampa_${CONFIG}_%a.out \
        --error=${LOG_DIR}/prep_grampa_${CONFIG}_%a.err \
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

# Submit ASTRAL job
if [ "$DO_ASTRAL" = true ]; then
    echo "Submitting ASTRAL job..."
    
    if [ -n "$PREP_JOB_ID" ] && [ "$PREP_JOB_ID" != "XXXXX" ]; then
        # Add dependency on prep job
        ASTRAL_CMD="sbatch --dependency=afterok:${PREP_JOB_ID} \
            --export=CONFIG=${CONFIG},NUM_REPLICATES=${NUM_REPLICATES} \
            --output=${LOG_DIR}/astral_grampa_${CONFIG}_%a.out \
            --error=${LOG_DIR}/astral_grampa_${CONFIG}_%a.err \
            ${ASTRAL_SCRIPT}"
    else
        ASTRAL_CMD="sbatch --export=CONFIG=${CONFIG},NUM_REPLICATES=${NUM_REPLICATES} \
            --output=${LOG_DIR}/astral_grampa_${CONFIG}_%a.out \
            --error=${LOG_DIR}/astral_grampa_${CONFIG}_%a.err \
            ${ASTRAL_SCRIPT}"
    fi
    
    if [ "$DRY_RUN" = true ]; then
        run_cmd "$ASTRAL_CMD"
        ASTRAL_JOB_ID="YYYYY"
    else
        ASTRAL_OUTPUT=$(run_cmd "$ASTRAL_CMD")
        ASTRAL_JOB_ID=$(echo "$ASTRAL_OUTPUT" | grep -oP '\d+$')
        echo "  ASTRAL job submitted: ${ASTRAL_JOB_ID}"
    fi
    echo ""
fi

# Submit GRANDMA_SPLIT run job
if [ "$DO_RUN" = true ]; then
    echo "Submitting GRANDMA_SPLIT run job..."
    
    EXPORT_VARS="CONFIG=${CONFIG},NUM_REPLICATES=${NUM_REPLICATES}"
    
    # Determine dependency
    DEPENDENCY=""
    if [ -n "$ASTRAL_JOB_ID" ] && [ "$ASTRAL_JOB_ID" != "YYYYY" ]; then
        DEPENDENCY="--dependency=afterok:${ASTRAL_JOB_ID}"
    elif [ -n "$PREP_JOB_ID" ] && [ "$PREP_JOB_ID" != "XXXXX" ]; then
        DEPENDENCY="--dependency=afterok:${PREP_JOB_ID}"
    fi
    
    RUN_CMD="sbatch --array=1-${ARRAY_SIZE} ${DEPENDENCY} \
        --export=${EXPORT_VARS} \
        --output=${LOG_DIR}/run_grandma_split_${CONFIG}_%a.out \
        --error=${LOG_DIR}/run_grandma_split_${CONFIG}_%a.err \
        ${RUN_SCRIPT}"
    
    if [ "$DRY_RUN" = true ]; then
        run_cmd "$RUN_CMD"
        RUN_JOB_ID="ZZZZZ"
    else
        RUN_OUTPUT=$(run_cmd "$RUN_CMD")
        RUN_JOB_ID=$(echo "$RUN_OUTPUT" | grep -oP '\d+$')
        echo "  GRANDMA_SPLIT run job submitted: ${RUN_JOB_ID}"
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
    echo "  - Output: clean_trees.tre, grampa_trees.tre, taxa_map.txt"
fi

if [ -n "$ASTRAL_JOB_ID" ]; then
    echo "ASTRAL job ID: ${ASTRAL_JOB_ID}"
    echo "  - 21 array tasks (one per network)"
    echo "  - Will process ${NUM_REPLICATES} replicates each"
    echo "  - Output: species.tre"
    if [ -n "$PREP_JOB_ID" ] && [ "$PREP_JOB_ID" != "XXXXX" ]; then
        echo "  - Depends on prep job ${PREP_JOB_ID}"
    fi
fi

if [ -n "$RUN_JOB_ID" ]; then
    echo "GRANDMA_SPLIT run job ID: ${RUN_JOB_ID}"
    echo "  - ${ARRAY_SIZE} array tasks (${NUM_NETWORKS} networks × ${NUM_REPLICATES} replicates)"
    if [ -n "$ASTRAL_JOB_ID" ] && [ "$ASTRAL_JOB_ID" != "YYYYY" ]; then
        echo "  - Depends on ASTRAL job ${ASTRAL_JOB_ID}"
    fi
fi

echo ""
echo "Monitor with:"
[ -n "$PREP_JOB_ID" ] && echo "  squeue -j ${PREP_JOB_ID}"
[ -n "$ASTRAL_JOB_ID" ] && echo "  squeue -j ${ASTRAL_JOB_ID}"
[ -n "$RUN_JOB_ID" ] && echo "  squeue -j ${RUN_JOB_ID}"
echo ""
echo "Check logs in:"
echo "  ${LOG_DIR}/"
echo "============================================================================"
