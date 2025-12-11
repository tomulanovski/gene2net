#!/bin/bash

# ============================================================================
# GRAMPA PIPELINE HELPER
# ============================================================================
# This script provides easy submission of the GRAMPA pipeline.
#
# Pipeline steps:
#   1. prep_grampa.sh    - Process gene trees (clean, fix substrings, add copy numbers)
#   2. run_astral_grampa.sh - Run ASTRAL to create species tree
#   3. run_grampa.sh     - Run GRAMPA
#
# Usage:
#   ./submit_grampa_pipeline.sh <config_name> [options]
#
# Examples:
#   # Run full pipeline
#   ./submit_grampa_pipeline.sh conf_ils_low_10M
#
#   # Run only prep step
#   ./submit_grampa_pipeline.sh conf_ils_low_10M --prep-only
#
#   # Run only ASTRAL step (assumes prep done)
#   ./submit_grampa_pipeline.sh conf_ils_low_10M --astral-only
#
#   # Run only GRAMPA step (assumes prep and ASTRAL done)
#   ./submit_grampa_pipeline.sh conf_ils_low_10M --run-only
#
#   # Skip prep (run ASTRAL + GRAMPA)
#   ./submit_grampa_pipeline.sh conf_ils_low_10M --skip-prep
#
#   # Dry run
#   ./submit_grampa_pipeline.sh conf_ils_low_10M --dry-run
# ============================================================================

set -euo pipefail

# ============================================================================
# DEFAULTS
# ============================================================================

NUM_REPLICATES=5
NUM_NETWORKS=21
PREP_ONLY=false
ASTRAL_ONLY=false
RUN_ONLY=false
SKIP_PREP=false
DRY_RUN=false

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_SCRIPT="${SCRIPT_DIR}/prep_grampa.sh"
ASTRAL_SCRIPT="${SCRIPT_DIR}/run_astral.sh"
RUN_SCRIPT="${SCRIPT_DIR}/run_grampa.sh"

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
  --astral-only     Only run ASTRAL (assumes prep is done)
  --run-only        Only run GRAMPA (assumes prep and ASTRAL done)
  --skip-prep       Skip prep, run ASTRAL + GRAMPA
  --dry-run         Show commands without executing
  -h, --help        Show this help message

Pipeline steps:
  1. prep_grampa.sh       - Process gene trees (21 jobs, one per network)
  2. run_astral_grampa.sh - Run ASTRAL (21 jobs, one per network)
  3. run_grampa.sh        - Run GRAMPA (105 jobs, one per networkªreplicate)

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
        --astral-only)
            ASTRAL_ONLY=true
            shift
            ;;
        --run-only)
            RUN_ONLY=true
            shift
            ;;
        --skip-prep)
            SKIP_PREP=true
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
conflicting=0
[ "$PREP_ONLY" = true ] && conflicting=$((conflicting + 1))
[ "$ASTRAL_ONLY" = true ] && conflicting=$((conflicting + 1))
[ "$RUN_ONLY" = true ] && conflicting=$((conflicting + 1))
[ "$SKIP_PREP" = true ] && conflicting=$((conflicting + 1))

if [ $conflicting -gt 1 ]; then
    echo "ERROR: Cannot combine --prep-only, --astral-only, --run-only, or --skip-prep"
    exit 1
fi

# Calculate array size for run step
ARRAY_SIZE=$((NUM_NETWORKS * NUM_REPLICATES))

# ============================================================================
# DISPLAY CONFIGURATION
# ============================================================================

echo "============================================================================"
echo "GRAMPA PIPELINE"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Replicates: ${NUM_REPLICATES}"
echo "Networks: ${NUM_NETWORKS}"
echo "Total GRAMPA jobs: ${ARRAY_SIZE}"
echo ""
echo "Steps to run:"

if [ "$PREP_ONLY" = true ]; then
    echo "  [x] Prep (process gene trees)"
    echo "  [ ] ASTRAL (skipped - --prep-only)"
    echo "  [ ] GRAMPA (skipped - --prep-only)"
elif [ "$ASTRAL_ONLY" = true ]; then
    echo "  [ ] Prep (skipped - --astral-only)"
    echo "  [x] ASTRAL"
    echo "  [ ] GRAMPA (skipped - --astral-only)"
elif [ "$RUN_ONLY" = true ]; then
    echo "  [ ] Prep (skipped - --run-only)"
    echo "  [ ] ASTRAL (skipped - --run-only)"
    echo "  [x] GRAMPA"
elif [ "$SKIP_PREP" = true ]; then
    echo "  [ ] Prep (skipped - --skip-prep)"
    echo "  [x] ASTRAL"
    echo "  [x] GRAMPA (with dependency)"
else
    echo "  [x] Prep (process gene trees)"
    echo "  [x] ASTRAL (with dependency)"
    echo "  [x] GRAMPA (with dependency)"
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
    DO_ASTRAL=false
    DO_RUN=false
elif [ "$ASTRAL_ONLY" = true ]; then
    DO_PREP=false
    DO_RUN=false
elif [ "$RUN_ONLY" = true ]; then
    DO_PREP=false
    DO_ASTRAL=false
elif [ "$SKIP_PREP" = true ]; then
    DO_PREP=false
fi

# Submit prep job
if [ "$DO_PREP" = true ]; then
    echo "Submitting prep job..."
    
    PREP_CMD="sbatch --export=CONFIG=${CONFIG},NUM_REPLICATES=${NUM_REPLICATES} \
        --output=${LOG_DIR}/prep_grampa_${CONFIG}_%A_%a.out \
        --error=${LOG_DIR}/prep_grampa_${CONFIG}_%A_%a.err \
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
            --output=${LOG_DIR}/astral_grampa_${CONFIG}_%A_%a.out \
            --error=${LOG_DIR}/astral_grampa_${CONFIG}_%A_%a.err \
            ${ASTRAL_SCRIPT}"
    else
        ASTRAL_CMD="sbatch --export=CONFIG=${CONFIG},NUM_REPLICATES=${NUM_REPLICATES} \
            --output=${LOG_DIR}/astral_grampa_${CONFIG}_%A_%a.out \
            --error=${LOG_DIR}/astral_grampa_${CONFIG}_%A_%a.err \
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

# Submit GRAMPA run job
if [ "$DO_RUN" = true ]; then
    echo "Submitting GRAMPA run job..."
    
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
        --output=${LOG_DIR}/run_grampa_${CONFIG}_%A_%a.out \
        --error=${LOG_DIR}/run_grampa_${CONFIG}_%A_%a.err \
        ${RUN_SCRIPT}"
    
    if [ "$DRY_RUN" = true ]; then
        run_cmd "$RUN_CMD"
        RUN_JOB_ID="ZZZZZ"
    else
        RUN_OUTPUT=$(run_cmd "$RUN_CMD")
        RUN_JOB_ID=$(echo "$RUN_OUTPUT" | grep -oP '\d+$')
        echo "  GRAMPA run job submitted: ${RUN_JOB_ID}"
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
    echo "GRAMPA run job ID: ${RUN_JOB_ID}"
    echo "  - ${ARRAY_SIZE} array tasks (${NUM_NETWORKS} networks ª ${NUM_REPLICATES} replicates)"
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