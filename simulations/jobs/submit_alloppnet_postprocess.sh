#!/bin/bash

# ============================================================================
# ALLOPPNET POST-PROCESSING SUBMISSION HELPER
# ============================================================================
# This script submits AlloppNET post-processing jobs for all replicates.
#
# AlloppNET post-processing runs on 8 networks with max 2 copies:
#   Bendiksby_2011, Ding_2023, Koenen_2020, Liu_2023,
#   Shahrestani_2015, Wisecaver_2023, Wu_2015, Zhao_2021
#
# Each replicate submission creates 8 array tasks (one per network).
# Each array task runs:
#   1. TreeAnnotator on sampledmultrees.txt (10% burnin)
#   2. Remove copy number suffixes (_0, _1)
#   3. Create alloppnet_result.tre
#
# Usage:
#   ./submit_alloppnet_postprocess.sh <config_name> [options]
#
# Examples:
#   # Submit all 5 replicates (40 total post-processing jobs)
#   ./submit_alloppnet_postprocess.sh conf_ils_low_10M
#
#   # Submit specific replicates
#   ./submit_alloppnet_postprocess.sh conf_ils_low_10M --replicates 1,3,5
#
#   # Dry run (show commands without submitting)
#   ./submit_alloppnet_postprocess.sh conf_ils_low_10M --dry-run
# ============================================================================

set -euo pipefail

# ============================================================================
# DEFAULTS
# ============================================================================

NUM_NETWORKS=8
ALL_REPLICATES="1 2 3 4 5"
REPLICATES=""
DRY_RUN=false

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARRAY_SCRIPT="${SCRIPT_DIR}/alloppnet_postprocess_array.sh"
LOG_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"

# ============================================================================
# USAGE
# ============================================================================

usage() {
    cat << EOF
Usage: $0 <config_name> [options]

Arguments:
  config_name       Configuration name (e.g., conf_ils_low_10M)

Options:
  --replicates R    Comma-separated replicate numbers (default: 1,2,3,4,5)
                    Examples: --replicates 1  or  --replicates 1,3,5
  --dry-run         Show commands without executing
  -h, --help        Show this help message

AlloppNET-compatible networks (8 networks):
  - Bendiksby_2011
  - Ding_2023
  - Koenen_2020
  - Liu_2023
  - Shahrestani_2015
  - Wisecaver_2023
  - Wu_2015
  - Zhao_2021

Post-processing steps (per job):
  1. Count trees in sampledmultrees.txt
  2. Calculate 10% burnin
  3. Run TreeAnnotator to create consensus tree
  4. Remove copy number suffixes (_0, _1)
  5. Create alloppnet_result.tre

Job structure:
  - One array job per replicate
  - Each array job has 8 tasks (one per network)
  - Total jobs = replicates × 8 networks

Examples:
  # Submit all replicates
  $0 conf_ils_low_10M

  # Submit just replicate 1
  $0 conf_ils_low_10M --replicates 1

  # Dry run
  $0 conf_ils_low_10M --dry-run

IMPORTANT:
  - Prerequisites: sampledmultrees.txt must exist (from alloppnet_run_array.sh)
  - Runtime: ~minutes per network-replicate
  - Each replicate submission launches 8 parallel jobs
  - Monitor with: squeue -u \$USER | grep alloppnet_postprocess
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
            # Parse comma-separated list
            REPLICATES=$(echo "$2" | tr ',' ' ')
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

# Use default replicates if not specified
if [ -z "$REPLICATES" ]; then
    REPLICATES="$ALL_REPLICATES"
fi

# Validate replicate numbers
for rep in $REPLICATES; do
    if ! [[ "$rep" =~ ^[1-5]$ ]]; then
        echo "ERROR: Invalid replicate number: $rep (must be 1-5)"
        exit 1
    fi
done

# Count replicates and total jobs
NUM_REPLICATES=$(echo "$REPLICATES" | wc -w)
TOTAL_JOBS=$((NUM_REPLICATES * NUM_NETWORKS))

# ============================================================================
# DISPLAY CONFIGURATION
# ============================================================================

echo "============================================================================"
echo "ALLOPPNET POST-PROCESSING SUBMISSION"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Replicates: ${REPLICATES// /, }"
echo "Networks per replicate: ${NUM_NETWORKS}"
echo "Total jobs: ${TOTAL_JOBS}"
echo "Estimated runtime: ~minutes per job"
echo ""
if [ "$DRY_RUN" = true ]; then
    echo "MODE: DRY RUN (commands will be shown but not executed)"
    echo ""
fi
echo "============================================================================"
echo ""

# ============================================================================
# SUBMIT JOBS
# ============================================================================

run_cmd() {
    if [ "$DRY_RUN" = true ]; then
        echo "[DRY RUN] Would execute: $*"
        return 0
    else
        eval "$@"
    fi
}

JOB_IDS=()

# Submit one array job per replicate
for rep in $REPLICATES; do
    echo "Submitting replicate ${rep} (8 networks)..."

    JOB_NAME="alloppnet_postprocess_${CONFIG}_rep${rep}"
    LOG_PREFIX="alloppnet_postprocess_${CONFIG}_rep${rep}"

    SUBMIT_CMD="sbatch \
        --job-name=${JOB_NAME} \
        --array=1-8 \
        --export=CONFIG=${CONFIG},REPLICATE=${rep} \
        --output=${LOG_DIR}/${LOG_PREFIX}_%A_%a.out \
        --error=${LOG_DIR}/${LOG_PREFIX}_%A_%a.err \
        ${ARRAY_SCRIPT}"

    if [ "$DRY_RUN" = true ]; then
        run_cmd "$SUBMIT_CMD"
        echo "  [DRY RUN] Job ID: XXXXX"
        JOB_IDS+=("XXXXX")
    else
        OUTPUT=$(run_cmd "$SUBMIT_CMD")
        JOB_ID=$(echo "$OUTPUT" | grep -oP '\d+$')
        echo "  ✓ Submitted: ${JOB_NAME} (Job ID: ${JOB_ID})"
        JOB_IDS+=("$JOB_ID")
    fi

    echo ""

    # Small delay to avoid overwhelming scheduler
    if [ "$DRY_RUN" = false ]; then
        sleep 1
    fi
done

# ============================================================================
# SUMMARY
# ============================================================================

echo "============================================================================"
echo "SUBMISSION SUMMARY"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Array jobs submitted: ${NUM_REPLICATES}"
echo "Total jobs: ${TOTAL_JOBS} (${NUM_NETWORKS} networks × ${NUM_REPLICATES} replicates)"
echo "Estimated runtime: ~minutes per job"
echo ""

if [ "$DRY_RUN" = false ]; then
    echo "Job IDs:"
    idx=0
    for rep in $REPLICATES; do
        echo "  Replicate ${rep}: ${JOB_IDS[$idx]} (8 array tasks)"
        idx=$((idx + 1))
    done
    echo ""
fi

echo "Networks (AlloppNET-compatible):"
echo "  - Bendiksby_2011"
echo "  - Ding_2023"
echo "  - Koenen_2020"
echo "  - Liu_2023"
echo "  - Shahrestani_2015"
echo "  - Wisecaver_2023"
echo "  - Wu_2015"
echo "  - Zhao_2021"
echo ""

echo "Post-processing steps (per job):"
echo "  1. Count trees in sampledmultrees.txt"
echo "  2. Calculate 10% burnin"
echo "  3. Run TreeAnnotator (consensus tree)"
echo "  4. Remove copy number suffixes (_0, _1)"
echo "  5. Create alloppnet_result.tre"
echo ""

if [ "$DRY_RUN" = false ]; then
    echo "Monitor jobs:"
    echo "  squeue -u \$USER | grep alloppnet_postprocess"
    echo "  squeue -u \$USER | grep ${CONFIG}"
    echo ""

    echo "Check specific job:"
    for job_id in "${JOB_IDS[@]}"; do
        if [ "$job_id" != "XXXXX" ]; then
            echo "  squeue -j ${job_id}"
            break
        fi
    done
    echo ""
fi

echo "Check logs:"
echo "  ${LOG_DIR}/alloppnet_postprocess_${CONFIG}_rep*"
echo ""

echo "Check status when complete:"
echo "  python simulations/scripts/check_pipeline_status.py ${CONFIG} --method alloppnet --step run --verbose"
echo "============================================================================"

exit 0

