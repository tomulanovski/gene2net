#!/bin/bash

# ============================================================================
# ALLOPPNET SUBMISSION HELPER
# ============================================================================
# This script submits AlloppNET jobs for all replicates.
#
# AlloppNET runs on 8 networks with max 2 copies:
#   Bendiksby_2011, Ding_2023, Koenen_2020, Liu_2023,
#   Shahrestani_2015, Wisecaver_2023, Wu_2015, Zhao_2021
#
# Each replicate submission creates 8 array tasks (one per network).
# Each array task runs the full pipeline:
#   1. Prepare input (PHY → NEXUS, ploidy, taxa table)
#   2. Generate BEAST XML
#   3. Run BEAST (100M iterations, ~5 days)
#   4. Summarize with TreeAnnotator
#   5. Post-process (remove copy suffixes)
#
# Usage:
#   ./submit_alloppnet.sh <config_name> [options]
#
# Examples:
#   # Submit all 5 replicates (40 total BEAST runs)
#   ./submit_alloppnet.sh conf_ils_low_10M
#
#   # Submit specific replicates
#   ./submit_alloppnet.sh conf_ils_low_10M --replicates 1,3,5
#
#   # Dry run (show commands without submitting)
#   ./submit_alloppnet.sh conf_ils_low_10M --dry-run
# ============================================================================

set -euo pipefail

# ============================================================================
# DEFAULTS
# ============================================================================

NUM_NETWORKS=8
ALL_REPLICATES="1 2 3 4 5"
REPLICATES=""
DRY_RUN=false
PREP_ONLY=false
RUN_ONLY=false

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_ARRAY_SCRIPT="${SCRIPT_DIR}/alloppnet_prep_array.sh"
RUN_ARRAY_SCRIPT="${SCRIPT_DIR}/alloppnet_run_array.sh"
OLD_ARRAY_SCRIPT="${SCRIPT_DIR}/alloppnet_array.sh"  # Deprecated but kept for reference
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
  --prep-only       Only run preparation step (fast, ~minutes)
  --run-only        Only run BEAST step (assumes prep done, ~5 days)
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

Pipeline steps:
  PREP (fast, ~minutes per network):
    1. Convert PHY → NEXUS (1000 alignments)
    2. Analyze copy distributions (kernel smoothing)
    3. Generate ploidy_level.json
    4. Create taxa_table.txt (homeolog pairing)
    5. Generate BEAST XML (AlloppDT scripts)

  RUN (slow, ~5 days per network):
    1. Run BEAST (100M iterations)
    2. Summarize with TreeAnnotator (10% burnin)
    3. Post-process (remove copy suffixes)

Job structure:
  - One array job per replicate
  - Each array job has 8 tasks (one per network)
  - Total jobs = replicates × 8 networks

Examples:
  # Submit full pipeline (prep + run)
  $0 conf_ils_low_10M

  # Two-step workflow (recommended)
  $0 conf_ils_low_10M --prep-only
  # ... validate prep outputs ...
  $0 conf_ils_low_10M --run-only

  # Submit just replicate 1
  $0 conf_ils_low_10M --replicates 1 --prep-only

  # Dry run
  $0 conf_ils_low_10M --dry-run

IMPORTANT:
  - PREP step: fast (~minutes), can validate before expensive BEAST run
  - RUN step: slow (~5 days per network-replicate)
  - Each replicate submission launches 8 parallel jobs
  - Monitor with: squeue -u \$USER | grep alloppnet
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

# Determine which step to run
if [ "$PREP_ONLY" = true ]; then
    STEP="PREP"
    ARRAY_SCRIPT="$PREP_ARRAY_SCRIPT"
    STEP_DESC="Preparation (NEXUS conversion + XML generation)"
    RUNTIME_EST="~minutes"
elif [ "$RUN_ONLY" = true ]; then
    STEP="RUN"
    ARRAY_SCRIPT="$RUN_ARRAY_SCRIPT"
    STEP_DESC="Run (BEAST + summarize)"
    RUNTIME_EST="~5 days per job"
else
    STEP="FULL"
    ARRAY_SCRIPT="$OLD_ARRAY_SCRIPT"
    STEP_DESC="Full pipeline (prep + run)"
    RUNTIME_EST="~5 days per job (BEAST dominates)"
    echo "WARNING: Full pipeline mode is deprecated. Consider using --prep-only then --run-only instead."
    echo ""
fi

# ============================================================================
# DISPLAY CONFIGURATION
# ============================================================================

echo "============================================================================"
echo "ALLOPPNET SUBMISSION"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Step: ${STEP} (${STEP_DESC})"
echo "Replicates: ${REPLICATES// /, }"
echo "Networks per replicate: ${NUM_NETWORKS}"
echo "Total jobs: ${TOTAL_JOBS}"
echo "Estimated runtime: ${RUNTIME_EST}"
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

    # Create job name with step
    if [ "$STEP" = "PREP" ]; then
        JOB_NAME="alloppnet_prep_${CONFIG}_rep${rep}"
        LOG_PREFIX="alloppnet_prep_${CONFIG}_rep${rep}"
    elif [ "$STEP" = "RUN" ]; then
        JOB_NAME="alloppnet_run_${CONFIG}_rep${rep}"
        LOG_PREFIX="alloppnet_run_${CONFIG}_rep${rep}"
    else
        JOB_NAME="alloppnet_${CONFIG}_rep${rep}"
        LOG_PREFIX="alloppnet_${CONFIG}_rep${rep}"
    fi

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
echo "Step: ${STEP} (${STEP_DESC})"
echo "Array jobs submitted: ${NUM_REPLICATES}"
echo "Total jobs: ${TOTAL_JOBS} (${NUM_NETWORKS} networks × ${NUM_REPLICATES} replicates)"
echo "Estimated runtime: ${RUNTIME_EST}"
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

if [ "$STEP" = "PREP" ]; then
    echo "Prep steps (per job):"
    echo "  1. Convert PHY → NEXUS (1000 alignments)"
    echo "  2. Analyze copy distributions (kernel smoothing)"
    echo "  3. Generate ploidy_level.json"
    echo "  4. Create taxa_table.txt"
    echo "  5. Generate BEAST XML"
elif [ "$STEP" = "RUN" ]; then
    echo "Run steps (per job):"
    echo "  1. Run BEAST (100M iterations, ~5 days)"
    echo "  2. Summarize with TreeAnnotator"
    echo "  3. Post-process (remove copy suffixes)"
else
    echo "Full pipeline steps (per job):"
    echo "  1. Prepare input (PHY → NEXUS, ploidy, taxa table)"
    echo "  2. Generate BEAST XML"
    echo "  3. Run BEAST (100M iterations, ~5 days)"
    echo "  4. Summarize with TreeAnnotator"
    echo "  5. Post-process (remove copy suffixes)"
fi
echo ""

if [ "$DRY_RUN" = false ]; then
    echo "Monitor jobs:"
    echo "  squeue -u \$USER | grep alloppnet"
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
if [ "$STEP" = "PREP" ]; then
    echo "  ${LOG_DIR}/alloppnet_prep_${CONFIG}_rep*"
elif [ "$STEP" = "RUN" ]; then
    echo "  ${LOG_DIR}/alloppnet_run_${CONFIG}_rep*"
else
    echo "  ${LOG_DIR}/alloppnet_${CONFIG}_rep*"
fi
echo ""

if [ "$STEP" = "PREP" ]; then
    echo "Next step (after prep completes):"
    echo "  1. Validate prep outputs:"
    echo "     python simulations/scripts/check_pipeline_status.py ${CONFIG} --method alloppnet --step prep"
    echo "  2. Submit run jobs:"
    echo "     $0 ${CONFIG} --run-only"
elif [ "$STEP" = "RUN" ]; then
    echo "Check status when complete:"
    echo "  python simulations/scripts/check_pipeline_status.py ${CONFIG} --method alloppnet --step run --verbose"
else
    echo "Check status when complete:"
    echo "  python simulations/scripts/check_pipeline_status.py ${CONFIG} --method alloppnet --step run"
fi
echo "============================================================================"

exit 0
