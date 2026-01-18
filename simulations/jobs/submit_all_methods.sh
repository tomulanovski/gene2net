#!/bin/bash

# ============================================================================
# MASTER PIPELINE ORCHESTRATOR
# ============================================================================
# This script orchestrates running multiple phylogenetic network inference
# methods across multiple configurations. It calls existing submit_*_pipeline.sh
# scripts in parallel while providing flexible filtering and parameter customization.
#
# Prerequisites:
#   - SimPhy simulations must already be completed
#   - validate with: python ../scripts/check_pipeline_status.py <CONFIG> --step simphy
#
# Usage:
#   ./submit_all_methods.sh [CONFIG1 CONFIG2 ...] [OPTIONS]
#
# Examples:
#   # Run all methods for three configs
#   ./submit_all_methods.sh conf_dup_loss_low conf_dup_loss_medium conf_dup_loss_high
#
#   # Prep only, then run separately
#   ./submit_all_methods.sh conf_dup_loss_low --prep-only
#   ./submit_all_methods.sh conf_dup_loss_low --run-only
#
#   # Specific methods only
#   ./submit_all_methods.sh conf_ils_low_10M --methods grampa,polyphest
#
#   # Custom method parameters
#   ./submit_all_methods.sh conf_ils_low_10M --methods polyphest --polyphest-percentile 50
#
#   # Dry run to preview
#   ./submit_all_methods.sh conf_dup_loss_low --dry-run
#
# Post-processing (run after SLURM jobs complete):
#   python ../scripts/check_pipeline_status.py <CONFIG> --step run
#   python ../scripts/postprocess_results.py <CONFIG>
#   python ../scripts/run_full_summary.py <CONFIG>
# ============================================================================

# Note: Not using 'set -e' to allow script to continue even if some submissions fail
set -uo pipefail

# ============================================================================
# DEFAULTS
# ============================================================================

# Default configurations
DEFAULT_CONFIGS=(
    "conf_ils_low_10M"
    "conf_ils_medium_10M"
    "conf_ils_high_10M"
)

# All available methods
ALL_METHODS=(grampa grandma_split polyphest padre mpsugar alloppnet)

# Common parameters
NUM_REPLICATES=5
PREP_ONLY=false
RUN_ONLY=false
DRY_RUN=false
VERBOSE=false

# Method-specific defaults
POLYPHEST_PERCENTILE=50
POLYPHEST_ISO_THRESHOLD=0.2
PADRE_JAVA_MEM="4g"
MPSUGAR_ITERATIONS=500
MPSUGAR_CHAINS=1

# GRAMPA-specific flags
GRAMPA_ASTRAL_ONLY=false
GRAMPA_SKIP_PREP=false

# Script paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

usage() {
    cat << EOF
Usage: $0 [CONFIG1 CONFIG2 ...] [OPTIONS]

Orchestrates method pipeline submissions across multiple configurations.
Assumes SimPhy simulations already completed.

Configuration selection:
  CONFIGS...                         One or more configuration names
                                     If omitted, defaults to: ${DEFAULT_CONFIGS[*]}

Method and step control:
  --methods METHOD1,METHOD2,...      Comma-separated list of methods to run
                                     Valid: ${ALL_METHODS[*]}
                                     Default: all methods
  --prep-only                        Only run preparation steps
  --run-only                         Only run method execution (assumes prep done)

Common parameters:
  --replicates N                     Number of replicates (default: 5)
  --dry-run                          Preview submissions without executing
  --verbose                          Show detailed output from each submission

Method-specific parameters:
  --polyphest-percentile N           Polyphest percentile (default: 60)
  --polyphest-iso-threshold N        Polyphest iso threshold (default: 0.2)
  --padre-java-mem SIZE              PADRE Java heap size (default: 4g)
  --mpsugar-iterations N             MP-SUGAR iterations (default: 500)
  --mpsugar-chains N                 MP-SUGAR chains (default: 1)

GRAMPA-specific (only with single config):
  --grampa-astral-only               Only run ASTRAL step
  --grampa-skip-prep                 Skip prep, run ASTRAL + GRAMPA

Other:
  -h, --help                         Show this help message

Examples:
  # Run all methods for three configs
  $0 conf_dup_loss_low conf_dup_loss_medium conf_dup_loss_high

  # Prep only
  $0 conf_dup_loss_low --prep-only

  # Specific methods
  $0 conf_ils_low_10M --methods grampa,polyphest

  # Custom parameters
  $0 conf_ils_low_10M --methods polyphest --polyphest-percentile 50

  # Dry run
  $0 conf_dup_loss_low --dry-run

Post-processing (after jobs complete):
  python ../scripts/check_pipeline_status.py <CONFIG> --step run
  python ../scripts/postprocess_results.py <CONFIG>
  python ../scripts/run_full_summary.py <CONFIG>
EOF
    exit 1
}

validate_method() {
    local method="$1"
    for valid_method in "${ALL_METHODS[@]}"; do
        if [ "$method" = "$valid_method" ]; then
            return 0
        fi
    done
    return 1
}

validate_config() {
    local config="$1"
    # Non-fatal validation - just warn if suspicious
    if [[ ! "$config" =~ ^conf_ ]]; then
        echo "WARNING: Configuration name doesn't start with 'conf_': $config" >&2
        echo "         Continuing anyway - you may know better" >&2
    fi
    return 0
}

log_info() {
    echo "  $1"
}

log_error() {
    echo "  ERROR: $1" >&2
}

# ============================================================================
# ARGUMENT PARSING
# ============================================================================

# Parse positional config arguments
CONFIGS=()
while [ $# -gt 0 ] && [[ ! "$1" =~ ^-- ]] && [[ ! "$1" =~ ^-h$ ]]; do
    CONFIGS+=("$1")
    shift
done

# If no configs provided, use defaults
if [ ${#CONFIGS[@]} -eq 0 ]; then
    CONFIGS=("${DEFAULT_CONFIGS[@]}")
fi

# Parse option flags
METHODS_FILTER=""
while [ $# -gt 0 ]; do
    case "$1" in
        --methods)
            METHODS_FILTER="$2"
            shift 2
            ;;
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
        --polyphest-percentile)
            POLYPHEST_PERCENTILE="$2"
            shift 2
            ;;
        --polyphest-iso-threshold)
            POLYPHEST_ISO_THRESHOLD="$2"
            shift 2
            ;;
        --padre-java-mem)
            PADRE_JAVA_MEM="$2"
            shift 2
            ;;
        --mpsugar-iterations)
            MPSUGAR_ITERATIONS="$2"
            shift 2
            ;;
        --mpsugar-chains)
            MPSUGAR_CHAINS="$2"
            shift 2
            ;;
        --grampa-astral-only)
            GRAMPA_ASTRAL_ONLY=true
            shift
            ;;
        --grampa-skip-prep)
            GRAMPA_SKIP_PREP=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            echo "Use -h or --help for usage information" >&2
            exit 1
            ;;
    esac
done

# Parse methods filter
if [ -n "$METHODS_FILTER" ]; then
    IFS=',' read -ra METHODS <<< "$METHODS_FILTER"
else
    METHODS=("${ALL_METHODS[@]}")
fi

# ============================================================================
# VALIDATION
# ============================================================================

# Validate methods
for method in "${METHODS[@]}"; do
    if ! validate_method "$method"; then
        echo "ERROR: Invalid method: $method" >&2
        echo "Valid methods: ${ALL_METHODS[*]}" >&2
        exit 1
    fi
done

# Validate configs
for config in "${CONFIGS[@]}"; do
    validate_config "$config"
done

# Validate conflicting options
if [ "$PREP_ONLY" = true ] && [ "$RUN_ONLY" = true ]; then
    echo "ERROR: Cannot specify both --prep-only and --run-only" >&2
    exit 1
fi

# GRAMPA-specific validation
if [ "$GRAMPA_ASTRAL_ONLY" = true ] || [ "$GRAMPA_SKIP_PREP" = true ]; then
    if [ ${#CONFIGS[@]} -gt 1 ]; then
        echo "ERROR: GRAMPA-specific flags only work with single config" >&2
        echo "       You specified ${#CONFIGS[@]} configs" >&2
        exit 1
    fi
    if [[ ! " ${METHODS[@]} " =~ " grampa " ]]; then
        echo "ERROR: GRAMPA-specific flags require --methods grampa" >&2
        exit 1
    fi
fi

# AlloppNET-specific validation (no longer needed - now supports prep/run split)
# AlloppNET now supports --prep-only and --run-only like other methods

# Validate that submit scripts exist
for method in "${METHODS[@]}"; do
    submit_script="${SCRIPT_DIR}/submit_${method}_pipeline.sh"
    if [ ! -f "$submit_script" ]; then
        echo "ERROR: Submit script not found: $submit_script" >&2
        exit 1
    fi
    if [ ! -x "$submit_script" ]; then
        echo "ERROR: Submit script not executable: $submit_script" >&2
        echo "       Run: chmod +x $submit_script" >&2
        exit 1
    fi
done

# ============================================================================
# DISPLAY CONFIGURATION SUMMARY
# ============================================================================

echo "============================================================================"
echo "MASTER PIPELINE ORCHESTRATOR"
echo "============================================================================"
echo ""
echo "Configurations (${#CONFIGS[@]}):"
for config in "${CONFIGS[@]}"; do
    echo "  - ${config}"
done
echo ""
echo "Methods (${#METHODS[@]}):"
for method in "${METHODS[@]}"; do
    echo "  - ${method}"
done
echo ""
echo "Common Parameters:"
echo "  Replicates: ${NUM_REPLICATES}"
echo "  Networks: 21 (AlloppNET: 8 compatible networks only)"
echo ""

# Show method-specific parameters only for selected methods
method_params_shown=false
if [[ " ${METHODS[@]} " =~ " polyphest " ]]; then
    if [ "$method_params_shown" = false ]; then
        echo "Method-Specific Parameters:"
        method_params_shown=true
    fi
    echo "  Polyphest:"
    echo "    Percentile: ${POLYPHEST_PERCENTILE}"
    echo "    Iso threshold: ${POLYPHEST_ISO_THRESHOLD}"
fi

if [[ " ${METHODS[@]} " =~ " padre " ]]; then
    if [ "$method_params_shown" = false ]; then
        echo "Method-Specific Parameters:"
        method_params_shown=true
    fi
    echo "  PADRE:"
    echo "    Java memory: ${PADRE_JAVA_MEM}"
fi

if [[ " ${METHODS[@]} " =~ " mpsugar " ]]; then
    if [ "$method_params_shown" = false ]; then
        echo "Method-Specific Parameters:"
        method_params_shown=true
    fi
    echo "  MP-SUGAR:"
    echo "    Iterations: ${MPSUGAR_ITERATIONS}"
    echo "    Chains: ${MPSUGAR_CHAINS}"
fi

if [[ " ${METHODS[@]} " =~ " alloppnet " ]]; then
    if [ "$method_params_shown" = false ]; then
        echo "Method-Specific Parameters:"
        method_params_shown=true
    fi
    echo "  AlloppNET:"
    echo "    Compatible networks: 8 (Bendiksby, Ding, Koenen, Liu,"
    echo "                            Shahrestani, Wisecaver, Wu, Zhao)"
    echo "    PREP: Fast (~minutes), generates BEAST XML with kernel smoothing"
    echo "    RUN: Slow (~5 days per network-replicate, 100M BEAST iterations)"
fi

if [ "$method_params_shown" = true ]; then
    echo ""
fi

echo "Steps to run:"
if [ "$PREP_ONLY" = true ]; then
    echo "  [x] Prep only"
elif [ "$RUN_ONLY" = true ]; then
    echo "  [x] Run only (assumes prep done)"
else
    echo "  [x] Prep + Run (full pipeline)"
fi
echo ""

if [ "$DRY_RUN" = true ]; then
    echo "MODE: DRY RUN (no jobs will be submitted)"
    echo ""
fi

echo "============================================================================"
echo ""

# ============================================================================
# COMMAND BUILDER FUNCTION
# ============================================================================

build_submit_command() {
    local method="$1"
    local config="$2"
    local submit_script="${SCRIPT_DIR}/submit_${method}_pipeline.sh"
    local cmd="${submit_script} ${config}"

    # Add common parameters
    # AlloppNET expects comma-separated list, others expect a single number
    if [ "$method" = "alloppnet" ]; then
        local replicate_list=$(seq -s ',' 1 ${NUM_REPLICATES})
        cmd="${cmd} --replicates ${replicate_list}"
    else
        cmd="${cmd} --replicates ${NUM_REPLICATES}"
    fi

    # Add step control flags
    if [ "$PREP_ONLY" = true ]; then
        cmd="${cmd} --prep-only"
    elif [ "$RUN_ONLY" = true ]; then
        cmd="${cmd} --run-only"
    fi

    # Add dry-run flag
    if [ "$DRY_RUN" = true ]; then
        cmd="${cmd} --dry-run"
    fi

    # Add method-specific parameters
    case "$method" in
        polyphest)
            cmd="${cmd} --percentile ${POLYPHEST_PERCENTILE}"
            cmd="${cmd} --iso-threshold ${POLYPHEST_ISO_THRESHOLD}"
            ;;
        padre)
            cmd="${cmd} --java-mem ${PADRE_JAVA_MEM}"
            ;;
        mpsugar)
            cmd="${cmd} --iterations ${MPSUGAR_ITERATIONS}"
            cmd="${cmd} --chains ${MPSUGAR_CHAINS}"
            ;;
        grampa)
            if [ "$GRAMPA_ASTRAL_ONLY" = true ]; then
                cmd="${cmd} --astral-only"
            elif [ "$GRAMPA_SKIP_PREP" = true ]; then
                cmd="${cmd} --skip-prep"
            fi
            ;;
    esac

    echo "$cmd"
}

# ============================================================================
# MAIN SUBMISSION LOOP
# ============================================================================

# Track all submissions
declare -A SUBMISSION_TRACKER
TOTAL_SUBMISSIONS=0
SUCCESSFUL_SUBMISSIONS=0
FAILED_SUBMISSIONS=0

# Iterate over configurations
for config in "${CONFIGS[@]}"; do
    echo "------------------------------------------------------------------------"
    echo "Configuration: ${config}"
    echo "------------------------------------------------------------------------"
    echo ""

    # Iterate over methods
    for method in "${METHODS[@]}"; do
        log_info "Submitting ${method} for ${config}..."

        # Build command
        cmd=$(build_submit_command "$method" "$config")

        # Show command if verbose
        if [ "$VERBOSE" = true ]; then
            echo "  Command: ${cmd}"
        fi

        if [ "$DRY_RUN" = true ]; then
            echo "  [DRY RUN] Would execute: ${cmd}"
            SUBMISSION_TRACKER["${config}_${method}"]="DRY_RUN"
            ((TOTAL_SUBMISSIONS++))
        else
            # Execute and capture output
            if output=$(eval "$cmd" 2>&1); then
                log_info "SUCCESS: ${method} submitted for ${config}"
                SUBMISSION_TRACKER["${config}_${method}"]="SUCCESS"
                ((SUCCESSFUL_SUBMISSIONS++))

                # Extract job IDs if verbose
                if [ "$VERBOSE" = true ]; then
                    echo "$output" | grep -E "job.*ID|Job submitted|job ID" || true
                fi
            else
                log_error "FAILED: ${method} failed for ${config}"
                SUBMISSION_TRACKER["${config}_${method}"]="FAILED"
                ((FAILED_SUBMISSIONS++))

                echo "  Error output:"
                echo "$output" | sed 's/^/    /'
            fi
            ((TOTAL_SUBMISSIONS++))
        fi
        echo ""
    done
    echo ""
done

# ============================================================================
# FINAL SUMMARY & MONITORING GUIDE
# ============================================================================

echo "============================================================================"
echo "SUBMISSION SUMMARY"
echo "============================================================================"
echo ""
echo "Total submissions: ${TOTAL_SUBMISSIONS}"
if [ "$DRY_RUN" = false ]; then
    echo "  Successful: ${SUCCESSFUL_SUBMISSIONS}"
    echo "  Failed: ${FAILED_SUBMISSIONS}"
fi
echo ""

# Display submission matrix
echo "Submission Matrix:"
echo ""
printf "%-25s" "Config/Method"
for method in "${METHODS[@]}"; do
    printf "%-15s" "${method}"
done
echo ""
printf '%.0s-' {1..80}
echo ""

for config in "${CONFIGS[@]}"; do
    printf "%-25s" "${config}"
    for method in "${METHODS[@]}"; do
        status="${SUBMISSION_TRACKER[${config}_${method}]:-}"
        case "$status" in
            SUCCESS)
                printf "%-15s" "✓"
                ;;
            FAILED)
                printf "%-15s" "✗"
                ;;
            DRY_RUN)
                printf "%-15s" "[dry]"
                ;;
            *)
                printf "%-15s" "-"
                ;;
        esac
    done
    echo ""
done
echo ""

# Monitoring commands
if [ "$DRY_RUN" = false ] && [ $SUCCESSFUL_SUBMISSIONS -gt 0 ]; then
    echo "Monitoring Commands:"
    echo ""
    echo "  # Check all jobs for your user"
    echo "  squeue -u \$USER"
    echo ""
    echo "  # Check specific method logs (tail recent output)"
    for method in "${METHODS[@]}"; do
        echo "  tail -f ${LOG_DIR}/*${method}*.out | head -50"
    done
    echo ""
    echo "  # Validate pipeline status (after prep completes)"
    for config in "${CONFIGS[@]}"; do
        echo "  python ../scripts/check_pipeline_status.py ${config} --step prep"
    done
    echo ""
    echo "  # Validate pipeline status (after run completes)"
    for config in "${CONFIGS[@]}"; do
        echo "  python ../scripts/check_pipeline_status.py ${config} --step run --verbose"
    done
    echo ""
fi

# Next steps
echo "Next Steps:"
echo ""
if [ "$PREP_ONLY" = true ]; then
    echo "  1. Wait for prep jobs to complete (monitor with squeue)"
    echo "  2. Validate prep output:"
    for config in "${CONFIGS[@]}"; do
        echo "     python ../scripts/check_pipeline_status.py ${config} --step prep"
    done
    echo "  3. Submit run jobs:"
    echo "     $0 ${CONFIGS[*]} --run-only"
elif [ "$RUN_ONLY" = false ] && [ "$DRY_RUN" = false ]; then
    echo "  1. Monitor job progress with squeue -u \$USER"
    echo "  2. After prep completes, validate with check_pipeline_status.py --step prep"
    echo "  3. After run completes, validate with check_pipeline_status.py --step run"
    echo "  4. Post-process results:"
    for config in "${CONFIGS[@]}"; do
        echo "     python ../scripts/postprocess_results.py ${config}"
    done
    echo "  5. Run full summary analysis:"
    for config in "${CONFIGS[@]}"; do
        echo "     python ../scripts/run_full_summary.py ${config}"
    done
elif [ "$DRY_RUN" = true ]; then
    echo "  This was a dry run. Remove --dry-run flag to actually submit jobs."
fi

echo ""
echo "Log directory: ${LOG_DIR}"
echo "============================================================================"
