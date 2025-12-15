#!/bin/bash
# ==============================================================================
# SCRIPT: submit_sequences.sh
# PURPOSE: Submit sequence simulation jobs using empirically-sampled GTR+Gamma parameters
#
# USAGE:
#   ./submit_sequences.sh <CONFIG_NAME> [NETWORK_SUBSET]
#
# ARGUMENTS:
#   1. CONFIG_NAME     : Configuration name (e.g., "ils_low_10M")
#   2. NETWORK_SUBSET  : (Optional) Specific network(s) to run.
#                        Can be a single network name or "all" (default: all)
#
# EXAMPLES:
#
#   --- Simulate sequences for all networks ---
#   ./submit_sequences.sh ils_low_10M
#   ./submit_sequences.sh ils_high_10M
#   ./submit_sequences.sh ils_med_10M
#
#   --- Simulate for a single network (for testing) ---
#   ./submit_sequences.sh ils_low_10M Ding_2023
#   ./submit_sequences.sh ils_low_10M Zhao_2021
#
# DESCRIPTION:
#   This script submits jobs to simulate DNA sequences along gene trees using
#   GTR+Gamma model parameters sampled from empirical distributions.
#
#   The empirical parameters come from 2,709 genes across three datasets:
#   - Zhao_2021 (982 genes)
#   - Ren_2024 (727 genes)
#   - Morales_Briones_2021 (1000 genes)
#
#   For each gene tree, the script samples:
#   - GTR rate parameters (AC, AG, AT, CG, CT, GT)
#   - Base frequencies (π_A, π_C, π_G, π_T)
#   - Gamma alpha parameter (rate heterogeneity)
#   - Alignment length
#
#   The same parameter set is used across all replicates for a given gene tree,
#   ensuring consistency while maintaining biological realism.
# ==============================================================================

# Parse arguments
CONFIG_NAME=$1
NETWORK_SUBSET=${2:-all}

if [ -z "$CONFIG_NAME" ]; then
    echo "ERROR: Configuration name is required"
    echo ""
    echo "Usage: $0 <CONFIG_NAME> [NETWORK_SUBSET]"
    echo ""
    echo "Examples:"
    echo "  $0 ils_low_10M              # All networks"
    echo "  $0 ils_low_10M Ding_2023    # Single network"
    echo ""
    exit 1
fi

# Paths
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations"
SCRIPT_PATH="${BASE_DIR}/jobs/simulate_sequences_1_dataset.sh"
GTR_PICKLE="${BASE_DIR}/distributions/gtr_parameters_all.pkl"

# Validation
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "ERROR: Sequence simulation script not found: $SCRIPT_PATH"
    exit 1
fi

if [ ! -f "$GTR_PICKLE" ]; then
    echo "ERROR: GTR parameters file not found: $GTR_PICKLE"
    echo "Please run the GTR parameter extraction script first:"
    echo "  python simulations/scripts/sequence_evolution/filter_and_extract_gtr.py"
    exit 1
fi

# All networks
ALL_NETWORKS=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# Determine which networks to process
if [ "$NETWORK_SUBSET" == "all" ]; then
    NETWORKS=("${ALL_NETWORKS[@]}")
else
    # Check if specified network is valid
    if [[ " ${ALL_NETWORKS[@]} " =~ " ${NETWORK_SUBSET} " ]]; then
        NETWORKS=("$NETWORK_SUBSET")
    else
        echo "ERROR: Unknown network: $NETWORK_SUBSET"
        echo ""
        echo "Available networks:"
        for net in "${ALL_NETWORKS[@]}"; do
            echo "  - $net"
        done
        exit 1
    fi
fi

# Display header
echo "============================================================================"
echo "Submitting Sequence Simulation Jobs"
echo "============================================================================"
echo "Configuration: ${CONFIG_NAME}"
echo "Networks: ${#NETWORKS[@]} (${NETWORK_SUBSET})"
echo "GTR parameters: ${GTR_PICKLE}"
echo "Date: $(date)"
echo "============================================================================"
echo ""
echo "The sequence simulation will automatically detect SimPhy batch mode:"
echo "  - Single batch (1000 trees in one file)"
echo "  - Batches of 10 (100 batches)"
echo "  - Batches of 1 (1000 batches)"
echo ""

# Submit jobs
job_count=0
for network in "${NETWORKS[@]}"; do
    echo "Submitting: ${network}"

    job_name="alisim_${network}_${CONFIG_NAME}"

    sbatch --job-name="${job_name}" \
        --output="${BASE_DIR}/logs/alisim_${CONFIG_NAME}_${network}_%a.out" \
        --error="${BASE_DIR}/logs/alisim_${CONFIG_NAME}_${network}_%a.err" \
        "$SCRIPT_PATH" \
        "${network}" \
        "${CONFIG_NAME}"

    if [ $? -eq 0 ]; then
        ((job_count++))
        echo "  ✓ Submitted: ${job_name}"
    else
        echo "  ✗ Failed: ${network}"
    fi

    # Small delay to avoid overwhelming scheduler
    sleep 1
done

# Summary
echo ""
echo "============================================================================"
echo "Submission Summary"
echo "============================================================================"
echo "Configuration: ${CONFIG_NAME}"
echo "Jobs submitted: ${job_count}/${#NETWORKS[@]}"
echo ""
echo "Each job will:"
echo "  - Process 1000 gene trees (array tasks)"
echo "  - Sample GTR+Gamma parameters from empirical distributions"
echo "  - Generate sequences across all replicates"
echo ""
echo "Monitoring commands:"
echo "  squeue -u \$USER | grep alisim_${CONFIG_NAME}"
echo "  sacct -X --format=JobName,State,ExitCode | grep alisim_${CONFIG_NAME}"
echo ""
echo "Logs directory:"
echo "  ${BASE_DIR}/logs/"
echo "============================================================================"
