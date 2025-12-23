#!/bin/bash
#SBATCH --job-name=astral_grampa
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/astral_grampa_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/astral_grampa_%A_%a.err
#SBATCH --time=5:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public 

# ============================================================================
# run_astral.sh
# ============================================================================
# Runs ASTRAL on clean trees to generate species tree for GRAMPA.
# Processes all replicates for a single network.
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M run_astral.sh
#   sbatch --export=CONFIG=conf_ils_low_10M,NUM_REPLICATES=3 run_astral.sh
#
# Environment variables:
#   CONFIG          - Configuration name (required, e.g., conf_ils_low_10M)
#   NUM_REPLICATES  - Number of replicates to process (default: 5)
# ============================================================================

set -eo pipefail

# Initialize LD_LIBRARY_PATH if not set (needed for conda activation)
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ============================================================================
# CONFIGURATION
# ============================================================================

# Required: configuration name
CONFIG="${CONFIG:?ERROR: CONFIG environment variable is required. Use --export=CONFIG=conf_name}"

# Optional: number of replicates (default: 5)
NUM_REPLICATES="${NUM_REPLICATES:-5}"

# ============================================================================
# PATHS
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"

# ASTRAL command (from conda environment)
ASTRAL_CMD="astral4"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# ============================================================================
# SETUP
# ============================================================================

# Get network for this array task
network_idx=$((SLURM_ARRAY_TASK_ID - 1))
network="${networks[$network_idx]}"

# Define paths
GRAMPA_INPUT_DIR="${BASE_DIR}/${network}/processed/${CONFIG}/grampa_input"

echo "============================================================================"
echo "RUN ASTRAL FOR GRAMPA - ${network}"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Replicates: ${NUM_REPLICATES}"
echo "Input directory: ${GRAMPA_INPUT_DIR}"
echo "Date: $(date)"
echo "============================================================================"
echo ""

# ============================================================================
# VALIDATION
# ============================================================================

if [ ! -d "$GRAMPA_INPUT_DIR" ]; then
    echo "ERROR: GRAMPA input directory not found: $GRAMPA_INPUT_DIR"
    echo "Did you run prep_grampa.sh first?"
    exit 1
fi

# ============================================================================
# ACTIVATE CONDA
# ============================================================================

source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

echo "? Conda environment activated: gene2net"

# Check if ASTRAL is available
if ! command -v $ASTRAL_CMD &> /dev/null; then
    echo "ERROR: $ASTRAL_CMD command not found. Check if it's installed in your conda environment."
    exit 1
fi

echo "? ASTRAL found: $(which $ASTRAL_CMD)"
echo ""

# ============================================================================
# MAIN PROCESSING
# ============================================================================

success_count=0
error_count=0

for replicate in $(seq 1 $NUM_REPLICATES); do
    echo "----------------------------------------------------------------------------"
    echo "REPLICATE ${replicate}/${NUM_REPLICATES}"
    echo "----------------------------------------------------------------------------"
    
    # Define paths for this replicate
    replicate_dir="${GRAMPA_INPUT_DIR}/replicate_${replicate}"
    input_file="${replicate_dir}/clean_trees.tre"
    output_file="${replicate_dir}/species.tre"
    log_file="${replicate_dir}/astral.log"
    
    if [ ! -d "$replicate_dir" ]; then
        echo "  ERROR: Replicate directory not found: $replicate_dir"
        error_count=$((error_count + 1))
        continue
    fi
    
    if [ ! -f "$input_file" ]; then
        echo "  ERROR: Clean trees file not found: $input_file"
        error_count=$((error_count + 1))
        continue
    fi
    
    echo "  Input: $input_file"
    echo "  Output: $output_file"
    echo ""
    
    # Run ASTRAL
    echo "  Running ASTRAL..."
    start_time=$(date +%s)
    
    $ASTRAL_CMD -o "$output_file" "$input_file" 2>"$log_file"
    
    exit_code=$?
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    
    if [ $exit_code -eq 0 ] && [ -f "$output_file" ]; then
        echo "  ? ASTRAL completed (${duration}s)"
        echo "  Species tree saved to: $output_file"
        success_count=$((success_count + 1))
    else
        echo "  ? ASTRAL failed (exit code: $exit_code)"
        echo "  Check log file: $log_file"
        error_count=$((error_count + 1))
    fi
    echo ""
done

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "============================================================================"
echo "SUMMARY - ${network}"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Successful replicates: ${success_count}/${NUM_REPLICATES}"
echo "Failed replicates: ${error_count}/${NUM_REPLICATES}"
echo ""

if [ $error_count -gt 0 ]; then
    echo "Status: ? COMPLETED WITH ERRORS"
    exit 1
else
    echo "Status: ? SUCCESS"
fi

echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"