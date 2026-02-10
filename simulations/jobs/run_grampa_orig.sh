#!/bin/bash
#SBATCH --job-name=run_grampa
#SBATCH --array=1-105
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_grampa_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_grampa_%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# RUN_GRAMPA.SH
# ============================================================================
# Runs GRAMPA on prepared gene trees and species tree.
# Uses a 2D array structure: 21 networks ª 5 replicates = 105 jobs
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M run_grampa.sh
#   sbatch --export=CONFIG=conf_ils_low_10M,NUM_REPLICATES=3 run_grampa.sh
#
# To run with different array size (e.g., 3 replicates = 63 jobs):
#   sbatch --array=1-63 --export=CONFIG=conf_ils_low_10M,NUM_REPLICATES=3 run_grampa.sh
#
# Environment variables:
#   CONFIG          - Configuration name (required, e.g., conf_ils_low_10M)
#   NUM_REPLICATES  - Number of replicates (default: 5)
# ============================================================================

set -eo pipefail

# Initialize LD_LIBRARY_PATH if not set (needed for conda activation)
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ============================================================================
# CONFIGURATION
# ============================================================================

# Required: configuration name
CONFIG="${CONFIG:?ERROR: CONFIG environment variable is required. Use --export=CONFIG=conf_name}"

# Optional parameters with defaults
NUM_REPLICATES="${NUM_REPLICATES:-5}"

# Number of networks (fixed)
NUM_NETWORKS=21

# ============================================================================
# PATHS
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# ============================================================================
# CALCULATE NETWORK AND REPLICATE FROM ARRAY TASK ID
# ============================================================================
# Array task IDs: 1-105 (for 21 networks ª 5 replicates)
# Mapping: task_id = (network_idx * NUM_REPLICATES) + replicate
# So: network_idx = (task_id - 1) / NUM_REPLICATES
#     replicate = ((task_id - 1) % NUM_REPLICATES) + 1

task_id=${SLURM_ARRAY_TASK_ID}
network_idx=$(( (task_id - 1) / NUM_REPLICATES ))
replicate=$(( ((task_id - 1) % NUM_REPLICATES) + 1 ))

# Validate network index
if [ $network_idx -ge $NUM_NETWORKS ]; then
    echo "ERROR: Invalid task ID ${task_id} - network index ${network_idx} exceeds ${NUM_NETWORKS}"
    exit 1
fi

network="${networks[$network_idx]}"

# ============================================================================
# PATHS FOR THIS JOB
# ============================================================================

INPUT_DIR="${BASE_DIR}/${network}/processed/${CONFIG}/grampa_input/replicate_${replicate}"
OUTPUT_DIR="${BASE_DIR}/${network}/results/${CONFIG}/grampa/replicate_${replicate}"

GENE_TREES="${INPUT_DIR}/grampa_trees.tre"
SPECIES_TREE="${INPUT_DIR}/species.tre"

# ============================================================================
# HEADER
# ============================================================================

echo "============================================================================"
echo "RUN GRAMPA"
echo "============================================================================"
echo "Network: ${network} (index: ${network_idx})"
echo "Replicate: ${replicate}"
echo "Task ID: ${task_id}"
echo "Configuration: ${CONFIG}"
echo ""
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Date: $(date)"
echo "============================================================================"
echo ""

# ============================================================================
# VALIDATION
# ============================================================================

if [ ! -f "$GENE_TREES" ]; then
    echo "ERROR: Gene trees file not found: $GENE_TREES"
    echo "Did you run prep_grampa.sh first?"
    exit 1
fi

if [ ! -f "$SPECIES_TREE" ]; then
    echo "ERROR: Species tree file not found: $SPECIES_TREE"
    echo "Did you run run_astral_grampa.sh first?"
    exit 1
fi

echo "? Input files validated"
echo "  Gene trees: $GENE_TREES"
echo "  Species tree: $SPECIES_TREE"
echo ""

# ============================================================================
# ACTIVATE CONDA
# ============================================================================

source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

echo "? Conda environment activated: gene2net"

# Check if GRAMPA is available
if ! command -v grampa &> /dev/null; then
    echo "ERROR: grampa command not found. Check if it's installed in your conda environment."
    exit 1
fi

echo "? GRAMPA found: $(which grampa)"
echo ""

# ============================================================================
# CREATE OUTPUT DIRECTORY
# ============================================================================

mkdir -p "$OUTPUT_DIR"
echo "? Output directory created: $OUTPUT_DIR"
echo ""

# ============================================================================
# RUN GRAMPA
# ============================================================================

echo "Starting GRAMPA..."
echo "============================================================================"
echo ""

start_time=$(date +%s)

grampa -g "$GENE_TREES" -s "$SPECIES_TREE" -o "$OUTPUT_DIR" --overwrite

exit_code=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "============================================================================"

if [ $exit_code -eq 0 ]; then
    echo "? GRAMPA COMPLETED SUCCESSFULLY"
    echo ""
    echo "Network: ${network}"
    echo "Replicate: ${replicate}"
    echo "Duration: ${duration} seconds ($(echo "scale=2; $duration/3600" | bc) hours)"
    echo "Results saved to: $OUTPUT_DIR"
    echo ""
    
    # List output files
    echo "Output files:"
    ls -la "$OUTPUT_DIR" 2>/dev/null || echo "  (unable to list)"
else
    echo "? GRAMPA FAILED"
    echo ""
    echo "Network: ${network}"
    echo "Replicate: ${replicate}"
    echo "Exit code: ${exit_code}"
    echo "Duration: ${duration} seconds"
fi

echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

exit $exit_code