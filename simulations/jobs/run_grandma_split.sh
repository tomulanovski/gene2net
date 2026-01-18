#!/bin/bash
#SBATCH --job-name=run_grandma_split
#SBATCH --array=1-105
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_grandma_split_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_grandma_split_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# RUN_GRANDMA_SPLIT.SH
# ============================================================================
# Runs GRANDMA_SPLIT on prepared gene trees and species tree.
# Uses the same inputs as GRAMPA (gene trees and species tree from grampa_input).
# Uses a 2D array structure: 21 networks × 5 replicates = 105 jobs
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M run_grandma_split.sh
#   sbatch --export=CONFIG=conf_ils_low_10M,NUM_REPLICATES=3 run_grandma_split.sh
#
# To run with different array size (e.g., 3 replicates = 63 jobs):
#   sbatch --array=1-63 --export=CONFIG=conf_ils_low_10M,NUM_REPLICATES=3 run_grandma_split.sh
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
GRANDMA_SPLIT_SCRIPT="/groups/itay_mayrose/ronenshtein/Grampa_revamp/run.py"

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
# Array task IDs: 1-105 (for 21 networks × 5 replicates)
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

# Use same input as GRAMPA
INPUT_DIR="${BASE_DIR}/${network}/processed/${CONFIG}/grampa_input/replicate_${replicate}"
# Output to grandma_split folder
OUTPUT_DIR="${BASE_DIR}/${network}/results/${CONFIG}/grandma_split/replicate_${replicate}"

GENE_TREES="${INPUT_DIR}/grampa_trees.tre"
SPECIES_TREE="${INPUT_DIR}/species.tre"

# ============================================================================
# HEADER
# ============================================================================

echo "============================================================================"
echo "RUN GRANDMA_SPLIT"
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
    echo "Did you run run_astral.sh first?"
    exit 1
fi

if [ ! -f "$GRANDMA_SPLIT_SCRIPT" ]; then
    echo "ERROR: GRANDMA_SPLIT script not found: $GRANDMA_SPLIT_SCRIPT"
    exit 1
fi

echo "✓ Input files validated"
echo "  Gene trees: $GENE_TREES"
echo "  Species tree: $SPECIES_TREE"
echo "  GRANDMA_SPLIT script: $GRANDMA_SPLIT_SCRIPT"
echo ""

# ============================================================================
# ACTIVATE CONDA
# ============================================================================

source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

echo "✓ Conda environment activated: gene2net"

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo "ERROR: Python command not found. Check conda environment."
    exit 1
fi

echo "✓ Python found: $(which python)"
echo ""

# ============================================================================
# CREATE OUTPUT DIRECTORY
# ============================================================================

mkdir -p "$OUTPUT_DIR"
echo "✓ Output directory created: $OUTPUT_DIR"
echo ""

# ============================================================================
# RUN GRANDMA_SPLIT
# ============================================================================

echo "Starting GRANDMA_SPLIT..."
echo "============================================================================"
echo ""

start_time=$(date +%s)

python "$GRANDMA_SPLIT_SCRIPT" \
    -g "$GENE_TREES" \
    -s "$SPECIES_TREE" \
    -o "$OUTPUT_DIR" \
    -m split

exit_code=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "============================================================================"

if [ $exit_code -eq 0 ]; then
    echo "✓ GRANDMA_SPLIT COMPLETED SUCCESSFULLY"
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
    echo "✗ GRANDMA_SPLIT FAILED"
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
