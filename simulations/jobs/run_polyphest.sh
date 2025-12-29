#!/bin/bash
#SBATCH --job-name=run_polyphest
#SBATCH --array=1-105
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_polyphest_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_polyphest_%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=50g
#SBATCH --cpus-per-task=16
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# RUN_POLYPHEST.SH
# ============================================================================
# Runs Polyphest on prepared gene trees and multi-set files.
# Uses a 2D array structure: 21 networks � 5 replicates = 105 jobs
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M run_polyphest.sh
#   sbatch --export=CONFIG=conf_ils_low_10M,NUM_REPLICATES=3 run_polyphest.sh
#
# To run with different array size (e.g., 3 replicates = 63 jobs):
#   sbatch --array=1-63 --export=CONFIG=conf_ils_low_10M,NUM_REPLICATES=3 run_polyphest.sh
#
# Environment variables:
#   CONFIG          - Configuration name (required, e.g., conf_ils_low_10M)
#   NUM_REPLICATES  - Number of replicates (default: 5)
#   PERCENTILE      - Polyphest percentile parameter (default: 60)
#   ISO_THRESHOLD   - Isomorphic threshold (default: 0.2)
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
PERCENTILE="${PERCENTILE:-60}"
ISO_THRESHOLD="${ISO_THRESHOLD:-0.2}"

# Number of networks (fixed)
NUM_NETWORKS=21

# ============================================================================
# PATHS
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"
POLYPHEST="/groups/itay_mayrose/tomulanovski/Polyphest/Polyphest.py"

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
# Array task IDs: 1-105 (for 21 networks � 5 replicates)
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

INPUT_DIR="${BASE_DIR}/${network}/processed/${CONFIG}/polyphest_input/replicate_${replicate}"
OUTPUT_DIR="${BASE_DIR}/${network}/results/${CONFIG}/polyphest_p${PERCENTILE}/replicate_${replicate}"

GENE_TREES="${INPUT_DIR}/polyphest_trees.tre"
MULTISET="${INPUT_DIR}/multi_set.txt"

# ============================================================================
# HEADER
# ============================================================================

echo "============================================================================"
echo "RUN POLYPHEST"
echo "============================================================================"
echo "Network: ${network} (index: ${network_idx})"
echo "Replicate: ${replicate}"
echo "Task ID: ${task_id}"
echo "Configuration: ${CONFIG}"
echo ""
echo "Parameters:"
echo "  Percentile: ${PERCENTILE}"
echo "  Isomorphic threshold: ${ISO_THRESHOLD}"
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
    echo "Did you run prep_polyphest.sh first?"
    exit 1
fi

if [ ! -f "$MULTISET" ]; then
    echo "ERROR: Multi-set file not found: $MULTISET"
    echo "Did you run prep_polyphest.sh first?"
    exit 1
fi

if [ ! -f "$POLYPHEST" ]; then
    echo "ERROR: Polyphest script not found: $POLYPHEST"
    exit 1
fi

echo "? Input files validated"
echo "  Gene trees: $GENE_TREES"
echo "  Multi-set: $MULTISET"
echo ""

# ============================================================================
# ACTIVATE CONDA
# ============================================================================

source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate polyphest || {
    echo "ERROR: Could not activate polyphest environment"
    exit 1
}

echo "? Conda environment activated: polyphest"
echo "  Python: $(which python)"
echo ""

# ============================================================================
# CREATE OUTPUT DIRECTORY
# ============================================================================

mkdir -p "$OUTPUT_DIR"
echo "? Output directory created: $OUTPUT_DIR"
echo ""

# ============================================================================
# RUN POLYPHEST
# ============================================================================

echo "Starting Polyphest..."
echo "============================================================================"
echo ""

start_time=$(date +%s)

python "$POLYPHEST" \
    --gene_tree_file "$GENE_TREES" \
    --consensus_multiset_file "$MULTISET" \
    --filter_strategy percentile \
    --output_dir "$OUTPUT_DIR" \
    --percentile "$PERCENTILE" \
    --use_near_isomorphic "True" \
    --isomorphic_threshold "$ISO_THRESHOLD"

exit_code=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "============================================================================"

if [ $exit_code -eq 0 ]; then
    echo "? POLYPHEST COMPLETED SUCCESSFULLY"
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
    echo "? POLYPHEST FAILED"
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