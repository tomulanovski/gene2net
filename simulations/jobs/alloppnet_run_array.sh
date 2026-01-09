#!/bin/bash
#SBATCH --job-name=alloppnet_run
#SBATCH --array=1-8
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alloppnet_run_%x_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alloppnet_run_%x_%A_%a.err
#SBATCH --time=120:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public 

# ============================================================================
# ALLOPPNET_RUN_ARRAY.SH
# ============================================================================
# Runs BEAST MCMC inference for AlloppNET.
# Array job for 8 AlloppNET-compatible networks.
#
# Run steps:
#   1. Run BEAST (100M iterations, ~5 days)
#
# Post-processing (TreeAnnotator + copy number removal) is handled separately
# by postprocess_results.py after BEAST completes or times out.
#
# Prerequisites:
#   - alloppnet_prep_array.sh must be completed first
#   - alloppnet.XML must exist in output directory
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M,REPLICATE=1 alloppnet_run_array.sh
#
# To run subset (e.g., just network 3):
#   sbatch --array=3 --export=CONFIG=conf_ils_low_10M,REPLICATE=1 alloppnet_run_array.sh
#
# Environment variables:
#   CONFIG    - Configuration name (required, e.g., conf_ils_low_10M)
#   REPLICATE - Replicate number (required, 1-5)
# ============================================================================

set -eo pipefail

# Initialize LD_LIBRARY_PATH if not set (needed for conda activation)
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ============================================================================
# CONFIGURATION
# ============================================================================

# Required parameters
CONFIG="${CONFIG:?ERROR: CONFIG environment variable is required. Use --export=CONFIG=conf_name}"
REPLICATE="${REPLICATE:?ERROR: REPLICATE environment variable is required. Use --export=REPLICATE=N}"

# Number of networks (fixed)
NUM_NETWORKS=8

# ============================================================================
# PATHS
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
SCRIPTS_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/alloppnet"
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"

# 8 AlloppNET-compatible networks (max 2 copies in network topology)
networks=(
  "Bendiksby_2011"
  "Ding_2023"
  "Koenen_2020"
  "Liu_2023"
  "Shahrestani_2015"
  "Wisecaver_2023"
  "Wu_2015"
  "Zhao_2021"
)

# ============================================================================
# CALCULATE NETWORK FROM ARRAY TASK ID
# ============================================================================
# Array task IDs: 1-8 (for 8 networks)
# Mapping: task_id = network_idx + 1

task_id=${SLURM_ARRAY_TASK_ID}
network_idx=$((task_id - 1))

# Validate network index
if [ $network_idx -ge $NUM_NETWORKS ]; then
    echo "ERROR: Invalid task ID ${task_id} - network index ${network_idx} exceeds ${NUM_NETWORKS}"
    exit 1
fi

network="${networks[$network_idx]}"

# ============================================================================
# PATHS FOR THIS JOB
# ============================================================================

OUTPUT_DIR="${BASE_DIR}/${network}/results/${CONFIG}/alloppnet/replicate_${REPLICATE}"

# ============================================================================
# HEADER
# ============================================================================

echo "============================================================================"
echo "AlloppNET RUN (BEAST Execution)"
echo "============================================================================"
echo "Network: ${network} (index: ${network_idx})"
echo "Replicate: ${REPLICATE}"
echo "Task ID: ${task_id}"
echo "Configuration: ${CONFIG}"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo "Date: $(date)"
echo "============================================================================"
echo ""

# ============================================================================
# VALIDATION
# ============================================================================

echo "Validating prerequisites..."

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "ERROR: Output directory not found: $OUTPUT_DIR"
    echo "Did you run alloppnet_prep_array.sh first?"
    exit 1
fi

if [ ! -f "${OUTPUT_DIR}/alloppnet.XML" ]; then
    echo "ERROR: BEAST XML file not found: ${OUTPUT_DIR}/alloppnet.XML"
    echo "Did you run alloppnet_prep_array.sh first?"
    exit 1
fi

echo "✓ Validation passed"
echo "  Found alloppnet.XML"
echo ""

# ============================================================================
# STEP 1: RUN BEAST
# ============================================================================

echo "============================================================================"
echo "[Step 1/1] Running BEAST"
echo "============================================================================"
echo "  - 100M iterations (~5 days)"
echo "  - Output: sampledmultrees.txt, sampledgtrees*.txt, sampledparams.txt"
echo ""

# Activate alloppnet environment (has BEAST)
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate alloppnet || {
    echo "ERROR: Could not activate alloppnet environment"
    echo "Please create environment: conda create -n alloppnet -c bioconda beast=1.8.4"
    exit 1
}

# Check if BEAST is available
if ! command -v beast &> /dev/null; then
    echo "ERROR: beast command not found. Check if it's installed in alloppnet environment."
    exit 1
fi

echo "✓ BEAST found: $(which beast)"
echo ""

start_time=$(date +%s)

# Run BEAST (must be in output directory)
cd "${OUTPUT_DIR}"
beast alloppnet.XML

exit_code=$?
end_time=$(date +%s)
duration=$((end_time - start_time))
hours=$(echo "scale=2; $duration/3600" | bc)

if [ $exit_code -ne 0 ]; then
    echo ""
    echo "ERROR: BEAST run failed (exit code: ${exit_code})"
    echo "Duration: ${duration}s (${hours} hours)"
    exit 1
fi

echo ""
echo "✓ BEAST run complete (${duration}s = ${hours} hours)"
echo ""
echo "Note: Post-processing (TreeAnnotator + copy number removal) will be"
echo "      handled by postprocess_results.py after BEAST completes."
echo ""

# ============================================================================
# VALIDATION OF OUTPUTS
# ============================================================================

echo "Validating outputs..."

# Check required output files (BEAST outputs only)
required_files=(
    "${OUTPUT_DIR}/sampledmultrees.txt"
    "${OUTPUT_DIR}/sampledparams.txt"
)

all_present=true
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file missing: $file"
        all_present=false
    fi
done

if [ "$all_present" = false ]; then
    echo "ERROR: Some required files are missing"
    exit 1
fi

# Count trees for summary
TOTAL_TREES=$(grep -c "tree STATE_" "${OUTPUT_DIR}/sampledmultrees.txt" || echo "0")

echo "✓ All required files present"
echo "  Found ${TOTAL_TREES} trees in sampledmultrees.txt"
echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo "============================================================================"
echo "ALLOPPNET RUN COMPLETE"
echo "============================================================================"
echo "Network: ${network}"
echo "Replicate: ${REPLICATE}"
echo "Configuration: ${CONFIG}"
echo ""
echo "Output files (${OUTPUT_DIR}):"
echo "  - alloppnet.XML (BEAST configuration)"
echo "  - sampledmultrees.txt (species networks, ${TOTAL_TREES} trees)"
echo "  - sampledgtrees*.txt (gene trees, 1000 files)"
echo "  - sampledparams.txt (MCMC parameters)"
echo "  - DBUGTUNE.txt (debug info)"
echo ""
echo "Next steps:"
echo "  1. Post-process results:"
echo "     python simulations/scripts/postprocess_results.py ${CONFIG} --methods alloppnet"
echo "  2. Validate final output:"
echo "     python simulations/scripts/check_pipeline_status.py ${CONFIG} --method alloppnet --step run"
echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

exit 0
