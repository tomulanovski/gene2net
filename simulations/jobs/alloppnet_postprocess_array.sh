#!/bin/bash
#SBATCH --job-name=alloppnet_postprocess
#SBATCH --array=1-8
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alloppnet_postprocess_%x_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alloppnet_postprocess_%x_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=power-general-public-pool
#SBATCH --qos=public

# ============================================================================
# ALLOPPNET_POSTPROCESS_ARRAY.SH
# ============================================================================
# Post-processes AlloppNET BEAST output (TreeAnnotator + copy number removal).
# Array job for 8 AlloppNET-compatible networks.
#
# Post-processing steps:
#   1. Run TreeAnnotator on sampledmultrees.txt (10% burnin, mean heights)
#   2. Remove copy number suffixes (_0, _1) using remove_copy_numbers.py
#   3. Create alloppnet_result.tre for summary pipeline
#
# Prerequisites:
#   - alloppnet_run_array.sh must be completed (or timed out)
#   - sampledmultrees.txt must exist in output directory
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M,REPLICATE=1 alloppnet_postprocess_array.sh
#
# To run subset (e.g., just network 3):
#   sbatch --array=3 --export=CONFIG=conf_ils_low_10M,REPLICATE=1 alloppnet_postprocess_array.sh
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
echo "AlloppNET Post-Processing"
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
    echo "Did you run alloppnet_run_array.sh first?"
    exit 1
fi

if [ ! -f "${OUTPUT_DIR}/sampledmultrees.txt" ]; then
    echo "ERROR: sampledmultrees.txt not found: ${OUTPUT_DIR}/sampledmultrees.txt"
    echo "Did you run alloppnet_run_array.sh first?"
    exit 1
fi

# Check if already processed
if [ -f "${OUTPUT_DIR}/alloppnet_result.tre" ]; then
    echo "WARNING: alloppnet_result.tre already exists. Skipping post-processing."
    echo "  To reprocess, delete: ${OUTPUT_DIR}/alloppnet_result.tre"
    exit 0
fi

echo "✓ Validation passed"
echo "  Found sampledmultrees.txt"
echo ""

# ============================================================================
# STEP 1: COUNT TREES AND CALCULATE BURNIN
# ============================================================================

echo "============================================================================"
echo "[Step 1/3] Counting trees and calculating burnin"
echo "============================================================================"

TOTAL_TREES=$(grep -c "tree STATE_" "${OUTPUT_DIR}/sampledmultrees.txt" || echo "0")
if [ $TOTAL_TREES -eq 0 ]; then
    echo "ERROR: No trees found in sampledmultrees.txt"
    exit 1
fi

USABLE_TREES=$((TOTAL_TREES - 1))  # Exclude STATE_0
BURNIN=$((USABLE_TREES / 10))

echo "  Total trees: ${TOTAL_TREES}"
echo "  Usable trees: ${USABLE_TREES}"
echo "  Burnin: ${BURNIN} trees (10%)"
echo ""

# ============================================================================
# STEP 2: RUN TREEANNOTATOR
# ============================================================================

echo "============================================================================"
echo "[Step 2/3] Running TreeAnnotator"
echo "============================================================================"
echo "  - 10% burnin (${BURNIN} trees)"
echo "  - Mean node heights"
echo ""

start_time=$(date +%s)

# Activate alloppnet environment (has TreeAnnotator)
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate alloppnet || {
    echo "ERROR: Could not activate alloppnet environment"
    exit 1
}

# Check if TreeAnnotator is available
if ! command -v treeannotator &> /dev/null; then
    echo "ERROR: treeannotator command not found. Check if it's installed in alloppnet environment."
    exit 1
fi

echo "✓ TreeAnnotator found: $(which treeannotator)"
echo ""

# Create consensus tree file
consensus_file="${OUTPUT_DIR}/alloppnet_consensus.tre"

# Run TreeAnnotator
treeannotator -burninTrees ${BURNIN} -heights mean \
    "${OUTPUT_DIR}/sampledmultrees.txt" \
    "${consensus_file}"

if [ $? -ne 0 ]; then
    echo "ERROR: TreeAnnotator failed"
    exit 1
fi

if [ ! -f "${consensus_file}" ]; then
    echo "ERROR: Consensus tree file not created: ${consensus_file}"
    exit 1
fi

end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "✓ Consensus tree generated (${duration}s)"
echo ""

# ============================================================================
# STEP 3: REMOVE COPY NUMBER SUFFIXES
# ============================================================================

echo "============================================================================"
echo "[Step 3/3] Removing copy number suffixes"
echo "============================================================================"
echo "  - Removing _0 and _1 suffixes from tip labels"
echo ""

start_time=$(date +%s)

# Activate gene2net environment for Python
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# Construct ploidy JSON path
# Path structure: ${BASE_DIR}/${network}/processed/${CONFIG}/alloppnet_input/replicate_${REPLICATE}/ploidy_level.json
PLOIDY_JSON="${BASE_DIR}/${network}/processed/${CONFIG}/alloppnet_input/replicate_${REPLICATE}/ploidy_level.json"

if [ ! -f "${PLOIDY_JSON}" ]; then
    echo "WARNING: Ploidy JSON not found: ${PLOIDY_JSON}"
    echo "         Will attempt auto-detection or use fallback method"
    # Run without --ploidy-json (will use auto-detection)
    python3 "${SCRIPTS_DIR}/remove_copy_numbers.py" \
        "${consensus_file}" \
        "${OUTPUT_DIR}/alloppnet_result.tre" \
        --verbose
else
    echo "  Using ploidy JSON: ${PLOIDY_JSON}"
    # Run with --ploidy-json
    python3 "${SCRIPTS_DIR}/remove_copy_numbers.py" \
        "${consensus_file}" \
        "${OUTPUT_DIR}/alloppnet_result.tre" \
        --ploidy-json "${PLOIDY_JSON}" \
        --verbose
fi

if [ $? -ne 0 ]; then
    echo "ERROR: Copy number removal failed"
    exit 1
fi

if [ ! -f "${OUTPUT_DIR}/alloppnet_result.tre" ]; then
    echo "ERROR: Final result file not created: ${OUTPUT_DIR}/alloppnet_result.tre"
    exit 1
fi

end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "✓ Copy number removal complete (${duration}s)"
echo ""

# ============================================================================
# VALIDATION OF OUTPUTS
# ============================================================================

echo "Validating outputs..."

if [ ! -f "${OUTPUT_DIR}/alloppnet_result.tre" ]; then
    echo "ERROR: Required file missing: alloppnet_result.tre"
    exit 1
fi

# Check file is non-empty
if [ ! -s "${OUTPUT_DIR}/alloppnet_result.tre" ]; then
    echo "ERROR: alloppnet_result.tre is empty"
    exit 1
fi

echo "✓ All required files present"
echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo "============================================================================"
echo "ALLOPPNET POST-PROCESSING COMPLETE"
echo "============================================================================"
echo "Network: ${network}"
echo "Replicate: ${REPLICATE}"
echo "Configuration: ${CONFIG}"
echo ""
echo "Output files (${OUTPUT_DIR}):"
echo "  - sampledmultrees.txt (${TOTAL_TREES} trees, input)"
echo "  - alloppnet_consensus.tre (TreeAnnotator consensus, intermediate)"
echo "  - alloppnet_result.tre (final MUL-tree for comparison)"
echo ""
echo "Final tree: ${OUTPUT_DIR}/alloppnet_result.tre"
echo ""
echo "Next step:"
echo "  Validate with: python simulations/scripts/check_pipeline_status.py ${CONFIG} --method alloppnet --step run"
echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

exit 0

