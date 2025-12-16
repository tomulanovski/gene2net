#!/bin/bash
#SBATCH --job-name=alloppnet
#SBATCH --array=1-8
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alloppnet_%x_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alloppnet_%x_%A_%a.err
#SBATCH --time=120:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# ALLOPPNET_ARRAY.SH
# ============================================================================
# Runs AlloppNET phylogenetic network inference on simulation data.
# Array job for 8 AlloppNET-compatible networks.
#
# Pipeline steps:
#   1. Prepare input (PHY → NEXUS, ploidy_level.json, taxa_table.txt)
#   2. Generate BEAST XML
#   3. Run BEAST (100M iterations, ~5 days)
#   4. Summarize with TreeAnnotator
#   5. Post-process (remove copy number suffixes)
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M,REPLICATE=1 alloppnet_array.sh
#
# To run subset (e.g., just network 3):
#   sbatch --array=3 --export=CONFIG=conf_ils_low_10M,REPLICATE=1 alloppnet_array.sh
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

ALIGNMENT_DIR="${BASE_DIR}/${network}/data/${CONFIG}/replicate_${REPLICATE}/1/alignments"
INPUT_DIR="${BASE_DIR}/${network}/processed/${CONFIG}/alloppnet_input/replicate_${REPLICATE}"
OUTPUT_DIR="${BASE_DIR}/${network}/results/${CONFIG}/alloppnet/replicate_${REPLICATE}"

# ============================================================================
# HEADER
# ============================================================================

echo "============================================================================"
echo "AlloppNET PHYLOGENETIC NETWORK INFERENCE"
echo "============================================================================"
echo "Network: ${network} (index: ${network_idx})"
echo "Replicate: ${REPLICATE}"
echo "Task ID: ${task_id}"
echo "Configuration: ${CONFIG}"
echo ""
echo "Alignment directory: ${ALIGNMENT_DIR}"
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Date: $(date)"
echo "============================================================================"
echo ""

# ============================================================================
# VALIDATION
# ============================================================================

echo "Validating paths..."

if [ ! -d "$ALIGNMENT_DIR" ]; then
    echo "ERROR: Alignment directory not found: $ALIGNMENT_DIR"
    echo "Did you run SimPhy sequence simulation first?"
    exit 1
fi

# Count alignment files
NUM_ALIGNMENTS=$(find "$ALIGNMENT_DIR" -name "alignment_*.phy" 2>/dev/null | wc -l)
if [ $NUM_ALIGNMENTS -eq 0 ]; then
    echo "ERROR: No alignment files found in $ALIGNMENT_DIR"
    exit 1
fi

echo "✓ Validation passed"
echo "  Found ${NUM_ALIGNMENTS} alignment files"
echo ""

# ============================================================================
# CREATE DIRECTORIES
# ============================================================================

mkdir -p "$INPUT_DIR"
mkdir -p "$OUTPUT_DIR"

echo "✓ Directories created"
echo ""

# ============================================================================
# STEP 1: PREPARE ALLOPPNET INPUT
# ============================================================================

echo "============================================================================"
echo "[Step 1/5] Preparing AlloppNET input"
echo "============================================================================"
echo "  - Converting PHY → NEXUS"
echo "  - Generating ploidy_level.json"
echo "  - Creating taxa_table.txt"
echo ""

start_time=$(date +%s)

# Activate gene2net environment for Python
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

python3 "${SCRIPTS_DIR}/prepare_alloppnet_input.py" \
    --network "${network}" \
    --config "${CONFIG}" \
    --replicate "${REPLICATE}" \
    --alignment-dir "${ALIGNMENT_DIR}" \
    --output-dir "${INPUT_DIR}" \
    --verbose

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to prepare AlloppNET input"
    exit 1
fi

end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "✓ Input preparation complete (${duration}s)"
echo ""

# ============================================================================
# STEP 2: GENERATE BEAST XML
# ============================================================================

echo "============================================================================"
echo "[Step 2/5] Generating BEAST XML"
echo "============================================================================"
echo "  - Using AlloppDT scripts"
echo "  - 100M BEAST iterations"
echo ""

start_time=$(date +%s)

Rscript "${SCRIPTS_DIR}/generate_beast_xml.r" \
    "${INPUT_DIR}" \
    "${OUTPUT_DIR}"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to generate BEAST XML"
    exit 1
fi

end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "✓ BEAST XML generated (${duration}s)"
echo ""

# ============================================================================
# STEP 3: RUN BEAST
# ============================================================================

echo "============================================================================"
echo "[Step 3/5] Running BEAST"
echo "============================================================================"
echo "  - 100M iterations (~5 days)"
echo "  - Output: sampledmultrees.txt, sampledgtrees*.txt, sampledparams.txt"
echo ""

# Activate alloppnet environment (has BEAST)
conda activate alloppnet || {
    echo "ERROR: Could not activate alloppnet environment"
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

# ============================================================================
# STEP 4: SUMMARIZE WITH TREEANNOTATOR
# ============================================================================

echo "============================================================================"
echo "[Step 4/5] Summarizing with TreeAnnotator"
echo "============================================================================"
echo "  - 10% burnin"
echo "  - Mean node heights"
echo ""

start_time=$(date +%s)

# Check if TreeAnnotator is available
if ! command -v treeannotator &> /dev/null; then
    echo "ERROR: treeannotator command not found. Check if it's installed in alloppnet environment."
    exit 1
fi

# Count trees and calculate 10% burnin
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

treeannotator -burninTrees ${BURNIN} -heights mean \
    "${OUTPUT_DIR}/sampledmultrees.txt" \
    "${OUTPUT_DIR}/alloppnet_consensus.tre"

if [ $? -ne 0 ]; then
    echo "ERROR: TreeAnnotator failed"
    exit 1
fi

end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "✓ Consensus tree generated (${duration}s)"
echo ""

# ============================================================================
# STEP 5: POST-PROCESSING
# ============================================================================

echo "============================================================================"
echo "[Step 5/5] Post-processing consensus tree"
echo "============================================================================"
echo "  - Removing copy number suffixes (_0, _1)"
echo ""

start_time=$(date +%s)

# Activate gene2net environment for Python
conda activate gene2net

python3 "${SCRIPTS_DIR}/remove_copy_numbers.py" \
    "${OUTPUT_DIR}/alloppnet_consensus.tre" \
    "${OUTPUT_DIR}/alloppnet_final.tre" \
    --verbose

if [ $? -ne 0 ]; then
    echo "ERROR: Post-processing failed"
    exit 1
fi

end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "✓ Post-processing complete (${duration}s)"
echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo "============================================================================"
echo "ALLOPPNET ANALYSIS COMPLETE"
echo "============================================================================"
echo "Network: ${network}"
echo "Replicate: ${REPLICATE}"
echo "Configuration: ${CONFIG}"
echo ""
echo "Output files:"
echo "  - alloppnet.XML (BEAST configuration)"
echo "  - sampledmultrees.txt (species networks)"
echo "  - sampledgtrees*.txt (gene trees)"
echo "  - sampledparams.txt (MCMC parameters)"
echo "  - DBUGTUNE.txt (debug info)"
echo "  - alloppnet_consensus.tre (TreeAnnotator consensus)"
echo "  - alloppnet_final.tre (final MUL-tree for comparison)"
echo ""
echo "Final tree: ${OUTPUT_DIR}/alloppnet_final.tre"
echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

exit 0
