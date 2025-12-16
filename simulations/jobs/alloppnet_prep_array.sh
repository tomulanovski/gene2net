#!/bin/bash
#SBATCH --job-name=alloppnet_prep
#SBATCH --array=1-8
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alloppnet_prep_%x_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alloppnet_prep_%x_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# ALLOPPNET_PREP_ARRAY.SH
# ============================================================================
# Prepares AlloppNET input files and generates BEAST XML.
# Array job for 8 AlloppNET-compatible networks.
#
# Prep steps:
#   1. Convert PHY alignments to NEXUS format
#   2. Generate ploidy_level.json (using kernel smoothing)
#   3. Create taxa_table.txt (homeolog pairing)
#   4. Generate BEAST XML using AlloppDT scripts
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M,REPLICATE=1 alloppnet_prep_array.sh
#
# To run subset (e.g., just network 3):
#   sbatch --array=3 --export=CONFIG=conf_ils_low_10M,REPLICATE=1 alloppnet_prep_array.sh
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
echo "AlloppNET PREPARATION (Prep Step)"
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
echo "[Step 1/2] Preparing AlloppNET input"
echo "============================================================================"
echo "  - Converting PHY → NEXUS (1000 alignments)"
echo "  - Analyzing copy number distributions (kernel smoothing)"
echo "  - Generating ploidy_level.json (robust ploidy detection)"
echo "  - Creating taxa_table.txt (homeolog pairing)"
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
    --kernel-width 2 \
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
echo "[Step 2/2] Generating BEAST XML"
echo "============================================================================"
echo "  - Using AlloppDT scripts"
echo "  - Reading NEXUS alignments and taxa_table.txt"
echo "  - Configuring 100M BEAST iterations"
echo ""

start_time=$(date +%s)

# Check that AlloppDT scripts exist
if [ ! -f "${SCRIPTS_DIR}/AlloppDT_5beastxml_toplevel.r" ]; then
    echo "ERROR: AlloppDT_5beastxml_toplevel.r not found in ${SCRIPTS_DIR}"
    echo "Please copy AlloppDT scripts first (see ALLOPPNET_GUIDE.md)"
    exit 1
fi

if [ ! -f "${SCRIPTS_DIR}/AlloppDT_6beastxml_bits.r" ]; then
    echo "ERROR: AlloppDT_6beastxml_bits.r not found in ${SCRIPTS_DIR}"
    echo "Please copy AlloppDT scripts first (see ALLOPPNET_GUIDE.md)"
    exit 1
fi

# Activate alloppnet environment for Rscript
conda activate alloppnet || {
    echo "ERROR: Could not activate alloppnet environment"
    echo "Please create environment: conda create -n alloppnet -c conda-forge r-base"
    exit 1
}

# Add trailing slash to INPUT_DIR for AlloppDT path construction
# AlloppDT concatenates paths without separator, so we need the trailing /
INPUT_DIR_WITH_SLASH="${INPUT_DIR}/"

Rscript "${SCRIPTS_DIR}/generate_beast_xml.r" \
    "${INPUT_DIR_WITH_SLASH}" \
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
# VALIDATION OF OUTPUTS
# ============================================================================

echo "Validating outputs..."

# Check NEXUS files
NUM_NEX=$(find "${INPUT_DIR}" -name "*.nex" 2>/dev/null | wc -l)
if [ $NUM_NEX -ne 1000 ]; then
    echo "WARNING: Expected 1000 NEXUS files, found ${NUM_NEX}"
fi

# Check required files
required_files=(
    "${INPUT_DIR}/ploidy_level.json"
    "${INPUT_DIR}/taxa_table.txt"
    "${OUTPUT_DIR}/alloppnet.XML"
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

echo "✓ All required files present"
echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo "============================================================================"
echo "ALLOPPNET PREP COMPLETE"
echo "============================================================================"
echo "Network: ${network}"
echo "Replicate: ${REPLICATE}"
echo "Configuration: ${CONFIG}"
echo ""
echo "Input files (${INPUT_DIR}):"
echo "  - ${NUM_NEX} NEXUS alignments"
echo "  - ploidy_level.json (kernel smoothing)"
echo "  - taxa_table.txt (homeolog pairing)"
echo ""
echo "Output files (${OUTPUT_DIR}):"
echo "  - alloppnet.XML (BEAST configuration, ready to run)"
echo ""
echo "Next step:"
echo "  Run BEAST with alloppnet_run_array.sh (takes ~5 days)"
echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

exit 0
