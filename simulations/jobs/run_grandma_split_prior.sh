#!/bin/bash
#SBATCH --job-name=run_grandma_split_prior
#SBATCH --array=1-105
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_grandma_split_prior_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_grandma_split_prior_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# RUN_GRANDMA_SPLIT_PRIOR.SH
# ============================================================================
# Runs GRANDMA_SPLIT with ploidy prior (-x) on simulated data.
# Uses the same gene trees / species tree as run_grandma_split.sh
# (from grampa_input), and derives the ploidy prior from the Polyphest
# multiset file (polyphest_input/replicate_N/multi_set.txt).
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M run_grandma_split_prior.sh
#
# Environment variables:
#   CONFIG          - Configuration name (required, e.g., conf_ils_low_10M)
#   NUM_REPLICATES  - Number of replicates (default: 5)
# ============================================================================

set -eo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ============================================================================
# CONFIGURATION
# ============================================================================

CONFIG="${CONFIG:?ERROR: CONFIG environment variable is required. Use --export=CONFIG=conf_name}"
NUM_REPLICATES="${NUM_REPLICATES:-5}"
NUM_NETWORKS=21

# ============================================================================
# PATHS
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"
GRANDMA_SPLIT_SCRIPT="/groups/itay_mayrose/ronenshtein/Grandma/gran.py"
GENERATE_PLOIDY_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/generate_ploidy_file.py"

networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# ============================================================================
# CALCULATE NETWORK AND REPLICATE FROM ARRAY TASK ID
# ============================================================================

task_id=${SLURM_ARRAY_TASK_ID}
network_idx=$(( (task_id - 1) / NUM_REPLICATES ))
replicate=$(( ((task_id - 1) % NUM_REPLICATES) + 1 ))

if [ $network_idx -ge $NUM_NETWORKS ]; then
    echo "ERROR: Invalid task ID ${task_id} - network index ${network_idx} exceeds ${NUM_NETWORKS}"
    exit 1
fi

network="${networks[$network_idx]}"

# ============================================================================
# PATHS FOR THIS JOB
# ============================================================================

GRAMPA_INPUT_DIR="${BASE_DIR}/${network}/processed/${CONFIG}/grampa_input/replicate_${replicate}"
POLYPHEST_INPUT_DIR="${BASE_DIR}/${network}/processed/${CONFIG}/polyphest_input/replicate_${replicate}"
OUTPUT_DIR="${BASE_DIR}/${network}/results/${CONFIG}/grandma_split_prior/replicate_${replicate}"

GENE_TREES="${GRAMPA_INPUT_DIR}/grampa_trees.tre"
SPECIES_TREE="${GRAMPA_INPUT_DIR}/species.tre"
MULTISET="${POLYPHEST_INPUT_DIR}/multi_set.txt"
PLOIDY_FILE="${GRAMPA_INPUT_DIR}/split_ploidies.txt"

# ============================================================================
# HEADER
# ============================================================================

echo "============================================================================"
echo "RUN GRANDMA_SPLIT WITH PLOIDY PRIOR"
echo "============================================================================"
echo "Network:          ${network} (index: ${network_idx})"
echo "Replicate:        ${replicate}"
echo "Task ID:          ${task_id}"
echo "Configuration:    ${CONFIG}"
echo ""
echo "Gene trees:       ${GENE_TREES}"
echo "Species tree:     ${SPECIES_TREE}"
echo "Multiset:         ${MULTISET}"
echo "Ploidy file:      ${PLOIDY_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Date:             $(date)"
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

if [ ! -f "$MULTISET" ]; then
    echo "ERROR: Multiset file not found: $MULTISET"
    echo "Did you run prep_polyphest.sh first?"
    exit 1
fi

if [ ! -f "$GRANDMA_SPLIT_SCRIPT" ]; then
    echo "ERROR: GRANDMA_SPLIT script not found: $GRANDMA_SPLIT_SCRIPT"
    exit 1
fi

echo "✓ Input files validated"
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
echo "✓ Python found: $(which python)"
echo ""

# ============================================================================
# GENERATE PLOIDY FILE FROM MULTISET
# ============================================================================

echo "Generating ploidy file from multiset..."
python "$GENERATE_PLOIDY_SCRIPT" "$MULTISET" "$PLOIDY_FILE"
echo "✓ Ploidy file written: ${PLOIDY_FILE}"
echo ""

# ============================================================================
# CREATE OUTPUT DIRECTORY
# ============================================================================

mkdir -p "$OUTPUT_DIR"
echo "✓ Output directory created: $OUTPUT_DIR"
echo ""

# ============================================================================
# RUN GRANDMA_SPLIT WITH PLOIDY PRIOR
# ============================================================================

echo "Starting GRANDMA_SPLIT with ploidy prior..."
echo "============================================================================"
echo ""

start_time=$(date +%s)

python "$GRANDMA_SPLIT_SCRIPT" \
    -g "$GENE_TREES" \
    -s "$SPECIES_TREE" \
    -o "$OUTPUT_DIR" \
    -m split \
    -x "$PLOIDY_FILE" \
    --strict_constraint \
    -c 15 \
    -i 50 \
    -p ${SLURM_CPUS_PER_TASK} \
    --plot \
    --debug \
    --v 3 \
    --overwrite

exit_code=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "============================================================================"

if [ $exit_code -eq 0 ]; then
    echo "✓ GRANDMA_SPLIT (PRIOR) COMPLETED SUCCESSFULLY"
    echo ""
    echo "Network:   ${network}"
    echo "Replicate: ${replicate}"
    echo "Duration:  ${duration} seconds ($(echo "scale=2; $duration/3600" | bc) hours)"
    echo "Results:   ${OUTPUT_DIR}"
    echo ""
    echo "Output files:"
    ls -la "$OUTPUT_DIR" 2>/dev/null || echo "  (unable to list)"
else
    echo "✗ GRANDMA_SPLIT (PRIOR) FAILED"
    echo ""
    echo "Network:   ${network}"
    echo "Replicate: ${replicate}"
    echo "Exit code: ${exit_code}"
    echo "Duration:  ${duration} seconds"
fi

echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

exit $exit_code
