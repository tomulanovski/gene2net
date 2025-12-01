#!/bin/bash
#SBATCH --job-name=polyphest_sims
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/polyphest_run_rep1_ILS_low_dup_low_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/polyphest_run_rep1_ILS_low_dup_low_job_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=50g
#SBATCH --cpus-per-task=16
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Check if conda is properly sourced
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda command not found. Check your conda installation path."
    exit 1
fi

# Activate polyphest environment
conda activate polyphest || {
    echo "ERROR: Could not activate polyphest environment"
    exit 1
}

# Verify the environment was activated
echo "Current conda environment: $CONDA_PREFIX"

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Polyphest script path
POLYPHEST="/groups/itay_mayrose/tomulanovski/Polyphest/Polyphest.py"

# Array of networks (must match the original array)
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# Get the network for this array task
network=${networks[$SLURM_ARRAY_TASK_ID-1]}

echo "=========================================="
echo "Running Polyphest for: ${network} - ILS_low_dup_low"
echo "Processing replicate: 1"
echo "=========================================="

# Define paths for replicate 1
INPUT_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/polyphest_input"
OUTPUT_DIR="${BASE_DIR}/${network}/results/ILS_low_dup_low/polyphest"

# Input files for replicate 1
GENE_TREES="${INPUT_DIR}/1_trees.tre"
MULTISET="${INPUT_DIR}/1_multi_set.txt"

# Check if input files exist
if [ ! -f "$GENE_TREES" ]; then
    echo "ERROR: Gene trees file not found: $GENE_TREES"
    exit 1
fi

if [ ! -f "$MULTISET" ]; then
    echo "ERROR: Multiset file not found: $MULTISET"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Input multiset: $MULTISET"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Run Polyphest
python $POLYPHEST \
    --gene_tree_file $GENE_TREES \
    --consensus_multiset_file $MULTISET \
    --filter_strategy percentile \
    --output_dir $OUTPUT_DIR \
    --percentile 60 \
    --use_near_isomorphic "True" \
    --isomorphic_threshold 0.2

# Check if Polyphest succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Polyphest completed successfully for ${network}"
    echo "Results saved to: $OUTPUT_DIR"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "ERROR: Polyphest failed for ${network}"
    echo "=========================================="
    exit 1
fi