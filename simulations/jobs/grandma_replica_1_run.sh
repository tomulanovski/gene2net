#!/bin/bash
#SBATCH --job-name=grandma_sim
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/grandma_run_rep1_ILS_low_dup_low_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/grandma_run_rep1_ILS_low_dup_low_job_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
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

# Activate gene2net environment
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# Verify the environment was activated
echo "Current conda environment: $CONDA_PREFIX"

# Check if GRAMPA is available
if ! command -v grampa &> /dev/null; then
    echo "ERROR: grampa command not found. Check if it's installed in your conda environment."
    exit 1
fi

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

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
echo "Running GRAMPA for: ${network} - ILS_low_dup_low"
echo "Processing replicate: 1"
echo "=========================================="

# Define paths for replicate 1
INPUT_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/grampa_input"
OUTPUT_DIR="${BASE_DIR}/${network}/results/ILS_low_dup_low/grandma"

# Input files for replicate 1
GENE_TREES="${INPUT_DIR}/1_trees.tre"
SPECIES_TREE="${INPUT_DIR}/1_species.tre"

# Check if input files exist
if [ ! -f "$GENE_TREES" ]; then
    echo "ERROR: Gene trees file not found: $GENE_TREES"
    exit 1
fi

if [ ! -f "$SPECIES_TREE" ]; then
    echo "ERROR: Species tree file not found: $SPECIES_TREE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Gene trees: $GENE_TREES"
echo "Species tree: $SPECIES_TREE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Run GRAMPA
grampa -g "$GENE_TREES" -s "$SPECIES_TREE" -o "$OUTPUT_DIR" --overwrite

python "/groups/itay_mayrose/ronenshtein/Grandma/grampa.py" \
    -s "$SPECIES_TREE" \
    -g "$GENE_TREES" \
    -o "$OUTPUT_DIR" \
    -i 0 -v 1

# Check if GRAMPA succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "GRANDMA completed successfully for ${network}"
    echo "Results saved to: $OUTPUT_DIR"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "ERROR: GRAMPA failed for ${network}"
    echo "=========================================="
    exit 1
fi

# Deactivate conda environment
conda deactivate