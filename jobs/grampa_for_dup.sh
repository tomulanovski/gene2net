#!/bin/bash
#SBATCH --job-name=grampa_array
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/grampa_dup/%x_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/grampa_dup/%x_%A_%a.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=power-general
#SBATCH --account=power-general-users
#SBATCH --array=0-10

# Define all datasets
DATASETS=("Koenen_2020" "Lawrence_2016" "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Popp_2005" "Wu_2015" "Ren_2024" "Zhao_2021" "Marcussen_2015" "Morales-Briones_2021")

# Create log directory if it doesn't exist
mkdir -p /groups/itay_mayrose/tomulanovski/gene2net/papers/grampa_dup

# Get current dataset based on array task ID
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}

echo "Processing dataset: $DATASET"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/papers"

# Set paths for current dataset
GENE_TREES="${BASE_DIR}/${DATASET}/gene_trees/grampa_trees.tre"
SPECIES_TREE="${BASE_DIR}/${DATASET}/gene_trees/grampa_species_tree.tre"
OUTPUT_DIR="${BASE_DIR}/${DATASET}/networks/grampa_for_duplication_analysis/"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Conda setup
CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Check if conda is properly sourced
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda command not found. Check your conda installation path."
    exit 1
fi

# Activate your specific environment
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# Verify the environment was activated
echo "Current conda environment: $CONDA_PREFIX"

# Check if input files exist
if [ ! -f "$GENE_TREES" ]; then
    echo "ERROR: Gene trees file not found: $GENE_TREES"
    exit 1
fi

if [ ! -f "$SPECIES_TREE" ]; then
    echo "ERROR: Species tree file not found: $SPECIES_TREE"
    exit 1
fi

# Check if GRAMPA is available
if ! command -v grampa &> /dev/null; then
    echo "ERROR: grampa command not found. Check if it's installed in your conda environment."
    exit 1
fi

echo "Running GRAMPA for $DATASET..."
echo "Gene trees: $GENE_TREES"
echo "Species tree: $SPECIES_TREE"
echo "Output directory: $OUTPUT_DIR"

# Run GRAMPA
export GRAMPA_SKIP_MULTREES=1
grampa -g $GENE_TREES -s $SPECIES_TREE --st-only -h1 "" -h2 "" --maps -o $OUTPUT_DIR --overwrite

# Check if GRAMPA completed successfully
if [ $? -eq 0 ]; then
    echo "GRAMPA completed successfully for $DATASET"
else
    echo "ERROR: GRAMPA failed for $DATASET"
    exit 1
fi

# Deactivate conda environment
conda deactivate

echo "Job completed for $DATASET"