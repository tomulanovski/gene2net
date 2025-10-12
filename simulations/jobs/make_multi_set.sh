#!/bin/bash
#SBATCH --job-name=copies_smooth
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/copies_smoothing_ILS_low_dup_low_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/copies_smoothing_ILS_low_dup_low_job_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=4g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Activate your specific environment
conda activate gene2net || {
    echo "ERROR: Could not activate polyphest environment"
    exit 1
}

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Path to the Python script
PYTHON_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/copies_smoothing_with_multiset.py"

# Array of networks (must match the original array)
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# Get the network for this array task
network=${networks[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing copy number smoothing for: ${network} - ILS_low_dup_low"

# Define paths - output to same directory as input. Currently just ILS_low_dup_low
INPUT_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/polyphest_input"

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Counter for processed files
processed_count=0
error_count=0

# Process each tree file in the input directory
for tree_file in "$INPUT_DIR"/*_trees.tre; do
    # Check if file exists
    if [ ! -f "$tree_file" ]; then
        echo "WARNING: No tree files found in $INPUT_DIR"
        continue
    fi
    
    # Extract replicate number from filename
    # Example: 1_trees.tre -> 1
    filename=$(basename "$tree_file")
    replicate="${filename%_trees.tre}"
    
    echo "Processing replicate: $replicate"
    
    # Define output file path
    output_multiset="${INPUT_DIR}/${replicate}_multi_set.txt"
    
    # Run the Python script
    python "$PYTHON_SCRIPT" \
        -i "$tree_file" \
        -m "$output_multiset" \
        -o "${INPUT_DIR}/${replicate}_distribution.tsv" \
        -e "full" \
        -k 2
    
    # Check if the script succeeded
    if [ $? -eq 0 ]; then
        echo "Successfully created multi_set for replicate $replicate"
        processed_count=$((processed_count + 1))
    else
        echo "ERROR: Failed to process replicate $replicate"
        error_count=$((error_count + 1))
    fi
    
    echo "---"
done

echo "=========================================="
echo "Processing completed for ${network}"
echo "Successfully processed: $processed_count replicates"
echo "Errors encountered: $error_count replicates"
echo "Output directory: $OUTPUT_DIR"
echo "=========================================="