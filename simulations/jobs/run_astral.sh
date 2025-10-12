#!/bin/bash
#SBATCH --job-name=astral4_batch
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/astral4_ILS_low_dup_low_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/astral4_ILS_low_dup_low_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Activate your specific environment
conda activate gene2net || {
    echo "ERROR: Could not activate polyphest environment"
    exit 1
}

# Path to ASTRAL-4 executable (now from conda)
ASTRAL_PATH="astral4"  # Will use the one from conda environment

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# Get the network for this array task
network=${networks[$SLURM_ARRAY_TASK_ID-1]}
echo "Running ASTRAL-4 for: ${network} - ILS_low_dup_low"

# Define paths
POLYPHEST_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/polyphest_input"
GRAMPA_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/grampa_input"

# Check if input directory exists
if [ ! -d "$POLYPHEST_DIR" ]; then
    echo "ERROR: Polyphest input directory not found: $POLYPHEST_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$GRAMPA_DIR"

# Count how many replicate tree files exist
replicate_count=0

# Process each replicate tree file
for input_file in "$POLYPHEST_DIR"/*_trees.tre; do
    # Check if file exists (in case no matches found)
    if [ ! -f "$input_file" ]; then
        echo "WARNING: No tree files found in $POLYPHEST_DIR"
        break
    fi
    
    # Extract replicate name (e.g., "1" from "1_trees.tre")
    filename=$(basename "$input_file")
    replicate="${filename%_trees.tre}"
    
    replicate_count=$((replicate_count + 1))
    
    echo ""
    echo "Processing replicate: $replicate"
    
    # Define output files
    output_file="${GRAMPA_DIR}/${replicate}_species.tre"
    log_file="${GRAMPA_DIR}/${replicate}_astral4.log"
    
    # Run ASTRAL-4
    echo "  Input:  $input_file"
    echo "  Output: $output_file"
    echo "  Running ASTRAL-4..."
    
    $ASTRAL_PATH -o "$output_file" "$input_file" 2>"$log_file"
    
    # Check if the run was successful
    if [ $? -eq 0 ]; then
        echo "  ? ASTRAL-4 completed successfully for replicate $replicate"
    else
        echo "  ? ERROR: ASTRAL-4 failed for replicate $replicate"
        echo "  Check log file: $log_file"
    fi
done

if [ $replicate_count -eq 0 ]; then
    echo "ERROR: No replicate tree files found in $POLYPHEST_DIR"
    exit 1
fi

echo ""
echo "=========================================="
echo "ASTRAL-4 processing completed for ${network}"
echo "Processed $replicate_count replicates"
echo "Species trees saved to: $GRAMPA_DIR"
echo "=========================================="