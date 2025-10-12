#!/bin/bash
#SBATCH --job-name=mpsugar_prep
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/mpsugar_ILS_low_dup_low_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/mpsugar_ILS_low_dup_low_job_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=2g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

# Set conda path
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"

# Activate conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Activate your specific environment
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

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

echo "Processing gene trees for MP-SUGAR: ${network} - ILS_low_dup_low"
echo "========================================================"

# Define paths
DATA_DIR="${BASE_DIR}/${network}/data/ILS_low_dup_low"
OUTPUT_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/mpallopp_input"

# Path to the Python processing script
SCRIPT_PATH="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/process_for_mpsugar.py"

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: Data directory not found: $DATA_DIR"
    exit 1
fi

# Check if Python script exists
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "ERROR: Python script not found: $SCRIPT_PATH"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each replicate directory
total_replicates=0
successful_replicates=0

for replicate_dir in "$DATA_DIR"/*/; do
    # Check if directory exists and is not empty
    if [ ! -d "$replicate_dir" ]; then
        continue
    fi
    
    # Get replicate number from directory name
    replicate=$(basename "$replicate_dir")
    
    echo "Processing replicate: $replicate"
    
    # Output files for this replicate
    output_file="${OUTPUT_DIR}/${replicate}_MPSUGAR_trees.nex"
    map_file="${OUTPUT_DIR}/${replicate}_taxon_map.json"
    
    # Run the Python script
    python "$SCRIPT_PATH" -i "$replicate_dir" -o "$output_file" -m "$map_file"
    
    if [ $? -eq 0 ]; then
        echo "? Successfully processed replicate $replicate"
        successful_replicates=$((successful_replicates + 1))
    else
        echo "? Failed to process replicate $replicate"
    fi
    
    total_replicates=$((total_replicates + 1))
    echo "----------------------------------------"
done

echo "========================================================"
echo "Processing completed for ${network}"
echo "Total replicates: $total_replicates"
echo "Successful: $successful_replicates"
echo "Failed: $((total_replicates - successful_replicates))"
echo "Output directory: $OUTPUT_DIR"