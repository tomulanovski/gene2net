#!/bin/bash
#SBATCH --job-name=process_grampa
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/process_ILS_low_dup_low_grampa_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/process_ILS_low_dup_low_grampa_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=2g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
SCRIPT_PATH="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/process_gene_trees_for_grampa.py"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# Get the network for this array task
network=${networks[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing gene trees for GRAMPA: ${network} - ILS_low_dup_low"

# Define paths
DATA_DIR="${BASE_DIR}/${network}/data/ILS_low_dup_low"
OUTPUT_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/grampa_input"

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: Data directory not found: $DATA_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each replicate directory
for replicate_dir in "$DATA_DIR"/*/; do
    if [ ! -d "$replicate_dir" ]; then
        continue
    fi
    
    replicate=$(basename "$replicate_dir")
    echo "Processing replicate: $replicate"
    
    output_file="${OUTPUT_DIR}/${replicate}_trees.tre"
    taxa_map_file="${OUTPUT_DIR}/${replicate}_taxa_map.txt"
    
    # Run the Python script to process all trees in this replicate
    # Now with taxa map output
    python "${SCRIPT_PATH}" "$replicate_dir" "$output_file" --taxa-map "$taxa_map_file" -v
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to process replicate $replicate"
    else
        echo "  - Trees: $output_file"
        echo "  - Taxa map: $taxa_map_file"
    fi
done

echo "Processing completed for ${network}"