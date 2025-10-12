#!/bin/bash
#SBATCH --job-name=process_simphy
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/process_ILS_low_dup_low_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/process_ILS_low_dup_low_job_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=2g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

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
echo "Processing gene trees for: ${network} - ILS_low_dup_low"

# Define paths. Need to change based on the combination. Currently makes gene trees for ILS_low_dup_low
DATA_DIR="${BASE_DIR}/${network}/data/ILS_low_dup_low"
OUTPUT_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/polyphest_input"

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: Data directory not found: $DATA_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each replicate directory
for replicate_dir in "$DATA_DIR"/*/; do
    # Check if directory exists and is not empty
    if [ ! -d "$replicate_dir" ]; then
        continue
    fi
    
    # Get replicate number from directory name
    replicate=$(basename "$replicate_dir")
    
    echo "Processing replicate: $replicate"
    
    # Output file for this replicate
    output_file="${OUTPUT_DIR}/${replicate}_trees.tre"
    
    # Remove output file if it exists
    rm -f "$output_file"
    
    # Find all gene tree files (starting with 'g')
    gene_tree_count=0
    for gene_tree in "$replicate_dir"/g_*; do
        # Check if file exists
        if [ ! -f "$gene_tree" ]; then
            continue
        fi
        
        # Process the gene tree: remove everything from first underscore onwards in taxa names
        # Using sed to replace pattern like "taxon_anything" with "taxon"
        sed 's/_[^,):]*\([,):]\)/\1/g' "$gene_tree" >> "$output_file"
        
        gene_tree_count=$((gene_tree_count + 1))
    done
    
    if [ $gene_tree_count -eq 0 ]; then
        echo "WARNING: No gene trees found in $replicate_dir"
        rm -f "$output_file"
    else
        echo "Processed $gene_tree_count gene trees for replicate $replicate"
        echo "Output saved to: $output_file"
    fi
done

echo "Processing completed for ${network}"