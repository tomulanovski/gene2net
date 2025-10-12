#!/bin/bash
#SBATCH --job-name=fix_species_trees_for_grampa
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/fix_species_ILS_low_dup_low_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/fix_species_ILS_low_dup_low_%a.err
#SBATCH --time=0:20:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
SCRIPT_PATH="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/fix_substrings_to_astral.py"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# Get the network for this array task
network=${networks[$SLURM_ARRAY_TASK_ID-1]}

echo "Fixing species trees with taxa maps: ${network} - ILS_low_dup_low"

# Define paths
GRAMPA_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/grampa_input"

# Check if directory exists
if [ ! -d "$GRAMPA_DIR" ]; then
    echo "ERROR: GRAMPA directory not found: $GRAMPA_DIR"
    exit 1
fi

# Counters
total_replicates=0
fixed_replicates=0
skipped_replicates=0
failed_replicates=0

# Process each replicate
for species_tree in "$GRAMPA_DIR"/*_species.tre; do
    # Check if file exists (in case no matches found)
    if [ ! -f "$species_tree" ]; then
        echo "WARNING: No species tree files found in $GRAMPA_DIR"
        break
    fi
    
    # Extract replicate name (e.g., "1" from "1_species.tre")
    filename=$(basename "$species_tree")
    replicate="${filename%_species.tre}"
    
    total_replicates=$((total_replicates + 1))
    
    # Check if taxa map exists for this replicate
    taxa_map="${GRAMPA_DIR}/${replicate}_taxa_map.txt"
    
    if [ ! -f "$taxa_map" ]; then
        echo "Replicate $replicate: No taxa map found - skipping (no changes needed)"
        skipped_replicates=$((skipped_replicates + 1))
        continue
    fi
    
    echo ""
    echo "Processing replicate: $replicate"
    echo "  Species tree: $species_tree"
    echo "  Taxa map: $taxa_map"
    
    # Apply the taxa map
    python "$SCRIPT_PATH" "$species_tree" "$taxa_map" -v
    
    # Check if successful
    if [ $? -eq 0 ]; then
        echo "  ? Successfully fixed species tree for replicate $replicate"
        fixed_replicates=$((fixed_replicates + 1))
    else
        echo "  ? ERROR: Failed to fix species tree for replicate $replicate"
        failed_replicates=$((failed_replicates + 1))
    fi
done

# Summary
echo ""
echo "=========================================="
echo "Species tree fixing completed for ${network}"
echo "Total replicates: $total_replicates"
echo "  Fixed: $fixed_replicates"
echo "  Skipped (no map): $skipped_replicates"
echo "  Failed: $failed_replicates"
echo "=========================================="

# Exit with error if any failures
if [ $failed_replicates -gt 0 ]; then
    exit 1
fi