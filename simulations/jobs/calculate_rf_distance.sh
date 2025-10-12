#!/bin/bash
#SBATCH --job-name=rf_all_networks
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/rf_all_networks.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/rf_all_networks.err
#SBATCH --time=8:00:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Activate your specific environment
conda activate polyphest || {
    echo "ERROR: Could not activate polyphest environment"
    exit 1
}

# Paths to scripts
RF_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/rf_ad.py"
AGGREGATE_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/aggregate_ad.py"

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

echo "=========================================="
echo "Starting RF Distance Calculation"
echo "=========================================="
echo "Total networks: ${#networks[@]}"
echo ""

# Process each network sequentially
for network in "${networks[@]}"; do
    echo "Processing: ${network}"
    echo "------------------------------------------"
    
    TREE_DIR="${BASE_DIR}/${network}/data/ILS_low_dup_low/1"
    OUTPUT_FILE="${BASE_DIR}/${network}/data/ILS_low_dup_low/rf_distance_results.txt"
    
    # Check if tree directory exists
    if [ ! -d "$TREE_DIR" ]; then
        echo "WARNING: Tree directory not found: $TREE_DIR"
        echo ""
        continue
    fi
    
    # Run the RF distance calculation
    python3 "$RF_SCRIPT" "$TREE_DIR" "$OUTPUT_FILE"
    
    if [ $? -eq 0 ]; then
        echo "? Successfully calculated RF distances for ${network}"
    else
        echo "? ERROR calculating RF distances for ${network}"
    fi
    echo ""
done

echo "=========================================="
echo "All RF Calculations Complete"
echo "=========================================="
echo ""
echo "Starting Aggregation..."
echo ""

# Run aggregation script
python3 "$AGGREGATE_SCRIPT"

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "All Processing Complete!"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "ERROR: Aggregation failed!"
    echo "=========================================="
    exit 1
fi