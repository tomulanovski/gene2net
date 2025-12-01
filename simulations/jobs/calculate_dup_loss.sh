#!/bin/bash
#SBATCH --job-name=extract_summarize_all
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/extract_dup_loss_100000_ne_50_million_height_all.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/extract_dup_loss_100000_ne_50_million_height_all.err
#SBATCH --time=8:00:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Activate your specific environment
conda activate polyphest || {
    echo "ERROR: Could not activate polyphest environment"
    exit 1
}

# Paths to scripts
EXTRACT_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/extract_sql.py"
SUMMARIZE_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/summarize_dup_loss.py"
AGGREGATE_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/aggregate_dup_loss.py"

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Data subdirectory,Database filename,final output file name (to be changed for different runs)
DATA_SUBDIR="data/high_dup_ne_100000_height_50_million_trees"
DB_FILENAME="high_dup_ne_100000_height_50_million_trees.db"
OUTPUT_FILE="dup_loss_summary_50_mil_high_dup.txt"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

echo "=========================================="
echo "Starting SQL Extraction and Summarization"
echo "=========================================="
echo "Total networks: ${#networks[@]}"
echo ""

# Step 1: Extract SQL databases to CSV and create summaries
for network in "${networks[@]}"; do
    echo "Processing: ${network}"
    echo "------------------------------------------"
    
    DATA_DIR="${BASE_DIR}/${network}/${DATA_SUBDIR}"
    DB_PATH="${DATA_DIR}/${DB_FILENAME}"
    
    # Check if database exists
    if [ ! -f "$DB_PATH" ]; then
        echo "WARNING: Database not found: $DB_PATH"
        echo ""
        continue
    fi
    
    # Extract SQL to CSV
    echo "Extracting database to CSV..."
    python3 "$EXTRACT_SCRIPT" "$DB_PATH"
    
    if [ $? -ne 0 ]; then
        echo "? ERROR extracting database for ${network}"
        echo ""
        continue
    fi
    
    echo "? Successfully extracted database"
    
    # Summarize duplications and losses
    echo "Summarizing duplications and losses..."
    python3 "$SUMMARIZE_SCRIPT" "$DATA_DIR"
    
    if [ $? -eq 0 ]; then
        echo "? Successfully created summary for ${network}"
    else
        echo "? ERROR creating summary for ${network}"
    fi
    
    echo ""
done

echo "=========================================="
echo "All Extractions and Summaries Complete"
echo "=========================================="
echo ""

# Step 2: Aggregate all network summaries
echo "Starting Aggregation Across All Networks..."
echo ""

python3 "$AGGREGATE_SCRIPT" "$DATA_SUBDIR" "$OUTPUT_FILE"

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