#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# CONFIGURATION
# ============================================================================

# Use the exported variable from the submitter, or default to a fallback
SIMPHY_CONFIG="${RF_CONFIG:-conf_ils_low_10M}"
NUM_REPLICATES=5

# ============================================================================
# SETUP
# ============================================================================

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate polyphest || {
    echo "ERROR: Could not activate polyphest environment"
    exit 1
}

RF_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/rf_ad.py"
AGGREGATE_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/aggregate_ad_reusable.py"
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

echo "=========================================="
echo "Starting RF Distance Calculation"
echo "=========================================="
echo "Configuration: ${SIMPHY_CONFIG}"
echo "Date: $(date)"
echo ""

# Process each network
for network in "${networks[@]}"; do
    echo "Processing: ${network}"
    echo "------------------------------------------"
    
    for rep in $(seq 1 $NUM_REPLICATES); do
        echo "  Replicate ${rep}/${NUM_REPLICATES}"
        
        REP_BASE="${BASE_DIR}/${network}/data/${SIMPHY_CONFIG}/replicate_${rep}"
        OUTPUT_FILE="${REP_BASE}/rf_distance_results.txt"
        
        # Check structure type
        if [ -d "${REP_BASE}/1" ]; then
            # Single batch structure
            echo "    Mode: Single batch (1000 trees)"
            TREE_DIR="${REP_BASE}/1"
            
            python3 "$RF_SCRIPT" "$TREE_DIR" "$OUTPUT_FILE" > /dev/null 2>&1
            
        elif [ -d "${REP_BASE}/batch_1" ]; then
            # Multi-batch structure - need to merge
            echo "    Mode: Multi-batch (10ª100)"
            
            # Create temporary merge directory
            MERGE_DIR="${REP_BASE}/merged_trees"
            mkdir -p "$MERGE_DIR"
            
            # Copy species tree from first batch
            cp "${REP_BASE}/batch_1/1/s_tree.trees" "$MERGE_DIR/"
            
            # Merge all gene trees with unique names
            batch_count=0
            for batch_dir in ${REP_BASE}/batch_*/1; do
                if [ -d "$batch_dir" ]; then
                    batch_num=$(basename $(dirname "$batch_dir") | sed 's/batch_//')
                    
                    # Copy gene trees with batch prefix
                    for tree_file in ${batch_dir}/g_trees*.trees; do
                        if [ -f "$tree_file" ]; then
                            basename=$(basename "$tree_file")
                            cp "$tree_file" "${MERGE_DIR}/batch${batch_num}_${basename}"
                        fi
                    done
                    ((batch_count++))
                fi
            done
            
            echo "    Merged trees from ${batch_count} batches"
            
            # Run RF analysis on merged directory
            python3 "$RF_SCRIPT" "$MERGE_DIR" "$OUTPUT_FILE" > /dev/null 2>&1
            
            # Clean up
            rm -rf "$MERGE_DIR"
            
        else
            echo "    WARNING: No tree directory found"
            continue
        fi
        
        if [ $? -eq 0 ]; then
            if [ -f "$OUTPUT_FILE" ]; then
                AD=$(grep "Average Distance (AD):" "$OUTPUT_FILE" | awk '{print $5}')
                echo "    ? AD = ${AD}"
            else
                 echo "    ? ERROR: Output file not created"
            fi
        else
            echo "    ? ERROR calculating RF distances"
        fi
    done
    echo ""
done

echo "=========================================="
echo "All RF Calculations Complete"
echo "=========================================="
echo ""
echo "Starting Aggregation..."
echo ""

# Run aggregation using the specific config
python3 "$AGGREGATE_SCRIPT" "${SIMPHY_CONFIG}"

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