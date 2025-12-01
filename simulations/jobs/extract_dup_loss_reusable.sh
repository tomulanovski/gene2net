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
SIMPHY_CONFIG="${EXTRACT_CONFIG:-conf_ils_low_10M}"
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

# Scripts
EXTRACT_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/extract_sql.py"
SUMMARIZE_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/summarize_dup_loss.py"
AGGREGATE_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/aggregate_dup_loss.py"

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
DATA_SUBDIR="data/${SIMPHY_CONFIG}"
OUTPUT_FILE="dup_loss_summary_${SIMPHY_CONFIG}.txt"

networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

echo "=========================================="
echo "Starting Extraction and Summarization"
echo "=========================================="
echo "Configuration: ${SIMPHY_CONFIG}"
echo ""

for network in "${networks[@]}"; do
    echo "Processing: ${network}"
    echo "------------------------------------------"
    
    for rep in $(seq 1 $NUM_REPLICATES); do
        REP_DIR="${BASE_DIR}/${network}/${DATA_SUBDIR}/replicate_${rep}"
        
        # Check if Batch Mode (looks for batch_1 folder)
        if [ -d "${REP_DIR}/batch_1" ]; then
            echo "  Replicate ${rep}: Detected Batch Mode (Merging batches...)"
            
            # 1. Reset Master Files
            # We will merge everything into these files inside replicate_X/
            > "${REP_DIR}/g_trees.trees"
            > "${REP_DIR}/l_trees.trees"
            
            extracted_count=0
            
            # 2. Loop through all batch folders
            for batch_dir in ${REP_DIR}/batch_*; do
                [ -d "$batch_dir" ] || continue
                
                # Find DB
                DB_PATH=$(find "${batch_dir}" -maxdepth 2 -name "*.db" | head -n 1)
                
                if [ -n "$DB_PATH" ]; then
                    # Extract this specific batch
                    python3 "$EXTRACT_SCRIPT" "$DB_PATH" > /dev/null 2>&1
                    
                    DB_DIR=$(dirname "$DB_PATH")
                    
                    # A. Handle Species Tree (Only take from the FIRST extracted batch)
                    if [ $extracted_count -eq 0 ] && [ -f "${DB_DIR}/s_tree.trees" ]; then
                        cp "${DB_DIR}/s_tree.trees" "${REP_DIR}/s_tree.trees"
                    fi

                    # B. Append Gene Trees (Take from ALL batches)
                    if [ -f "${DB_DIR}/g_trees.trees" ]; then
                        cat "${DB_DIR}/g_trees.trees" >> "${REP_DIR}/g_trees.trees"
                    fi
                    
                    # C. Append Locus Trees (Take from ALL batches)
                    if [ -f "${DB_DIR}/l_trees.trees" ]; then
                        cat "${DB_DIR}/l_trees.trees" >> "${REP_DIR}/l_trees.trees"
                    fi
                    
                    ((extracted_count++))
                fi
            done
            echo "    Merged data from ${extracted_count} batches."
            
            # 3. VERIFY MERGE AND RUN SUMMARY
            if [ -s "${REP_DIR}/g_trees.trees" ]; then
                echo "    Running Summary on Merged Files..."
                # Run summary on the main replicate folder where we just put the big files
                python3 "$SUMMARIZE_SCRIPT" "$REP_DIR"
            else
                echo "    ERROR: Merged gene tree file is empty! Check extraction."
            fi

        else
            # Single Mode (Standard)
            DB_PATH=$(find "${REP_DIR}" -maxdepth 2 -name "*.db" | head -n 1)
            
            if [ -n "$DB_PATH" ]; then
                echo "  Replicate ${rep}: Single Mode. DB: $(basename "$DB_PATH")"
                python3 "$EXTRACT_SCRIPT" "$DB_PATH"
                
                ACTUAL_DB_DIR=$(dirname "$DB_PATH")
                python3 "$SUMMARIZE_SCRIPT" "$ACTUAL_DB_DIR"
                
                # Move summary up if needed
                if [ "$ACTUAL_DB_DIR" != "$REP_DIR" ] && [ -f "${ACTUAL_DB_DIR}/dup_loss_summary.txt" ]; then
                     cp "${ACTUAL_DB_DIR}/dup_loss_summary.txt" "${REP_DIR}/"
                fi
            else
                echo "    WARNING: No .db file found in Replicate ${rep}"
            fi
        fi
    done
    echo ""
done

echo "=========================================="
echo "Starting Aggregation"
echo "=========================================="
python3 "$AGGREGATE_SCRIPT" "$DATA_SUBDIR" "$OUTPUT_FILE"

if [ $? -eq 0 ]; then
    echo "SUCCESS: Aggregated results saved to ${OUTPUT_FILE}"
else
    echo "ERROR: Aggregation failed"
    exit 1
fi