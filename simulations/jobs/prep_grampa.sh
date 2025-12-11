#!/bin/bash
#SBATCH --job-name=prep_grampa
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/prep_grampa_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/prep_grampa_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=4g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# PREP_GRAMPA.SH
# ============================================================================
# Prepares SimPhy output for GRAMPA by:
#   1. Concatenating gene trees from all batches (if applicable) into temp dir
#   2. Running Python script to:
#      - Clean taxa names (remove suffixes)
#      - Fix substring issues
#      - Output clean_trees.tre (for ASTRAL)
#      - Output grampa_trees.tre (with copy numbers for GRAMPA)
#      - Output taxa_map.txt (if any fixes needed)
#
# Usage:
#   sbatch --export=CONFIG=conf_ils_low_10M prep_grampa.sh
#   sbatch --export=CONFIG=conf_ils_low_10M,NUM_REPLICATES=3 prep_grampa.sh
#
# Environment variables:
#   CONFIG          - Configuration name (required, e.g., conf_ils_low_10M)
#   NUM_REPLICATES  - Number of replicates to process (default: 5)
# ============================================================================

set -eo pipefail

# Initialize LD_LIBRARY_PATH if not set (needed for conda activation)
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ============================================================================
# CONFIGURATION
# ============================================================================

# Required: configuration name
CONFIG="${CONFIG:?ERROR: CONFIG environment variable is required. Use --export=CONFIG=conf_name}"

# Optional: number of replicates (default: 5)
NUM_REPLICATES="${NUM_REPLICATES:-5}"

# ============================================================================
# PATHS
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"
PYTHON_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/process_gene_trees_for_grampa.py"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# ============================================================================
# SETUP
# ============================================================================

# Get network for this array task
network_idx=$((SLURM_ARRAY_TASK_ID - 1))
network="${networks[$network_idx]}"

# Define paths
DATA_DIR="${BASE_DIR}/${network}/data/${CONFIG}"
OUTPUT_BASE="${BASE_DIR}/${network}/processed/${CONFIG}/grampa_input"

echo "============================================================================"
echo "PREP GRAMPA - ${network}"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Replicates: ${NUM_REPLICATES}"
echo "Data directory: ${DATA_DIR}"
echo "Output directory: ${OUTPUT_BASE}"
echo "Date: $(date)"
echo "============================================================================"
echo ""

# ============================================================================
# VALIDATION
# ============================================================================

if [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: Data directory not found: $DATA_DIR"
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# ============================================================================
# ACTIVATE CONDA
# ============================================================================

source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

echo "? Conda environment activated: gene2net"
echo ""

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

collect_gene_trees() {
    # Collect all gene trees into a temporary directory
    # Args: replicate_data_dir, temp_dir
    local replicate_data_dir=$1
    local temp_dir=$2
    
    local gene_tree_count=0
    
    # Check if this is a batched structure or single batch
    if [ -d "${replicate_data_dir}/batch_1" ]; then
        # Multiple batches structure: replicate_X/batch_Y/1/g_*
        echo "    Detected batched structure"
        
        for batch_dir in "$replicate_data_dir"/batch_*/; do
            if [ ! -d "$batch_dir" ]; then
                continue
            fi
            
            local inner_dir="${batch_dir}1"
            if [ ! -d "$inner_dir" ]; then
                echo "    WARNING: No '1' subdirectory in $(basename $batch_dir)"
                continue
            fi
            
            for gene_tree in "$inner_dir"/g_*; do
                if [ ! -f "$gene_tree" ]; then
                    continue
                fi
                
                # Copy to temp dir with unique name
                gene_tree_count=$((gene_tree_count + 1))
                cp "$gene_tree" "${temp_dir}/g_tree_${gene_tree_count}"
            done
        done
    else
        # Single batch structure: replicate_X/1/g_*
        echo "    Detected single batch structure"
        
        local inner_dir="${replicate_data_dir}/1"
        if [ ! -d "$inner_dir" ]; then
            echo "    ERROR: No '1' subdirectory found in ${replicate_data_dir}"
            return 1
        fi
        
        for gene_tree in "$inner_dir"/g_*; do
            if [ ! -f "$gene_tree" ]; then
                continue
            fi
            
            # Copy to temp dir with unique name
            gene_tree_count=$((gene_tree_count + 1))
            cp "$gene_tree" "${temp_dir}/g_tree_${gene_tree_count}"
        done
    fi
    
    if [ $gene_tree_count -eq 0 ]; then
        echo "    ERROR: No gene trees found"
        return 1
    fi
    
    echo "    Collected $gene_tree_count gene trees"
    return 0
}

# ============================================================================
# MAIN PROCESSING
# ============================================================================

success_count=0
error_count=0

for replicate in $(seq 1 $NUM_REPLICATES); do
    echo "----------------------------------------------------------------------------"
    echo "REPLICATE ${replicate}/${NUM_REPLICATES}"
    echo "----------------------------------------------------------------------------"
    
    # Find the replicate data directory
    replicate_data_dir="${DATA_DIR}/replicate_${replicate}"
    
    if [ ! -d "$replicate_data_dir" ]; then
        echo "  ERROR: Replicate directory not found: $replicate_data_dir"
        error_count=$((error_count + 1))
        continue
    fi
    
    # Create output directory for this replicate
    output_dir="${OUTPUT_BASE}/replicate_${replicate}"
    mkdir -p "$output_dir"
    
    # Create temporary directory for collecting gene trees
    temp_dir=$(mktemp -d)
    
    echo "  Input: $replicate_data_dir"
    echo "  Output: $output_dir"
    echo "  Temp: $temp_dir"
    echo ""
    
    # Step 1: Collect gene trees into temp directory
    echo "  Step 1: Collecting gene trees..."
    if ! collect_gene_trees "$replicate_data_dir" "$temp_dir"; then
        echo "  ? Failed to collect gene trees"
        rm -rf "$temp_dir"
        error_count=$((error_count + 1))
        continue
    fi
    echo "  ? Gene trees collected"
    echo ""
    
    # Step 2: Run Python script to process trees
    echo "  Step 2: Processing trees (clean, fix substrings, add copy numbers)..."
    python "$PYTHON_SCRIPT" "$temp_dir" "$output_dir" -v
    
    if [ $? -ne 0 ]; then
        echo "  ? Failed to process gene trees"
        rm -rf "$temp_dir"
        error_count=$((error_count + 1))
        continue
    fi
    echo "  ? Trees processed"
    echo ""
    
    # Cleanup temp directory
    rm -rf "$temp_dir"
    
    echo "  ? Replicate ${replicate} completed successfully"
    success_count=$((success_count + 1))
done

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "============================================================================"
echo "SUMMARY - ${network}"
echo "============================================================================"
echo "Configuration: ${CONFIG}"
echo "Successful replicates: ${success_count}/${NUM_REPLICATES}"
echo "Failed replicates: ${error_count}/${NUM_REPLICATES}"
echo "Output directory: ${OUTPUT_BASE}"
echo ""
echo "Output files per replicate:"
echo "  - clean_trees.tre   (for ASTRAL)"
echo "  - grampa_trees.tre  (for GRAMPA)"
echo "  - taxa_map.txt      (if substring fixes needed)"
echo ""

if [ $error_count -gt 0 ]; then
    echo "Status: ? COMPLETED WITH ERRORS"
    exit 1
else
    echo "Status: ? SUCCESS"
fi

echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"