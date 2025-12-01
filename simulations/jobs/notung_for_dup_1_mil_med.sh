#!/bin/bash
#SBATCH --job-name=notung_simulated
#SBATCH --array=0-20
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/notung_simulated_1_mil_med_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/notung_simulated_1_mil_med_%a.err
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# Get current network based on array task ID
NETWORK=${networks[$SLURM_ARRAY_TASK_ID]}

echo "Processing simulated dataset: $NETWORK"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "========================================"

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Set paths for current dataset
GENE_TREES_DIR="${BASE_DIR}/${NETWORK}/data/med_dup_ne_100000_height_1_million_trees/1"
GENE_TREES_NHX_DIR="${BASE_DIR}/${NETWORK}/data/med_dup_ne_100000_height_1_million_trees/gene_trees_nhx"
SPECIES_TREE_NEWICK="${BASE_DIR}/${NETWORK}/species_tree_1_million.tre"
OUTPUT_DIR="${BASE_DIR}/${NETWORK}/data/med_dup_ne_100000_height_1_million_trees/reconciled_for_dup"
BATCH_FILE="${OUTPUT_DIR}/notung_batch.txt"

# Script paths
CONVERT_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/add_nhx_tags.py"

# Notung jar path
NOTUNG_JAR="/groups/itay_mayrose/tomulanovski/Notung-2.9.1.5/Notung-2.9.1.5.jar"

# Conda setup
CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# Check if input files exist
if [ ! -d "$GENE_TREES_DIR" ]; then
    echo "ERROR: Gene trees directory not found: $GENE_TREES_DIR"
    exit 1
fi

if [ ! -f "$SPECIES_TREE_NEWICK" ]; then
    echo "ERROR: Species tree Newick file not found: $SPECIES_TREE_NEWICK"
    exit 1
fi

# Create output directories
mkdir -p "$GENE_TREES_NHX_DIR"
mkdir -p "$OUTPUT_DIR"

echo ""
echo "Gene trees directory: $GENE_TREES_DIR"
echo "Species tree Newick: $SPECIES_TREE_NEWICK"
echo "Output directory: $OUTPUT_DIR"

# Step 1: (Removed conversion from NEXUS to Newick, using existing Newick directly)

# Step 2: Convert gene trees to NHX format
echo ""
echo "Step 2: Converting gene trees to NHX format..."

tree_count=0
for gene_tree in ${GENE_TREES_DIR}/g_trees*.trees; do
    if [ -f "$gene_tree" ]; then
        filename=$(basename "$gene_tree")
        output_nhx="${GENE_TREES_NHX_DIR}/${filename}"
        
        python "$CONVERT_SCRIPT" "$SPECIES_TREE_NEWICK" "$gene_tree" "$output_nhx"
        
        if [ $? -ne 0 ]; then
            echo "ERROR: Failed to convert $gene_tree to NHX format"
            exit 1
        fi
        
        ((tree_count++))
        
        # Progress update every 100 trees
        if [ $((tree_count % 100)) -eq 0 ]; then
            echo "Converted $tree_count trees..."
        fi
    fi
done

echo "Successfully converted $tree_count gene trees to NHX format"

# Step 3: Create Notung batch file
echo ""
echo "Step 3: Creating Notung batch file..."

# Write species tree contents as first line of batch file
echo "$(<"$SPECIES_TREE_NEWICK")" > "$BATCH_FILE"

# Add all NHX gene tree files
for gene_tree in ${GENE_TREES_NHX_DIR}/g_trees*.trees; do
    if [ -f "$gene_tree" ]; then
        echo "$gene_tree" >> "$BATCH_FILE"
    fi
done

# Verify batch file was created
if [ ! -f "$BATCH_FILE" ]; then
    echo "ERROR: Batch file not created: $BATCH_FILE"
    exit 1
fi

echo "Batch file created at: $BATCH_FILE"
echo "Number of gene trees: $(wc -l < "$BATCH_FILE" | awk '{print $1-1}')"

# Step 4: Run Notung reconciliation
echo ""
echo "Step 4: Running Notung reconciliation..."

# Run in headless mode (no GUI)
export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"

java -Xmx25G -jar ${NOTUNG_JAR} \
    -b ${BATCH_FILE} \
    --absfilenames \
    --reconcile \
    --treestats \
    --speciestag nhx \
    --outputdir ${OUTPUT_DIR} \
    --log
    
if [ $? -eq 0 ]; then
    echo ""
    echo "SUCCESS: Notung completed successfully for $NETWORK"
    echo "Output directory: $OUTPUT_DIR"
else
    echo ""
    echo "ERROR: Notung failed for $NETWORK"
    exit 1
fi

# Deactivate conda
conda deactivate

echo ""
echo "========================================"
echo "Job completed for $NETWORK"
echo "========================================"
