#!/bin/bash
#SBATCH --job-name=notung_reconcile_only
#SBATCH --array=0-20
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/notung_reconcile_only_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/notung_reconcile_only_%a.err
#SBATCH --time=50:00:00
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
GENE_TREES_NHX_DIR="${BASE_DIR}/${NETWORK}/data/low_dup_ne_100000_height_1_million_trees/gene_trees_nhx"
SPECIES_TREE_NEWICK="${BASE_DIR}/${NETWORK}/species_tree_1_million.tre"
OUTPUT_DIR="${BASE_DIR}/${NETWORK}/data/low_dup_ne_100000_height_1_million_trees/reconciled_for_dup"
BATCH_FILE="${OUTPUT_DIR}/notung_batch.txt"

# Notung jar path
NOTUNG_JAR="/groups/itay_mayrose/tomulanovski/Notung-2.9.1.5/Notung-2.9.1.5.jar"

# Check if input files exist
if [ ! -d "$GENE_TREES_NHX_DIR" ]; then
    echo "ERROR: NHX gene trees directory not found: $GENE_TREES_NHX_DIR"
    exit 1
fi

if [ ! -f "$SPECIES_TREE_NEWICK" ]; then
    echo "ERROR: Species tree not found: $SPECIES_TREE_NEWICK"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo ""
echo "NHX gene trees directory: $GENE_TREES_NHX_DIR"
echo "Species tree: $SPECIES_TREE_NEWICK"
echo "Output directory: $OUTPUT_DIR"

# Step 1: Create Notung batch file
echo ""
echo "Step 1: Creating Notung batch file..."

# Write species tree as first line
echo "$SPECIES_TREE_NEWICK" > "$BATCH_FILE"

# Add all NHX gene tree files
gene_count=0
for gene_tree in ${GENE_TREES_NHX_DIR}/g_trees*.trees; do
    if [ -f "$gene_tree" ]; then
        echo "$gene_tree" >> "$BATCH_FILE"
        ((gene_count++))
    fi
done

# Verify batch file was created
if [ ! -f "$BATCH_FILE" ]; then
    echo "ERROR: Batch file not created: $BATCH_FILE"
    exit 1
fi

echo "Batch file created at: $BATCH_FILE"
echo "Number of gene trees: $gene_count"

if [ $gene_count -eq 0 ]; then
    echo "ERROR: No gene trees found in $GENE_TREES_NHX_DIR"
    exit 1
fi

# Step 2: Run Notung reconciliation
echo ""
echo "Step 2: Running Notung reconciliation..."

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
    
    # Count output files
    reconciled_count=$(ls ${OUTPUT_DIR}/*.reconciled 2>/dev/null | wc -l)
    echo "Reconciled trees created: $reconciled_count"
else
    echo ""
    echo "ERROR: Notung failed for $NETWORK"
    exit 1
fi

echo ""
echo "========================================"
echo "Job completed for $NETWORK"
echo "========================================"
