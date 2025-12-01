#!/bin/bash
#SBATCH --job-name=notung_array
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/notung_dup/%x_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/notung_dup/%x_%a.err
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner
#SBATCH --array=0-10

# Define all datasets
DATASETS=("Koenen_2020" "Lawrence_2016" "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Popp_2005" "Wu_2015" "Ren_2024" "Zhao_2021" "Marcussen_2015" "Morales-Briones_2021")

# Create log directory if it doesn't exist
mkdir -p /groups/itay_mayrose/tomulanovski/gene2net/papers/notung_dup

# Get current dataset based on array task ID
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}
echo "Processing dataset: $DATASET"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "========================================"

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/papers"

# Set paths for current dataset
GENE_TREES="${BASE_DIR}/${DATASET}/gene_trees/grampa_trees.tre"
SPECIES_TREE="${BASE_DIR}/${DATASET}/gene_trees/grampa_species_tree.tre"
TREES_DIR="${BASE_DIR}/${DATASET}/trees_for_dup"
OUTPUT_DIR="${BASE_DIR}/${DATASET}/reconciled_for_dup"
BATCH_FILE="${TREES_DIR}/notung_batch.txt"

# Notung jar path
NOTUNG_JAR="/groups/itay_mayrose/tomulanovski/Notung-2.9.1.5/Notung-2.9.1.5.jar"

# Python scripts
SPLIT_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/scripts/split_file_to_trees.py"
BATCH_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/scripts/create_notung_batch.py"

# Conda setup
CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Check if conda is properly sourced
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda command not found. Check your conda installation path."
    exit 1
fi

# Activate environment
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

echo "Current conda environment: $CONDA_PREFIX"

# Check if input files exist
if [ ! -f "$GENE_TREES" ]; then
    echo "ERROR: Gene trees file not found: $GENE_TREES"
    exit 1
fi

if [ ! -f "$SPECIES_TREE" ]; then
    echo "ERROR: Species tree file not found: $SPECIES_TREE"
    exit 1
fi

# Step 1: Split gene trees into individual files
echo ""
echo "Step 1: Splitting gene trees into individual files..."
mkdir -p "$TREES_DIR"
python "$SPLIT_SCRIPT" "$GENE_TREES" "$TREES_DIR"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to split trees for $DATASET"
    exit 1
fi

# Step 2: Create Notung batch file with species tree as first line
echo ""
echo "Step 2: Creating Notung batch file..."
python "$BATCH_SCRIPT" "$TREES_DIR" "$SPECIES_TREE"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create batch file for $DATASET"
    exit 1
fi

# Verify batch file was created
if [ ! -f "$BATCH_FILE" ]; then
    echo "ERROR: Batch file not created: $BATCH_FILE"
    exit 1
fi

echo "Batch file created at: $BATCH_FILE"
echo "First few lines of batch file:"
head -n 3 "$BATCH_FILE"

# Step 3: Run Notung reconciliation
echo ""
echo "Step 3: Running Notung reconciliation..."
mkdir -p "$OUTPUT_DIR"

# Run in headless mode (no GUI)
export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"

java -Xmx25G -jar ${NOTUNG_JAR} \
    -b ${BATCH_FILE} \
    --absfilenames \
    --reconcile \
    --treestats \
    --speciestag postfix \
    --outputdir ${OUTPUT_DIR} \
    --log
    
if [ $? -eq 0 ]; then
    echo ""
    echo "SUCCESS: Notung completed successfully for $DATASET"
    echo "Output directory: $OUTPUT_DIR"
else
    echo ""
    echo "ERROR: Notung failed for $DATASET"
    exit 1
fi

# Deactivate conda environment
conda deactivate

echo ""
echo "========================================"
echo "Job completed for $DATASET"
echo "========================================"