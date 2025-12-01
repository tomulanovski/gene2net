#!/bin/bash
#SBATCH --job-name=alisim
#SBATCH --array=1-1000
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alisim_1_mil_low_dup_%x_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alisim_1_mil_low_dup_%x_job_%a.err
#SBATCH --time=0:10:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2

# Usage: sbatch --job-name=alisim_Ding_2023 simulate_sequences.sh Ding_2023

# currently it looks just for gene trees that are in this format: g_treesxxxx.trees where xxxx are digits. the first gene tree name is g_trees0001.trees

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Activate your specific environment
conda activate gene2net

# Get network name from command line argument or job name
if [ -z "$1" ]; then
    # Extract network name from job name (after alisim_)
    NETWORK=$(echo $SLURM_JOB_NAME | sed 's/alisim_//')
else
    NETWORK=$1
fi

echo "Processing network: ${NETWORK}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
INPUT_DIR="${BASE_DIR}/${NETWORK}/data/low_dup_ne_100000_height_1_million_trees/1"
OUTPUT_DIR="${INPUT_DIR}/alignments"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Format the gene tree number with leading zeros (4 digits)
GENE_NUM=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

# Define input and output files
GENE_TREE="${INPUT_DIR}/g_trees${GENE_NUM}.trees"
OUTPUT_PREFIX="${OUTPUT_DIR}/alignment_${GENE_NUM}"

# Check if gene tree file exists
if [ ! -f "$GENE_TREE" ]; then
    echo "WARNING: Gene tree file not found: $GENE_TREE"
    echo "This might be beyond the number of gene trees generated."
    exit 0
fi

echo "Input gene tree: $GENE_TREE"
echo "Output alignment: ${OUTPUT_PREFIX}.phy"

# Run AliSim
iqtree --alisim "$OUTPUT_PREFIX" \
    -m "GTR+G" \
    -t "$GENE_TREE" \
    --length 100 \
    -seed $SLURM_ARRAY_TASK_ID

# Check if simulation was successful
if [ $? -eq 0 ]; then
    echo "Successfully simulated alignment: ${OUTPUT_PREFIX}.phy"
else
    echo "ERROR: AliSim failed for gene tree: $GENE_TREE"
    exit 1
fi

echo "Job completed for ${NETWORK}, gene tree ${GENE_NUM}"