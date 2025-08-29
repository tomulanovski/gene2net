#!/bin/bash
#SBATCH --job-name=MPALLOP_reformat
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/Ding_2023/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/Ding_2023/%x.e%j
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=power-general
#SBATCH --account=power-general-users

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Check if conda is properly sourced
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda command not found. Check your conda installation path."
    exit 1
fi


# Activate your specific environment
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

INPUT_FILE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Ding_2023/gene_trees/polyphest_input_gene_trees.tre"
OUTPUT_FILE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Ding_2023/gene_trees/mpallop_input_gene_trees.tre"


java -jar "/groups/itay_mayrose/tomulanovski/phylonet/TreeConversion.jar" $INPUT_FILE $OUTPUT_FILE 
