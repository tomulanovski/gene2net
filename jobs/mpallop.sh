#!/bin/bash
#SBATCH --job-name=MPALLOP_analysis
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/Ren_2024/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/Ren_2024/%x.e%j
#SBATCH --mem=32G
#SBATCH --cpus-per-task=12
#SBATCH --partition=itaym
#SBATCH --account=itaym-users
#SBATCH --time=24:00:00

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

MPALLOP="/groups/itay_mayrose/tomulanovski/phylonet/PhyloNet_3.8.2.jar"

INPUT_FILE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Ren_2024/gene_trees/mapaalop_test.nex"

java -jar $MPALLOP $INPUT_FILE


conda deactivate