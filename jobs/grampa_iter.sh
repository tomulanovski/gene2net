#!/bin/bash
#SBATCH --job-name=grampa_iter
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/Lawrence_2016/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/Lawrence_2016/%x.e%j
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym
#SBATCH --account=itaym-users


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

# Verify the environment was activated
echo "Current conda environment: $CONDA_PREFIX"

# Change paths according to specific run!!

GENE_TREES="/groups/itay_mayrose/tomulanovski/gene2net/papers/Lawrence_2016/gene_trees/grampa_trees.tre"
SPECIES_TREE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Lawrence_2016/gene_trees/astral_4.tre"
OUTPUT_DIR="/groups/itay_mayrose/tomulanovski/gene2net/papers/Lawrence_2016/networks/grampa/analysis/iter"


# Check if GRAMPA is available
if ! command -v grampa &> /dev/null; then
    echo "ERROR: grampa command not found. Check if it's installed in your conda environment."
    exit 1
fi

# Run GRAMPA
python "/groups/itay_mayrose/ronenshtein/Grandma/grampa.py" -g $GENE_TREES -s $SPECIES_TREE -o $OUTPUT_DIR --overwrite
# Deactivate conda environment
conda deactivate