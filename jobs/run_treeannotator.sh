#!/bin/bash
#SBATCH --job-name=treeannotator
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/%x.e%j
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2

# Usage: sbatch run_treeannotator.sh /path/to/input.trees /path/to/output.tree

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.trees> <output.tree>"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=$2

# Count the number of trees (excluding the STATE_0 tree)
TOTAL_TREES=$(grep -c "tree STATE_" "$INPUT_FILE")
USABLE_TREES=$((TOTAL_TREES - 1))

# Calculate 10% burn-in (round down)
BURNIN=$((USABLE_TREES / 10))

echo "Total trees found: $TOTAL_TREES"
echo "Usable trees (excluding STATE_0): $USABLE_TREES"
echo "Using $BURNIN trees as burn-in (~10%)"

# Load conda
CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate alloppnet

# Run TreeAnnotator
treeannotator -burninTrees $BURNIN -heights mean "$INPUT_FILE" "$OUTPUT_FILE"

echo "Done. Summary tree written to: $OUTPUT_FILE"
