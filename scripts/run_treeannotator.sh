#!/bin/bash

# Usage: /groups/itay_mayrose/tomulanovski/gene2net/scripts/run_treeannotator.sh /path/to/input.trees /path/to/output.tree

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

# Load conda if not already available
source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

# Activate environment with TreeAnnotator
conda activate alloppnet

# Run TreeAnnotator
treeannotator -burninTrees $BURNIN -heights mean "$INPUT_FILE" "$OUTPUT_FILE"
