#!/bin/bash
#SBATCH --job-name=iqtree_genetrees
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2

# Usage:
#   sbatch run_iqtree.sh <input_dir> <output_dir>
#
# Runs IQ-TREE with GTR+G on all .fasta/.fa files in input_dir,
# writing gene trees (.treefile) to output_dir.
#
# Example:
#   sbatch run_iqtree.sh /path/to/alignments /path/to/gene_trees

set -euo pipefail

# --- Arguments ---
INPUT_DIR="${1:?ERROR: Please provide input directory as first argument}"
OUTPUT_DIR="${2:?ERROR: Please provide output directory as second argument}"

# --- Conda setup ---
CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"

if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda command not found."
    exit 1
fi

conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

echo "Current conda environment: $CONDA_PREFIX"

# --- Validate inputs ---
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# --- Run IQ-TREE on each alignment ---
COUNT=0
for ALIGNMENT in "$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fa; do
    # Skip if no files match the glob
    [ -e "$ALIGNMENT" ] || continue

    BASENAME=$(basename "$ALIGNMENT")
    NAME="${BASENAME%.*}"
    PREFIX="$OUTPUT_DIR/$NAME"

    echo "Processing: $BASENAME"
    iqtree2 -s "$ALIGNMENT" -m GTR+G -pre "$PREFIX" -nt "$SLURM_CPUS_PER_TASK" -quiet

    COUNT=$((COUNT + 1))
done

echo "Done. Processed $COUNT alignment(s). Trees written to: $OUTPUT_DIR"

conda deactivate
