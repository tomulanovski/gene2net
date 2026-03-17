#!/bin/bash
#SBATCH --job-name=grandma_split_empirical
#SBATCH --array=1-12
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/logs/grandma_split_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/logs/grandma_split_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# RUN_GRANDMA_SPLIT_EMPIRICAL.SH
# ============================================================================
# Runs GRANDMA_SPLIT on empirical datasets.
# Input gene trees and species tree are taken from the same locations as GRAMPA.
# Uses a SLURM array: one job per dataset (12 datasets total).
#
# Usage:
#   sbatch run_grandma_split_empirical.sh
# ============================================================================

set -eo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ============================================================================
# CONFIGURATION
# ============================================================================

CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"
GRANDMA_SPLIT_SCRIPT="/groups/itay_mayrose/ronenshtein/Grandma/gran.py"
PAPERS_DIR="/groups/itay_mayrose/tomulanovski/gene2net/papers"

datasets=(
  "Ren_2024"
  "Zhao_2021"
  "Ding_2023"
  "Koenen_2020"
  "Marcussen_2012"
  "Wu_2015"
  "Lawrence_2016"
  "Popp_2005"
  "Diaz-Perez_2018"
  "Marcussen_2015"
  "Morales-Briones_2021"
  "Wisecaver_2023"
)

# ============================================================================
# DETERMINE DATASET FROM ARRAY TASK ID
# ============================================================================

dataset_idx=$(( SLURM_ARRAY_TASK_ID - 1 ))
dataset="${datasets[$dataset_idx]}"

# ============================================================================
# PATHS FOR THIS JOB
# ============================================================================

# Datasets with stripped internal node labels (rooting added labels that
# interfere with GRAMPA/GRANDMA_SPLIT; use grampa_trees_striped.tre for these)
stripped_datasets=("Ren_2024" "Koenen_2020" "Wu_2015" "Popp_2005" "Marcussen_2015" "Wisecaver_2023")

gene_trees_file="grampa_trees.tre"
for d in "${stripped_datasets[@]}"; do
    if [ "$dataset" = "$d" ]; then
        gene_trees_file="grampa_trees_striped.tre"
        break
    fi
done

GENE_TREES="${PAPERS_DIR}/${dataset}/gene_trees/${gene_trees_file}"
SPECIES_TREE="${PAPERS_DIR}/${dataset}/gene_trees/grampa_species_tree.tre"
OUTPUT_DIR="${PAPERS_DIR}/${dataset}/networks/grandma_split"

# ============================================================================
# HEADER
# ============================================================================

echo "============================================================================"
echo "RUN GRANDMA_SPLIT (EMPIRICAL)"
echo "============================================================================"
echo "Dataset:          ${dataset}"
echo "Task ID:          ${SLURM_ARRAY_TASK_ID}"
echo "Gene trees:       ${GENE_TREES} (${gene_trees_file})"
echo "Species tree:     ${SPECIES_TREE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Date:             $(date)"
echo "============================================================================"
echo ""

# ============================================================================
# VALIDATION
# ============================================================================

if [ ! -f "$GENE_TREES" ]; then
    echo "ERROR: Gene trees file not found: $GENE_TREES"
    exit 1
fi

if [ ! -f "$SPECIES_TREE" ]; then
    echo "ERROR: Species tree file not found: $SPECIES_TREE"
    exit 1
fi

if [ ! -f "$GRANDMA_SPLIT_SCRIPT" ]; then
    echo "ERROR: GRANDMA_SPLIT script not found: $GRANDMA_SPLIT_SCRIPT"
    exit 1
fi

echo "✓ Input files validated"
echo ""

# ============================================================================
# ACTIVATE CONDA
# ============================================================================

source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

echo "✓ Conda environment activated: gene2net"

if ! command -v python &> /dev/null; then
    echo "ERROR: Python command not found. Check conda environment."
    exit 1
fi

echo "✓ Python found: $(which python)"
echo ""

# ============================================================================
# CREATE OUTPUT DIRECTORY
# ============================================================================

mkdir -p "$OUTPUT_DIR"
echo "✓ Output directory created: $OUTPUT_DIR"
echo ""

# ============================================================================
# RUN GRANDMA_SPLIT
# ============================================================================

echo "Starting GRANDMA_SPLIT..."
echo "============================================================================"
echo ""

start_time=$(date +%s)

python "$GRANDMA_SPLIT_SCRIPT" \
    -g "$GENE_TREES" \
    -s "$SPECIES_TREE" \
    -o "$OUTPUT_DIR" \
    -m split \
    -c 15   \
    -i 50 \
    -p ${SLURM_CPUS_PER_TASK} \
    --plot \
    --debug \
    --v 3 \
    --overwrite

exit_code=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "============================================================================"

if [ $exit_code -eq 0 ]; then
    echo "✓ GRANDMA_SPLIT COMPLETED SUCCESSFULLY"
    echo ""
    echo "Dataset:  ${dataset}"
    echo "Duration: ${duration} seconds ($(echo "scale=2; $duration/3600" | bc) hours)"
    echo "Results:  ${OUTPUT_DIR}"
    echo ""
    echo "Output files:"
    ls -la "$OUTPUT_DIR" 2>/dev/null || echo "  (unable to list)"
else
    echo "✗ GRANDMA_SPLIT FAILED"
    echo ""
    echo "Dataset:    ${dataset}"
    echo "Exit code:  ${exit_code}"
    echo "Duration:   ${duration} seconds"
fi

echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

exit $exit_code
