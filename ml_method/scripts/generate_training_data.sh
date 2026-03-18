#!/bin/bash
#SBATCH --job-name=gen_training
#SBATCH --output=logs/generate_%A_%a.out
#SBATCH --error=logs/generate_%A_%a.err
#SBATCH --array=0-4499
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=itaym

# Training data generation for Gene2Net-GNN
# Generates 500 base examples × 9 augmentations (3 ILS × 3 GDL) = 4500 total
#
# Usage:
#   sbatch generate_training_data.sh /path/to/output_dir
#
#   # Or for a smaller test run:
#   sbatch --array=0-9 generate_training_data.sh /path/to/output_dir

set -euo pipefail

OUTPUT_BASE="${1:-/groups/itay_mayrose/tomulanovski/gene2net/ml_method/data/training}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Configuration
N_BASE_EXAMPLES=500
ILS_LEVELS=("low" "medium" "high")
N_ILS=${#ILS_LEVELS[@]}

# Compute which base example and ILS level this job handles
BASE_IDX=$(( SLURM_ARRAY_TASK_ID / N_ILS ))
ILS_IDX=$(( SLURM_ARRAY_TASK_ID % N_ILS ))
ILS_LEVEL="${ILS_LEVELS[$ILS_IDX]}"

# Vary number of species and events
# Use base_idx to deterministically vary parameters
N_SPECIES=$(( 5 + (BASE_IDX % 16) * 5 ))  # 5, 10, 15, ..., 80
if [ $N_SPECIES -gt 100 ]; then N_SPECIES=100; fi

N_EVENTS=$(( (BASE_IDX / 16) % 16 ))  # 0, 1, 2, ..., 15
if [ $N_EVENTS -gt 15 ]; then N_EVENTS=15; fi

SEED=$((BASE_IDX * 1000 + ILS_IDX))

EXAMPLE_DIR="${OUTPUT_BASE}/example_$(printf '%04d' $SLURM_ARRAY_TASK_ID)"

echo "=== Generating example ${SLURM_ARRAY_TASK_ID} ==="
echo "  Species: ${N_SPECIES}, Events: ${N_EVENTS}, ILS: ${ILS_LEVEL}, Seed: ${SEED}"

# Activate conda
source activate gene2net

# Create log dir
mkdir -p "$(dirname "$0")/../logs"

# Run generation
python "${SCRIPT_DIR}/generate_one_example.py" \
    --output-dir "${EXAMPLE_DIR}" \
    --n-species "${N_SPECIES}" \
    --n-events "${N_EVENTS}" \
    --ils-level "${ILS_LEVEL}" \
    --n-gene-trees 1000 \
    --seed "${SEED}"

echo "=== Done ==="
