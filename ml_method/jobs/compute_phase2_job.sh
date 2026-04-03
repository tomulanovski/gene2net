#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

# Compute Phase 2 gene tree summary features for training samples.
# CPU-only — no GPU needed.

set -euo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate final_project

export PYTHONPATH="${PYTHONPATH:-}:${BASE_DIR}"
export PYTHONUNBUFFERED=1

echo "============================================================================"
echo "Phase 2: Computing Gene Tree Summary Features"
echo "============================================================================"
echo "Input: ${INPUT_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo "Max gene trees: ${MAX_GENE_TREES:-0}"
echo "Date: $(date)"
echo "============================================================================"

python "${BASE_DIR}/scripts/compute_phase2_features.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --max-gene-trees "${MAX_GENE_TREES:-0}"
