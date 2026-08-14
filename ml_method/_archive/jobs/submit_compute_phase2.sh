#!/bin/bash
# ==============================================================================
# Submit Phase 2 feature extraction: compute gene tree summary features
#
# USAGE:
#   ./submit_compute_phase2.sh [--input-dir DIR] [--output-dir DIR] [--general] [--max-gene-trees N]
#
# EXAMPLES:
#   ./submit_compute_phase2.sh                              # defaults (ils_low)
#   ./submit_compute_phase2.sh --general                    # use general partition
#   ./submit_compute_phase2.sh --max-gene-trees 300         # limit gene trees
# ==============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
DEFAULT_INPUT="${BASE_DIR}/data/mul_trees_2k/training/ils_low"
DEFAULT_OUTPUT="${BASE_DIR}/data/mul_trees_2k/training_phase2/ils_low"

# Parse arguments
INPUT_DIR=""
OUTPUT_DIR=""
USE_GENERAL=false
MAX_GENE_TREES=0

while [ $# -gt 0 ]; do
    case "$1" in
        --general) USE_GENERAL=true; shift;;
        --input-dir) INPUT_DIR="$2"; shift 2;;
        --output-dir) OUTPUT_DIR="$2"; shift 2;;
        --max-gene-trees) MAX_GENE_TREES="$2"; shift 2;;
        *) shift;;
    esac
done

INPUT_DIR="${INPUT_DIR:-$DEFAULT_INPUT}"
OUTPUT_DIR="${OUTPUT_DIR:-$DEFAULT_OUTPUT}"

# Set partition
if [ "$USE_GENERAL" = true ]; then
    PARTITION="gpu-general-pool"
    ACCOUNT=""
    QOS="public"
else
    PARTITION="gpu-rotemhsh-pool"
    ACCOUNT="--account=itaym-users_v2"
    QOS="owner"
fi

LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

echo "================================================"
echo "Phase 2: Gene Tree Feature Extraction"
echo "================================================"
echo "Input: ${INPUT_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo "Max gene trees: ${MAX_GENE_TREES}"
echo "Partition: ${PARTITION}"
echo "================================================"

sbatch \
    --job-name="gene2net_phase2_feat" \
    --output="${LOG_DIR}/phase2_feat_%j.out" \
    --error="${LOG_DIR}/phase2_feat_%j.err" \
    --time=02:00:00 \
    --mem=16G \
    --cpus-per-task=2 \
    --partition=${PARTITION} \
    ${ACCOUNT} \
    --qos=${QOS} \
    --export=ALL,INPUT_DIR="${INPUT_DIR}",OUTPUT_DIR="${OUTPUT_DIR}",BASE_DIR="${BASE_DIR}",MAX_GENE_TREES="${MAX_GENE_TREES}" \
    "${BASE_DIR}/jobs/compute_phase2_job.sh"
