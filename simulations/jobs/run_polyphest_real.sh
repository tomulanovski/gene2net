#!/bin/bash
#SBATCH --job-name=run_polyphest_real
#SBATCH --array=1-105
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_polyphest_real_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_polyphest_real_%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=50g
#SBATCH --cpus-per-task=16
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner
#SBATCH --constraint=localscratch

# ============================================================================
# RUN_POLYPHEST_REAL.SH
# ============================================================================
# Runs Polyphest with the REAL (ground-truth) ploidy prior instead of the
# inferred multiset. The real multiset is generated in-job from the network's
# MUL-tree (simulations/networks/<network>.tre), where each species appears
# exactly its true-copy-number times.
#
# Same gene trees as run_polyphest.sh (from polyphest_input). Only the multiset
# source differs, and output goes to a separate polyphest_real_p<PCT>/ dir so
# the inferred-prior results are left intact.
#
# Usage:
#   sbatch --export=CONFIG=conf_dup_loss_low_10M run_polyphest_real.sh
#   sbatch --export=CONFIG=conf_dup_loss_low_10M,PERCENTILE=70 run_polyphest_real.sh
#
# Environment variables:
#   CONFIG          - Configuration name (required, e.g., conf_dup_loss_low_10M)
#   NUM_REPLICATES  - Number of replicates (default: 5)
#   PERCENTILE      - Polyphest percentile parameter (default: 60)
#   ISO_THRESHOLD   - Isomorphic threshold (default: 0.2)
# ============================================================================

set -eo pipefail

# Initialize LD_LIBRARY_PATH if not set (needed for conda activation)
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ============================================================================
# CONFIGURATION
# ============================================================================

CONFIG="${CONFIG:?ERROR: CONFIG environment variable is required. Use --export=CONFIG=conf_name}"

NUM_REPLICATES="${NUM_REPLICATES:-5}"
PERCENTILE="${PERCENTILE:-60}"
ISO_THRESHOLD="${ISO_THRESHOLD:-0.2}"

NUM_NETWORKS=21

# ============================================================================
# PATHS
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
NETWORKS_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks"
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"
POLYPHEST="/groups/itay_mayrose/tomulanovski/Polyphest/Polyphest.py"
MULTISET_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/make_real_multiset.py"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# ============================================================================
# CALCULATE NETWORK AND REPLICATE FROM ARRAY TASK ID
# ============================================================================

task_id=${SLURM_ARRAY_TASK_ID}
network_idx=$(( (task_id - 1) / NUM_REPLICATES ))
replicate=$(( ((task_id - 1) % NUM_REPLICATES) + 1 ))

if [ $network_idx -ge $NUM_NETWORKS ]; then
    echo "ERROR: Invalid task ID ${task_id} - network index ${network_idx} exceeds ${NUM_NETWORKS}"
    exit 1
fi

network="${networks[$network_idx]}"

# ============================================================================
# PATHS FOR THIS JOB
# ============================================================================

INPUT_DIR="${BASE_DIR}/${network}/processed/${CONFIG}/polyphest_input/replicate_${replicate}"
OUTPUT_DIR="${BASE_DIR}/${network}/results/${CONFIG}/polyphest_real_p${PERCENTILE}/replicate_${replicate}"

GENE_TREES="${INPUT_DIR}/polyphest_trees.tre"
MULTISET="${INPUT_DIR}/multi_set_real.txt"
NETWORK_TRE="${NETWORKS_DIR}/${network}.tre"

# ============================================================================
# HEADER
# ============================================================================

echo "============================================================================"
echo "RUN POLYPHEST (REAL PLOIDY PRIOR)"
echo "============================================================================"
echo "Network: ${network} (index: ${network_idx})"
echo "Replicate: ${replicate}"
echo "Task ID: ${task_id}"
echo "Configuration: ${CONFIG}"
echo ""
echo "Parameters:"
echo "  Percentile: ${PERCENTILE}"
echo "  Isomorphic threshold: ${ISO_THRESHOLD}"
echo ""
echo "Network MUL-tree: ${NETWORK_TRE}"
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Date: $(date)"
echo "============================================================================"
echo ""

# ============================================================================
# VALIDATION
# ============================================================================

if [ ! -f "$GENE_TREES" ]; then
    echo "ERROR: Gene trees file not found: $GENE_TREES"
    echo "Did you run prep_polyphest.sh first?"
    exit 1
fi

if [ ! -f "$NETWORK_TRE" ]; then
    echo "ERROR: Network MUL-tree not found: $NETWORK_TRE"
    exit 1
fi

if [ ! -f "$MULTISET_SCRIPT" ]; then
    echo "ERROR: Multiset script not found: $MULTISET_SCRIPT"
    exit 1
fi

if [ ! -f "$POLYPHEST" ]; then
    echo "ERROR: Polyphest script not found: $POLYPHEST"
    exit 1
fi

echo "✓ Input files validated"
echo ""

# ============================================================================
# CREATE OUTPUT DIRECTORY
# ============================================================================

mkdir -p "$OUTPUT_DIR"
echo "✓ Output directory created: $OUTPUT_DIR"
echo ""

# ============================================================================
# GENERATE REAL MULTISET FROM MUL-TREE (needs ete3 -> gene2net env)
# ============================================================================

source "$CONDA_PATH/etc/profile.d/conda.sh"

conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

echo "Generating real multiset from MUL-tree..."
python "$MULTISET_SCRIPT" -i "$NETWORK_TRE" -m "$MULTISET"
echo "✓ Real multiset written: ${MULTISET}"
echo ""

# ============================================================================
# ACTIVATE POLYPHEST ENV
# ============================================================================

conda activate polyphest || {
    echo "ERROR: Could not activate polyphest environment"
    exit 1
}

echo "✓ Conda environment activated: polyphest"
echo "  Python: $(which python)"
echo ""

# ============================================================================
# RUN POLYPHEST
# ============================================================================

echo "Starting Polyphest..."
echo "============================================================================"
echo ""

# Send PuLP/CBC's multi-GB <uuid>-pulp.mps solver temp files to node-local
# scratch instead of /tmp, and delete them when the job ends (even if killed).
# Does not affect the Polyphest run itself (output goes to --output_dir below).
if [ -d /localscratch ]; then
    export TMPDIR="$(mktemp -d -p /localscratch polyphest.XXXXXX)"
else
    export TMPDIR="$(mktemp -d)"
fi
trap 'rm -rf "$TMPDIR"' EXIT
echo "Solver temp dir (TMPDIR): $TMPDIR"

start_time=$(date +%s)

python "$POLYPHEST" \
    --gene_tree_file "$GENE_TREES" \
    --consensus_multiset_file "$MULTISET" \
    --filter_strategy percentile \
    --output_dir "$OUTPUT_DIR" \
    --percentile "$PERCENTILE" \
    --use_near_isomorphic "True" \
    --isomorphic_threshold "$ISO_THRESHOLD"

exit_code=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

echo ""
echo "============================================================================"

if [ $exit_code -eq 0 ]; then
    echo "✓ POLYPHEST (REAL PRIOR) COMPLETED SUCCESSFULLY"
    echo ""
    echo "Network: ${network}"
    echo "Replicate: ${replicate}"
    echo "Duration: ${duration} seconds ($(echo "scale=2; $duration/3600" | bc) hours)"
    echo "Results saved to: $OUTPUT_DIR"
    echo ""
    echo "Output files:"
    ls -la "$OUTPUT_DIR" 2>/dev/null || echo "  (unable to list)"
else
    echo "✗ POLYPHEST (REAL PRIOR) FAILED"
    echo ""
    echo "Network: ${network}"
    echo "Replicate: ${replicate}"
    echo "Exit code: ${exit_code}"
    echo "Duration: ${duration} seconds"
fi

echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

exit $exit_code
