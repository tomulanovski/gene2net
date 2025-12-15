#!/bin/bash
#SBATCH --job-name=alisim
#SBATCH --array=1-1000
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alisim_%x_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alisim_%x_job_%a.err
#SBATCH --time=0:30:00
#SBATCH --mem=2g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2

# Usage: sbatch --job-name=alisim_Ding_2023_ils_low_10M simulate_sequences_1_dataset.sh Ding_2023 ils_low_10M
#
# This script simulates DNA sequences along gene trees using GTR+Gamma parameters
# sampled from empirical distributions (2,709 genes from Zhao_2021, Ren_2024, Morales_Briones_2021).
#
# Arguments:
#   $1: Network name (e.g., "Ding_2023")
#   $2: Configuration name (e.g., "ils_low_10M")

# ============================================================================
# SETUP AND CONFIGURATION
# ============================================================================

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Activate environment
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# Get network name and configuration from command line arguments or job name
if [ -z "$1" ]; then
    # Extract network name from job name (format: alisim_NetworkName_ConfigName)
    NETWORK=$(echo $SLURM_JOB_NAME | sed 's/alisim_//' | cut -d'_' -f1-2)
else
    NETWORK=$1
fi

if [ -z "$2" ]; then
    # Try to extract configuration from job name
    CONFIGURATION=$(echo $SLURM_JOB_NAME | sed 's/alisim_//' | cut -d'_' -f3-)
    if [ -z "$CONFIGURATION" ]; then
        echo "ERROR: Configuration not specified"
        echo "Usage: sbatch --job-name=alisim_NetworkName_ConfigName $0 NetworkName ConfigName"
        exit 1
    fi
else
    CONFIGURATION=$2
fi

echo "============================================================================"
echo "AliSim Sequence Simulation"
echo "============================================================================"
echo "Network: ${NETWORK}"
echo "Configuration: ${CONFIGURATION}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "============================================================================"
echo ""

# ============================================================================
# PATHS
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
GTR_PICKLE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions/gtr_parameters_all.pkl"
SAMPLE_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/sequence_evolution/sample_gtr_parameters.py"

# Determine number of replicates by checking how many exist
REPLICATE_DIRS=("${BASE_DIR}/${NETWORK}/data/${CONFIGURATION}"/replicate_*)
NUM_REPLICATES=${#REPLICATE_DIRS[@]}

if [ $NUM_REPLICATES -eq 0 ]; then
    echo "ERROR: No replicate directories found in ${BASE_DIR}/${NETWORK}/data/${CONFIGURATION}"
    exit 1
fi

echo "Found ${NUM_REPLICATES} replicates"
echo ""

# ============================================================================
# VALIDATION
# ============================================================================

if [ ! -f "$GTR_PICKLE" ]; then
    echo "ERROR: GTR parameters file not found: $GTR_PICKLE"
    exit 1
fi

if [ ! -f "$SAMPLE_SCRIPT" ]; then
    echo "ERROR: Sampling script not found: $SAMPLE_SCRIPT"
    exit 1
fi

# ============================================================================
# PROCESS EACH REPLICATE
# ============================================================================

GENE_NUM=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
echo "Processing gene tree: g_trees${GENE_NUM}.trees"
echo ""

# Create a unique seed for parameter sampling (based on array task ID)
PARAM_SEED=$((SLURM_ARRAY_TASK_ID * 1000 + RANDOM))

# Sample GTR+Gamma parameters and alignment length from empirical distribution
echo "Sampling GTR+Gamma parameters from empirical distribution..."
echo "  Pickle file: ${GTR_PICKLE}"
echo "  Parameter seed: ${PARAM_SEED}"

# Sample parameters (output as bash variables)
PARAMS=$(python3 "$SAMPLE_SCRIPT" "$GTR_PICKLE" --seed $PARAM_SEED --output-format bash)

if [ $? -ne 0 ] || [ -z "$PARAMS" ]; then
    echo "ERROR: Failed to sample GTR parameters"
    exit 1
fi

# Load the sampled parameters into environment
eval "$PARAMS"

echo "  Model: ${IQTREE_MODEL}"
echo "  Alignment length: ${ALIGNMENT_LENGTH} bp"
echo "  Alpha: ${ALPHA}"
echo ""

# Process each replicate
for replicate in $(seq 1 $NUM_REPLICATES); do
    echo "----------------------------------------------------------------------------"
    echo "REPLICATE ${replicate}/${NUM_REPLICATES}"
    echo "----------------------------------------------------------------------------"

    REPLICATE_DIR="${BASE_DIR}/${NETWORK}/data/${CONFIGURATION}/replicate_${replicate}"
    OUTPUT_DIR="${REPLICATE_DIR}/1/alignments"

    # Create output directory if it doesn't exist
    mkdir -p "$OUTPUT_DIR"

    OUTPUT_PREFIX="${OUTPUT_DIR}/alignment_${GENE_NUM}"

    # ========================================================================
    # SMART GENE TREE DETECTION
    # SimPhy creates tree files numbered WITHIN each batch:
    # - Single batch: replicate_N/1/g_trees0001.trees ... g_trees1000.trees
    # - Batches of 10: Each batch has g_trees1.trees ... g_trees10.trees
    # - Batches of 1: Each batch has g_trees1.trees
    # ========================================================================

    GENE_TREE_FILE=""

    # Try single batch first: replicate_N/1/g_treesXXXX.trees
    SINGLE_BATCH_FILE="${REPLICATE_DIR}/1/g_trees${GENE_NUM}.trees"

    if [ -f "$SINGLE_BATCH_FILE" ]; then
        # Single batch mode - tree file exists directly
        GENE_TREE_FILE="$SINGLE_BATCH_FILE"
        echo "  Mode: Single batch"
        echo "  File: replicate_${replicate}/1/g_trees${GENE_NUM}.trees"

    else
        # Multi-batch mode - calculate which batch and tree number within batch
        # Count batches to determine trees per batch
        NUM_BATCHES=$(ls -d ${REPLICATE_DIR}/batch_* 2>/dev/null | wc -l)

        if [ $NUM_BATCHES -eq 0 ]; then
            echo "  ERROR: No gene trees found in ${REPLICATE_DIR}"
            exit 1
        fi

        # Calculate trees per batch (1000 total / number of batches)
        TREES_PER_BATCH=$((1000 / NUM_BATCHES))

        # Calculate which batch contains this gene tree (1-indexed)
        # Gene tree 42 with batches of 10: batch = ceiling(42/10) = 5
        BATCH_NUM=$(( (SLURM_ARRAY_TASK_ID - 1) / TREES_PER_BATCH + 1 ))

        # Calculate tree number WITHIN that batch (1-indexed)
        # Gene tree 42 with batches of 10: tree_in_batch = ((42-1) % 10) + 1 = 2
        TREE_IN_BATCH=$(( (SLURM_ARRAY_TASK_ID - 1) % TREES_PER_BATCH + 1 ))

        # Format tree number (no leading zeros for batch tree numbering)
        GENE_TREE_FILE="${REPLICATE_DIR}/batch_${BATCH_NUM}/1/g_trees${TREE_IN_BATCH}.trees"

        if [ ! -f "$GENE_TREE_FILE" ]; then
            echo "  ERROR: Gene tree file not found: $GENE_TREE_FILE"
            echo "  Calculated: Gene tree ${SLURM_ARRAY_TASK_ID} → batch ${BATCH_NUM}, tree ${TREE_IN_BATCH}"
            exit 1
        fi

        echo "  Mode: ${NUM_BATCHES} batches of ${TREES_PER_BATCH} trees each"
        echo "  Gene tree ${SLURM_ARRAY_TASK_ID} → batch_${BATCH_NUM}/1/g_trees${TREE_IN_BATCH}.trees"
    fi

    echo "  Output: ${OUTPUT_PREFIX}.phy"

    # Create unique seed for this replicate's simulation
    SIM_SEED=$((SLURM_ARRAY_TASK_ID * 10000 + replicate * 1000 + RANDOM % 1000))
    echo "  Simulation seed: ${SIM_SEED}"

    # Run AliSim with sampled empirical parameters
    # The gene tree file contains exactly ONE tree
    iqtree --alisim "$OUTPUT_PREFIX" \
        -m "$IQTREE_MODEL" \
        -t "$GENE_TREE_FILE" \
        --length "$ALIGNMENT_LENGTH" \
        -seed "$SIM_SEED" \
        --quiet

    # Check if simulation was successful
    if [ $? -eq 0 ] && [ -f "${OUTPUT_PREFIX}.phy" ]; then
        echo "  ✓ Success: ${OUTPUT_PREFIX}.phy"
    else
        echo "  ✗ ERROR: AliSim failed for gene tree: $GENE_TREE"
        exit 1
    fi
    echo ""
done

echo "============================================================================"
echo "Completed gene tree ${GENE_NUM} across all ${NUM_REPLICATES} replicates"
echo "============================================================================"
echo "Sampled parameters:"
echo "  Model: ${IQTREE_MODEL}"
echo "  Length: ${ALIGNMENT_LENGTH} bp"
echo "  Alpha: ${ALPHA}"
echo "============================================================================"
