#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ==============================================================================
# Run ASTRAL on SimPhy gene trees for one MUL-tree (ML training data).
# Called by submit_astral_training.sh — do not run directly.
#
# For each MUL-tree index:
#   1. Collects individual gene tree files from SimPhy output
#   2. Concatenates into a single gene_trees.tre file
#   3. Runs ASTRAL 4 to infer the species tree
#
# Our species names (sp_0, sp_1, ...) are clean — no need for the
# substring fixing / GRAMPA preprocessing that simulations/ pipeline does.
# ==============================================================================

set -eo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

CONFIGURATION="${CONFIG_NAME:?ERROR: CONFIG_NAME required}"
MUL_TREES_DIR="${MUL_TREES_DIR:?ERROR: MUL_TREES_DIR required}"
NUM_REPS="${NUM_REPLICATES:-1}"

# ============================================================================
# PATHS
# ============================================================================

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || { echo "ERROR: Could not activate gene2net"; exit 1; }

ASTRAL_CMD="astral4"
if ! command -v $ASTRAL_CMD &> /dev/null; then
    echo "ERROR: $ASTRAL_CMD not found in conda environment"
    exit 1
fi

IDX=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")
SIMPHY_DIR="${MUL_TREES_DIR}/simphy/${CONFIGURATION}/${IDX}"

echo "============================================================================"
echo "ASTRAL Training - MUL-tree ${IDX}"
echo "============================================================================"
echo "Config: ${CONFIGURATION}"
echo "SimPhy dir: ${SIMPHY_DIR}"
echo "Date: $(date)"
echo "============================================================================"
echo ""

if [ ! -d "$SIMPHY_DIR" ]; then
    echo "ERROR: SimPhy output not found: $SIMPHY_DIR"
    echo "Did you run submit_simphy_training.sh first?"
    exit 1
fi

# ============================================================================
# PROCESS EACH REPLICATE
# ============================================================================

success_count=0

for replicate in $(seq 1 $NUM_REPS); do
    echo "--- Replicate ${replicate}/${NUM_REPS} ---"

    REP_DIR="${SIMPHY_DIR}/replicate_${replicate}"
    INNER_DIR="${REP_DIR}/1"

    if [ ! -d "$INNER_DIR" ]; then
        echo "  ERROR: No SimPhy output at ${INNER_DIR}"
        continue
    fi

    # Step 1: Collect and concatenate gene trees
    GENE_TREES_FILE="${REP_DIR}/gene_trees.tre"
    echo "  Collecting gene trees..."

    > "$GENE_TREES_FILE"  # empty the file
    tree_count=0

    for gt_file in "$INNER_DIR"/g_trees*.trees; do
        if [ ! -f "$gt_file" ]; then
            continue
        fi
        # Each g_treesN.trees file contains one gene tree per line
        cat "$gt_file" >> "$GENE_TREES_FILE"
        tree_count=$((tree_count + 1))
    done

    if [ $tree_count -eq 0 ]; then
        echo "  ERROR: No gene tree files found in ${INNER_DIR}"
        continue
    fi

    total_trees=$(wc -l < "$GENE_TREES_FILE")
    echo "  Collected ${tree_count} files, ${total_trees} gene trees total"

    # Step 2: Run ASTRAL
    SPECIES_TREE_FILE="${REP_DIR}/astral_species.tre"
    ASTRAL_LOG="${REP_DIR}/astral.log"

    echo "  Running ASTRAL..."
    start_time=$(date +%s)

    $ASTRAL_CMD -o "$SPECIES_TREE_FILE" "$GENE_TREES_FILE" 2>"$ASTRAL_LOG"
    exit_code=$?

    end_time=$(date +%s)
    duration=$((end_time - start_time))

    if [ $exit_code -eq 0 ] && [ -f "$SPECIES_TREE_FILE" ]; then
        echo "  ASTRAL completed (${duration}s)"
        echo "  Species tree: $SPECIES_TREE_FILE"
        success_count=$((success_count + 1))
    else
        echo "  ASTRAL failed (exit code: $exit_code)"
        echo "  Log: $ASTRAL_LOG"
    fi
    echo ""
done

echo "============================================================================"
echo "MUL-tree ${IDX}: ${success_count}/${NUM_REPS} replicates completed"
echo "============================================================================"

[ $success_count -eq $NUM_REPS ] || exit 1
