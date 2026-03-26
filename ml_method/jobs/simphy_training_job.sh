#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ==============================================================================
# SimPhy gene tree simulation for one MUL-tree (ML training data).
# Called by submit_simphy_training.sh — do not run directly.
#
# Mirrors simulations/jobs/simphy_reusable.sh but adapted for generated MUL-trees.
# ==============================================================================

set -euo pipefail

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# ============================================================================
# CONFIGURATION
# ============================================================================

CONFIGURATION="${CONFIG_NAME:-ils_low}"
NE="${NE_VAL:-200000}"
NUM_TREES="${NUM_GENE_TREES:-500}"
NUM_REPS="${NUM_REPLICATES:-1}"

# Gene dup/loss — default 0 for training (pure ILS)
DUPLICATION_RATE=0
LOSS_RATE=0
TRANSFER_RATE=0
GC_RATE=0
GENERATION_TIME=1

# Timeout per SimPhy run
TIMEOUT_SECONDS=1800

# ============================================================================
# LOAD CONDA
# ============================================================================

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || { echo "ERROR: Could not activate gene2net"; exit 1; }

# ============================================================================
# PATHS
# ============================================================================

SIMPHY_BIN="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulators/simphy/SimPhy_1.0.2/bin/simphy_lnx64"
PICKLE_FILE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions/substitution_rates.pkl"

TREE_INDEX=$(( SLURM_ARRAY_TASK_ID + ${INDEX_OFFSET:-0} ))
IDX=$(printf "%04d" "$TREE_INDEX")
SPECIES_TREE="${MUL_TREES_DIR}/mul_tree_${IDX}.nex"
OUTPUT_BASE="${MUL_TREES_DIR}/simphy/${CONFIGURATION}/${IDX}"

# ============================================================================
# VALIDATION
# ============================================================================

if [ ! -f "$SPECIES_TREE" ]; then
    echo "ERROR: MUL-tree not found: $SPECIES_TREE"
    exit 1
fi

if [ ! -x "$SIMPHY_BIN" ]; then
    echo "ERROR: SimPhy binary not found: $SIMPHY_BIN"
    exit 1
fi

echo "============================================================================"
echo "SimPhy Training - MUL-tree ${IDX}"
echo "============================================================================"
echo "Species tree: ${SPECIES_TREE}"
echo "Config: ${CONFIGURATION} | Ne: ${NE}"
echo "Gene trees: ${NUM_TREES} | Replicates: ${NUM_REPS}"
echo "Output: ${OUTPUT_BASE}"
echo "Date: $(date)"
echo "============================================================================"

# ============================================================================
# SAMPLE SUBSTITUTION RATE (same as simphy_reusable.sh)
# ============================================================================

sample_substitution_rate() {
    python3 - "$1" << 'PYEOF'
import pickle
import random
import sys

pickle_file = sys.argv[1]
try:
    with open(pickle_file, 'rb') as f:
        rates_dict = pickle.load(f)
    all_rates = []
    for dataset_name, rates_data in rates_dict.items():
        all_rates.extend([r['rate'] for r in rates_data])
    sampled_rate = random.choice(all_rates)
    print(f"{sampled_rate:.15e}")
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF
}

SUB_RATE=$(sample_substitution_rate "$PICKLE_FILE")
if [ $? -ne 0 ] || [ -z "$SUB_RATE" ]; then
    echo "ERROR: Failed to sample substitution rate"
    exit 1
fi
echo "Substitution rate: ${SUB_RATE}"

# ============================================================================
# RUN SIMPHY
# ============================================================================

run_simphy() {
    local output_dir=$1
    local num_trees=$2
    local sub_rate=$3
    local seed=$4
    local species_tree=$5
    local timeout=$6

    mkdir -p "$output_dir"

    timeout ${timeout} $SIMPHY_BIN \
        -rs 1 \
        -rl f:${num_trees} \
        -rg 1 \
        -SR "$species_tree" \
        -sp f:${NE} \
        -su f:${sub_rate} \
        -sg f:${GENERATION_TIME} \
        -lb f:${DUPLICATION_RATE} \
        -ld f:${LOSS_RATE} \
        -lt f:${TRANSFER_RATE} \
        -lg f:${GC_RATE} \
        -cs ${seed} \
        -o "$output_dir" \
        -v 1 \
        -od 1 \
        -op 1 \
        -oc 1

    local exit_code=$?

    if [ $exit_code -eq 124 ]; then
        echo "  TIMEOUT: Exceeded ${timeout} seconds"
        return 1
    elif [ $exit_code -ne 0 ]; then
        echo "  FAILED: Exit code ${exit_code}"
        return 1
    fi

    if [ ! -f "${output_dir}/1/l_trees.trees" ]; then
        echo "  FAILED: No locus trees file created"
        return 1
    fi

    return 0
}

BASE_SEED=$((TREE_INDEX * 100000 + 42))
MAX_RETRIES=10

# Clean previous output
if [ -d "$OUTPUT_BASE" ]; then
    rm -rf "$OUTPUT_BASE"
fi

success=0

for replicate in $(seq 1 $NUM_REPS); do
    echo ""
    echo "--- Replicate ${replicate}/${NUM_REPS} ---"

    RUN_OUTPUT="${OUTPUT_BASE}/replicate_${replicate}"
    batch_succeeded=0

    for retry in $(seq 1 $MAX_RETRIES); do
        SEED=$((BASE_SEED + replicate * 10000 + retry))

        if [ $retry -gt 1 ]; then
            echo "  Retry ${retry}/${MAX_RETRIES}"
            [ -d "$RUN_OUTPUT" ] && rm -rf "$RUN_OUTPUT"
        fi

        run_simphy "$RUN_OUTPUT" "$NUM_TREES" "$SUB_RATE" "$SEED" "$SPECIES_TREE" "$TIMEOUT_SECONDS"

        if [ $? -eq 0 ]; then
            echo "  Replicate ${replicate} completed (attempt ${retry})"
            batch_succeeded=1
            break
        fi
    done

    if [ $batch_succeeded -eq 0 ]; then
        echo "  Replicate ${replicate} FAILED after ${MAX_RETRIES} attempts"
        exit 1
    fi
done

# Save config summary
cat > "${OUTPUT_BASE}/simulation_config.txt" << EOF
MUL-tree: ${SPECIES_TREE}
Configuration: ${CONFIGURATION}
Ne: ${NE}
Gene trees: ${NUM_TREES}
Replicates: ${NUM_REPS}
Substitution rate: ${SUB_RATE}
Base seed: ${BASE_SEED}
Date: $(date)
EOF

echo ""
echo "SUCCESS - MUL-tree ${IDX} complete"
echo "Output: ${OUTPUT_BASE}"
