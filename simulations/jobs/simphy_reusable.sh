#!/bin/bash
#SBATCH --array=1-21
#SBATCH --time=24:00:00
#SBATCH --mem=64g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# ============================================================================
# SIMULATION CONFIGURATION - EDIT THESE FOR DIFFERENT RUNS
# ============================================================================

# Configuration name (used in output directories and job tracking)
CONFIGURATION="${SIMPHY_CONFIG:-ils_low_10M}"

# Tree height (in generations)
TREE_HEIGHT="${SIMPHY_TREE_HEIGHT:-10M}"

# ILS parameters
NE="${SIMPHY_NE:-200000}"  # Effective population size

# Gene duplication/loss parameters (per gene per million years)
# Note: SimPhy uses rates per generation, so we'll convert
DUPLICATION_RATE_PER_MY="${SIMPHY_DUP:-0}"  # e.g., 0.001 for low
LOSS_RATE_PER_MY="${SIMPHY_LOSS:-0}"  # e.g., 0.001 for low
TRANSFER_RATE_PER_MY="${SIMPHY_TRANSFER:-0}"  # e.g., 0 for no HGT
GC_RATE_PER_MY="${SIMPHY_GC:-0}"  # e.g., 0 for no GC

# Generation time (years)
GENERATION_TIME=1

# Number of replicates
NUM_REPLICATES=5

# Timeout per SimPhy run (seconds)
TIMEOUT_SECONDS=3600

# ============================================================================
# CONVERT RATES: per MY ? per generation
# ============================================================================
# Formula: rate_per_gen = rate_per_MY / (1,000,000 / gen_time)
# For gen_time=1: rate_per_gen = rate_per_MY / 1,000,000

DUPLICATION_RATE=$(awk "BEGIN {print ${DUPLICATION_RATE_PER_MY} / 1000000}")
LOSS_RATE=$(awk "BEGIN {print ${LOSS_RATE_PER_MY} / 1000000}")
TRANSFER_RATE=$(awk "BEGIN {print ${TRANSFER_RATE_PER_MY} / 1000000}")
GC_RATE=$(awk "BEGIN {print ${GC_RATE_PER_MY} / 1000000}")

# ============================================================================
# LOAD CONDA ENVIRONMENT
# ============================================================================

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

echo "? Conda environment activated: gene2net"
echo "? Python: $(which python3)"
echo ""

# ============================================================================
# PATHS AND CONFIGURATION
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
PICKLE_FILE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions/substitution_rates.pkl"
SIMPHY_BIN="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulators/simphy/SimPhy_1.0.2/bin/simphy_lnx64"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

sample_substitution_rate() {
    python3 - "$1" << 'EOF'
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
EOF
}

run_simphy_single_replicate() {
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
        echo "  ? TIMEOUT: Exceeded ${timeout} seconds"
        return 1
    elif [ $exit_code -ne 0 ]; then
        echo "  ? FAILED: Exit code ${exit_code}"
        return 1
    fi
    
    if [ ! -d "${output_dir}/1" ]; then
        echo "  ? FAILED: Output directory not created"
        return 1
    fi
    
    if [ ! -f "${output_dir}/1/l_trees.trees" ]; then
        echo "  ? FAILED: No locus trees file created"
        return 1
    fi
    
    return 0
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

network_idx=$((SLURM_ARRAY_TASK_ID - 1))
network=${networks[$network_idx]}

echo "============================================================================"
echo "SimPhy Simulation"
echo "============================================================================"
echo "Network: ${network}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Configuration: ${CONFIGURATION}"
if [ "$TREE_HEIGHT" == "0" ]; then
    echo "Tree Height: single-label (10M height, no hybridization)"
else
    echo "Tree Height: ${TREE_HEIGHT}"
fi
echo ""
echo "Parameters:"
echo "  Effective Population Size (Ne): ${NE}"
echo "  Duplication rate: ${DUPLICATION_RATE_PER_MY} per MY (${DUPLICATION_RATE} per gen)"
echo "  Loss rate: ${LOSS_RATE_PER_MY} per MY (${LOSS_RATE} per gen)"
echo "  Transfer rate: ${TRANSFER_RATE_PER_MY} per MY (${TRANSFER_RATE} per gen)"
echo "  GC rate: ${GC_RATE_PER_MY} per MY (${GC_RATE} per gen)"
echo "  Generation time: ${GENERATION_TIME} year(s)"
echo "  Number of Replicates: ${NUM_REPLICATES}"
echo "  Timeout: ${TIMEOUT_SECONDS} seconds"
echo "Date: $(date)"
echo "============================================================================"
echo ""

# Define paths
if [ "$TREE_HEIGHT" == "0" ]; then
    # Single-label tree (no hybridization events, 10M height)
    SPECIES_TREE="${BASE_DIR}/${network}/single_label.nex"
elif [ "$TREE_HEIGHT" == "50M" ]; then
    SPECIES_TREE="${BASE_DIR}/${network}/tree_50_mil.nex"
else
    SPECIES_TREE="${BASE_DIR}/${network}/tree_10_mil.nex"
fi
OUTPUT_BASE="${BASE_DIR}/${network}/data/${CONFIGURATION}"

# Validation
if [ ! -f "$SPECIES_TREE" ]; then
    echo "ERROR: Species tree not found: $SPECIES_TREE"
    exit 1
fi

if [ ! -f "$PICKLE_FILE" ]; then
    echo "ERROR: Substitution rates pickle file not found: $PICKLE_FILE"
    exit 1
fi

if [ ! -x "$SIMPHY_BIN" ]; then
    echo "ERROR: SimPhy binary not found or not executable: $SIMPHY_BIN"
    exit 1
fi

echo "? Species tree: $SPECIES_TREE"
echo "? Pickle file: $PICKLE_FILE"
echo "? SimPhy binary: $SIMPHY_BIN"
echo ""

# Sample substitution rate
echo "Sampling substitution rate for all ${NUM_REPLICATES} replicates..."
SUB_RATE=$(sample_substitution_rate "$PICKLE_FILE")

if [ $? -ne 0 ] || [ -z "$SUB_RATE" ]; then
    echo "ERROR: Failed to sample substitution rate"
    exit 1
fi

echo "? Sampled substitution rate: ${SUB_RATE}"
echo ""

# Generate base seed
BASE_SEED=$((SLURM_ARRAY_TASK_ID * 100000 + RANDOM))
echo "Base seed: ${BASE_SEED}"
echo ""

# Batch configurations to try
BATCH_CONFIGS=(
    "1000:1"
    "1:1000"
)

# Maximum retries for single batch (1000 trees)
MAX_SINGLE_BATCH_RETRIES=10

# Maximum retries for small batch approach
MAX_BATCH_RETRIES=100

success=0

for config in "${BATCH_CONFIGS[@]}"; do
    IFS=':' read -r trees_per_run num_runs <<< "$config"
    total_trees=$((trees_per_run * num_runs))
    
    echo "============================================================================"
    echo "ATTEMPTING CONFIGURATION:"
    echo "  ${trees_per_run} trees ? ${num_runs} runs ? ${NUM_REPLICATES} replicates"
    echo "  = ${total_trees} trees per replicate"
    if [ $num_runs -gt 1 ]; then
        echo "  (with up to ${MAX_BATCH_RETRIES} retries per batch)"
    else
        echo "  (with up to ${MAX_SINGLE_BATCH_RETRIES} retries)"
    fi
    echo "============================================================================"
    echo ""
    
    if [ -d "$OUTPUT_BASE" ]; then
        echo "Cleaning previous attempt..."
        rm -rf "$OUTPUT_BASE"
    fi
    mkdir -p "$OUTPUT_BASE"
    
    all_runs_succeeded=1
    
    for replicate in $(seq 1 $NUM_REPLICATES); do
        echo "----------------------------------------------------------------------------"
        echo "REPLICATE ${replicate}/${NUM_REPLICATES}"
        echo "----------------------------------------------------------------------------"
        
        replicate_success=1
        
        for run in $(seq 1 $num_runs); do
            if [ $num_runs -eq 1 ]; then
                # Single batch mode - with retries
                RUN_OUTPUT="${OUTPUT_BASE}/replicate_${replicate}"
                if [ $run -eq 1 ]; then
                    echo "Mode: Single batch (${MAX_SINGLE_BATCH_RETRIES} retries)"
                    echo "Output: ${RUN_OUTPUT}"
                fi
                
                batch_succeeded=0
                
                for retry in $(seq 1 $MAX_SINGLE_BATCH_RETRIES); do
                    # Use different seed for each retry attempt
                    SEED=$((BASE_SEED + replicate * 10000 + run * 100 + retry))
                    
                    if [ $retry -eq 1 ]; then
                        echo "  Attempt ${retry}/${MAX_SINGLE_BATCH_RETRIES}"
                    else
                        echo "  Retry ${retry}/${MAX_SINGLE_BATCH_RETRIES}"
                        # Clean up failed attempt before retrying
                        if [ -d "$RUN_OUTPUT" ]; then
                            rm -rf "$RUN_OUTPUT"
                        fi
                    fi
                    
                    echo "    Seed: ${SEED}"
                    echo "    Trees: ${trees_per_run}"
                    echo "    Starting SimPhy (timeout: ${TIMEOUT_SECONDS}s)..."
                    
                    start_time=$(date +%s)
                    
                    run_simphy_single_replicate "$RUN_OUTPUT" "$trees_per_run" "$SUB_RATE" "$SEED" "$SPECIES_TREE" "$TIMEOUT_SECONDS"
                    
                    run_exit_code=$?
                    end_time=$(date +%s)
                    duration=$((end_time - start_time))
                    
                    if [ $run_exit_code -eq 0 ]; then
                        echo "  ? Completed on attempt ${retry} (duration: ${duration}s)"
                        batch_succeeded=1
                        break
                    else
                        echo "  ? Attempt ${retry} FAILED (duration: ${duration}s)"
                        if [ $retry -lt $MAX_SINGLE_BATCH_RETRIES ]; then
                            echo "  ? Retrying..."
                        fi
                    fi
                done
                
                if [ $batch_succeeded -eq 0 ]; then
                    echo "? Single batch FAILED after ${MAX_SINGLE_BATCH_RETRIES} attempts"
                    replicate_success=0
                    all_runs_succeeded=0
                    break
                fi
            else
                # Multiple batch mode - with retries
                RUN_OUTPUT="${OUTPUT_BASE}/replicate_${replicate}/batch_${run}"
                if [ $run -eq 1 ]; then
                    echo "Mode: Multiple batches (${num_runs} batches per replicate, ${MAX_BATCH_RETRIES} retries each)"
                fi
                echo "  Batch ${run}/${num_runs}"
                echo "  Output: ${RUN_OUTPUT}"
                
                batch_succeeded=0
                
                for retry in $(seq 1 $MAX_BATCH_RETRIES); do
                    # Use different seed for each retry attempt
                    SEED=$((BASE_SEED + replicate * 10000 + run * 100 + retry))
                    
                    if [ $retry -eq 1 ]; then
                        echo "    Attempt ${retry}/${MAX_BATCH_RETRIES}"
                    else
                        echo "    Retry ${retry}/${MAX_BATCH_RETRIES}"
                        # Clean up failed attempt before retrying
                        if [ -d "$RUN_OUTPUT" ]; then
                            rm -rf "$RUN_OUTPUT"
                        fi
                    fi
                    
                    echo "    Seed: ${SEED}"
                    echo "    Trees: ${trees_per_run}"
                    echo "    Starting SimPhy (timeout: ${TIMEOUT_SECONDS}s)..."
                    
                    start_time=$(date +%s)
                    
                    run_simphy_single_replicate "$RUN_OUTPUT" "$trees_per_run" "$SUB_RATE" "$SEED" "$SPECIES_TREE" "$TIMEOUT_SECONDS"
                    
                    run_exit_code=$?
                    end_time=$(date +%s)
                    duration=$((end_time - start_time))
                    
                    if [ $run_exit_code -eq 0 ]; then
                        echo "    ? Batch ${run} completed on attempt ${retry} (duration: ${duration}s)"
                        batch_succeeded=1
                        break
                    else
                        echo "    ? Attempt ${retry} FAILED (duration: ${duration}s)"
                        if [ $retry -lt $MAX_BATCH_RETRIES ]; then
                            echo "    ? Retrying..."
                        fi
                    fi
                done
                
                if [ $batch_succeeded -eq 0 ]; then
                    echo "  ? Batch ${run} FAILED after ${MAX_BATCH_RETRIES} attempts"
                    replicate_success=0
                    all_runs_succeeded=0
                    break
                fi
            fi
        done
        
        if [ $replicate_success -eq 1 ]; then
            echo "? Replicate ${replicate} completed successfully"
        else
            echo "? Replicate ${replicate} FAILED"
            break
        fi
        echo ""
    done
    
    if [ $all_runs_succeeded -eq 1 ]; then
        echo ""
        echo "============================================================================"
        echo "??? SUCCESS ???"
        echo "============================================================================"
        
        # Save configuration
        cat > "${OUTPUT_BASE}/simulation_config.txt" << EOF
================================================================================
Simulation Configuration Summary
================================================================================

Network: ${network}
Configuration: ${CONFIGURATION}
Tree Height: ${TREE_HEIGHT}
Date: $(date)

Parameters:
-----------
Effective Population Size (Ne): ${NE}
Duplication rate: ${DUPLICATION_RATE_PER_MY} per MY (${DUPLICATION_RATE} per gen)
Loss rate: ${LOSS_RATE_PER_MY} per MY (${LOSS_RATE} per gen)
Transfer rate: ${TRANSFER_RATE_PER_MY} per MY (${TRANSFER_RATE} per gen)
GC rate: ${GC_RATE_PER_MY} per MY (${GC_RATE} per gen)
Generation time: ${GENERATION_TIME} year(s)
Number of Replicates: ${NUM_REPLICATES}
Substitution rate (sampled): ${SUB_RATE}

Batch Configuration:
--------------------
Trees per run: ${trees_per_run}
Number of runs per replicate: ${num_runs}
Total trees per replicate: ${total_trees}
Max retries per batch: ${MAX_BATCH_RETRIES}

Base seed: ${BASE_SEED}
================================================================================
EOF
        
        success=1
        break
    else
        echo ""
        echo "============================================================================"
        echo "??? CONFIGURATION FAILED ???"
        echo "============================================================================"
        
        if [ "$config" == "1000:1" ]; then
            echo "? Falling back to smaller batches (1 tree ? 1000 runs with ${MAX_BATCH_RETRIES} retries each)"
        fi
        echo ""
    fi
done

# Final summary
echo ""
echo "============================================================================"
echo "FINAL SUMMARY - ${network}"
echo "============================================================================"

if [ $success -eq 1 ]; then
    echo "Status: ? SUCCESS"
    echo "Configuration: ${trees_per_run} trees ? ${num_runs} runs ? ${NUM_REPLICATES} replicates"
    echo "Output: ${OUTPUT_BASE}"
else
    echo "Status: ? FAILURE"
fi

echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

exit $((1 - success))