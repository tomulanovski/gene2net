#!/bin/bash
#SBATCH --job-name=simphy_ils_low
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/simphy_ils_low_10M_net_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/simphy_ils_low_10M_net_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=64g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

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
echo "? Python version: $(python3 --version)"
echo ""

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
PICKLE_FILE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions/substitution_rates.pkl"
SIMPHY_BIN="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulators/simphy/SimPhy_1.0.2/bin/simphy_lnx64"

# Simulation parameters
TREE_HEIGHT="10M"
CONFIGURATION="ils_low"
NUM_REPLICATES=5
NE=10000  # REMEMBER: This gives very low ILS! You may want 500000-1000000

# Timeout for each SimPhy run (in seconds)
TIMEOUT_SECONDS=3600  # 1 hour

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
    
    # Combine all rates from all datasets
    all_rates = []
    for dataset_name, rates_data in rates_dict.items():
        all_rates.extend([r['rate'] for r in rates_data])
    
    # Sample one rate
    sampled_rate = random.choice(all_rates)
    
    # Print with high precision (15 decimal places in scientific notation)
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
    
    # Run SimPhy for SINGLE replicate (note: -rs 1, not -rs 5)
    timeout ${timeout} $SIMPHY_BIN \
        -rs 1 \
        -rl f:${num_trees} \
        -rg 1 \
        -SR "$species_tree" \
        -sp f:${NE} \
        -su f:${sub_rate} \
        -sg f:1 \
        -lb f:0 \
        -ld f:0 \
        -lt f:0 \
        -lg f:0 \
        -cs ${seed} \
        -o "$output_dir" \
        -v 1 \
        -od 1 \
        -op 1 \
        -oc 1
    
    local exit_code=$?
    
    # Check exit status
    if [ $exit_code -eq 124 ]; then
        echo "  ? TIMEOUT: Exceeded ${timeout} seconds"
        return 1
    elif [ $exit_code -ne 0 ]; then
        echo "  ? FAILED: Exit code ${exit_code}"
        return 1
    fi
    
    # Verify output was created (should be folder "1")
    if [ ! -d "${output_dir}/1" ]; then
        echo "  ? FAILED: Output directory not created"
        return 1
    fi
    
    # Check if trees were created
    if [ ! -f "${output_dir}/1/l_trees.trees" ]; then
        echo "  ? FAILED: No locus trees file created"
        return 1
    fi
    
    return 0
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Get network for this array task
network_idx=$((SLURM_ARRAY_TASK_ID - 1))
network=${networks[$network_idx]}

echo "============================================================================"
echo "SimPhy Simulation - ILS Only (Low)"
echo "============================================================================"
echo "Network: ${network}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Configuration: ${CONFIGURATION}"
echo "Tree Height: ${TREE_HEIGHT}"
echo "Effective Population Size (Ne): ${NE}"
echo "Number of Replicates: ${NUM_REPLICATES}"
echo "Timeout per attempt: ${TIMEOUT_SECONDS} seconds"
echo "Date: $(date)"
echo "============================================================================"
echo ""

# Define paths
SPECIES_TREE="${BASE_DIR}/${network}/tree_10_mil.nex"
OUTPUT_BASE="${BASE_DIR}/${network}/data/conf_${CONFIGURATION}_${TREE_HEIGHT}"

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

# ============================================================================
# SAMPLE SUBSTITUTION RATE (once for all replicates)
# ============================================================================

echo "Sampling substitution rate for all ${NUM_REPLICATES} replicates..."
SUB_RATE=$(sample_substitution_rate "$PICKLE_FILE")

if [ $? -ne 0 ] || [ -z "$SUB_RATE" ]; then
    echo "ERROR: Failed to sample substitution rate"
    exit 1
fi

echo "? Sampled substitution rate: ${SUB_RATE}"
echo "  (This rate will be used for all 5 replicates of ${network})"
echo ""

# Generate base seed for this network
BASE_SEED=$((SLURM_ARRAY_TASK_ID * 100000 + RANDOM))
echo "Base seed: ${BASE_SEED}"
echo ""

# ============================================================================
# ADAPTIVE BATCH STRATEGY
# ============================================================================

# Try different batch sizes (trees per replicate : number of runs)
BATCH_CONFIGS=(
    "1000:1"    # Try 1000 trees in 1 run first
    "10:100"    # If that fails, do 10 trees ª 100 runs
)

success=0

for config in "${BATCH_CONFIGS[@]}"; do
    IFS=':' read -r trees_per_run num_runs <<< "$config"
    total_trees=$((trees_per_run * num_runs))
    
    echo "============================================================================"
    echo "ATTEMPTING CONFIGURATION:"
    echo "  ${trees_per_run} trees ª ${num_runs} runs ª ${NUM_REPLICATES} replicates"
    echo "  = ${total_trees} trees per replicate"
    echo "============================================================================"
    echo ""
    
    # Clean output directory if it exists from previous attempt
    if [ -d "$OUTPUT_BASE" ]; then
        echo "Cleaning previous attempt..."
        rm -rf "$OUTPUT_BASE"
    fi
    mkdir -p "$OUTPUT_BASE"
    
    all_runs_succeeded=1
    
    # Loop through replicates (5 replicates)
    for replicate in $(seq 1 $NUM_REPLICATES); do
        echo "----------------------------------------------------------------------------"
        echo "REPLICATE ${replicate}/${NUM_REPLICATES}"
        echo "----------------------------------------------------------------------------"
        
        replicate_success=1
        
        # Loop through batches (for each replicate)
        for run in $(seq 1 $num_runs); do
            
            # Generate unique seed for this replicate and batch
            SEED=$((BASE_SEED + replicate * 10000 + run * 100))
            
            # Determine output directory
            if [ $num_runs -eq 1 ]; then
                # Single run - output directly to replicate folder
                RUN_OUTPUT="${OUTPUT_BASE}/replicate_${replicate}"
                if [ $run -eq 1 ]; then
                    echo "Mode: Single batch"
                    echo "Output: ${RUN_OUTPUT}"
                fi
            else
                # Multiple runs - use batch subdirectories
                RUN_OUTPUT="${OUTPUT_BASE}/replicate_${replicate}/batch_${run}"
                if [ $run -eq 1 ]; then
                    echo "Mode: Multiple batches (${num_runs} batches per replicate)"
                fi
                echo "  Batch ${run}/${num_runs}"
                echo "  Output: ${RUN_OUTPUT}"
            fi
            
            echo "  Seed: ${SEED}"
            echo "  Trees: ${trees_per_run}"
            echo "  Starting SimPhy (timeout: ${TIMEOUT_SECONDS}s)..."
            
            start_time=$(date +%s)
            
            # Run SimPhy for this batch
            run_simphy_single_replicate "$RUN_OUTPUT" "$trees_per_run" "$SUB_RATE" "$SEED" "$SPECIES_TREE" "$TIMEOUT_SECONDS"
            
            run_exit_code=$?
            end_time=$(date +%s)
            duration=$((end_time - start_time))
            
            if [ $run_exit_code -ne 0 ]; then
                echo "  ? Batch ${run} FAILED (duration: ${duration}s)"
                replicate_success=0
                all_runs_succeeded=0
                break
            else
                echo "  ? Batch ${run} completed (duration: ${duration}s)"
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
    
    # Check if all replicates succeeded
    if [ $all_runs_succeeded -eq 1 ]; then
        echo ""
        echo "============================================================================"
        echo "??? SUCCESS ???"
        echo "============================================================================"
        echo "Configuration: ${trees_per_run} trees ª ${num_runs} runs ª ${NUM_REPLICATES} replicates"
        echo "Total trees per replicate: ${total_trees}"
        echo "============================================================================"
        echo ""
        
        # Save detailed configuration
        cat > "${OUTPUT_BASE}/simulation_config.txt" << EOF
================================================================================
Simulation Configuration Summary
================================================================================

Network Information:
--------------------
Network Name: ${network}
Array Task ID: ${SLURM_ARRAY_TASK_ID}
Date: $(date)

Simulation Parameters:
----------------------
Configuration: ${CONFIGURATION}
Tree Height: ${TREE_HEIGHT}
Effective Population Size (Ne): ${NE}
Number of Replicates: ${NUM_REPLICATES}
Generation Time: 1 year

Substitution Rate:
------------------
Sampled Rate: ${SUB_RATE} (substitutions per site per generation)
Note: This rate is shared across all ${NUM_REPLICATES} replicates

Batch Configuration:
--------------------
Trees per run: ${trees_per_run}
Number of runs per replicate: ${num_runs}
Total trees per replicate: ${total_trees}

Seeds Used:
-----------
Base seed: ${BASE_SEED}
Replicate 1: Seeds ${BASE_SEED}+10000+100, +10000+200, ... (+10000+$((num_runs*100)))
Replicate 2: Seeds ${BASE_SEED}+20000+100, +20000+200, ... (+20000+$((num_runs*100)))
...
Replicate 5: Seeds ${BASE_SEED}+50000+100, +50000+200, ... (+50000+$((num_runs*100)))

Evolutionary Processes:
-----------------------
Gene Duplication Rate: 0 (ILS only)
Gene Loss Rate: 0 (ILS only)
Horizontal Transfer Rate: 0 (ILS only)
Gene Conversion Rate: 0 (ILS only)

Output Structure:
-----------------
${OUTPUT_BASE}/
  replicate_1/
EOF

        if [ $num_runs -eq 1 ]; then
            cat >> "${OUTPUT_BASE}/simulation_config.txt" << EOF
    1/           (${total_trees} gene trees)
  replicate_2/
    1/           (${total_trees} gene trees)
  ...
  replicate_5/
    1/           (${total_trees} gene trees)
EOF
        else
            cat >> "${OUTPUT_BASE}/simulation_config.txt" << EOF
    batch_1/
      1/         (${trees_per_run} gene trees)
    batch_2/
      1/         (${trees_per_run} gene trees)
    ...
    batch_${num_runs}/
      1/         (${trees_per_run} gene trees)
  replicate_2/
    batch_1/ through batch_${num_runs}/
  ...
  replicate_5/
    batch_1/ through batch_${num_runs}/

Total per replicate: ${num_runs} batches ª ${trees_per_run} trees = ${total_trees} trees
EOF
        fi

        cat >> "${OUTPUT_BASE}/simulation_config.txt" << EOF

================================================================================
EOF
        
        success=1
        break
    else
        echo ""
        echo "============================================================================"
        echo "??? CONFIGURATION FAILED ???"
        echo "============================================================================"
        echo "Configuration: ${trees_per_run} trees ª ${num_runs} runs"
        
        if [ "$config" == "1000:1" ]; then
            echo ""
            echo "? Falling back to smaller batches (10 trees ª 100 runs)"
            echo "  This will take longer but should work for complex networks"
        fi
        echo "============================================================================"
        echo ""
    fi
done

# ============================================================================
# FINAL SUMMARY
# ============================================================================

echo ""
echo "============================================================================"
echo "FINAL SUMMARY - ${network}"
echo "============================================================================"

if [ $success -eq 1 ]; then
    echo "Status: ? SUCCESS"
    echo ""
    echo "Final Configuration:"
    echo "  ${trees_per_run} trees ª ${num_runs} runs ª ${NUM_REPLICATES} replicates"
    echo "  Total: ${total_trees} trees per replicate"
    echo ""
    echo "Output directory: ${OUTPUT_BASE}"
    echo ""
    
    # Verify structure
    echo "Verifying output structure..."
    echo ""
    
    for rep in $(seq 1 $NUM_REPLICATES); do
        if [ $num_runs -eq 1 ]; then
            rep_dir="${OUTPUT_BASE}/replicate_${rep}/1"
        else
            rep_dir="${OUTPUT_BASE}/replicate_${rep}/batch_1/1"
        fi
        
        if [ -d "$rep_dir" ]; then
            echo "  ? Replicate ${rep}: Found"
        else
            echo "  ? Replicate ${rep}: MISSING!"
        fi
    done
    
    echo ""
    echo "Configuration saved to: ${OUTPUT_BASE}/simulation_config.txt"
else
    echo "Status: ? FAILURE"
    echo ""
    echo "ERROR: All batch configurations failed for ${network}"
    echo ""
    echo "Attempted configurations:"
    echo "  1. 1000 trees ª 1 run ª 5 replicates"
    echo "  2. 10 trees ª 100 runs ª 5 replicates"
fi

echo "============================================================================"
echo "Job completed: $(date)"
echo "============================================================================"

# Exit with appropriate code
if [ $success -eq 1 ]; then
    exit 0
else
    exit 1
fi
```

## Key Changes:

1. **Changed function**: `run_simphy_single_replicate` now uses `-rs 1` (not `-rs 5`)
2. **Added replicate loop**: Calls SimPhy 5 separate times (once per replicate)
3. **Better seed management**: Each replicate gets unique seeds
4. **Correct verification**: Now expects folder structure per replicate

## Output Structure Will Be:
```
conf_ils_low_10M/
??? replicate_1/
?   ??? 1/              # SimPhy output (1000 trees)
??? replicate_2/
?   ??? 1/              # SimPhy output (1000 trees)
??? replicate_3/
?   ??? 1/
??? replicate_4/
?   ??? 1/
??? replicate_5/
?   ??? 1/
??? simulation_config.txt
```

Or if it needs small batches:
```
conf_ils_low_10M/
??? replicate_1/
?   ??? batch_1/1/      # 10 trees
?   ??? batch_2/1/      # 10 trees
?   ??? ...
??? replicate_2/
?   ??? batch_1/1/
?   ??? ...
??? ...