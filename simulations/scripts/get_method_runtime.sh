#!/bin/bash
# ==============================================================================
# SCRIPT: get_method_runtime.sh
# PURPOSE: Calculate average runtime for methods across all datasets
#
# USAGE:
#   ./get_method_runtime.sh <CONFIG> [METHOD]
#
# EXAMPLES:
#   ./get_method_runtime.sh conf_ils_low_10M              # All methods
#   ./get_method_runtime.sh conf_ils_low_10M polyphest    # Single method
#   ./get_method_runtime.sh conf_ils_low_10M --verbose    # All methods, verbose
#
# OUTPUT:
#   - Prints to terminal and saves to simulations/analysis/runtime_<CONFIG>.txt
#   - Average runtime in seconds and hours
#   - Statistics: how many jobs completed / total tasks
#   - Min/max runtimes
#   - Uses highest job ID for each task (most recent run)
# ==============================================================================

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Arguments
CONFIG=$1
METHOD=$2
VERBOSE=0

if [ -z "$CONFIG" ]; then
    echo "Usage: $0 <CONFIG> [METHOD]"
    echo ""
    echo "Examples:"
    echo "  $0 conf_ils_low_10M              # All methods"
    echo "  $0 conf_ils_low_10M polyphest    # Single method"
    echo "  $0 conf_ils_low_10M --verbose    # All methods, verbose"
    exit 1
fi

# Check for verbose flag
if [ "$METHOD" == "--verbose" ] || [ "$METHOD" == "-v" ]; then
    VERBOSE=1
    METHOD=""
fi
if [ "$3" == "--verbose" ] || [ "$3" == "-v" ]; then
    VERBOSE=1
fi

# Paths
LOG_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"
ANALYSIS_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis"
OUTPUT_FILE="${ANALYSIS_DIR}/runtime_${CONFIG}.txt"
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Networks (21 total)
NETWORKS=(
    "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
    "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019"
    "Popp_2005" "Wu_2015" "Liu_2023" "Ren_2024" "Marcussen_2011"
    "Marcussen_2012" "Sessa_2012b" "Zhao_2021" "Hori_2014"
    "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# All methods
ALL_METHODS=("grampa" "polyphest" "padre" "mpsugar")

# If method specified, use only that one
if [ -n "$METHOD" ]; then
    METHODS=("$METHOD")
else
    METHODS=("${ALL_METHODS[@]}")
fi

# Create analysis directory if it doesn't exist
mkdir -p "$ANALYSIS_DIR"

# Function to map task_id (1-105) to network and replicate
# task_id = (network_index * 5) + replicate, where network_index is 0-based
get_network_replicate() {
    local task_id=$1
    local network_idx=$(( (task_id - 1) / 5 ))
    local replicate=$(( ((task_id - 1) % 5) + 1 ))
    echo "${NETWORKS[$network_idx]} $replicate"
}

# Function to check if output file exists for a method/task
check_output_file_exists() {
    local method=$1
    local task_id=$2
    local network_replicate=$(get_network_replicate "$task_id")
    local network=$(echo "$network_replicate" | awk '{print $1}')
    local replicate=$(echo "$network_replicate" | awk '{print $2}')
    
    local output_file=""
    
    case "$method" in
        "padre")
            # PADRE writes to processed directory
            output_file="${BASE_DIR}/${network}/processed/${CONFIG}/padre_input/replicate_${replicate}/padre_trees-result.tre"
            ;;
        "grampa")
            output_file="${BASE_DIR}/${network}/results/${CONFIG}/grampa/replicate_${replicate}/grampa-scores.txt"
            ;;
        "polyphest")
            # Polyphest uses percentile directories - check common ones (50, 70, 90)
            for percentile in 50 70 90; do
                local polyphest_file="${BASE_DIR}/${network}/results/${CONFIG}/polyphest_p${percentile}/replicate_${replicate}/polyphest_trees-polyphest.txt"
                if [ -f "$polyphest_file" ]; then
                    echo "1"
                    return 0
                fi
            done
            echo "0"
            return 1
            ;;
        "mpsugar")
            output_file="${BASE_DIR}/${network}/results/${CONFIG}/mpsugar/replicate_${replicate}/mpsugar_results.txt"
            ;;
        *)
            echo "0"
            return 1
            ;;
    esac
    
    if [ -f "$output_file" ]; then
        echo "1"
        return 0
    else
        echo "0"
        return 1
    fi
}

# Function to analyze a single method
analyze_method() {
    local method=$1
    local pattern="run_${method}_${CONFIG}_*_*.out"
    
    # Total expected tasks (21 networks × 5 replicates = 105)
    local total_tasks=105
    
    # Arrays to store results
    declare -A task_durations
    declare -A task_job_ids
    
    # Find all matching log files and extract durations
    for log_file in ${LOG_DIR}/${pattern}; do
        # Check if file exists (in case no matches)
        if [ ! -f "$log_file" ]; then
            continue
        fi
        
        # Extract job_id and task_id from filename
        # Format: run_<method>_<config>_<job_id>_<task_id>.out
        local basename=$(basename "$log_file")
        
        # Remove the prefix and suffix to get job_id_task_id
        local temp="${basename#run_${method}_${CONFIG}_}"
        temp="${temp%.out}"
        
        # Split by underscore to get job_id and task_id
        local job_id=$(echo "$temp" | cut -d'_' -f1)
        local task_id=$(echo "$temp" | cut -d'_' -f2)
        
        # Check if job completed successfully (case-insensitive)
        if ! grep -qi "COMPLETED SUCCESSFULLY" "$log_file"; then
            continue  # Skip failed jobs
        fi
        
        # Extract duration from log file
        local duration=$(grep "Duration:" "$log_file" | grep -oP '\d+(?= seconds)' | head -1)
        
        # Filter out 0-second durations (suspicious/failed jobs)
        if [ -n "$duration" ] && [ "$duration" -gt 0 ]; then
            # Check if this is the highest job_id for this task
            if [ -z "${task_job_ids[$task_id]}" ] || [ "$job_id" -gt "${task_job_ids[$task_id]}" ]; then
                # Validate that output file exists before counting
                if [ "$(check_output_file_exists "$method" "$task_id")" == "1" ]; then
                    task_job_ids[$task_id]=$job_id
                    task_durations[$task_id]=$duration
                fi
            fi
        fi
    done
    
    # Calculate statistics
    local total_found=0
    local sum_duration=0
    local min_duration=""
    local max_duration=""
    
    for task_id in $(seq 1 $total_tasks); do
        if [ -n "${task_durations[$task_id]}" ]; then
            local duration=${task_durations[$task_id]}
            ((total_found++))
            sum_duration=$((sum_duration + duration))
            
            # Track min/max
            if [ -z "$min_duration" ] || [ "$duration" -lt "$min_duration" ]; then
                min_duration=$duration
            fi
            if [ -z "$max_duration" ] || [ "$duration" -gt "$max_duration" ]; then
                max_duration=$duration
            fi
        fi
    done
    
    # Return results as space-separated string
    echo "$total_found $sum_duration $min_duration $max_duration"
}

# Start output
{
    echo "============================================================================"
    echo "METHOD RUNTIME ANALYSIS"
    echo "============================================================================"
    echo "Configuration: ${CONFIG}"
    echo "Date: $(date)"
    echo "Log directory: ${LOG_DIR}"
    echo ""
    
    if [ ${#METHODS[@]} -eq 1 ]; then
        echo "Method: ${METHODS[0]}"
    else
        echo "Methods: All (${ALL_METHODS[@]})"
    fi
    echo ""
    echo "============================================================================"
    echo ""
    
    # Analyze each method
    for method in "${METHODS[@]}"; do
        echo "----------------------------------------------------------------------------"
        echo "METHOD: ${method^^}"
        echo "----------------------------------------------------------------------------"
        echo ""
        
        # Check if logs exist
        pattern="run_${method}_${CONFIG}_*_*.out"
        if ! ls ${LOG_DIR}/${pattern} 1> /dev/null 2>&1; then
            echo "⚠ No log files found for ${method}"
            echo "  Pattern: ${LOG_DIR}/${pattern}"
            echo ""
            continue
        fi
        
        # Analyze method
        result=$(analyze_method "$method")
        total_found=$(echo "$result" | awk '{print $1}')
        sum_duration=$(echo "$result" | awk '{print $2}')
        min_duration=$(echo "$result" | awk '{print $3}')
        max_duration=$(echo "$result" | awk '{print $4}')
        
        if [ "$total_found" -eq 0 ]; then
            echo "✗ No completed jobs found with duration information"
            echo ""
            echo "Possible reasons:"
            echo "  - Jobs are still running"
            echo "  - Jobs failed before completion"
            echo "  - Jobs hit time limit (check SLURM sacct for elapsed time)"
            echo ""
            continue
        fi
        
        # Calculate average
        avg_duration=$((sum_duration / total_found))
        avg_hours=$(echo "scale=2; $avg_duration / 3600" | bc)
        
        # Calculate coverage percentage
        coverage=$(echo "scale=1; ($total_found * 100) / 105" | bc)
        
        # Display statistics
        echo "✓ Completed jobs: ${total_found} / 105 (${coverage}%)"
        echo ""
        echo "Average runtime:"
        echo "  ${avg_duration} seconds"
        echo "  ${avg_hours} hours"
        echo ""
        echo "Min runtime: ${min_duration} seconds ($(echo "scale=2; $min_duration/3600" | bc) hours)"
        echo "Max runtime: ${max_duration} seconds ($(echo "scale=2; $max_duration/3600" | bc) hours)"
        echo ""
        
        # Suggest time limit (max * 2)
        suggested_limit=$((max_duration * 2))
        suggested_hours=$(echo "scale=0; ($suggested_limit / 3600) + 1" | bc)
        
        echo "RECOMMENDED TIME LIMIT:"
        echo "  ${suggested_hours}:00:00 (2x max runtime)"
        echo ""
        
        # Warn if coverage is low
        if [ $total_found -lt 52 ]; then
            echo "⚠ WARNING: Low coverage (<50%). Some jobs may be incomplete."
            echo ""
        fi
    done
    
    echo "============================================================================"
    echo "SUMMARY"
    echo "============================================================================"
    echo ""
    echo "Report saved to: ${OUTPUT_FILE}"
    echo ""
    echo "To analyze with verbose per-task breakdown:"
    echo "  $0 ${CONFIG} [METHOD] --verbose"
    echo ""
    echo "============================================================================"
    
} | tee "$OUTPUT_FILE"

echo ""
echo -e "${GREEN}✓ Analysis complete${NC}"
echo -e "${CYAN}Report saved to: ${OUTPUT_FILE}${NC}"
