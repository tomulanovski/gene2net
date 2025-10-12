#!/bin/bash
# ==============================================================================
# SLURM Error Log Checker
# ==============================================================================
# This script checks if SLURM job error files are empty (successful jobs) or
# contain error messages (failed jobs).
#
# USAGE:
#   ./check_error_logs.sh                    # Check all error files automatically
#   ./check_error_logs.sh [pattern]          # Check specific pattern
#
# EXAMPLES:
#
#   1. Check all jobs automatically (finds all patterns):
#      ./check_error_logs.sh
#
#   2. Check only processing jobs:
#      ./check_error_logs.sh "process_ILS_low_dup_low_job_*.err"
#
#   3. Check only SimPhy jobs:
#      ./check_error_logs.sh "ILS_low_dup_low_job_*.err"
#
#   4. Check copies smoothing jobs:
#      ./check_error_logs.sh "copies_smoothing_ILS_low_dup_low_job_*.err"
#
#   5. Check specific job number (e.g., job 5):
#      ./check_error_logs.sh "*_job_5.err"
#
# OUTPUT:
#   - Green text: Empty error files (successful jobs)
#   - Red text: Non-empty error files with preview of error content
#   - Summary: Count of empty vs non-empty files
#
# NOTE: Pattern argument should NOT include the full path, just the filename pattern
# ==============================================================================

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to check error files
check_errors() {
    local pattern=$1
    local job_name=$2
    
    echo "=========================================="
    echo -e "${BLUE}Checking: $job_name${NC}"
    echo "Pattern: $pattern"
    echo "=========================================="
    
    # Find all matching error files
    error_files=($(ls $pattern 2>/dev/null))
    
    if [ ${#error_files[@]} -eq 0 ]; then
        echo -e "${YELLOW}WARNING: No error files found matching pattern${NC}"
        echo ""
        return
    fi
    
    local total_files=${#error_files[@]}
    local empty_files=0
    local non_empty_files=0
    
    echo "Found $total_files error files"
    echo ""
    
    # Check each file
    for file in "${error_files[@]}"; do
        if [ ! -s "$file" ]; then
            # File is empty
            empty_files=$((empty_files + 1))
        else
            # File has content
            non_empty_files=$((non_empty_files + 1))
            echo -e "${RED}NON-EMPTY:${NC} $file ($(wc -c < "$file") bytes)"
            echo "  First few lines:"
            head -n 5 "$file" | sed 's/^/    /'
            echo ""
        fi
    done
    
    # Summary
    echo "----------------------------------------"
    echo -e "${GREEN}Empty files: $empty_files / $total_files${NC}"
    if [ $non_empty_files -gt 0 ]; then
        echo -e "${RED}Non-empty files: $non_empty_files / $total_files${NC}"
    fi
    echo ""
}

# Function to extract job type pattern from filename
# Example: "process_ILS_low_dup_low_job_1.err" -> "process_ILS_low_dup_low_job"
get_job_pattern() {
    local filename=$(basename "$1")
    # Remove _NUMBER.err from the end to get the base pattern
    echo "$filename" | sed 's/_[0-9]*\.err$//'
}

# Function to get a nice display name from pattern
get_display_name() {
    local pattern=$1
    # Remove _job suffix and replace underscores with spaces, capitalize
    echo "$pattern" | sed 's/_job$//' | sed 's/_/ /g' | awk '{for(i=1;i<=NF;i++) $i=toupper(substr($i,1,1)) tolower(substr($i,2));}1'
}

# Main execution
LOG_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"

# Check if specific pattern provided as argument
if [ $# -gt 0 ]; then
    # User provided specific pattern
    full_pattern="${LOG_DIR}/$1"
    check_errors "$full_pattern" "Custom pattern: $1"
else
    # Automatically find all unique job patterns
    echo "=========================================="
    echo "Auto-detecting all error file patterns..."
    echo "=========================================="
    echo ""
    
    # Find all .err files and extract unique patterns
    declare -A patterns
    
    for err_file in "${LOG_DIR}"/*.err; do
        if [ -f "$err_file" ]; then
            base_pattern=$(get_job_pattern "$err_file")
            patterns["$base_pattern"]=1
        fi
    done
    
    # Check if any patterns found
    if [ ${#patterns[@]} -eq 0 ]; then
        echo -e "${YELLOW}No error files found in ${LOG_DIR}${NC}"
        exit 0
    fi
    
    echo -e "${GREEN}Found ${#patterns[@]} different job types${NC}"
    echo ""
    
    # Check each unique pattern
    for pattern in "${!patterns[@]}"; do
        display_name=$(get_display_name "$pattern")
        check_errors "${LOG_DIR}/${pattern}_*.err" "$display_name"
    done
fi

echo "=========================================="
echo "Check complete!"
echo "=========================================="