#!/bin/bash

# ============================================================================
# NETWORK INVESTIGATION SCRIPT
# ============================================================================
# Investigates why a specific network failed or has missing results
#
# Usage:
#   ./investigate_network.sh <network_name> <config> [method]
#
# Examples:
#   ./investigate_network.sh Diaz-Perez_2018 conf_ils_low_10M
#   ./investigate_network.sh Diaz-Perez_2018 conf_ils_low_10M polyphest
#   ./investigate_network.sh Marcussen_2011 conf_ils_medium_10M grampa
# ============================================================================

set -uo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# ============================================================================
# CONFIGURATION
# ============================================================================

LOGS_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Network list (0-indexed)
NETWORKS=(
  "Bendiksby_2011"       # 0
  "Koenen_2020"          # 1
  "Brysting_2007"        # 2
  "Lawrence_2016"        # 3
  "Diaz-Perez_2018"      # 4
  "Wisecaver_2023"       # 5
  "Ding_2023"            # 6
  "Liang_2019"           # 7
  "Popp_2005"            # 8
  "Wu_2015"              # 9
  "Liu_2023"             # 10
  "Ren_2024"             # 11
  "Marcussen_2011"       # 12
  "Marcussen_2012"       # 13
  "Sessa_2012b"          # 14
  "Zhao_2021"            # 15
  "Hori_2014"            # 16
  "Marcussen_2015"       # 17
  "Shahrestani_2015"     # 18
  "Morales-Briones_2021" # 19
  "Soza_2014"            # 20
)

ALL_METHODS=(grampa polyphest padre mpsugar)

# ============================================================================
# USAGE
# ============================================================================

usage() {
  cat << EOF
Usage: $0 <network_name> <config> [method]

Arguments:
  network_name    Network to investigate (e.g., Diaz-Perez_2018)
  config          Configuration name (e.g., conf_ils_low_10M)
  method          Optional: specific method to check (grampa, polyphest, padre, mpsugar)

Examples:
  $0 Diaz-Perez_2018 conf_ils_low_10M
  $0 Diaz-Perez_2018 conf_ils_low_10M polyphest
  $0 Marcussen_2011 conf_ils_medium_10M grampa

Available networks:
$(for i in "${!NETWORKS[@]}"; do echo "  - ${NETWORKS[$i]} (index $i, tasks $((i*5+1))-$((i*5+5)))"; done)

EOF
  exit 1
}

# ============================================================================
# ARGUMENT PARSING
# ============================================================================

if [ $# -lt 2 ]; then
  usage
fi

NETWORK="$1"
CONFIG="$2"
FILTER_METHOD="${3:-}"

# Find network index
NETWORK_INDEX=-1
for i in "${!NETWORKS[@]}"; do
  if [ "${NETWORKS[$i]}" = "$NETWORK" ]; then
    NETWORK_INDEX=$i
    break
  fi
done

if [ $NETWORK_INDEX -eq -1 ]; then
  echo "ERROR: Network not found: $NETWORK"
  echo ""
  echo "Available networks:"
  for net in "${NETWORKS[@]}"; do
    echo "  - $net"
  done
  exit 1
fi

# Calculate task IDs
TASK_START=$((NETWORK_INDEX * 5 + 1))
TASK_END=$((NETWORK_INDEX * 5 + 5))
TASKS=($(seq $TASK_START $TASK_END))

# Determine which methods to check
if [ -n "$FILTER_METHOD" ]; then
  METHODS=("$FILTER_METHOD")
else
  METHODS=("${ALL_METHODS[@]}")
fi

# ============================================================================
# HEADER
# ============================================================================

echo "========================================================================"
echo "INVESTIGATING: $NETWORK"
echo "========================================================================"
echo "Network index: $NETWORK_INDEX"
echo "Tasks: ${TASKS[*]} (replicates 1-5)"
echo "Configuration: $CONFIG"
echo "Checking methods: ${METHODS[*]}"
echo "========================================================================"
echo ""

# ============================================================================
# CHECK LOGS
# ============================================================================

cd "$LOGS_DIR" || { echo "ERROR: Cannot access $LOGS_DIR"; exit 1; }

for method in "${METHODS[@]}"; do
  echo "--- $method ---"
  echo ""

  found_any=false
  success_count=0
  failed_count=0
  missing_count=0

  for task in "${TASKS[@]}"; do
    rep=$((task - TASK_START + 1))

    # Find run log
    run_log=$(ls run_${method}_${CONFIG}_*_${task}.out 2>/dev/null | head -1)

    if [ -z "$run_log" ]; then
      echo -e "  Rep $rep (task $task): ${RED}NO LOG FOUND${NC} (never ran?)"
      ((missing_count++))
      continue
    fi

    found_any=true

    # Check status
    if grep -q "COMPLETED SUCCESSFULLY\|completed successfully" "$run_log" 2>/dev/null; then
      echo -e "  Rep $rep (task $task): ${GREEN}SUCCESS${NC}"
      ((success_count++))
    elif grep -q "FAILED" "$run_log" 2>/dev/null; then
      echo -e "  Rep $rep (task $task): ${RED}FAILED${NC}"
      ((failed_count++))

      # Show why it failed
      echo "    Last 10 lines of log:"
      tail -10 "$run_log" | sed 's/^/      /'
      echo ""

      # Check error log
      err_log="${run_log%.out}.err"
      if [ -s "$err_log" ]; then
        echo "    Error log ($err_log):"
        cat "$err_log" | sed 's/^/      /'
        echo ""
      fi
    else
      echo -e "  Rep $rep (task $task): ${YELLOW}UNKNOWN${NC} (incomplete or running?)"

      # Show last few lines
      echo "    Last 5 lines:"
      tail -5 "$run_log" | sed 's/^/      /'
      echo ""
    fi
  done

  if [ "$found_any" = false ]; then
    echo "  No logs found - job was never submitted!"
    echo ""
    echo "  To submit:"
    echo "    ./submit_${method}_pipeline.sh $CONFIG"
  fi

  echo ""
  echo "  Summary: $success_count success, $failed_count failed, $missing_count missing"
  echo ""
done

# ============================================================================
# CHECK OUTPUT FILES
# ============================================================================

echo "========================================================================"
echo "OUTPUT FILES"
echo "========================================================================"
echo ""

# Expand methods to include polyphest variants
CHECK_METHODS=()
for method in "${METHODS[@]}"; do
  if [ "$method" = "polyphest" ]; then
    CHECK_METHODS+=(polyphest_p50 polyphest_p70 polyphest_p90)
  else
    CHECK_METHODS+=("$method")
  fi
done

for method in "${CHECK_METHODS[@]}"; do
  echo "$method:"

  for rep in {1..5}; do
    # Determine result file based on method
    if [ "$method" = "grampa" ]; then
      result_file="${BASE_DIR}/${NETWORK}/results/${CONFIG}/${method}/replicate_${rep}/grampa_result.tre"
    elif [[ "$method" == polyphest* ]]; then
      result_file="${BASE_DIR}/${NETWORK}/results/${CONFIG}/${method}/replicate_${rep}/polyphest_result.tre"
    elif [ "$method" = "padre" ]; then
      result_file="${BASE_DIR}/${NETWORK}/results/${CONFIG}/${method}/replicate_${rep}/padre_result.tre"
    elif [ "$method" = "mpsugar" ]; then
      result_file="${BASE_DIR}/${NETWORK}/results/${CONFIG}/${method}/replicate_${rep}/mpsugar_result.tre"
    fi

    if [ -f "$result_file" ]; then
      size=$(stat -c%s "$result_file" 2>/dev/null || stat -f%z "$result_file" 2>/dev/null)
      echo -e "  Rep $rep: ${GREEN}EXISTS${NC} ($size bytes)"

      # Show preview if small
      if [ "$size" -lt 500 ]; then
        echo "    Content:"
        cat "$result_file" | sed 's/^/      /'
      fi
    else
      echo -e "  Rep $rep: ${RED}NOT FOUND${NC}"
      echo "    Expected: $result_file"
    fi
  done
  echo ""
done

# ============================================================================
# CHECK PREP STEP (for methods that need it)
# ============================================================================

if [ "$FILTER_METHOD" = "polyphest" ] || [ -z "$FILTER_METHOD" ]; then
  echo "========================================================================"
  echo "PREP STEP CHECK (Polyphest)"
  echo "========================================================================"
  echo ""

  cd "$LOGS_DIR"

  for task in "${TASKS[@]}"; do
    rep=$((task - TASK_START + 1))

    prep_log=$(ls prep_polyphest_${CONFIG}_*_${task}.out 2>/dev/null | head -1)

    if [ -z "$prep_log" ]; then
      echo -e "Rep $rep (task $task): ${RED}NO PREP LOG${NC}"
    elif grep -q "completed successfully\|COMPLETED" "$prep_log" 2>/dev/null; then
      echo -e "Rep $rep (task $task): ${GREEN}PREP OK${NC}"
    else
      echo -e "Rep $rep (task $task): ${RED}PREP FAILED${NC}"
      echo "  Last 5 lines:"
      tail -5 "$prep_log" | sed 's/^/    /'
    fi
  done
  echo ""
fi

# ============================================================================
# RECOMMENDATIONS
# ============================================================================

echo "========================================================================"
echo "RECOMMENDATIONS"
echo "========================================================================"
echo ""

# Count overall status
total_expected=$((${#METHODS[@]} * 5))
if [ "$FILTER_METHOD" = "polyphest" ]; then
  total_expected=$((3 * 5))  # p50, p70, p90
fi

echo "Next steps:"
echo ""

# Check if any logs are missing
if grep -q "NO LOG FOUND" <(echo "$output"); then
  echo "1. Some jobs never ran. Submit them:"
  for method in "${METHODS[@]}"; do
    echo "   ./submit_${method}_pipeline.sh $CONFIG"
  done
  echo ""
fi

# Check if any failed
if grep -q "FAILED" <(echo "$output"); then
  echo "2. Some jobs failed. Check the error messages above."
  echo "   Common issues:"
  echo "   - Missing input files (check prep step)"
  echo "   - Timeout (increase --time in job script)"
  echo "   - Out of memory (increase --mem in job script)"
  echo "   - Software error (check error logs)"
  echo ""
fi

# Check if outputs missing
echo "3. After jobs complete, run post-processing:"
echo "   python simulations/scripts/postprocess_results.py $CONFIG"
echo ""

echo "4. Verify with pipeline status checker:"
echo "   python simulations/scripts/check_pipeline_status.py $CONFIG --step run --verbose"
echo ""

echo "========================================================================"
