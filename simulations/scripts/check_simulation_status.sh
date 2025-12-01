#!/bin/bash
# check_simulation_status.sh - Check status of SimPhy simulations
#
# USAGE:
#   Default (checks conf_ils_low_10M):
#     ./check_simulation_status.sh
#
#
#   Or set and run:
#     export CHECK_CONFIG="conf_ils_medium_10M"
#     ./check_simulation_status.sh
#
# CONFIGURATION
# ============================================================================
SIMPHY_CONFIG="${CHECK_CONFIG:-conf_ils_low_10M}"

# ============================================================================
# PATHS
# ============================================================================
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations"
LOG_DIR="${BASE_DIR}/logs"
SIM_DIR="${BASE_DIR}/simulations"

networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

echo "================================================================================"
echo "SimPhy Simulation Status Summary"
echo "================================================================================"
echo "Looking for Configuration: ${SIMPHY_CONFIG}"
echo ""

success_count=0
fail_count=0
pending_count=0

success_1000=0
success_10x100=0

declare -a successful_networks
declare -a failed_networks
declare -a pending_networks

for i in "${!networks[@]}"; do
    network="${networks[$i]}"
    array_id=$((i + 1))
    
    # ------------------------------------------------------------------------
    # INTELLIGENT LOG SELECTION
    # ------------------------------------------------------------------------
    # 1. Define the specific log name (Best Practice)
    #    Expected: simphy_conf_ils_low_10M_1.out
    specific_log="${LOG_DIR}/simphy_${SIMPHY_CONFIG}_${array_id}.out"
    
    # 2. Define the generic log name (Old/Default SLURM behavior)
    #    Expected: simphy_simphy_1.out
    generic_log="${LOG_DIR}/simphy_simphy_${array_id}.out"
    
    log_file=""

    if [ -f "$specific_log" ]; then
        # Priority 1: Found the specifically named file. Use it.
        log_file="$specific_log"
    elif [ -f "$generic_log" ]; then
        # Priority 2: specific file missing, but generic exists.
        # MUST VALIDATE: Does this generic file actually belong to the config we want?
        # We grep for the config name in the file header.
        if grep -q "Configuration: ${SIMPHY_CONFIG}" "$generic_log"; then
            log_file="$generic_log"
        fi
    fi
    
    # ------------------------------------------------------------------------
    # STATUS CHECK
    # ------------------------------------------------------------------------
    
    printf "%-30s" "${network}:"
    
    # If no valid log file was found
    if [ -z "$log_file" ]; then
        echo " [PENDING] Log not found for ${SIMPHY_CONFIG}"
        pending_networks+=("$network")
        ((pending_count++))
        continue
    fi
    
    # Check for success (look for "SUCCESS" text)
    if grep -q "SUCCESS" "$log_file" 2>/dev/null; then
        # Check which batch configuration worked
        # We use 'strings' to avoid binary issues if file encoding is weird
        config_line=$(strings "$log_file" 2>/dev/null | grep "Trees per run:" | head -1)
        
        if echo "$config_line" | grep -q "1000"; then
            echo " [SUCCESS] 1000×1 (fast)"
            ((success_1000++))
        elif echo "$config_line" | grep -q "10"; then
            echo " [SUCCESS] 10×100 (slow)"
            ((success_10x100++))
        else
            echo " [SUCCESS] (config details unclear)"
        fi
        successful_networks+=("$network")
        ((success_count++))
        
    # Check for failure
    elif grep -q "FAILURE" "$log_file" 2>/dev/null; then
        echo " [FAILED] All batches failed"
        failed_networks+=("$network")
        ((fail_count++))
        
    # Check if still running
    elif grep -q "Starting SimPhy" "$log_file" 2>/dev/null; then
        echo " [RUNNING] In progress"
        pending_networks+=("$network")
        ((pending_count++))
    else
        echo " [UNKNOWN] Status unclear"
        pending_networks+=("$network")
        ((pending_count++))
    fi
done

echo ""
echo "================================================================================"
echo "SUMMARY"
echo "================================================================================"
echo "Configuration: ${SIMPHY_CONFIG}"
echo "Total networks:        21"
echo "Successful:            ${success_count}"
echo "  - 1000×1 (fast):     ${success_1000}"
echo "  - 10×100 (slow):     ${success_10x100}"
echo "Failed:                ${fail_count}"
echo "Pending/Running:       ${pending_count}"
echo ""

if [ ${#successful_networks[@]} -gt 0 ]; then
    echo "Successful networks:"
    for net in "${successful_networks[@]}"; do
        echo "  ✓ $net"
    done
    echo ""
fi

if [ ${#failed_networks[@]} -gt 0 ]; then
    echo "Failed networks:"
    for net in "${failed_networks[@]}"; do
        echo "  ✗ $net"
    done
    echo ""
fi

# Detailed breakdown
echo "================================================================================"
echo "CONFIGURATION BREAKDOWN (Sampled Rates)"
echo "================================================================================"
echo ""

for i in "${!networks[@]}"; do
    network="${networks[$i]}"
    array_id=$((i + 1))
    
    # Re-logic to find the file for the breakdown section
    specific_log="${LOG_DIR}/simphy_${SIMPHY_CONFIG}_${array_id}.out"
    generic_log="${LOG_DIR}/simphy_simphy_${array_id}.out"
    log_file=""

    if [ -f "$specific_log" ]; then
        log_file="$specific_log"
    elif [ -f "$generic_log" ] && grep -q "Configuration: ${SIMPHY_CONFIG}" "$generic_log"; then
        log_file="$generic_log"
    fi
    
    if [ -n "$log_file" ] && grep -q "SUCCESS" "$log_file"; then
        sub_rate=$(strings "$log_file" 2>/dev/null | grep "Sampled substitution rate:" | head -1 | awk '{print $NF}')
        
        printf "%-25s: " "$network"
        echo "Rate: $sub_rate"
    fi
done

echo ""
echo "================================================================================"