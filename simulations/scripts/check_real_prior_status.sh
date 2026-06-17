#!/bin/bash

# ============================================================================
# CHECK_REAL_PRIOR_STATUS.SH
# ============================================================================
# Verifies the real-prior runs (run_polyphest_real.sh and
# run_grandma_split_prior_real.sh) completed for all 9 dup_loss configs.
#
# Checks the RAW method outputs (before postprocess_results.py):
#   Polyphest: polyphest_real_p<PCT>/replicate_N/polyphest_trees-polyphest.txt
#   GRANDMA:   grandma_split_prior_real/replicate_N/final_multree.tre
# and scans the job logs for failure markers.
#
# Usage:
#   ./check_real_prior_status.sh            # summary + count of missing
#   ./check_real_prior_status.sh --list     # also list every missing output
# ============================================================================

set -uo pipefail

LIST_MISSING=false
if [ "${1:-}" = "--list" ]; then
    LIST_MISSING=true
fi

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
LOG_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"

CONFIGS=(
    "conf_dup_loss_low_10M"
    "conf_dup_loss_medium_10M"
    "conf_dup_loss_high_10M"
    "conf_dup_loss_low_10M_ne1M"
    "conf_dup_loss_medium_10M_ne1M"
    "conf_dup_loss_high_10M_ne1M"
    "conf_dup_loss_low_10M_ne2M"
    "conf_dup_loss_medium_10M_ne2M"
    "conf_dup_loss_high_10M_ne2M"
)

PERCENTILES=(50 70 90)

networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

NUM_REPLICATES=5

# Returns 0 if file exists and is non-empty
exists_nonempty() { [ -s "$1" ]; }

# --------------------------------------------------------------------------
# POLYPHEST REAL
# --------------------------------------------------------------------------
echo "============================================================================"
echo "POLYPHEST (REAL PRIOR)  -  output: polyphest_trees-polyphest.txt"
echo "============================================================================"
poly_found=0; poly_total=0
declare -a poly_missing=()
for config in "${CONFIGS[@]}"; do
    for pct in "${PERCENTILES[@]}"; do
        c_found=0; c_total=0
        for network in "${networks[@]}"; do
            for rep in $(seq 1 $NUM_REPLICATES); do
                f="${BASE_DIR}/${network}/results/${config}/polyphest_real_p${pct}/replicate_${rep}/polyphest_trees-polyphest.txt"
                c_total=$((c_total+1))
                if exists_nonempty "$f"; then
                    c_found=$((c_found+1))
                else
                    poly_missing+=("$f")
                fi
            done
        done
        printf "  %-32s p%-3s %3d/%-3d\n" "$config" "$pct" "$c_found" "$c_total"
        poly_found=$((poly_found+c_found)); poly_total=$((poly_total+c_total))
    done
done
echo "  --------------------------------------------------------------------------"
printf "  POLYPHEST TOTAL: %d/%d\n" "$poly_found" "$poly_total"
echo ""

# --------------------------------------------------------------------------
# GRANDMA_SPLIT_PRIOR REAL
# --------------------------------------------------------------------------
echo "============================================================================"
echo "GRANDMA_SPLIT (REAL PRIOR)  -  output: final_multree.tre"
echo "============================================================================"
gr_found=0; gr_total=0
declare -a gr_missing=()
for config in "${CONFIGS[@]}"; do
    c_found=0; c_total=0
    for network in "${networks[@]}"; do
        for rep in $(seq 1 $NUM_REPLICATES); do
            f="${BASE_DIR}/${network}/results/${config}/grandma_split_prior_real/replicate_${rep}/final_multree.tre"
            c_total=$((c_total+1))
            if exists_nonempty "$f"; then
                c_found=$((c_found+1))
            else
                gr_missing+=("$f")
            fi
        done
    done
    printf "  %-32s %3d/%-3d\n" "$config" "$c_found" "$c_total"
    gr_found=$((gr_found+c_found)); gr_total=$((gr_total+c_total))
done
echo "  --------------------------------------------------------------------------"
printf "  GRANDMA TOTAL: %d/%d\n" "$gr_found" "$gr_total"
echo ""

# --------------------------------------------------------------------------
# LOG FAILURE SCAN
# --------------------------------------------------------------------------
echo "============================================================================"
echo "LOG FAILURE SCAN (markers written by the run scripts)"
echo "============================================================================"
poly_fail=$(grep -l "POLYPHEST (REAL PRIOR) FAILED" "${LOG_DIR}"/run_polyphest_real_*.out 2>/dev/null | wc -l)
gr_fail=$(grep -l "GRANDMA_SPLIT (REAL PRIOR) FAILED" "${LOG_DIR}"/run_grandma_split_prior_real_*.out 2>/dev/null | wc -l)
echo "  Polyphest logs with FAILED marker: ${poly_fail}"
echo "  GRANDMA   logs with FAILED marker: ${gr_fail}"
echo ""
echo "  (errors in .err files worth a look:)"
echo "    grep -l . ${LOG_DIR}/run_polyphest_real_*.err 2>/dev/null | head"
echo "    grep -l . ${LOG_DIR}/run_grandma_split_prior_real_*.err 2>/dev/null | head"
echo ""

# --------------------------------------------------------------------------
# MISSING LIST (optional)
# --------------------------------------------------------------------------
if [ "$LIST_MISSING" = true ]; then
    echo "============================================================================"
    echo "MISSING OUTPUTS"
    echo "============================================================================"
    if [ ${#poly_missing[@]} -gt 0 ]; then
        echo "Polyphest (${#poly_missing[@]}):"
        printf '  %s\n' "${poly_missing[@]}"
    fi
    if [ ${#gr_missing[@]} -gt 0 ]; then
        echo "GRANDMA (${#gr_missing[@]}):"
        printf '  %s\n' "${gr_missing[@]}"
    fi
fi

echo "Expected totals: Polyphest 2835 (9x3x21x5), GRANDMA 945 (9x21x5)."
