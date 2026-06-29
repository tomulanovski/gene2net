#!/bin/bash
# ============================================================================
# CHECK_DUP_LOSS_RATES.SH
# ============================================================================
# Audits the dup/loss rate (and Ne) actually used by SimPhy for each
# conf_dup_loss_* config, reading the simulation_config.txt that simphy_reusable.sh
# writes per network/data/config.
#
# Intended (per MY):  low = 0.001   medium = 0.01   high = 0.1
# Flags any config whose recorded rate != the intended rate for its tier.
#
# Usage:
#   ./check_dup_loss_rates.sh
#   ./check_dup_loss_rates.sh --all-networks   # check every network, not just one
# ============================================================================

set -uo pipefail

BASE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

CONFIGS=(
    conf_dup_loss_low_10M       conf_dup_loss_medium_10M       conf_dup_loss_high_10M
    conf_dup_loss_low_10M_ne1M  conf_dup_loss_medium_10M_ne1M  conf_dup_loss_high_10M_ne1M
    conf_dup_loss_low_10M_ne2M  conf_dup_loss_medium_10M_ne2M  conf_dup_loss_high_10M_ne2M
)

networks=(
  Bendiksby_2011 Koenen_2020 Brysting_2007 Lawrence_2016
  Diaz-Perez_2018 Wisecaver_2023 Ding_2023 Liang_2019 Popp_2005 Wu_2015
  Liu_2023 Ren_2024 Marcussen_2011 Marcussen_2012 Sessa_2012b Zhao_2021
  Hori_2014 Marcussen_2015 Shahrestani_2015 Morales-Briones_2021 Soza_2014
)

ALL_NETWORKS=false
[ "${1:-}" = "--all-networks" ] && ALL_NETWORKS=true

# intended rate (per MY) for a tier
intended_rate() {
    case "$1" in
        *low*)    echo "0.001" ;;
        *medium*) echo "0.01"  ;;
        *high*)   echo "0.1"   ;;
        *)        echo "?"     ;;
    esac
}

# pull a numeric field (per MY) from a simulation_config.txt
get_field() {  # $1 = file, $2 = label regex
    grep -E "$2" "$1" 2>/dev/null | head -1 | grep -oE '[0-9.eE+-]+ per MY' | grep -oE '^[0-9.eE+-]+'
}
get_ne() { grep -E "Population Size" "$1" 2>/dev/null | head -1 | grep -oE '[0-9]+'; }

printf "%-32s %-10s %-10s %-10s %-10s %s\n" "CONFIG" "Ne" "DUP/MY" "LOSS/MY" "INTENDED" "STATUS"
printf '%.0s-' {1..90}; echo

for C in "${CONFIGS[@]}"; do
    want=$(intended_rate "$C")

    if [ "$ALL_NETWORKS" = true ]; then
        net_list=("${networks[@]}")
    else
        # use the first network that has a config file
        net_list=()
        for N in "${networks[@]}"; do
            [ -f "$BASE/$N/data/$C/simulation_config.txt" ] && { net_list=("$N"); break; }
        done
    fi

    if [ ${#net_list[@]} -eq 0 ]; then
        printf "%-32s %-10s %-10s %-10s %-10s %s\n" "$C" "-" "-" "-" "$want" "NO CONFIG FILE"
        continue
    fi

    for N in "${net_list[@]}"; do
        f="$BASE/$N/data/$C/simulation_config.txt"
        [ -f "$f" ] || { [ "$ALL_NETWORKS" = true ] && printf "%-32s %-10s %-10s %-10s %-10s %s\n" "$C" "?" "?" "?" "$want" "MISSING ($N)"; continue; }
        ne=$(get_ne "$f")
        dup=$(get_field "$f" "Duplication rate")
        loss=$(get_field "$f" "Loss rate")
        status="OK"
        [ "$dup" != "$want" ] && status="MISMATCH (want $want)"
        label="$C"
        [ "$ALL_NETWORKS" = true ] && label="$C [$N]"
        printf "%-32s %-10s %-10s %-10s %-10s %s\n" "$label" "${ne:-?}" "${dup:-?}" "${loss:-?}" "$want" "$status"
    done
done
