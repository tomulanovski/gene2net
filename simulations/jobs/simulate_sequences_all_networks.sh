#!/bin/bash

# Usage: ./simulate_sequences_all_networks.sh <configuration>
#
# Example: ./simulate_sequences_all_networks.sh ils_low_10M
#
# This script submits sequence simulation jobs for all 21 networks using the specified configuration.
# Each job will sample GTR+Gamma parameters from empirical distributions for realistic sequence evolution.

if [ -z "$1" ]; then
    echo "ERROR: Configuration not specified"
    echo ""
    echo "Usage: $0 <configuration>"
    echo ""
    echo "Examples:"
    echo "  $0 ils_low_10M"
    echo "  $0 ils_high_10M"
    echo "  $0 ils_med_10M"
    echo ""
    exit 1
fi

CONFIGURATION=$1

echo "============================================================================"
echo "Submitting Sequence Simulation Jobs for All Networks"
echo "============================================================================"
echo "Configuration: ${CONFIGURATION}"
echo "Date: $(date)"
echo "============================================================================"
echo ""
echo "This will submit jobs for 21 networks, each processing 1000 gene trees"
echo "using empirically-sampled GTR+Gamma parameters."
echo ""

networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

job_count=0

for network in "${networks[@]}"; do
    echo "Submitting: ${network}"

    job_name="alisim_${network}_${CONFIGURATION}"

    sbatch --job-name="${job_name}" \
        "/groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs/simulate_sequences_1_dataset.sh" \
        "${network}" \
        "${CONFIGURATION}"

    if [ $? -eq 0 ]; then
        ((job_count++))
        echo "  ✓ Job submitted: ${job_name}"
    else
        echo "  ✗ Failed to submit job for ${network}"
    fi

    # Small delay to avoid overwhelming scheduler
    sleep 1
done

echo ""
echo "============================================================================"
echo "Summary"
echo "============================================================================"
echo "Configuration: ${CONFIGURATION}"
echo "Jobs submitted: ${job_count}/${#networks[@]}"
echo ""
echo "To monitor jobs:"
echo "  squeue -u \$USER | grep alisim"
echo ""
echo "To check job status:"
echo "  sacct -X --format=JobName,State,ExitCode | grep alisim"
echo "============================================================================"
