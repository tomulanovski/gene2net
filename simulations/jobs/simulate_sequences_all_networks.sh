#!/bin/bash

# usage: "/groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs/simulate_sequences_all_networks.sh"

networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

for network in "${networks[@]}"; do
    echo "Submitting jobs for: ${network}"
    sbatch --job-name=alisim_${network} "/groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs/simulate_sequences_1_dataset.sh" ${network}
    sleep 1  # Small delay to avoid overwhelming scheduler
done