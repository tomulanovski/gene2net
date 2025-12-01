#!/bin/bash
#SBATCH --job-name=simphy_normalized
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/simphy_height_1_million_trees_100000_ne_low_dup_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/simphy_height_1_million_trees_100000_ne_low_dup_job_%a.err
#SBATCH --time=1:00:00
#SBATCH --mem=4g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Array of networks
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# Get the network for this array task
network=${networks[$SLURM_ARRAY_TASK_ID-1]}
echo "Running SimPhy for: ${network} - no dup (height 1 million)"

# Define paths - CHANGED TO USE NORMALIZED TREES
SPECIES_TREE="${BASE_DIR}/${network}/species_tree_1_million.nex"
OUTPUT_DIR="${BASE_DIR}/${network}/data/low_dup_ne_100000_height_1_million_trees"

# Check if species tree exists
if [ ! -f "$SPECIES_TREE" ]; then
    echo "ERROR: Species tree not found: $SPECIES_TREE"
    exit 1
fi

echo "Using height 1 million species tree: $SPECIES_TREE"

# Run SimPhy with specified duplication and loss rates (1e-9 per generation)
# effective population size (Ne) = 1,000,000
/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulators/simphy/SimPhy_1.0.2/bin/simphy_lnx64 -rs 1 \
    -rl f:1000 \
    -rg 1 \
    -SR "$SPECIES_TREE" \
    -sp f:100000 \
    -su f:0.0001 \
    -sg f:1 \
    -lb f:1e-9 \
    -ld f:1e-9 \
    -lt f:0 \
    -lg f:0 \
    -cs 12345 \
    -o "$OUTPUT_DIR" \
    -v 1 \
    -od 1 \
    -op 1

echo "SimPhy completed for ${network} - ILS_low_dup_low (height 1 million trees)"
