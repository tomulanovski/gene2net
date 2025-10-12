#!/bin/bash
#SBATCH --job-name=padre_sims
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/padre_run_rep1_ILS_low_dup_low_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/padre_run_rep1_ILS_low_dup_low_job_%a.err
#SBATCH --time=100:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# PADRE jar path
PADRE_JAR="/groups/itay_mayrose/tomulanovski/padre/padre-cli.jar"

# Array of networks (must match the original array)
networks=(
  "Bendiksby_2011" "Koenen_2020" "Brysting_2007" "Lawrence_2016"
  "Diaz-Perez_2018" "Wisecaver_2023" "Ding_2023" "Liang_2019" "Popp_2005" "Wu_2015"
  "Liu_2023" "Ren_2024" "Marcussen_2011" "Marcussen_2012" "Sessa_2012b" "Zhao_2021"
  "Hori_2014" "Marcussen_2015" "Shahrestani_2015" "Morales-Briones_2021" "Soza_2014"
)

# Get the network for this array task
network=${networks[$SLURM_ARRAY_TASK_ID-1]}

echo "=========================================="
echo "Running PADRE for: ${network} - ILS_low_dup_low"
echo "Processing replicate: 1"
echo "=========================================="

# Define paths for replicate 1
INPUT_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/padre_input"
OUTPUT_DIR="${BASE_DIR}/${network}/results/ILS_low_dup_low/padre"

# Input file for replicate 1
INPUT_TREES="${INPUT_DIR}/1_trees.tre"

# Check if input file exists
if [ ! -f "$INPUT_TREES" ]; then
    echo "ERROR: Gene trees file not found: $INPUT_TREES"
    exit 1
fi

# Check if PADRE jar exists
if [ ! -f "$PADRE_JAR" ]; then
    echo "ERROR: PADRE jar file not found: $PADRE_JAR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Change to output directory so PADRE writes results there
cd "$OUTPUT_DIR"

echo "Input trees: $INPUT_TREES"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Run PADRE
java -Djava.awt.headless=true -Xmx4g -jar "$PADRE_JAR" -i "$INPUT_TREES"

# Check if PADRE succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "PADRE completed successfully for ${network}"
    echo "Results saved to: $OUTPUT_DIR"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "ERROR: PADRE failed for ${network}"
    echo "=========================================="
    exit 1
fi