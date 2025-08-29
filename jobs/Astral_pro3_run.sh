#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=power-general
#SBATCH --account=power-general-users
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/%x.e%j

# Define input and output files

INPUT_FILE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/gene_trees/concatenated_trees.tre"
OUTPUT_FILE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/gene_trees/species_astral_pro3_tre.tre"
LOG_FILE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/gene_trees/astral_pro3_run.log"

# Path to ASTRAL-3 executable
# Change this to match your installation path
ASTRAL_PATH="/groups/itay_mayrose/tomulanovski/ASTER-Linux/bin/astral-pro3"

# Run ASTRAL-4
echo "Starting ASTRAL pro analysis..."
echo "Input file: $INPUT_FILE"
echo "Output will be written to: $OUTPUT_FILE"
echo "Log will be written to: $LOG_FILE"

cd /groups/itay_mayrose/tomulanovski/gene2net/papers/Koenen_2020/gene_trees/

$ASTRAL_PATH -o $OUTPUT_FILE $INPUT_FILE 2>$LOG_FILE

# Check if the run was successful
if [ $? -eq 0 ]; then
    echo "ASTRAL-4 analysis completed successfully."
    echo "Species tree saved to: $OUTPUT_FILE"
else
    echo "Error: ASTRAL-4 analysis failed. Check the log file for details."
fi