#!/bin/bash
#SBATCH --job-name=mpsugar_sims
#SBATCH --array=1-21
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/mpsugar_run_rep1_ILS_low_dup_low_job_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/mpsugar_run_rep1_ILS_low_dup_low_job_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=50g
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

# Set conda path
CONDA_PATH="/groups/itay_mayrose/tomulanovski/miniconda3"

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Activate gene2net environment (assuming MP-SUGAR/PhyNetPy is installed here)
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# Verify the environment was activated
echo "Current conda environment: $CONDA_PREFIX"

# Base directory
BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# MP-SUGAR script path
MPSUGAR_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/mpsugr.py"

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
echo "Running MP-SUGAR for: ${network} - ILS_low_dup_low"
echo "Processing replicate: 1"
echo "=========================================="

# Define paths for replicate 1
INPUT_DIR="${BASE_DIR}/${network}/processed/ILS_low_dup_low/mpallopp_input"
OUTPUT_DIR="${BASE_DIR}/${network}/results/ILS_low_dup_low/mpallopp"

# Input files for replicate 1
GENE_TREES="${INPUT_DIR}/1_MPSUGAR_trees.nex"
TAXON_MAP="${INPUT_DIR}/1_taxon_map.json"

# Check if input files exist
if [ ! -f "$GENE_TREES" ]; then
    echo "ERROR: Gene trees file not found: $GENE_TREES"
    exit 1
fi

if [ ! -f "$TAXON_MAP" ]; then
    echo "ERROR: Taxon map file not found: $TAXON_MAP"
    exit 1
fi

# Check if MP-SUGAR script exists
if [ ! -f "$MPSUGAR_SCRIPT" ]; then
    echo "ERROR: MP-SUGAR script not found: $MPSUGAR_SCRIPT"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Output file
OUTPUT_FILE="${OUTPUT_DIR}/1_mpsugar_results.txt"

echo "Input gene trees: $GENE_TREES"
echo "Input taxon map: $TAXON_MAP"
echo "Output file: $OUTPUT_FILE"
echo ""

# Run MP-SUGAR
# Adjust --iterations and --chains as needed
python "$MPSUGAR_SCRIPT" \
    --trees "$GENE_TREES" \
    --map "$TAXON_MAP" \
    --output "$OUTPUT_FILE" \
    --iterations 500 \
    --chains 1

# Check if MP-SUGAR succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "MP-SUGAR completed successfully for ${network}"
    echo "Results saved to: $OUTPUT_FILE"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "ERROR: MP-SUGAR failed for ${network}"
    echo "=========================================="
    exit 1
fi