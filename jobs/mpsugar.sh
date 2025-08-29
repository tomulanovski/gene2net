#!/bin/bash
#SBATCH --job-name=mpsugar_run
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/MPSUGAR_TEST/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/MPSUGAR_TEST/%x.e%j
#SBATCH --time=12:00:00
#SBATCH --mem=100g
#SBATCH --cpus-per-task=10
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

# job path : "/groups/itay_mayrose/tomulanovski/gene2net/jobs/mpsugar.sh"

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Check if conda is properly sourced
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda command not found. Check your conda installation path."
    exit 1
fi

# Activate environment
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# Verify the environment was activated
echo "Current conda environment: $CONDA_PREFIX"

# Force Python to flush output immediately
export PYTHONUNBUFFERED=1

# Define file paths
PYTHON_SCRIPT="/groups/itay_mayrose/tomulanovski/gene2net/scripts/MPSUGAR.py"

echo "Starting MP-SUGAR analysis..."
echo "Python script: $PYTHON_SCRIPT"
echo "Start time: $(date)"

# Run the MP-SUGAR analysis
python -u $PYTHON_SCRIPT

echo "End time: $(date)"
echo "MP-SUGAR analysis bash job ended"