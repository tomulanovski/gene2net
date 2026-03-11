#!/bin/bash
#SBATCH --job-name=jones_xml
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/%x.e%j
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate alloppnet

NEX_DIR="${1:?ERROR: Please provide path to jones_alloppnet directory as first argument}"

Rscript /groups/itay_mayrose/tomulanovski/gene2net/scripts/jones_script.r "$NEX_DIR"
