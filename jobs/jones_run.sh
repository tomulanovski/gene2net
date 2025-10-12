#!/bin/bash
#SBATCH --job-name=beast_pgic
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/Sessa_2012b/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/Sessa_2012b/%x.e%j
#SBATCH --time=120:00:00
#SBATCH --mem=32G
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate alloppnet

cd /groups/itay_mayrose/tomulanovski/gene2net/papers/Sessa_2012b/

INPUT_FILE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Sessa_2012b/jones_analysis/sessa_2012_b_10mil.XML"

beast $INPUT_FILE