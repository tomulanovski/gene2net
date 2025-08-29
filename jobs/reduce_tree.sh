#!/bin/bash
#SBATCH --job-name=tree_pruning
#SBATCH --output=pruning_%j.out
#SBATCH --error=pruning_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=power-general
#SBATCH --account=power-general-users

# Set your parameters here - easy to modify for different runs
INPUT_TREE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Hori_2014/gene_trees/output_tree.nwk"
OUTPUT_TREE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Hori_2014/gene_trees/output_tree_like_paper.nwk"

# Define taxa in the correct format for the script (species:k,species:k)
TAXA_STRING="Dryopteris_varia:1,Dryopteris_bissetiana:2,Dryopteris_sacrosancta:3,Dryopteris_saxifragivaria:2,Dryopteris_pacifica_a:2,Dryopteris_pacifica_b:3,Dryopteris_pacifica_c:3,Dryopteris_protobissetiana:1,Dryopteris_saxifraga:1,Dryopteris_kobayashii:3,Dryopteris_chinensis:1,Dryopteris_gymnophylla:1,Dryopteris_sordidipes:1,Dryopteris_caudipinna:1,Dryopteris_expansa:1,Dryopteris_sieboldii:1"

# Conda setup
CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Activate your environment
conda activate gene2net

# Show the active environment (helpful for debugging)
echo "Using conda environment: $CONDA_PREFIX"

# Run the pruning script with correct parameters
python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/reduce_tree_most_diverse_k.py" -i "$INPUT_TREE" -taxa "$TAXA_STRING" -o "$OUTPUT_TREE" -v

# Print completion message
echo "Job completed"