#!/bin/bash
#SBATCH --job-name=polyphest_run
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/%x.e%j
#SBATCH --time=3-00:00:00
#SBATCH --mem=50g
#SBATCH --cpus-per-task=32
#SBATCH --partition=itaym
#SBATCH --account=itaym-users


# job path : "/groups/itay_mayrose/tomulanovski/gene2net/jobs/polyphest_run.sh"


CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Check if conda is properly sourced
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda command not found. Check your conda installation path."
    exit 1
fi


# Activate your specific environment
conda activate polyphest || {
    echo "ERROR: Could not activate polyphest environment"
    exit 1
}

# Verify the environment was activated
echo "Current conda environment: $CONDA_PREFIX"


POLYPHEST="/groups/itay_mayrose/tomulanovski/Polyphest/Polyphest.py" 


GENE_TREES="/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/gene_trees/concatenated_trees.tre"
MULTISET="/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/multi_set_paper.txt"
OUTPUT_DIR="/groups/itay_mayrose/tomulanovski/gene2net/papers/Marcussen_2012/networks/polyphest"

python $POLYPHEST --gene_tree_file $GENE_TREES --consensus_multiset_file $MULTISET --filter_strategy percentile --output_dir $OUTPUT_DIR --percentile 95 --use_near_isomorphic "True" --isomorphic_threshold 0.2



