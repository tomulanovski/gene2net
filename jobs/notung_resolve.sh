#!/bin/bash
#SBATCH --job-name=notung_rooting
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/papers/Wu_2015/%x.o%j
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/papers/Wu_2015/%x.e%j
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=power-general
#SBATCH --account=power-general-users


CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3

# Source conda
source "$CONDA_PATH/etc/profile.d/conda.sh"

# Check if conda is properly sourced
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda command not found. Check your conda installation path."
    exit 1
fi


# Activate your specific environment
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# Paths to files
BATCH_FILE="/groups/itay_mayrose/tomulanovski/gene2net/papers/Wu_2015/notung_batch.txt"
NOTUNG_JAR="/groups/itay_mayrose/tomulanovski/Notung-2.9.1.5/Notung-2.9.1.5.jar"

# Create output directory
OUTPUT_DIR="/groups/itay_mayrose/tomulanovski/gene2net/papers/Wu_2015/rooted_trees"
mkdir -p ${OUTPUT_DIR}

# Run Notung directly on the gene tree file to root it
echo "Running Notung to root trees..."
java -Xmx25G -jar ${NOTUNG_JAR} -b ${BATCH_FILE} --absfilenames --root --resolve --speciestag postfix --nolosses --treeoutput newick --outputdir ${OUTPUT_DIR} --log

echo "Processing complete!"