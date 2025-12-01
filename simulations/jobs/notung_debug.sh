#!/bin/bash
#SBATCH --job-name=notung_debug
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/notung_debug.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/notung_debug.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Create output directory
OUTPUT_DIR="/groups/itay_mayrose/tomulanovski/gene2net/simulations/notung_debug"
mkdir -p "$OUTPUT_DIR"

# Paths
GENE_TREE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/notung_debug/gene_tree_1_no_dup_diaz_nhx.tre"
SPECIES_TREE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations/Diaz-Perez_2018/species_tree_1_million.tre"
NOTUNG_JAR="/groups/itay_mayrose/tomulanovski/Notung-2.9.1.5/Notung-2.9.1.5.jar"

# Conda setup
CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate gene2net

echo "Running Notung reconciliation for a single tree..."
export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"

java -Xmx25G -jar ${NOTUNG_JAR} \
    -g ${GENE_TREE} \
    -s ${SPECIES_TREE} \
    --reconcile \
    --treestats \
    --speciestag nhx \
    --outputdir ${OUTPUT_DIR} \
    --log

conda deactivate
