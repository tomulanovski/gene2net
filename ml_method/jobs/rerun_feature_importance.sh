#!/bin/bash
#SBATCH --job-name=feat_importance
#SBATCH --partition=itaym-pool
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --output=feat_importance_%j.out

# Recompute permutation importance on the SHIPPED model (reconstruct_final), because the
# numbers currently in the feature-importance chapter section were run on an earlier model
# (reconstruct_op_away_clade), i.e. likely stale. Permutation importance is forward passes on
# the validation split, no training, so this is cheap. Two scripts:
#   permutation_importance.py          -> detection F1 drops (node + edge features)
#   partner_permutation_importance.py  -> allopolyploid partner-accuracy drops (node + edge + pairwise)
# The scripts live under _archive/experiments/scripts. Their sibling import
# (partner imports permutation_importance) works because Python puts the script's own dir on the
# path; gene2net_gnn is found via PYTHONPATH=the ml_method root.

set -o pipefail
CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate final_project
cd /groups/itay_mayrose/tomulanovski/gene2net/ml_method
export PYTHONPATH="$PWD"

# The six configs the shipped model was trained/selected on (rooted, away-labelled).
DATA_DIR=""
for c in ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
  DATA_DIR="$DATA_DIR data/mul_trees_2k/training_rooted/$c"
done

MODEL=output/reconstruct_final
CFG=configs/reconstruct_final.yaml

echo "================ detection F1 permutation importance ($(date +%H:%M:%S)) ================"
python _archive/experiments/scripts/permutation_importance.py \
    --data-dir $DATA_DIR --model-dir "$MODEL" --config "$CFG" \
    || echo "  !! detection permutation importance FAILED"

echo "================ partner-accuracy permutation importance ($(date +%H:%M:%S)) ================"
python _archive/experiments/scripts/partner_permutation_importance.py \
    --data-dir $DATA_DIR --model-dir "$MODEL" --config "$CFG" \
    || echo "  !! partner permutation importance FAILED"

echo "================ DONE ($(date +%H:%M:%S)) ================"
