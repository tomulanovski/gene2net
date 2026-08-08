#!/bin/bash
#SBATCH --job-name=g2n_partner_perm
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/partner_perm_%j.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/partner_perm_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner

# Partner-accuracy permutation importance on the shipped model. CPU only (an eval),
# but memory-hungry because it holds the validation samples plus their pairwise
# features. Submit with:  sbatch jobs/partner_perm_job.sh
set -euo pipefail
source /groups/itay_mayrose/tomulanovski/miniconda3/etc/profile.d/conda.sh
conda activate final_project

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
cd "$BASE_DIR"

DATA_DIR=""
for c in ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
    DATA_DIR="$DATA_DIR ${BASE_DIR}/data/mul_trees_2k/training_rooted/$c"
done

python scripts/partner_permutation_importance.py \
    --data-dir ${DATA_DIR} \
    --model-dir "${BASE_DIR}/output/reconstruct_op_away_clade" \
    --config "${BASE_DIR}/configs/reconstruct_oneparter.yaml" \
    --max-samples 400
