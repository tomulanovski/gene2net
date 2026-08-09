#!/bin/bash
#SBATCH --job-name=g2n_optuna
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/optuna_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/optuna_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --array=0-7
#SBATCH --partition=gpu-rotemhsh-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner
# 8 array tasks x 13 trials = 104 trials against ONE shared SQLite study.
# Prereq (run ONCE before sbatch, to create the DB and avoid a create race):
#   python scripts/optuna_reconstruct.py --data-dir <one dir> --config configs/reconstruct_oneparter.yaml \
#     --study-name g2n_reconstruct --storage sqlite:///$PWD/optuna/g2n_reconstruct.db \
#     --n-trials 0 --out-root output/optuna/g2n_reconstruct
set -euo pipefail
source /groups/itay_mayrose/tomulanovski/miniconda3/etc/profile.d/conda.sh
conda activate final_project

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
export PYTHONPATH="${PYTHONPATH:-}:${BASE_DIR}"
export PYTHONUNBUFFERED=1
cd "$BASE_DIR"
mkdir -p optuna logs

DATA_DIR=""
for c in ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
    DATA_DIR="$DATA_DIR ${BASE_DIR}/data/mul_trees_2k/training_rooted/$c"
done

STUDY=g2n_reconstruct
STORAGE="sqlite:///${BASE_DIR}/optuna/${STUDY}.db"

echo "Array task ${SLURM_ARRAY_TASK_ID:-?} | study $STUDY | storage $STORAGE"

python scripts/optuna_reconstruct.py \
    --data-dir ${DATA_DIR} \
    --config "${BASE_DIR}/configs/reconstruct_oneparter.yaml" \
    --study-name "$STUDY" --storage "$STORAGE" \
    --n-trials 13 --out-root "${BASE_DIR}/output/optuna/${STUDY}"
