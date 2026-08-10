#!/bin/bash
#SBATCH --job-name=g2n_optuna
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/optuna_%A_%a.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/ml_method/logs/optuna_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --array=0-103%4
#SBATCH --partition=gpu-rotemhsh-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner
# JOB-PER-TRIAL: 104 array tasks, each runs ONE Optuna trial (~3h) against the ONE
# shared SQLite study, then exits and releases its GPU. Short jobs backfill onto the
# busy single-node GPU pool far better than long persistent workers, and %4 caps
# concurrency so we do not flood the queue. Total ~104 trials, up to 4 at a time.
# Trial output dirs use the study-assigned trial.number, so there is no collision
# between array tasks. Prereq (run ONCE first, creates the DB):
#   python scripts/optuna_reconstruct.py --config configs/reconstruct_oneparter.yaml \
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

STUDY=g2n_reconstruct_ext          # extended-space study (hidden up to 512, depth up to 5)
STORAGE="sqlite:///${BASE_DIR}/optuna/${STUDY}.db"

echo "Array task ${SLURM_ARRAY_TASK_ID:-?} | one trial | study $STUDY"

python scripts/optuna_reconstruct.py \
    --data-dir ${DATA_DIR} \
    --config "${BASE_DIR}/configs/reconstruct_oneparter.yaml" \
    --study-name "$STUDY" --storage "$STORAGE" \
    --n-trials 1 --out-root "${BASE_DIR}/output/optuna/${STUDY}"
