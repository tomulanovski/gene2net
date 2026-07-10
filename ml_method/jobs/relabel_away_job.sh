#!/bin/bash
# Re-relabel all rooted training configs with the AWAY-parent target (fixes the
# ~55% partner==home labelling bug). One-time preprocessing; ~6 min/config.
# Submit with: sbatch jobs/relabel_away_job.sh
#SBATCH --job-name=relabel_away
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=itaym-pool
#SBATCH --account=itaym-users_v2
#SBATCH --qos=owner
#SBATCH --output=logs/relabel_away_%j.out
#SBATCH --error=logs/relabel_away_%j.err

set -euo pipefail

CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate final_project

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
cd "$BASE_DIR"
mkdir -p logs

for cfg in ils_low ils_medium ils_high \
           dup_loss_low dup_loss_medium dup_loss_high \
           dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
    echo "============ $cfg ============"
    python scripts/relabel_from_metadata.py --training-subdir "training_rooted/$cfg" --away-parent
done

echo "ALL DONE"
