#!/bin/bash
#SBATCH --job-name=val_oracle_mu
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=val_oracle_mu_%j.out

# Recompute the diagnostic oracle under the mu-distance, so "backbone is the wall" can be
# re-examined under the new headline metric. Three points on the held-out validation split:
#   true_bb   : ground-truth metadata events on the TRUE backbone   -> pure build-faithfulness
#               floor (expected ~0, isomorphic).
#   astral_bb : ground-truth metadata events on the ASTRAL backbone -> ceiling that isolates
#               ASTRAL backbone error from event-prediction error.
#   model     : the model's own events on the ASTRAL backbone (bound_driven) -> the actual model.
# If model ~ astral_bb under mu, event prediction is saturated and the backbone is the wall under
# mu too. If model >> astral_bb, event prediction still has room and the edit-era conclusion does
# not carry over.
#
# STEP 1 (final_project) reconstructs each point; STEP 2 (gene2net) scores them under mu.

set -o pipefail   # not -u: the gene2net conda activation references an unset LD_LIBRARY_PATH
CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
cd /groups/itay_mayrose/tomulanovski/gene2net/ml_method

# The SAME six training dirs, in the SAME order used to train reconstruct_final, so the val split
# is reproduced identically (see submit_finalize.sh / score_val_reconstruction.py).
DATA_DIR=""
for c in ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
  DATA_DIR="$DATA_DIR data/mul_trees_2k/training_rooted/$c"
done

COMMON="--data-dir $DATA_DIR --model-dir output/reconstruct_final \
        --model-config configs/reconstruct_final.yaml --away-labels"

echo "================ STEP 1: reconstruct oracle points (final_project) ================"
conda activate final_project
python scripts/score_val_reconstruction.py $COMMON --events metadata --backbone true \
    --out-base output/oracle_mu/true_bb   || echo "  !! true_bb FAILED"
python scripts/score_val_reconstruction.py $COMMON --events metadata --backbone sample \
    --out-base output/oracle_mu/astral_bb || echo "  !! astral_bb FAILED"
python scripts/score_val_reconstruction.py $COMMON --events model --backbone sample \
    --strategy bound_driven --out-base output/oracle_mu/model || echo "  !! model FAILED"

echo "================ STEP 2: score under mu (gene2net) ================"
conda activate gene2net
for P in true_bb astral_bb model; do
  echo ">>> score $P"
  python scripts/score_reconstructions.py --recon-dir "output/oracle_mu/$P/bound_driven" --workers 8 \
      || echo "  !! score FAILED: $P"
done

echo "================ DONE ================"
