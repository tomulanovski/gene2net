#!/bin/bash
#SBATCH --job-name=finalize_score
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --output=finalize_score_%j.out

# Scoring half of the finalization sweep, for when the benchmark (STEP 1) already ran.
# Scores every config/replicate/decode, aggregates over replicates, and runs head-to-head,
# all in the gene2net env. No `set -u`: the gene2net conda activation script references an
# unset LD_LIBRARY_PATH, which would otherwise abort under -u. Per-command failures are logged.

set -o pipefail

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate gene2net
cd /groups/itay_mayrose/tomulanovski/gene2net/ml_method

CONFIGS="conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
conf_dup_loss_low_10M_ne1M conf_dup_loss_medium_10M_ne1M conf_dup_loss_high_10M_ne1M \
conf_dup_loss_medium_10M_ne1M_fix025 conf_dup_loss_medium_10M_ne1M_fix050 conf_dup_loss_medium_10M_ne1M_fix075"

echo "================ score ($(date +%H:%M:%S)) ================"
for CFG in $CONFIGS; do
  for R in 1 2 3 4 5; do
    for S in bound_driven detect_only; do
      python scripts/score_reconstructions.py \
          --recon-dir "output/reconstruct_final/final/$CFG/rep$R/$S" --workers 8 \
          || echo "  !! score FAILED: $CFG rep$R $S"
    done
  done
done

echo "================ aggregate + head-to-head ================"
for CFG in $CONFIGS; do
  for S in bound_driven detect_only; do
    python scripts/aggregate_replicates.py \
        --rep-dirs "output/reconstruct_final/final/$CFG/rep1/$S" \
                   "output/reconstruct_final/final/$CFG/rep2/$S" \
                   "output/reconstruct_final/final/$CFG/rep3/$S" \
                   "output/reconstruct_final/final/$CFG/rep4/$S" \
                   "output/reconstruct_final/final/$CFG/rep5/$S" \
        --out "output/reconstruct_final/final/$CFG/agg/$S/reconstruction_scores.csv" \
        || echo "  !! aggregate FAILED: $CFG $S"
    python scripts/head_to_head.py --config "$CFG" \
        --gnn-scores "output/reconstruct_final/final/$CFG/agg/$S" \
        --out "output/reconstruct_final/final/$CFG/h2h_$S.csv" \
        || echo "  !! head_to_head FAILED (competitor comparison may be missing): $CFG $S"
  done
done

echo "================ DONE  ($(date +%H:%M:%S)) ================"
