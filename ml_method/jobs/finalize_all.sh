#!/bin/bash
#SBATCH --job-name=finalize_all
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=finalize_all_%j.out

# Full finalization sweep: benchmark (both decodes, kernel ploidy, 5 replicates) then
# score, aggregate over replicates, and head-to-head, across all 9 configs. Runs the
# benchmark in the final_project env and the scoring/aggregation in the gene2net env,
# switching conda envs mid-script. Per-command failures are tolerated (logged, not fatal)
# so one bad config or a missing competitor comparison does not abort the whole sweep.

set -uo pipefail

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
cd /groups/itay_mayrose/tomulanovski/gene2net/ml_method

CONFIGS="conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
conf_dup_loss_low_10M_ne1M conf_dup_loss_medium_10M_ne1M conf_dup_loss_high_10M_ne1M \
conf_dup_loss_medium_10M_ne1M_fix025 conf_dup_loss_medium_10M_ne1M_fix050 conf_dup_loss_medium_10M_ne1M_fix075"

echo "================ STEP 1: benchmark (final_project) ================"
conda activate final_project
for CFG in $CONFIGS; do
  for R in 1 2 3 4 5; do
    echo ">>> benchmark $CFG rep$R  ($(date +%H:%M:%S))"
    python scripts/benchmark_networks.py --model-dir output/reconstruct_final --config "$CFG" \
        --model-config configs/reconstruct_final.yaml --replicate "$R" \
        --strategies bound_driven,detect_only --threshold 0.5 --copy-bound kernel \
        --out-base "output/reconstruct_final/final/$CFG/rep$R" \
        || echo "  !! benchmark FAILED: $CFG rep$R"
  done
done

echo "================ STEP 2: score (gene2net) ================"
conda activate gene2net
for CFG in $CONFIGS; do
  for R in 1 2 3 4 5; do
    for S in bound_driven detect_only; do
      python scripts/score_reconstructions.py \
          --recon-dir "output/reconstruct_final/final/$CFG/rep$R/$S" --workers 8 \
          || echo "  !! score FAILED: $CFG rep$R $S"
    done
  done
done

echo "================ STEP 3-4: aggregate + head-to-head ================"
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
