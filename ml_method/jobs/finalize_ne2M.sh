#!/bin/bash
#SBATCH --job-name=finalize_ne2M
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=finalize_ne2M_%j.out

# Full pipeline for the Ne=2M dup/loss arm, now that its rates are fixed
# (check_dup_loss_rates.sh all OK, 2026-08-19) and its gene trees + competitor
# inferences are current. Self-contained so it can run alongside resummary_mu.sh
# (which covers the other 9 configs) without conflict. Three steps, two envs:
#   STEP 1 (final_project): GNN benchmark on the ne2M configs (both decodes, 5 reps, kernel ploidy)
#   STEP 2 (gene2net):      competitor comparisons re-scored under mu (--force-recompute, MANDATORY:
#                           the compute_comparisons cache is keyed on file hashes, not metric code)
#   STEP 3 (gene2net):      GNN score + aggregate over reps + head-to-head
#
# detect_only threshold is 0.5 here to match the rest of the benchmark. The threshold
# is applied at edge-selection time (not scoring time), so if the val sweep picks a
# different mu-calibrated default, the detect_only benchmark must be RE-RUN at that
# threshold for ALL configs, ne2M included -- not just re-scored.

set -o pipefail   # not -u: the gene2net conda activation references an unset LD_LIBRARY_PATH

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"

ROOT=/groups/itay_mayrose/tomulanovski/gene2net
ML=$ROOT/ml_method

CONFIGS="conf_dup_loss_low_10M_ne2M conf_dup_loss_medium_10M_ne2M conf_dup_loss_high_10M_ne2M"

echo "================ STEP 1: GNN benchmark (final_project) ================"
conda activate final_project
cd "$ML"
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

echo "================ STEP 2: competitor re-score under mu (gene2net) ================"
conda activate gene2net
cd "$ROOT"
for CFG in $CONFIGS; do
  echo ">>> run_full_summary $CFG  ($(date +%H:%M:%S))"
  python simulations/scripts/run_full_summary.py "$CFG" --force-recompute \
      || echo "  !! run_full_summary FAILED: $CFG"
done

echo "================ STEP 3: GNN score + aggregate + head-to-head (gene2net) ================"
cd "$ML"
for CFG in $CONFIGS; do
  for R in 1 2 3 4 5; do
    for S in bound_driven detect_only; do
      python scripts/score_reconstructions.py \
          --recon-dir "output/reconstruct_final/final/$CFG/rep$R/$S" --workers 8 \
          || echo "  !! score FAILED: $CFG rep$R $S"
    done
  done
done
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
