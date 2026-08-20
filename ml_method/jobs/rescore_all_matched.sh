#!/bin/bash
#SBATCH --job-name=rescore_all_matched
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=rescore_all_matched_%j.out

# Re-score ALL 15 configs with the current scoring code so the partial-match jaccards
# (ret_leaf_jaccard_matched / ret_sisters_jaccard_matched) are populated everywhere, on both
# the competitor and GNN sides. mu_distance is unaffected by the change -- only the jaccards
# and the new _matched columns move -- so this is a re-score, not re-inference. The *_result.tre
# postprocessing already exists for every config (finalize_ne2M / finalize_10M / pp_summary_ne2M),
# so run_full_summary finds its inputs.
#
# After this, python scripts/make_results_tables.py shows mu | num_rets | ret_leaf | ret_leaf_m |
# ret_sis | ret_sis_m | n for all 15 configs.

set -o pipefail   # not -u: the gene2net conda activation references an unset LD_LIBRARY_PATH
CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate gene2net

ROOT=/groups/itay_mayrose/tomulanovski/gene2net
ML=$ROOT/ml_method

CONFIGS="conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
conf_dup_loss_low_10M conf_dup_loss_medium_10M conf_dup_loss_high_10M \
conf_dup_loss_low_10M_ne1M conf_dup_loss_medium_10M_ne1M conf_dup_loss_high_10M_ne1M \
conf_dup_loss_low_10M_ne2M conf_dup_loss_medium_10M_ne2M conf_dup_loss_high_10M_ne2M \
conf_dup_loss_medium_10M_ne1M_fix025 conf_dup_loss_medium_10M_ne1M_fix050 conf_dup_loss_medium_10M_ne1M_fix075"

echo "================ competitor re-score, all 15 (mu + matched jaccards) ================"
cd "$ROOT"
for CFG in $CONFIGS; do
  echo ">>> run_full_summary $CFG  ($(date +%H:%M:%S))"
  python simulations/scripts/run_full_summary.py "$CFG" --force-recompute \
      || echo "  !! run_full_summary FAILED: $CFG"
done

echo "================ GNN re-score + aggregate, all 15 ================"
cd "$ML"
for CFG in $CONFIGS; do
  for R in 1 2 3 4 5; do
    for S in bound_driven detect_only; do
      python scripts/score_reconstructions.py \
          --recon-dir "output/reconstruct_final/final/$CFG/rep$R/$S" --workers 8 \
          || echo "  !! score FAILED: $CFG rep$R $S"
    done
  done
  for S in bound_driven detect_only; do
    python scripts/aggregate_replicates.py \
        --rep-dirs "output/reconstruct_final/final/$CFG/rep1/$S" \
                   "output/reconstruct_final/final/$CFG/rep2/$S" \
                   "output/reconstruct_final/final/$CFG/rep3/$S" \
                   "output/reconstruct_final/final/$CFG/rep4/$S" \
                   "output/reconstruct_final/final/$CFG/rep5/$S" \
        --out "output/reconstruct_final/final/$CFG/agg/$S/reconstruction_scores.csv" \
        || echo "  !! aggregate FAILED: $CFG $S"
  done
done

echo "================ DONE  ($(date +%H:%M:%S)) ================"
