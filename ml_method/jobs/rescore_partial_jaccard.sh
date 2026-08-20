#!/bin/bash
#SBATCH --job-name=rescore_partial
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --output=rescore_partial_%j.out

# Quick look at the partial_match=True reticulation jaccards on three representative configs,
# to see whether isolating placement (no penalty for unmatched reticulations) changes the
# picture. Re-scores competitors (run_full_summary --force-recompute) and the GNN reconstructions
# (score_reconstructions) with the new scoring code, then prints the per-method ret_leaf/ret_sisters
# means via head_to_head. Compare against docs/generated_tables_mu.md (standard jaccards).
# mu_distance is unaffected by the change; only the two jaccards move.
#
# NOTE: this re-scores ONLY these 3 configs, so the other 12 configs' CSVs still hold the standard
# jaccards -- do NOT run the full make_results_tables against a mix. Read the printed summaries here.

set -o pipefail   # not -u: gene2net conda activation references an unset LD_LIBRARY_PATH
CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate gene2net

ROOT=/groups/itay_mayrose/tomulanovski/gene2net
ML=$ROOT/ml_method
CONFIGS="conf_ils_low_10M conf_dup_loss_high_10M_ne1M conf_dup_loss_medium_10M_ne1M_fix050"

echo "================ competitor re-score (partial jaccards) ================"
cd "$ROOT"
for CFG in $CONFIGS; do
  echo ">>> run_full_summary $CFG"
  python simulations/scripts/run_full_summary.py "$CFG" --force-recompute \
      || echo "  !! run_full_summary FAILED: $CFG"
done

echo "================ GNN re-score + aggregate + head-to-head ================"
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
    echo "================ head-to-head: $CFG / $S ================"
    python scripts/head_to_head.py --config "$CFG" \
        --gnn-scores "output/reconstruct_final/final/$CFG/agg/$S" \
        --out "output/reconstruct_final/final/$CFG/h2h_partial_$S.csv" \
        || echo "  !! head_to_head FAILED: $CFG $S"
  done
done

echo "================ DONE ================"
