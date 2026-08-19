#!/bin/bash
#SBATCH --job-name=resummary_mu
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=resummary_mu_%j.out

# Re-score the COMPETITOR comparisons under the new headline metric (normalized
# mu-distance), for all nine configs. The mu metric lives in the shared
# compare_reticulations.pairwise_compare, so this re-run picks it up on the
# competitor side, exactly as score_reconstructions does on the GNN side.
#
# --force-recompute is MANDATORY: the comparison cache is keyed on gt/inf file
# hashes only, NOT on the metric-code version, so a plain re-run would reload
# pre-mu cached results that have no mu_distance column. --force-recompute
# ignores the cache and recomputes every pair, writing a fresh comparisons_raw.csv.
#
# Runs in the gene2net env (the comparison core / ReticulateTree; final_project
# errors with a GLIBCXX ImportError). mu-distance is far cheaper than the folded
# or MUL-tree GED it replaces, so the fix075 OOM/timeout risk is much reduced;
# 64G is kept as headroom. Per-config failures are logged, not fatal.

set -o pipefail   # not -u: the gene2net conda activation references an unset LD_LIBRARY_PATH

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate gene2net
# run_full_summary resolves simulations/summary_config.yaml relative to the repo ROOT.
cd /groups/itay_mayrose/tomulanovski/gene2net

CONFIGS="conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
conf_dup_loss_low_10M_ne1M conf_dup_loss_medium_10M_ne1M conf_dup_loss_high_10M_ne1M \
conf_dup_loss_medium_10M_ne1M_fix025 conf_dup_loss_medium_10M_ne1M_fix050 conf_dup_loss_medium_10M_ne1M_fix075"

for CFG in $CONFIGS; do
  echo "================ re-score competitors (mu): $CFG  ($(date +%H:%M:%S)) ================"
  python simulations/scripts/run_full_summary.py "$CFG" --force-recompute \
      || echo "  !! run_full_summary FAILED: $CFG"
done

echo "================ DONE  ($(date +%H:%M:%S)) ================"
