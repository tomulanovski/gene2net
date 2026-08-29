#!/bin/bash
#SBATCH --job-name=resummary_empirical_mu
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=resummary_empirical_mu_%j.out

# Re-score the EMPIRICAL pairwise comparisons under the headline metric
# (normalized mu-distance), replacing the approximate MUL-tree edit distance.
#
# The metric itself lives in the shared compare_reticulations.pairwise_compare,
# which scripts/compute_comparisons.py imports from simulations/scripts, so this
# re-run picks up mu_distance and mu_scored with no change to the engine. This
# is the empirical counterpart of ml_method/jobs/resummary_mu.sh, which does the
# same for the simulated configs.
#
# --force-recompute is MANDATORY: the comparison cache under
# analysis/real_data_summary/cache is keyed on the network file hashes only, NOT
# on the metric-code version, so a plain re-run would reload pre-mu cached
# results that carry an edit_distance_multree column and no mu_distance at all.
#
# Runs in the gene2net env (the comparison core / ReticulateTree / phylonetwork).
# mu-distance is far cheaper than the graph edit distance it replaces, so the
# 64G here is headroom rather than a requirement.

set -o pipefail   # not -u: the gene2net conda activation references an unset LD_LIBRARY_PATH

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# run_analysis.py resolves scripts/papers_config.yaml relative to the repo ROOT.
cd /groups/itay_mayrose/tomulanovski/gene2net || exit 1

echo "================ re-score empirical (mu)  ($(date +%H:%M:%S)) ================"
python scripts/run_analysis.py --force-recompute

echo "================ DONE  ($(date +%H:%M:%S)) ================"
echo
echo "Numbers for manuscript Section 3.3 land in analysis/real_data_summary/:"
echo "  mu_scoring_report.csv       how many pairs the mu-distance could score"
echo "  pairwise_summary.csv        mean/std/min/max per method pair and metric"
echo "  per_network_comparisons.csv per-dataset values, wide format"
echo "  method_rankings.csv         avg_mu_distance per method"
