#!/bin/bash
#SBATCH --job-name=pp_summary_ne2M
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=pp_summary_ne2M_%j.out

# Recover the ne2M medium/high competitor mu scores. The methods RAN (final_multree.tre,
# grandma.log, etc. present under replicate_N/) but postprocess_results.py was never run
# for these two configs, so the standardized <method>_result.tre files the inventory keys
# on (collect_results.py:78) are missing -> run_full_summary saw "0 valid combinations".
# Fix: postprocess (creates the *_result.tre), then re-score competitors under mu
# (--force-recompute). No method re-runs. The GNN side for ne2M is already scored.

set -o pipefail   # not -u: the gene2net conda activation references an unset LD_LIBRARY_PATH
CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate gene2net
cd /groups/itay_mayrose/tomulanovski/gene2net

CONFIGS="conf_dup_loss_medium_10M_ne2M conf_dup_loss_high_10M_ne2M"

echo "================ postprocess ($(date +%H:%M:%S)) ================"
for CFG in $CONFIGS; do
  echo ">>> postprocess $CFG"
  python simulations/scripts/postprocess_results.py "$CFG" \
      || echo "  !! postprocess FAILED: $CFG"
done

echo "================ competitor mu re-score ($(date +%H:%M:%S)) ================"
for CFG in $CONFIGS; do
  echo ">>> run_full_summary $CFG"
  python simulations/scripts/run_full_summary.py "$CFG" --force-recompute \
      || echo "  !! run_full_summary FAILED: $CFG"
done

echo "================ DONE ($(date +%H:%M:%S)) ================"
