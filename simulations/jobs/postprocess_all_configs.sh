#!/bin/bash
#SBATCH --job-name=postprocess_all
#SBATCH --partition=itaym-pool
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --output=postprocess_all_%j.out

# Post-process program outputs into the standardized {program}_result.tre files
# that the summary pipeline reads, for all 15 configs.
#
# Why this has to run BEFORE run_full_summary: collect_results decides
# inferred_exists purely by whether {method_dir}/replicate_N/{output_file}
# exists. A config that was never post-processed therefore looks exactly like
# one where every method failed. That corrupts the completion rates in 3.2.1 AND
# biases every accuracy mean, because the mean is then taken over whichever runs
# happened to have been post-processed. Nothing in the summary output would look
# wrong.
#
# Safe to re-run. Post-processing an already-processed config just rewrites the
# same result files.
#
# ONE CAVEAT, read before assuming this job finishing means the work is done:
# the alloppnet branch of postprocess_results.py SUBMITS ITS OWN sbatch jobs
# rather than converting in place. So when this script exits, alloppnet results
# may still be pending. Check with squeue before starting the summary, or the
# summary will record alloppnet as incomplete.

set -o pipefail   # not -u: the gene2net conda activation references an unset LD_LIBRARY_PATH

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
conda activate gene2net || {
    echo "ERROR: Could not activate gene2net environment"
    exit 1
}

# postprocess_results.py resolves simulations/summary_config.yaml from the repo ROOT.
cd /groups/itay_mayrose/tomulanovski/gene2net || exit 1

CONFIGS="conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
conf_dup_loss_low_10M conf_dup_loss_medium_10M conf_dup_loss_high_10M \
conf_dup_loss_low_10M_ne1M conf_dup_loss_medium_10M_ne1M conf_dup_loss_high_10M_ne1M \
conf_dup_loss_low_10M_ne2M conf_dup_loss_medium_10M_ne2M conf_dup_loss_high_10M_ne2M \
conf_dup_loss_medium_10M_ne1M_fix025 conf_dup_loss_medium_10M_ne1M_fix050 conf_dup_loss_medium_10M_ne1M_fix075"

FAILED=""
for CFG in $CONFIGS; do
  echo "================ postprocess: $CFG  ($(date +%H:%M:%S)) ================"
  if ! python simulations/scripts/postprocess_results.py "$CFG"; then
      echo "  !! postprocess FAILED: $CFG"
      FAILED="$FAILED $CFG"
  fi
done

echo
echo "================ DONE  ($(date +%H:%M:%S)) ================"
if [ -n "$FAILED" ]; then
  echo "FAILED configs:$FAILED"
  echo "Do NOT start the summary until these are resolved. A config that failed"
  echo "here is indistinguishable from a config where the methods failed."
  exit 1
fi
echo "All 15 configs post-processed."
echo
echo "Next: confirm no alloppnet jobs are still queued"
echo "    squeue -u \$USER"
echo "then re-check coverage, then run:"
echo "    sbatch ml_method/jobs/resummary_mu.sh"
