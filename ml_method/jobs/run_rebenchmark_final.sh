#!/bin/bash
# Re-benchmark the locked final model on the six configs (bound_driven), capturing
# runtime (--profile) and copy-number under-shoot (--ploidy-check) for each, then
# scoring. Run in final_project on a compute node (light: ~5 s/net, CPU is fine).
#   bash jobs/run_rebenchmark_final.sh
# Afterward: head_to_head against the baselines (needs reconstruction_scores.csv,
# which score_reconstructions writes).
set -euo pipefail
cd /groups/itay_mayrose/tomulanovski/gene2net/ml_method
mkdir -p logs

CONFIGS="conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
         conf_dup_loss_low_10M_ne1M conf_dup_loss_medium_10M_ne1M conf_dup_loss_high_10M_ne1M"

for c in $CONFIGS; do
    echo "==================== $c ===================="
    python scripts/benchmark_networks.py \
        --model-dir output/reconstruct_final \
        --model-config configs/reconstruct_final.yaml \
        --config "$c" \
        --strategies bound_driven \
        --parents head \
        --profile --ploidy-check \
        --out-base "output/reconstruct_final/benchmark/$c" 2>&1 | tee "logs/rebench_${c}.out"

    python scripts/score_reconstructions.py \
        --recon-dir "output/reconstruct_final/benchmark/$c/bound_driven" \
        --workers 8 2>&1 | tee "logs/score_${c}.out"
done

echo "Done. Per-config profile + ploidy-check in logs/rebench_*.out; scores in logs/score_*.out"
