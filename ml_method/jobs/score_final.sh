#!/bin/bash
# Score the re-benchmark outputs. MUST run in the gene2net env (the comparison core
# / ReticulateTree lives there; final_project errors with a GLIBCXX ImportError).
#   conda activate gene2net
#   bash jobs/score_final.sh
set -euo pipefail
cd /groups/itay_mayrose/tomulanovski/gene2net/ml_method
mkdir -p logs

CONFIGS="conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
         conf_dup_loss_low_10M_ne1M conf_dup_loss_medium_10M_ne1M conf_dup_loss_high_10M_ne1M"

for c in $CONFIGS; do
    echo "==================== score $c ===================="
    python scripts/score_reconstructions.py \
        --recon-dir "output/reconstruct_final/benchmark/$c/bound_driven" \
        --workers 8 2>&1 | tee "logs/score_${c}.out"
done

echo "Done. Per-config scores in logs/score_*.out and each recon-dir's reconstruction_scores.csv"
