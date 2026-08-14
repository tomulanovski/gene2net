#!/bin/bash
# ==============================================================================
# The away-parent label experiment. Tests whether fixing the ~56% partner==home
# labelling bug (label_audit) raises one-parter PARTNER ACCURACY.
#
# label_audit found ASTRAL places X next to a TRUE parent ~94% of the time, so
# "the away parent = the parent that is NOT X's home" is well-defined ~94% of the
# time. The away relabel makes the training target consistent (always the non-home
# parent). This run trains the one-parter twice, differing ONLY in the labels:
#   orig : labels_clade.pkl  (original clade labels, partner == metadata partner)
#   away : labels_away.pkl   (partner retargeted to the non-home parent)
# and we compare best allo accuracy. Same one-parter config, from scratch, seed=42.
#
# Relabels write SEPARATE sidecars (labels_clade.pkl vs labels_away.pkl) so the two
# label sets coexist and neither clobbers the other.
#
# USAGE:  ./jobs/submit_away_experiment.sh
# ==============================================================================
set -euo pipefail

BASE_DIR="/groups/itay_mayrose/tomulanovski/gene2net/ml_method"
CONDA_PATH=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA_PATH/etc/profile.d/conda.sh"
conda activate final_project
cd "$BASE_DIR"

CONFIGS="ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M"

echo "### Relabeling (writes labels_clade.pkl AND labels_away.pkl per sample) ###"
for cfg in $CONFIGS; do
    echo "--- $cfg : original clade labels ---"
    python scripts/relabel_from_metadata.py --training-subdir "training_rooted/$cfg"
    echo "--- $cfg : away-parent labels ---"
    python scripts/relabel_from_metadata.py --training-subdir "training_rooted/$cfg" --away-parent
done

DATA_DIR=""
for cfg in $CONFIGS; do
    DATA_DIR="$DATA_DIR ${BASE_DIR}/data/mul_trees_2k/training_rooted/$cfg"
done

PARTITION="gpu-rotemhsh-pool"; ACCOUNT="--account=itaym-users_v2"; QOS="owner"
LOG_DIR="${BASE_DIR}/logs"; mkdir -p "$LOG_DIR"

submit () {  # $1 = tag, $2 = extra export (CLADE_LABELS=1 or AWAY_LABELS=1)
    local tag="$1"; local flag="$2"
    echo "Submitting one-parter ($tag)"
    sbatch \
        --job-name="g2n_away_${tag}" \
        --output="${LOG_DIR}/away_${tag}_%j.out" \
        --error="${LOG_DIR}/away_${tag}_%j.err" \
        --time=1-00:00:00 --mem=128G --cpus-per-task=4 \
        --partition=${PARTITION} ${ACCOUNT} --qos=${QOS} --gres=gpu:1 \
        --export=ALL,DATA_DIR="${DATA_DIR}",BASE_DIR="${BASE_DIR}",CONFIG="${BASE_DIR}/configs/reconstruct_oneparter.yaml",OUTPUT_DIR="${BASE_DIR}/output/reconstruct_op_${tag}",INIT_FROM="",${flag} \
        "${BASE_DIR}/jobs/train_reconstruct_job.sh"
}

echo "### Submitting two one-parter trainings (orig vs away) ###"
submit "orig" "CLADE_LABELS=1"
submit "away" "AWAY_LABELS=1"

echo "Done. Watch: squeue -u \$USER"
echo "Read: compare best allo in logs/away_orig_*.out vs logs/away_away_*.out"
