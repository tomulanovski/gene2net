#!/bin/bash
#SBATCH --job-name=val_thresh
#SBATCH --partition=itaym-pool
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=val_thresh_%j.out

# Calibrate the ploidy-free (detect_only) decode threshold on the HELD-OUT
# validation split, leak-free. For each candidate threshold we reconstruct the
# same subsample of validation networks with detect_only and score them on the
# real reconstruction metric (edit distance / reticulation Jaccards). The
# threshold that minimizes validation edit distance is the published default;
# the whole curve is the sensitivity analysis. The subsample is seeded and
# identical across thresholds (SUBSAMPLE_SEED in score_val_reconstruction.py),
# so the comparison is apples-to-apples.
#
# Reconstruction (STEP 1) needs torch -> final_project env.
# Scoring (STEP 2) needs the comparison core / ReticulateTree -> gene2net env.
# No `set -u`: the gene2net conda activation references an unset LD_LIBRARY_PATH.

set -o pipefail

CONDA=/groups/itay_mayrose/tomulanovski/miniconda3
source "$CONDA/etc/profile.d/conda.sh"
cd /groups/itay_mayrose/tomulanovski/gene2net/ml_method

# The SAME six training dirs, in the SAME order used to train reconstruct_final
# (submit_finalize.sh). Order is load-bearing: the val split is reproduced from it.
DATA_DIR=""
for c in ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M; do
  DATA_DIR="$DATA_DIR data/mul_trees_2k/training_rooted/$c"
done

THRESHOLDS="0.3 0.4 0.5 0.6 0.7 0.8 0.9"

echo "================ STEP 1: val reconstruction (final_project) ================"
conda activate final_project
for TH in $THRESHOLDS; do
  TAG="t${TH/./}"
  echo ">>> detect_only threshold=$TH -> $TAG  ($(date +%H:%M:%S))"
  python scripts/score_val_reconstruction.py \
      --data-dir $DATA_DIR \
      --model-dir output/reconstruct_final \
      --model-config configs/reconstruct_final.yaml \
      --away-labels \
      --strategy detect_only --threshold "$TH" \
      --out-base "output/reconstruct_final/val_thresh/$TAG" \
      || echo "  !! val recon FAILED: threshold $TH"
done

echo "================ STEP 2: score (gene2net) ================"
conda activate gene2net
for TH in $THRESHOLDS; do
  TAG="t${TH/./}"
  python scripts/score_reconstructions.py \
      --recon-dir "output/reconstruct_final/val_thresh/$TAG/detect_only" --workers 8 \
      || echo "  !! score FAILED: threshold $TH"
done

echo "================ STEP 3: threshold curve (mean over val) ================"
python - <<'PY'
import glob, os
import pandas as pd
rows = []
for tag_dir in sorted(glob.glob("output/reconstruct_final/val_thresh/t*")):
    tag = os.path.basename(tag_dir)
    th = float("0." + tag[2:])       # t05 -> 0.5, t03 -> 0.3
    csv = os.path.join(tag_dir, "detect_only", "reconstruction_scores.csv")
    if not os.path.exists(csv):
        print(f"  MISSING scores for {tag} ({csv})")
        continue
    df = pd.read_csv(csv)
    m = df.mean(numeric_only=True)
    m["threshold"] = th
    m["n"] = len(df)
    rows.append(m)
if not rows:
    raise SystemExit("No threshold scores found under output/reconstruct_final/val_thresh/")
curve = pd.DataFrame(rows).sort_values("threshold")
front = ["threshold", "n"]
cols = front + [c for c in curve.columns if c not in front]
curve = curve[cols]
out = "output/reconstruct_final/val_thresh/threshold_curve.csv"
curve.to_csv(out, index=False)
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 50)
print(curve.to_string(index=False))
print(f"\nWrote {out}")
PY

echo "================ DONE  ($(date +%H:%M:%S)) ================"
