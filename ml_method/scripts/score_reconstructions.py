"""Score inferred MUL-trees against ground truth, reusing the existing
comparison core (reticulate_tree + compare_reticulations) so the numbers are
directly comparable to how GRAMPA/Polyphest were scored.

Expects the layout produced by reconstruct_mul_tree.py:
    <recon-dir>/sample_NNNN/output.tre        (inferred MUL-tree)
    <recon-dir>/sample_NNNN/ground_truth.nex  (true MUL-tree)

For each sample it builds ReticulateTree objects (with the same is_multree
conversion used in compare_nets.run_comparison_analysis), computes the pairwise
metrics, extracts the ground_truth-vs-gene2net values, and aggregates across
samples (mean/median). Metrics are distances — lower is better.

Run in the env where the comparison code imports work (the `gene2net` env).

Usage:
    python scripts/score_reconstructions.py \
        --recon-dir output/reconstruct_aligned/mul_trees/ils_low \
        --sim-scripts /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts
"""
import argparse
import os
import sys

import pandas as pd

METRICS = ["edit_distance", "ret_leaf_jaccard", "ret_sisters_jaccard",
           "num_rets_diff", "ploidy_diff"]


def newick_from_file(path):
    """Read a Newick string, extracting it from NEXUS if needed."""
    txt = open(path).read().strip()
    low = txt.lower()
    if low.startswith("#nexus") or "begin trees" in low:
        for line in txt.splitlines():
            s = line.strip()
            if s.lower().startswith("tree") and "=" in s:
                return s.split("=", 1)[1].strip()
    return txt


def build_rt(tree_str, ReticulateTree):
    """Mirror run_comparison_analysis: apply MUL-tree conversion when duplicated."""
    rt_temp = ReticulateTree(tree_str)
    if rt_temp.check_duplicated():
        return ReticulateTree(tree_str, is_multree=True)
    return rt_temp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--recon-dir", required=True,
                        help="Dir with sample_NNNN/{output.tre,ground_truth.nex}")
    parser.add_argument("--sim-scripts",
                        default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts",
                        help="Path to simulations/scripts for the comparison imports")
    parser.add_argument("--out-csv", default=None)
    args = parser.parse_args()

    sys.path.insert(0, args.sim_scripts)
    from reticulate_tree import ReticulateTree
    from compare_reticulations import pairwise_comparison

    rows = []
    sample_dirs = sorted(
        d for d in os.listdir(args.recon_dir)
        if os.path.isdir(os.path.join(args.recon_dir, d)) and d.startswith("sample_")
    )

    for name in sample_dirs:
        sdir = os.path.join(args.recon_dir, name)
        inf_path = os.path.join(sdir, "output.tre")
        gt_path = os.path.join(sdir, "ground_truth.nex")
        if not (os.path.exists(inf_path) and os.path.exists(gt_path)):
            continue
        try:
            inf = newick_from_file(inf_path)
            gt = newick_from_file(gt_path)
            rt_gt = build_rt(gt, ReticulateTree)
            rt_inf = build_rt(inf, ReticulateTree)
            df = pd.DataFrame([
                {"name": "ground_truth", "object": rt_gt, **rt_gt.measure(printout=False)},
                {"name": "gene2net", "object": rt_inf, **rt_inf.measure(printout=False)},
            ]).set_index("name")
            comp = pairwise_comparison(df, debug=False)
            row = {"sample": name}
            for m in METRICS:
                try:
                    row[m] = float(comp[m].loc["ground_truth", "gene2net"])
                except Exception:
                    row[m] = float("nan")
            rows.append(row)
        except Exception as e:
            print(f"  skip {name}: {type(e).__name__}: {e}")

    if not rows:
        print("No samples scored — check the recon-dir layout.")
        return

    df = pd.DataFrame(rows)
    print(f"\nScored {len(df)} samples from {args.recon_dir}\n")
    print("Per-metric summary (distances — lower is better):")
    print(df[METRICS].describe().loc[["mean", "50%", "min", "max"]].to_string())

    out_csv = args.out_csv or os.path.join(args.recon_dir, "reconstruction_scores.csv")
    df.to_csv(out_csv, index=False)
    print(f"\nSaved per-sample scores to {out_csv}")


if __name__ == "__main__":
    main()
