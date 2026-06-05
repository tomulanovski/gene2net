"""Trivial single-feature baseline: how far does one raw edge feature get,
with no model at all?

Sweeps a threshold directly on one edge-feature column over the validation
edges and reports the best F1. Used to check how much the GNN adds over just
reading the dominant hand-crafted feature (e.g. frac_clade_dup).

Usage:
    python scripts/feature_baseline.py --data-dir <9 dirs> --feature frac_clade_dup
"""
import argparse
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset

EDGE_COLS = [
    "concordance", "branch_length", "clade_size", "depth",
    "dup_synchrony", "mirrored_sister_frac",
    "copy_pair_div_mean", "copy_pair_div_cv", "frac_clade_dup",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, nargs="+")
    parser.add_argument("--feature", default="frac_clade_dup", choices=EDGE_COLS)
    parser.add_argument("--expected-edge-dim", type=int, default=9)
    parser.add_argument("--max-samples", type=int, default=1500)
    args = parser.parse_args()

    col = EDGE_COLS.index(args.feature)

    datasets = [Gene2NetDataset(d) for d in args.data_dir]
    all_pairs = [(ds, i) for ds in datasets for i in range(len(ds))]
    indices = list(range(len(all_pairs)))
    random.seed(42)
    random.shuffle(indices)
    n_val = int(len(all_pairs) * 0.2)
    val_indices = indices[:min(n_val, args.max_samples)]

    vals, targets = [], []
    for idx in val_indices:
        ds, li = all_pairs[idx]
        try:
            s = ds[li]
        except Exception:
            continue
        if s.labels is None:
            continue
        ef = s.species_tree_edge_features
        if ef is None or ef.shape[1] != args.expected_edge_dim:
            continue
        wgd = s.labels.wgd_counts
        n_edges = ef.shape[0]
        for e in range(n_edges):
            vals.append(float(ef[e, col]))
            targets.append(1 if (wgd is not None and e < len(wgd) and wgd[e] > 0) else 0)

    vals = np.array(vals)
    targets = np.array(targets)
    print(f"Feature: {args.feature}  |  edges: {len(targets)}  positive: {targets.sum()} "
          f"({100*targets.mean():.1f}%)\n")

    print(f"{'thresh':>8} | {'Prec':>6} | {'Rec':>6} | {'F1':>6}")
    print("-" * 36)
    lo, hi = float(vals.min()), float(vals.max())
    best_f1, best_t = 0.0, lo
    for t in np.linspace(lo, hi, 40):
        preds = (vals >= t).astype(int)
        tp = int(((preds == 1) & (targets == 1)).sum())
        fp = int(((preds == 1) & (targets == 0)).sum())
        fn = int(((preds == 0) & (targets == 1)).sum())
        prec = tp / max(tp + fp, 1)
        rec = tp / max(tp + fn, 1)
        f1 = 2 * prec * rec / max(prec + rec, 1e-8)
        if f1 > best_f1:
            best_f1, best_t = f1, t
    # Print a few rows around the best for context
    for t in np.linspace(lo, hi, 11):
        preds = (vals >= t).astype(int)
        tp = int(((preds == 1) & (targets == 1)).sum())
        fp = int(((preds == 1) & (targets == 0)).sum())
        fn = int(((preds == 0) & (targets == 1)).sum())
        prec = tp / max(tp + fp, 1)
        rec = tp / max(tp + fn, 1)
        f1 = 2 * prec * rec / max(prec + rec, 1e-8)
        print(f"{t:>8.3f} | {prec:>6.3f} | {rec:>6.3f} | {f1:>6.3f}")

    print(f"\nBest single-feature F1 ({args.feature}): {best_f1:.3f} at threshold {best_t:.3f}")
    print("Compare to the full GNN's tuned F1 (~0.80). A large gap = the GNN")
    print("adds real value beyond this feature; a small gap = the feature is doing the work.")


if __name__ == "__main__":
    main()
