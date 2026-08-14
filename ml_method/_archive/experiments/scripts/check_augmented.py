#!/usr/bin/env python
"""Sanity-check augmented samples before training.

Verifies that every sample in a training dir has 9-dim edge features, that the
5 new detection features are populated (not all zero), and — most importantly —
that they are higher on true WGD edges than on non-WGD edges (i.e. they carry
signal). Run on one config dir after augmenting it.

Usage:
    python check_augmented.py --data-dir data/mul_trees_2k/training/ils_low
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset

# Column layout of species_tree_edge_features (see dataset.py).
COLS = [
    "concordance", "branch_length", "clade_size", "depth",        # base (0-3)
    "dup_synchrony", "mirrored_sister_frac",                       # new (4-5)
    "copy_pair_div_mean", "copy_pair_div_cv", "frac_clade_dup",    # new (6-8)
]
NEW_COLS = list(range(4, 9))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, help="A training/<config> directory")
    parser.add_argument("--max-samples", type=int, default=0, help="0 = all")
    args = parser.parse_args()

    ds = Gene2NetDataset(args.data_dir)
    n = len(ds)
    if args.max_samples:
        n = min(n, args.max_samples)
    print(f"Samples found: {len(ds)} (checking {n})\n")

    dim_counts = {}
    bad_dim = 0
    # Accumulators for label-split feature means.
    sums = {c: [0.0, 0.0] for c in NEW_COLS}   # col -> [sum_pos, sum_neg]
    counts = [0, 0]                            # [n_pos_edges, n_neg_edges]
    nonzero = {c: 0 for c in NEW_COLS}
    total_edges = 0

    for i in range(n):
        s = ds[i]
        ef = s.species_tree_edge_features
        d = ef.shape[1] if ef is not None else 0
        dim_counts[d] = dim_counts.get(d, 0) + 1
        if d != 9:
            bad_dim += 1
            continue

        labels = s.labels
        wgd_counts = labels.wgd_counts if labels is not None else None
        n_edges = ef.shape[0]
        for e in range(n_edges):
            total_edges += 1
            is_wgd = 1 if (wgd_counts is not None and e < len(wgd_counts) and wgd_counts[e] > 0) else 0
            counts[0 if is_wgd else 1] += 1
            for c in NEW_COLS:
                v = float(ef[e, c])
                sums[c][0 if is_wgd else 1] += v
                if abs(v) > 1e-9:
                    nonzero[c] += 1

    print("Edge-feature dim histogram:", dim_counts)
    if bad_dim:
        print(f"  WARNING: {bad_dim} samples are NOT 9-dim — augment incomplete!\n")
    else:
        print("  OK: all checked samples are 9-dim.\n")

    print(f"Total edges: {total_edges}  (WGD: {counts[0]}, non-WGD: {counts[1]})\n")
    print(f"{'feature':>20} | {'% nonzero':>9} | {'mean WGD':>9} | {'mean non':>9} | {'ratio':>6}")
    print("-" * 70)
    for c in NEW_COLS:
        pct_nz = 100.0 * nonzero[c] / total_edges if total_edges else 0.0
        mean_pos = sums[c][0] / counts[0] if counts[0] else 0.0
        mean_neg = sums[c][1] / counts[1] if counts[1] else 0.0
        ratio = (mean_pos / mean_neg) if mean_neg > 1e-9 else float("inf")
        print(f"{COLS[c]:>20} | {pct_nz:>8.1f}% | {mean_pos:>9.4f} | {mean_neg:>9.4f} | {ratio:>6.2f}")

    print("\nInterpretation: for a useful detection feature, 'mean WGD' should be")
    print("clearly higher than 'mean non' (ratio > 1). dup_synchrony, mirrored_sister")
    print("and frac_clade_dup are the ones expected to separate most.")


if __name__ == "__main__":
    main()
