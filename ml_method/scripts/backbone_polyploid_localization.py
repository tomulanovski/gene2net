"""Localize ASTRAL's topology error: is it concentrated at the polyploid species?

We feed ASTRAL multi-copy (polyploid) gene trees, which it isn't built for. This
checks whether that corrupts the WHOLE backbone or just where the polyploids sit.

For each sample, compare the ASTRAL species tree to the TRUE diploid species tree
two ways (unrooted Robinson-Foulds, so this is pure topology error, not rooting):
  - full:   all species
  - pruned: polyploid species removed (single-copy backbone only)

If the full RF is high but the pruned RF collapses to ~0, ASTRAL's error is
concentrated at the polyploids -> a bounded fix (build the backbone from the
single-copy species where ASTRAL is reliable, then LEARN to attach the polyploids)
directly attacks the edit-distance gap. If the pruned RF stays high, the error is
diffuse and we'd need a cluster-based backbone instead.

Read-only. Run in any env with ete3.

Usage:
    python scripts/backbone_polyploid_localization.py --config ils_low --max-samples 300
"""
import argparse
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree


def load_tree_any(path):
    txt = open(path).read()
    if "#nexus" in txt.lower() or "begin trees" in txt.lower():
        for line in txt.splitlines():
            s = line.strip()
            if s.lower().startswith("tree") and "=" in s:
                return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(txt.strip(), format=1)


def norm_rf(t1, t2):
    """Unrooted normalized RF, or None if incomparable (ete3 returns 'NaN' string
    when the trees are too small to have informative splits)."""
    try:
        v = float(t1.compare(t2, unrooted=True)["norm_rf"])
    except (TypeError, ValueError, Exception):
        return None
    if v != v:  # NaN
        return None
    return v


def prune_to(tree, keep_names):
    t = tree.copy()
    keep = [l for l in t.get_leaf_names() if l in keep_names]
    if len(keep) < 3:
        return None
    t.prune(keep)
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--config", default="ils_low")
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-samples", type=int, default=300)
    args = ap.parse_args()

    sim_dir = os.path.join(args.data_root, "simphy", args.config)
    idxs = sorted(d for d in os.listdir(sim_dir)
                  if os.path.isdir(os.path.join(sim_dir, d)))[:args.max_samples]

    full_rfs, pruned_rfs = [], []
    n_full_ok = n_pruned_ok = 0
    n = skipped = 0
    for idx in idxs:
        astral_path = os.path.join(sim_dir, idx, f"replicate_{args.replicate}", "astral_species.tre")
        true_path = os.path.join(args.data_root, f"species_tree_{idx}.nex")
        md_path = os.path.join(args.data_root, f"metadata_{idx}.json")
        if not (os.path.exists(astral_path) and os.path.exists(true_path) and os.path.exists(md_path)):
            skipped += 1
            continue

        astral = load_tree_any(astral_path)
        true = load_tree_any(true_path)
        poly = set(json.load(open(md_path)).get("polyploid_species", {}))

        fr = norm_rf(astral, true)
        diploids = (set(true.get_leaf_names()) & set(astral.get_leaf_names())) - poly
        a_p = prune_to(astral, diploids)
        t_p = prune_to(true, diploids)
        if fr is None or a_p is None or t_p is None:
            skipped += 1
            continue
        pr = norm_rf(a_p, t_p)
        if pr is None:
            skipped += 1
            continue

        full_rfs.append(fr)
        pruned_rfs.append(pr)
        n_full_ok += (fr == 0)
        n_pruned_ok += (pr == 0)
        n += 1

    if n == 0:
        print("No comparable samples found — check paths.")
        return

    print(f"\nBackbone polyploid-localization on {n} samples (config {args.config}), {skipped} skipped\n")
    print("Full tree (all species), unrooted topology error vs the true tree:")
    print(f"  mean norm_rf = {statistics.mean(full_rfs):.3f}   median = {statistics.median(full_rfs):.3f}")
    print(f"  exactly correct topology (rf=0): {n_full_ok}/{n} ({100 * n_full_ok / n:.0f}%)")
    print("\nPruned tree (polyploids removed -> single-copy backbone only):")
    print(f"  mean norm_rf = {statistics.mean(pruned_rfs):.3f}   median = {statistics.median(pruned_rfs):.3f}")
    print(f"  exactly correct topology (rf=0): {n_pruned_ok}/{n} ({100 * n_pruned_ok / n:.0f}%)")
    print("\nReading:")
    print("  pruned RF << full RF  -> error concentrated at polyploids -> bounded fix viable")
    print("  pruned RF ~ full RF   -> error is diffuse -> need a cluster-based backbone")


if __name__ == "__main__":
    main()
