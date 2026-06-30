"""Does ASTRAL place each polyploid next to ONE of its two true parents, or neither?

An allopolyploid X is a hybrid of two lineages (symmetric — no real "home" vs
"partner"). In the true MUL-tree X appears twice, so it has two true neighbor
sets. ASTRAL puts X once. The question: is ASTRAL's single placement next to one
of those two true neighbors (then the current pipeline already handles X fine —
it keeps that placement and predicts the other parent), or next to neither (then
ASTRAL's placement is genuinely wrong and re-placing it has headroom)?

For each polyploid X:
  - true neighbor sets  = species in each copy's sibling subtree, from mul_tree_{idx}.nex
  - astral neighbor set = species in X's sibling subtree, from astral_species.tre
  - best Jaccard of the astral set vs the two true sets -> "hits a true parent"
    if >= threshold, else "neither".

Read-only. Run in any env with ete3.

Usage:
    python scripts/backbone_polyploid_placement_accuracy.py --config ils_low --max-samples 300
"""
import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree
from gene2net_gnn.data.mul_tree_generator import get_polyploid_species


def load_tree_any(path):
    txt = open(path).read()
    if "#nexus" in txt.lower() or "begin trees" in txt.lower():
        for line in txt.splitlines():
            s = line.strip()
            if s.lower().startswith("tree") and "=" in s:
                return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(txt.strip(), format=1)


def sibling_species_sets(tree, x):
    """Species-name set in the sibling subtree of every leaf named x (self-copies removed)."""
    out = []
    for leaf in tree.get_leaves():
        if leaf.name == x and leaf.up is not None:
            sibs = set()
            for child in leaf.up.children:
                if child is leaf:
                    continue
                sibs |= set(child.get_leaf_names())
            sibs.discard(x)
            if sibs:
                out.append(sibs)
    return out


def jaccard(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--config", default="ils_low")
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-samples", type=int, default=300)
    ap.add_argument("--threshold", type=float, default=0.5)
    args = ap.parse_args()

    sim_dir = os.path.join(args.data_root, "simphy", args.config)
    idxs = sorted(d for d in os.listdir(sim_dir)
                  if os.path.isdir(os.path.join(sim_dir, d)))[:args.max_samples]

    best_jaccards = []
    n_hit = n_neither = 0
    n_poly = n_samples = 0
    for idx in idxs:
        astral_path = os.path.join(sim_dir, idx, f"replicate_{args.replicate}", "astral_species.tre")
        mul_path = os.path.join(args.data_root, f"mul_tree_{idx}.nex")
        if not (os.path.exists(astral_path) and os.path.exists(mul_path)):
            continue
        astral = load_tree_any(astral_path)
        mul = load_tree_any(mul_path)
        polys = get_polyploid_species(mul)
        if not polys:
            continue
        n_samples += 1
        for x in polys:
            true_sets = sibling_species_sets(mul, x)      # two (or more) true neighbors
            astral_sets = sibling_species_sets(astral, x)  # one (ASTRAL is single-copy)
            if not true_sets or not astral_sets:
                continue
            na = astral_sets[0]
            # strip the polyploids themselves from comparison so the match is to
            # the diploid context, not to other (also-displaced) polyploids
            best = max(jaccard(na, ts) for ts in true_sets)
            best_jaccards.append(best)
            if best >= args.threshold:
                n_hit += 1
            else:
                n_neither += 1
            n_poly += 1

    if n_poly == 0:
        print("No polyploids measured — check paths.")
        return

    print(f"\nPolyploid placement accuracy: {n_poly} polyploids over {n_samples} samples "
          f"(config {args.config})\n")
    print(f"ASTRAL placement next to ONE of the two true parents (Jaccard >= {args.threshold}):")
    print(f"  hits a true parent: {n_hit}/{n_poly} ({100 * n_hit / n_poly:.0f}%)")
    print(f"  neither (wrong):    {n_neither}/{n_poly} ({100 * n_neither / n_poly:.0f}%)")
    print(f"  mean best-Jaccard:  {statistics.mean(best_jaccards):.3f}   "
          f"median {statistics.median(best_jaccards):.3f}")
    print("\nReading:")
    print("  mostly 'hits a true parent' -> current pipeline already places polyploids OK;")
    print("    the edit gap is rooting/prediction, NOT polyploid placement -> backbone fix has little headroom.")
    print("  lots of 'neither' -> ASTRAL mis-places polyploids -> re-placing them (the backbone fix) has headroom.")


if __name__ == "__main__":
    main()
