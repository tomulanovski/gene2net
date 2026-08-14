"""Verify training_rooted/ samples are actually rooted (at the true root), vs the
old unrooted training/ samples.

For each index, compares the basal split (root's two children's leaf-name sets) of:
  - the rooted sample      (training_rooted/<config>/sample_NNNN)
  - the unrooted sample     (training/<config>/sample_NNNN)
against the true tree's basal split (species_tree_NNNN.nex).

If the rooted set matches the true root far more often than the unrooted set, the
re-packaging did its job (and used the fixed hybrid_root). Run in final_project.

Usage:
    python scripts/check_rooted_packaging.py --config ils_low --max-samples 100
"""
import argparse
import os
import pickle
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree


def load_nexus(path):
    for line in open(path).read().split("\n"):
        s = line.strip()
        if s.lower().startswith("tree") and "=" in s:
            return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(open(path).read().strip(), format=1)


def true_basal_split(path):
    t = load_nexus(path)
    return frozenset(frozenset(c.get_leaf_names()) for c in t.children)


def sample_basal_split(sample_pkl):
    with open(sample_pkl, "rb") as f:
        d = pickle.load(f)
    ei = d["species_tree_edge_index"]
    names = d["species_tree_node_names"]
    is_leaf = d["species_tree_is_leaf"]
    children = {}
    for k in range(0, ei.shape[1], 2):
        p = int(ei[0, k]); c = int(ei[1, k])
        children.setdefault(p, []).append(c)
    memo = {}

    def leaves(n):
        if n in memo:
            return memo[n]
        if bool(is_leaf[n]):
            r = {names[n]}
        else:
            r = set()
            for c in children.get(n, []):
                r |= leaves(c)
        memo[n] = r
        return r

    all_children = {c for cs in children.values() for c in cs}
    root = next(n for n in range(len(names)) if n not in all_children)
    return frozenset(frozenset(leaves(c)) for c in children.get(root, []))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--config", default="ils_low")
    ap.add_argument("--max-samples", type=int, default=100)
    args = ap.parse_args()

    rooted_dir = os.path.join(args.data_root, "training_rooted", args.config)
    unrooted_dir = os.path.join(args.data_root, "training", args.config)

    n = rooted_match = unrooted_match = rooted_vs_unrooted_differ = 0
    for name in sorted(os.listdir(rooted_dir))[:args.max_samples]:
        if not name.startswith("sample_"):
            continue
        idx = name.replace("sample_", "")
        true_path = os.path.join(args.data_root, f"species_tree_{idx}.nex")
        rooted_pkl = os.path.join(rooted_dir, name, "sample.pkl")
        unrooted_pkl = os.path.join(unrooted_dir, name, "sample.pkl")
        if not (os.path.exists(true_path) and os.path.exists(rooted_pkl) and os.path.exists(unrooted_pkl)):
            continue
        try:
            true_bs = true_basal_split(true_path)
            r_bs = sample_basal_split(rooted_pkl)
            u_bs = sample_basal_split(unrooted_pkl)
        except Exception:
            continue
        n += 1
        rooted_match += (r_bs == true_bs)
        unrooted_match += (u_bs == true_bs)
        rooted_vs_unrooted_differ += (r_bs != u_bs)

    if n == 0:
        print("No comparable samples found.")
        return
    print(f"\nRooted-packaging check on {n} samples ({args.config})\n")
    print(f"  rooted sample root == true root:   {rooted_match}/{n} ({100*rooted_match/n:.0f}%)")
    print(f"  unrooted sample root == true root: {unrooted_match}/{n} ({100*unrooted_match/n:.0f}%)")
    print(f"  rooted differs from unrooted:      {rooted_vs_unrooted_differ}/{n} ({100*rooted_vs_unrooted_differ/n:.0f}%)")
    print("\nReading: rooted% >> unrooted% and 'differs' high -> the rooted repackage is real and working.")


if __name__ == "__main__":
    main()
