"""Micro-profile the feature computation to find the hot function before optimizing.

Times each feature block separately on a few networks: per-species copy counts,
per-species co-clustering summary, per-edge base features, and per-edge WGD-detection
features. The runtime profile showed feature computation is ~97% of our per-network
compute (7.3 s), so this says which block to optimize.

Read-only. gene2net env. Uses the sim data (same feature code as inference).
Usage:
  python scripts/profile_features.py --config ils_low --n 5 --max-gene-trees 500
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ete3 import Tree

from gene2net_gnn.data.features import (
    compute_copy_count_features,
    compute_clustering_summary,
    compute_species_tree_edge_features,
    compute_species_tree_edge_detection_features,
)


def load_trees(path, maxt):
    out = []
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        try:
            out.append(Tree(s, format=1))
        except Exception:
            pass
        if len(out) >= maxt:
            break
    return out


def load_nexus_or_newick(path):
    for line in open(path).read().splitlines():
        s = line.strip()
        if s.lower().startswith("tree") and "=" in s:
            return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(open(path).read().strip(), format=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mul-trees-dir", default="data/mul_trees_2k")
    ap.add_argument("--config", default="ils_low")
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-gene-trees", type=int, default=500)
    args = ap.parse_args()

    sim = os.path.join(args.mul_trees_dir, "simphy", args.config)
    t = {"copy": 0.0, "cluster": 0.0, "edge_base": 0.0, "edge_det": 0.0, "n": 0, "n_gt": 0}
    for idx in sorted(os.listdir(sim))[:args.n]:
        gt = os.path.join(sim, idx, f"replicate_{args.replicate}", "gene_trees.tre")
        sp = os.path.join(sim, idx, f"replicate_{args.replicate}", "astral_species.tre")
        if not (os.path.exists(gt) and os.path.exists(sp)):
            continue
        trees = load_trees(gt, args.max_gene_trees)
        astral = load_nexus_or_newick(sp)
        species = set(astral.get_leaf_names())
        if not trees:
            continue
        t["n"] += 1
        t["n_gt"] += len(trees)

        s = time.perf_counter()
        for name in species:
            compute_copy_count_features(trees, name)
        t["copy"] += time.perf_counter() - s

        s = time.perf_counter()
        for name in species:
            compute_clustering_summary(trees, name, species)
        t["cluster"] += time.perf_counter() - s

        s = time.perf_counter()
        compute_species_tree_edge_features(astral, trees)
        t["edge_base"] += time.perf_counter() - s

        s = time.perf_counter()
        compute_species_tree_edge_detection_features(astral, trees)
        t["edge_det"] += time.perf_counter() - s

    n = max(t["n"], 1)
    print(f"\nFeature micro-profile ({args.config}) — {t['n']} networks, "
          f"mean {t['n_gt']/n:.0f} gene trees, {len(species)} species\n")
    for k, label in [("copy", "per-species copy counts"),
                     ("cluster", "per-species co-clustering summary"),
                     ("edge_base", "per-edge base (concordance, bl, size, depth)"),
                     ("edge_det", "per-edge WGD-detection features")]:
        print(f"  {label:<45}{t[k]/n:>8.3f} s/net")
    total = (t["copy"] + t["cluster"] + t["edge_base"] + t["edge_det"]) / n
    print(f"  {'TOTAL':<45}{total:>8.3f} s/net")


if __name__ == "__main__":
    main()
