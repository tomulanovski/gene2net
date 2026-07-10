"""Stage-0 phaser baseline: can gene-tree co-clustering recover an allopolyploid's
SECOND parent better than the current partner head's 0.44?

No training. For each allo event (single-species target X), rank other species by how
often a copy of X sits sister to them across gene trees, and check whether the true
partner_clade appears in the top-1/2/3. See docs/phaser_prototype_spec.md.

Run in the gene2net env (ete3).
Usage:
  python scripts/phaser_baseline.py --mul-trees-dir data/mul_trees_2k --config ils_low --n 200
  # the config where you compete:
  python scripts/phaser_baseline.py --mul-trees-dir data/mul_trees_2k --config dup_loss_high_ne1M --n 200
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree

from gene2net_gnn.data.features import compute_clustering_profile


def load_gene_trees(path, max_trees=500):
    trees = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                trees.append(Tree(line, format=1))
            except Exception:
                pass
            if len(trees) >= max_trees:
                break
    return trees


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mul-trees-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--max-gene-trees", type=int, default=500)
    args = ap.parse_args()

    n_allo = n_single = 0
    top1 = top2 = top3 = 0
    skipped = 0

    for idx in range(args.start, args.start + args.n):
        s = f"{idx:04d}"
        gt_path = os.path.join(args.mul_trees_dir, "simphy", args.config, s,
                               f"replicate_{args.replicate}", "gene_trees.tre")
        md_path = os.path.join(args.mul_trees_dir, f"metadata_{idx}.json")
        if not (os.path.exists(gt_path) and os.path.exists(md_path)):
            skipped += 1
            continue
        gene_trees = load_gene_trees(gt_path, args.max_gene_trees)
        if not gene_trees:
            skipped += 1
            continue
        all_species = set()
        for t in gene_trees:
            all_species.update(t.get_leaf_names())
        with open(md_path) as f:
            events = json.load(f).get("events", [])

        for ev in events:
            if ev.get("event_type") != "allo":
                continue
            n_allo += 1
            target = ev.get("target_clade") or []
            partner = set(ev.get("partner_clade") or [])
            if len(target) != 1 or not partner:
                continue          # focus on single-species targets (the common case)
            X = target[0]
            if X not in all_species:
                continue
            n_single += 1
            profile = compute_clustering_profile(gene_trees, X, all_species)
            ranked = sorted(profile, key=profile.get, reverse=True)
            top1 += int(any(sp in partner for sp in ranked[:1]))
            top2 += int(any(sp in partner for sp in ranked[:2]))
            top3 += int(any(sp in partner for sp in ranked[:3]))

    print(f"\nPhaser baseline (co-clustering top-k recovers the true partner) — config {args.config}")
    print(f"  allo events: {n_allo}   single-species (scored): {n_single}   samples skipped: {skipped}")
    if n_single:
        print(f"  partner in top-1: {top1}/{n_single} ({100*top1/n_single:.1f}%)")
        print(f"  partner in top-2: {top2}/{n_single} ({100*top2/n_single:.1f}%)   <- compare to 0.44")
        print(f"  partner in top-3: {top3}/{n_single} ({100*top3/n_single:.1f}%)")
        print("\nReading: top-2 >> 44% -> signal is rich, current head under-uses it -> build the phaser.")
        print("         top-2 ~ 44%  -> signal is the wall -> phasing won't beat it.")


if __name__ == "__main__":
    main()
