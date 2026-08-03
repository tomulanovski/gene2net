"""Signal check: does branch-length distance to a candidate parent discriminate the
TRUE partner, BEFORE we build a Ks-to-parent feature?

For each single-species allopolyploid event (target X, true partner P), rank candidate
species by the mean-over-gene-trees of the MIN branch-length distance from X's copies to
that species. Report how often the true partner P lands top-1/2/3 by that ranking, next
to the co-clustering baseline (which the phaser check put at ~top-1 0.42 / top-2 0.79 at
ils_low). If Ks-distance beats or complements co-clustering, the feature is worth building;
if it does not, hand-crafted partner is exhausted and the answer is phasing.

Read-only. Run in the gene2net env (ete3). Small defaults; raise --n / --max-gene-trees
if it's fast enough.

Usage:
    python scripts/partner_ks_signal.py --config ils_low --n 60 --max-gene-trees 60
"""
import argparse
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ete3 import Tree


def load_gene_trees(path, max_trees):
    trees = []
    for line in open(path):
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
    ap.add_argument("--mul-trees-dir", default="data/mul_trees_2k")
    ap.add_argument("--config", default="ils_low")
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-gene-trees", type=int, default=60)
    args = ap.parse_args()

    ks_rank, cc_rank, both_rank = [], [], []   # rank of true partner (1 = best)
    n_events = 0
    sim = os.path.join(args.mul_trees_dir, "simphy", args.config)
    for idx in sorted(os.listdir(sim))[:args.n]:
        gt = os.path.join(sim, idx, f"replicate_{args.replicate}", "gene_trees.tre")
        md = os.path.join(args.mul_trees_dir, f"metadata_{idx}.json")
        if not (os.path.exists(gt) and os.path.exists(md)):
            continue
        events = json.load(open(md)).get("events", [])
        allo = [e for e in events if e.get("event_type") == "allo"
                and len(e.get("target_clade") or []) == 1]
        if not allo:
            continue
        trees = load_gene_trees(gt, args.max_gene_trees)
        if not trees:
            continue
        all_sp = set()
        for t in trees:
            all_sp.update(t.get_leaf_names())

        for e in allo:
            X = next(iter(e["target_clade"]))
            P = next(iter(e.get("partner_clade") or []), None)
            if P is None or X not in all_sp:
                continue
            cands = [s for s in all_sp if s != X]
            ks_d = defaultdict(list)     # candidate -> list of min-distances per tree
            cc_c = defaultdict(int)      # candidate -> cherry-sister count
            used = 0
            for t in trees:
                xl = [l for l in t.get_leaves() if l.name == X]
                if not xl:
                    continue
                used += 1
                for l in xl:                      # co-cluster (cherry) tally
                    if l.up:
                        for c in l.up.children:
                            if c is not l and c.is_leaf() and c.name != X:
                                cc_c[c.name] += 1
                leaves_by_sp = defaultdict(list)  # branch-length distance tally
                for l in t.get_leaves():
                    if l.name in all_sp and l.name != X:
                        leaves_by_sp[l.name].append(l)
                for s, sl in leaves_by_sp.items():
                    ks_d[s].append(min(t.get_distance(a, b) for a in xl for b in sl))
            if used == 0:
                continue
            ks_list = sorted(ks_d, key=lambda s: sum(ks_d[s]) / len(ks_d[s]))       # ascending dist
            cc_list = sorted(cc_c, key=lambda s: -cc_c[s])                          # descending count
            # simple combination: average of the two ranks
            rankpos = {}
            for r, s in enumerate(ks_list):
                rankpos[s] = rankpos.get(s, 0) + r
            for r, s in enumerate(cc_list):
                rankpos[s] = rankpos.get(s, 0) + r
            both_list = sorted(rankpos, key=lambda s: rankpos[s])
            if P in ks_list:
                ks_rank.append(ks_list.index(P) + 1)
            if P in cc_list:
                cc_rank.append(cc_list.index(P) + 1)
            if P in both_list:
                both_rank.append(both_list.index(P) + 1)
            n_events += 1

    def topk(ranks, k):
        return sum(1 for r in ranks if r <= k) / len(ranks) if ranks else 0.0

    print(f"\nPartner signal check ({args.config}) — {n_events} single-species allo events\n")
    print(f"{'method':<14}{'top-1':>8}{'top-2':>8}{'top-3':>8}{'n':>6}")
    for name, ranks in [("Ks-distance", ks_rank), ("co-cluster", cc_rank), ("combined", both_rank)]:
        print(f"{name:<14}{topk(ranks,1):>8.2f}{topk(ranks,2):>8.2f}{topk(ranks,3):>8.2f}{len(ranks):>6}")
    print("\nReading: if Ks-distance top-2 clearly beats co-cluster (~0.79 at ils_low), or the")
    print("combined beats both, the Ks-to-parent feature is worth building. If it's <= co-cluster,")
    print("hand-crafted partner is exhausted -> the lever is phasing, not another feature.")


if __name__ == "__main__":
    main()
