"""Who beats the true partner for top-1 co-clustering?

For each single-species allopolyploid event (target X, true partner P), rank
candidate species by cherry co-clustering with X across gene trees. When the true
partner P is NOT the top candidate, record WHO is (the confuser C) and classify C:
  - is C a diploid (mean copies ~1) or itself a polyploid (mean copies > 1.5)?
  - the theory says the confuser is X's HOME parent -- usually a diploid neighbor.
Reports the breakdown plus concrete X -> C (true P at rank r) examples so we can
see whether the confuser is systematically a diploid home lineage or something else.

Read-only. gene2net env. Usage:
    python scripts/partner_confuser.py --config ils_low --n 80 --max-gene-trees 60
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
    ap.add_argument("--n", type=int, default=80)
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-gene-trees", type=int, default=60)
    ap.add_argument("--examples", type=int, default=15)
    args = ap.parse_args()

    n_events = 0
    partner_is_top1 = 0
    confuser_diploid = 0       # top-1 confuser has mean copies < 1.5
    confuser_polyploid = 0     # top-1 confuser has mean copies >= 1.5
    confuser_is_a_partner = 0  # confuser is the TRUE partner of some OTHER event in the sample
    examples = []

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
        # mean copy count per species across these trees
        copies = defaultdict(list)
        for t in trees:
            c = defaultdict(int)
            for l in t.get_leaf_names():
                c[l] += 1
            for s in all_sp:
                copies[s].append(c.get(s, 0))
        mean_copies = {s: (sum(v) / len(v) if v else 0.0) for s, v in copies.items()}
        partners_here = {next(iter(e.get("partner_clade") or []), None) for e in allo}

        for e in allo:
            X = next(iter(e["target_clade"]))
            P = next(iter(e.get("partner_clade") or []), None)
            if P is None or X not in all_sp:
                continue
            cc = defaultdict(int)
            used = 0
            for t in trees:
                xl = [l for l in t.get_leaves() if l.name == X]
                if not xl:
                    continue
                used += 1
                for l in xl:
                    if l.up:
                        for c in l.up.children:
                            if c is not l and c.is_leaf() and c.name != X:
                                cc[c.name] += 1
            if used == 0 or not cc:
                continue
            ranked = sorted(cc, key=lambda s: -cc[s])
            n_events += 1
            top1 = ranked[0]
            if top1 == P:
                partner_is_top1 += 1
                continue
            # partner not top-1: characterize the confuser
            mc = mean_copies.get(top1, 0.0)
            if mc >= 1.5:
                confuser_polyploid += 1
            else:
                confuser_diploid += 1
            if top1 in partners_here and top1 != P:
                confuser_is_a_partner += 1
            if len(examples) < args.examples:
                prank = ranked.index(P) + 1 if P in ranked else -1
                examples.append((idx, X, top1, round(mc, 2), P, prank))

    print(f"\nConfuser analysis ({args.config}) — {n_events} single-species allo events\n")
    if n_events == 0:
        print("no events"); return
    wrong = n_events - partner_is_top1
    print(f"true partner IS top-1        : {partner_is_top1}/{n_events} = {partner_is_top1/n_events:.2f}")
    print(f"true partner NOT top-1       : {wrong}/{n_events} = {wrong/n_events:.2f}")
    if wrong:
        print(f"  of those, top-1 confuser is:")
        print(f"    a DIPLOID (mean copies <1.5) : {confuser_diploid}/{wrong} = {confuser_diploid/wrong:.2f}")
        print(f"    a POLYPLOID (>=1.5)          : {confuser_polyploid}/{wrong} = {confuser_polyploid/wrong:.2f}")
        print(f"    the true partner of ANOTHER event : {confuser_is_a_partner}/{wrong} = {confuser_is_a_partner/wrong:.2f}")
    print(f"\nExamples (idx: X  -> confuser C [meanCopies]   true partner P @rank):")
    for idx, X, C, mc, P, pr in examples:
        print(f"  {idx}: {X:>6} -> {C:>6} [{mc}]   partner {P:>6} @{pr}")


if __name__ == "__main__":
    main()
