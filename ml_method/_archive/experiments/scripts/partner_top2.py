"""Who is the #1 confuser relative to the true partner?

For each single-species allo event where the true partner P is NOT top-1 (it is
usually #2), take the #1 confuser C and ask: does C co-cluster with P in the gene
trees? Two cases with different implications:
  (a) C and P do NOT co-cluster  -> C is in a different clade = the HOME parent.
      The top-2 are the two parents. Pure phasing story.
  (b) C and P DO co-cluster       -> C is a RELATIVE of P (same clade). Co-clustering
      found the right region, argmax picked the wrong tip. A clade-vs-tip resolution
      problem, potentially fixable without full phasing.

Reports the split (mean coclust(C,P), fraction that are P-relatives) plus examples.
Read-only. gene2net env. Usage:
    python scripts/partner_top2.py --config ils_low --n 120 --max-gene-trees 80
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


def cherry_counts(trees):
    """species-pair -> #trees they are direct cherry-sisters; and per-species X
    cherry tallies are derivable from it. Returns dict[(a,b)] with a<b and the
    per-tree count usable for both X-ranking and C-P co-clustering."""
    pair = defaultdict(int)
    for t in trees:
        seen = set()
        # group leaf children by parent
        for node in t.traverse():
            kids = [c for c in node.children if c.is_leaf()]
            names = sorted({k.name for k in kids})
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    seen.add((names[i], names[j]))
        for a, b in seen:
            pair[(a, b)] += 1
    return pair


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mul-trees-dir", default="data/mul_trees_2k")
    ap.add_argument("--config", default="ils_low")
    ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-gene-trees", type=int, default=80)
    ap.add_argument("--relative-thresh", type=float, default=0.20,
                    help="coclust(C,P) >= this -> count C as a relative of P")
    args = ap.parse_args()

    cp_vals = []          # coclust(C,P) when partner is not top-1
    relatives = 0
    n_confused = 0
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
        n_trees = len(trees)
        pair = cherry_counts(trees)

        def cc(a, b):
            return pair[(a, b)] if a < b else pair[(b, a)]

        all_sp = set()
        for t in trees:
            all_sp.update(t.get_leaf_names())

        for e in allo:
            X = next(iter(e["target_clade"]))
            P = next(iter(e.get("partner_clade") or []), None)
            if P is None or X not in all_sp:
                continue
            cands = [s for s in all_sp if s != X and cc(X, s) > 0]
            if not cands:
                continue
            ranked = sorted(cands, key=lambda s: -cc(X, s))
            if ranked[0] == P:
                continue  # partner already top-1
            C = ranked[0]
            n_confused += 1
            v = cc(C, P) / n_trees
            cp_vals.append(v)
            if v >= args.relative_thresh:
                relatives += 1
            if len(examples) < 15:
                prank = ranked.index(P) + 1 if P in ranked else -1
                examples.append((idx, X, C, P, prank, round(v, 2)))

    print(f"\nTop-2 identity ({args.config}) — {n_confused} confused allo events "
          f"(partner not top-1)\n")
    if n_confused:
        mean_cp = sum(cp_vals) / len(cp_vals)
        print(f"mean coclust(confuser C, true partner P) : {mean_cp:.3f}")
        print(f"C is a RELATIVE of P (coclust >= {args.relative_thresh}) : "
              f"{relatives}/{n_confused} = {relatives/n_confused:.2f}   (case b: clade-vs-tip)")
        print(f"C is ELSEWHERE / the home (coclust < {args.relative_thresh}) : "
              f"{n_confused-relatives}/{n_confused} = {(n_confused-relatives)/n_confused:.2f}   (case a: two parents)")
    print(f"\nExamples (idx: X  confuser C  partner P @rank   coclust(C,P)):")
    for idx, X, C, P, pr, v in examples:
        print(f"  {idx}: {X:>6}  C={C:>6}  P={P:>6} @{pr}   cc(C,P)={v}")


if __name__ == "__main__":
    main()
