"""Home/away discriminator signal-check: does the away partner separate from the
home confuser by WHEN it co-clusters with X?

Theory: X's home subgenome sisters its (diploid) home parent even in single-copy
gene trees, while the away partner only appears when X is DUPLICATED (both copies
present). So ranking candidates by the CONTRAST
    contrast[s] = coclust(s | X duplicated)  -  coclust(s | X single-copy)
should push the away partner above the home, where raw co-clustering ranks the
home first. This is distinct from ablation B (which used only the duplicated-tree
term, not the contrast). Only works if there is a DOMINANT subgenome (home copy
retained more often) -- this check tells us whether that holds in the sim.

Reports how often the true partner is top-1/2/3 under RAW co-clustering vs the
CONTRAST. If contrast top-1 clearly beats raw (~0.42), build the feature; if flat,
partner is closed and phasing is the answer.

Read-only. gene2net env. Usage:
    python scripts/partner_contrast.py --config ils_low --n 120 --max-gene-trees 80
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


def cherry_partners(tree, X):
    """Species that are a direct cherry-sister of some copy of X in this tree."""
    out = set()
    for l in tree.get_leaves():
        if l.name != X or l.up is None:
            continue
        for c in l.up.children:
            if c is not l and c.is_leaf() and c.name != X:
                out.add(c.name)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mul-trees-dir", default="data/mul_trees_2k")
    ap.add_argument("--config", default="ils_low")
    ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-gene-trees", type=int, default=80)
    args = ap.parse_args()

    raw_rank, con_rank = [], []
    n_events = 0
    n_no_single = 0
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

        for e in allo:
            X = next(iter(e["target_clade"]))
            P = next(iter(e.get("partner_clade") or []), None)
            if P is None:
                continue
            dup_cc = defaultdict(int); single_cc = defaultdict(int)
            n_dup = n_single = 0
            for t in trees:
                names = t.get_leaf_names()
                xc = names.count(X)
                if xc == 0:
                    continue
                partners = cherry_partners(t, X)
                if xc >= 2:
                    n_dup += 1
                    for s in partners:
                        dup_cc[s] += 1
                else:
                    n_single += 1
                    for s in partners:
                        single_cc[s] += 1
            if n_dup == 0:
                continue
            if n_single == 0:
                n_no_single += 1
            cands = set(dup_cc) | set(single_cc)
            if not cands:
                continue
            raw = {s: (dup_cc[s] + single_cc[s]) / (n_dup + n_single) for s in cands}
            dfrac = {s: dup_cc[s] / n_dup for s in cands}
            sfrac = {s: (single_cc[s] / n_single if n_single else 0.0) for s in cands}
            contrast = {s: dfrac[s] - sfrac[s] for s in cands}

            raw_list = sorted(cands, key=lambda s: -raw[s])
            con_list = sorted(cands, key=lambda s: -contrast[s])
            if P in raw_list:
                raw_rank.append(raw_list.index(P) + 1)
            if P in con_list:
                con_rank.append(con_list.index(P) + 1)
            n_events += 1

    def topk(ranks, k):
        return sum(1 for r in ranks if r <= k) / len(ranks) if ranks else 0.0

    print(f"\nHome/away contrast check ({args.config}) — {n_events} allo events "
          f"({n_no_single} had no single-copy trees)\n")
    print(f"{'ranking':<12}{'top-1':>8}{'top-2':>8}{'top-3':>8}{'n':>6}")
    for name, ranks in [("raw coclust", raw_rank), ("contrast", con_rank)]:
        print(f"{name:<12}{topk(ranks,1):>8.2f}{topk(ranks,2):>8.2f}{topk(ranks,3):>8.2f}{len(ranks):>6}")
    print("\nRead: if contrast top-1 clearly beats raw (~0.42), the away partner separates")
    print("from the home by duplication timing -> build the feature. If contrast <= raw,")
    print("there is no dominant-subgenome signal to exploit -> partner is closed, phasing is it.")


if __name__ == "__main__":
    main()
