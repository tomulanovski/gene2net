"""Verify the partner-label bug directly from ground truth (no model).

For each allopolyploid event, compare X's ASTRAL home (its sibling in the rooted
species tree the model uses) to X's TWO true parents:
  B = partner_clade (metadata)
  A = X's sibling in the TRUE species tree (species_tree_<idx>.nex)  [the other parent]

Reports how often the ASTRAL home is:
  the labelled partner (B)  -> the labelling bug (target == home; build collapses)
  the other parent  (A)     -> label is fine (partner is the away parent)
  neither                   -> ASTRAL misplaced X (the ~25% ambiguous case)

Run in gene2net (ete3).
Usage:
  python scripts/label_audit.py --data-root data/mul_trees_2k --configs ils_low dup_loss_high_ne1M --max 200
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.data.metadata_labels import sample_edge_bipartitions


def load_nexus(path):
    for line in open(path).read().splitlines():
        s = line.strip()
        if s.lower().startswith("tree") and "=" in s:
            return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(open(path).read().strip(), format=1)


def sibling_clade(tree, sp):
    nodes = tree.search_nodes(name=sp)
    if not nodes or nodes[0].up is None:
        return None
    x = nodes[0]
    sib = set()
    for ch in x.up.get_children():
        if ch is not x:
            sib |= set(ch.get_leaf_names())
    sib.discard(sp)
    return sib


def astral_home(sample, sp):
    sd = {
        "species_tree_edge_index": sample.species_tree_edge_index,
        "species_list": sample.species_list,
        "species_tree_node_names": sample.species_tree_node_names,
        "species_tree_is_leaf": sample.species_tree_is_leaf,
    }
    clades = [c for _, c in sample_edge_bipartitions(sd)]
    x = frozenset({sp})
    supersets = [c for c in clades if x < c]
    if not supersets:
        return None
    return set(min(supersets, key=len) - x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--subdir", default="training_rooted")
    ap.add_argument("--configs", nargs="+", required=True)
    ap.add_argument("--max", type=int, default=200)
    args = ap.parse_args()

    for cfg in args.configs:
        train_dir = os.path.join(args.data_root, args.subdir, cfg)
        names = sorted(d for d in os.listdir(train_dir) if d.startswith("sample_"))[:args.max]
        n_allo = n_single = n_home_B = n_home_A = n_neither = n_A_missing = 0
        for name in names:
            idx = name.replace("sample_", "")
            md_path = os.path.join(args.data_root, f"metadata_{idx}.json")
            tt_path = os.path.join(args.data_root, f"species_tree_{idx}.nex")
            if not (os.path.exists(md_path) and os.path.exists(tt_path)):
                continue
            try:
                sample = Gene2NetSample.load(os.path.join(train_dir, name))
                true_tree = load_nexus(tt_path)
                events = json.load(open(md_path)).get("events", [])
            except Exception:
                continue
            for ev in events:
                if ev.get("event_type") != "allo":
                    continue
                n_allo += 1
                target = ev.get("target_clade") or []
                B = set(ev.get("partner_clade") or [])
                if len(target) != 1 or not B:
                    continue
                X = target[0]
                H = astral_home(sample, X)
                if H is None:
                    continue
                n_single += 1
                A = sibling_clade(true_tree, X)   # the other (true home) parent
                if A is None:
                    n_A_missing += 1
                home_is_B = bool(H & B)
                home_is_A = A is not None and bool(H & A)
                if home_is_B:
                    n_home_B += 1
                elif home_is_A:
                    n_home_A += 1
                else:
                    n_neither += 1

        print(f"\n=== {cfg} — {n_single} single-species allo events scored ({n_allo} allo total) ===")
        if n_single:
            print(f"  ASTRAL home == labelled partner B (THE BUG):  {n_home_B}/{n_single} ({100*n_home_B/n_single:.1f}%)")
            print(f"  ASTRAL home == other parent A (label fine):    {n_home_A}/{n_single} ({100*n_home_A/n_single:.1f}%)")
            print(f"  ASTRAL home == neither (X misplaced):          {n_neither}/{n_single} ({100*n_neither/n_single:.1f}%)")
            if n_A_missing:
                print(f"  (note: X absent from true species tree in {n_A_missing} cases -> A unknown)")


if __name__ == "__main__":
    main()
