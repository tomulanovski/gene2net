"""Viability check for extending the away-relabel to CLADE-level allo events.

The single-species relabel skips allo events whose target is a clade (|target|>1) --
about a third of allo events. This audits those clade-level allo events to see whether
they are fixable the same way:
  - is the target clade MONOPHYLETIC in the true tree?  (its sibling A is well-defined)
  - is it MONOPHYLETIC in the ASTRAL tree the model uses? (its home H is well-defined)
  - for the ones monophyletic in BOTH (the fixable set), does the same partner==home
    bug occur, i.e. is the clade's ASTRAL home the labelled partner B?

Read like label_audit. If most clade allo are monophyletic in both and show the bug,
extending the relabel is a clean win. If many are non-monophyletic in ASTRAL, "the
clade's home" is fuzzy and the extension helps only the monophyletic subset.

gene2net env. Usage:
  python scripts/clade_allo_audit.py --configs ils_low dup_loss_high_ne1M --max 300
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


def clade_mrca(tree, C):
    """(mrca_node, is_monophyletic) for species set C in tree; (None, False) if a
    species is missing."""
    leaves = [l for l in tree.get_leaves() if l.name in C]
    if len({l.name for l in leaves}) < len(C):
        return None, False
    mrca = tree.get_common_ancestor(leaves) if len(leaves) > 1 else leaves[0]
    return mrca, set(mrca.get_leaf_names()) == set(C)


def clade_sibling(tree, C):
    """Sibling leaf set of C's MRCA in the true tree (the other parent A); None if C is
    not monophyletic or has no sibling."""
    mrca, mono = clade_mrca(tree, C)
    if mrca is None or not mono or mrca.up is None:
        return None
    sib = set()
    for ch in mrca.up.get_children():
        if ch is not mrca:
            sib |= set(ch.get_leaf_names())
    return sib - set(C)


def clade_astral_home(sample, C):
    """(home_leafset, is_monophyletic_in_astral) for clade C. Monophyletic iff C is
    itself an ASTRAL edge-clade. Home = smallest ASTRAL clade strictly containing C,
    minus C."""
    sd = {
        "species_tree_edge_index": sample.species_tree_edge_index,
        "species_list": sample.species_list,
        "species_tree_node_names": sample.species_tree_node_names,
        "species_tree_is_leaf": sample.species_tree_is_leaf,
    }
    clades = [c for _, c in sample_edge_bipartitions(sd)]
    x = frozenset(C)
    mono = x in set(clades)
    supersets = [c for c in clades if x < c]
    if not supersets:
        return None, mono
    return set(min(supersets, key=len) - x), mono


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--subdir", default="training_rooted")
    ap.add_argument("--configs", nargs="+", required=True)
    ap.add_argument("--max", type=int, default=300)
    args = ap.parse_args()

    for cfg in args.configs:
        train_dir = os.path.join(args.data_root, args.subdir, cfg)
        names = sorted(d for d in os.listdir(train_dir) if d.startswith("sample_"))[:args.max]
        n_clade = 0
        mono_true = mono_astral = mono_both = 0
        fixable = home_B = home_A = neither = 0
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
                target = ev.get("target_clade") or []
                B = set(ev.get("partner_clade") or [])
                if len(target) <= 1 or not B:
                    continue
                n_clade += 1
                C = set(target)
                _, tmono = clade_mrca(true_tree, C)
                H, amono = clade_astral_home(sample, C)
                if tmono:
                    mono_true += 1
                if amono:
                    mono_astral += 1
                if tmono and amono and H is not None:
                    mono_both += 1
                    A = clade_sibling(true_tree, C)
                    fixable += 1
                    if bool(H & B):
                        home_B += 1
                    elif A is not None and bool(H & A):
                        home_A += 1
                    else:
                        neither += 1

        print(f"\n=== {cfg} — {n_clade} clade-level allo events ===")
        if n_clade:
            print(f"  monophyletic in TRUE tree   : {mono_true}/{n_clade} ({100*mono_true/n_clade:.1f}%)")
            print(f"  monophyletic in ASTRAL      : {mono_astral}/{n_clade} ({100*mono_astral/n_clade:.1f}%)")
            print(f"  monophyletic in BOTH (fixable): {mono_both}/{n_clade} ({100*mono_both/n_clade:.1f}%)")
            if fixable:
                print(f"    of fixable, home == B (BUG): {home_B}/{fixable} ({100*home_B/fixable:.1f}%)")
                print(f"    of fixable, home == A (fine): {home_A}/{fixable} ({100*home_A/fixable:.1f}%)")
                print(f"    of fixable, home == neither : {neither}/{fixable} ({100*neither/fixable:.1f}%)")


if __name__ == "__main__":
    main()
