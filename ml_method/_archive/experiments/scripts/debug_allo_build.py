"""Dump the allo decompose->build round-trip for ONE sample + species, to pinpoint
why an allopolyploid's two homeologs land next to the wrong parents.

Shows, for the target polyploid species:
  - where each of its copies sits in the ground-truth MUL-tree (sibling leaf sets),
  - the events decompose_mul_tree extracts for it (wgd_clade / partner_clade),
  - where it sits in the backbone used for building,
  - where each copy lands in the REBUILT tree (sibling leaf sets).

Compare the GT sibling sets to the rebuilt sibling sets: that's the ret_sisters
mismatch, at the source.

Run in gene2net env (ete3).
Usage:
  python scripts/debug_allo_build.py --mul-trees-dir data/mul_trees_2k --config ils_low --idx 27 --species sp39
  # backbone: true species tree (default) or astral
  python scripts/debug_allo_build.py ... --backbone astral
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree

from gene2net_gnn.data.label_extractor import decompose_mul_tree, extract_backbone
from gene2net_gnn.inference.mul_tree_builder import build_mul_tree


def load_nexus_tree(path):
    for line in open(path).read().splitlines():
        s = line.strip()
        if s.lower().startswith("tree") and "=" in s:
            return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(open(path).read().strip(), format=1)


def copy_siblings(tree, species):
    """For each copy of `species`, the leaf-name set of its siblings (minus itself)."""
    out = []
    for leaf in tree.get_leaves():
        if leaf.name != species or leaf.up is None:
            continue
        sib = set()
        for ch in leaf.up.get_children():
            if ch is not leaf:
                sib |= set(ch.get_leaf_names())
        sib.discard(species)
        out.append(sorted(sib))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mul-trees-dir", required=True)
    ap.add_argument("--config", default="ils_low")
    ap.add_argument("--idx", type=int, required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--backbone", choices=["true", "astral"], default="true")
    args = ap.parse_args()

    s = f"{args.idx:04d}"
    gt = load_nexus_tree(os.path.join(args.mul_trees_dir, f"mul_tree_{s}.nex"))
    if args.backbone == "true":
        bb = load_nexus_tree(os.path.join(args.mul_trees_dir, f"species_tree_{s}.nex"))
    else:
        bb = Tree(open(os.path.join(args.mul_trees_dir, "simphy", args.config, s,
                  f"replicate_{args.replicate}", "astral_species.tre")).read().strip(), format=1)

    sp = args.species
    print(f"\n=== sample {s}, species {sp}, backbone={args.backbone} ===")

    print(f"\nGT MUL-tree: {sp} appears {sum(1 for l in gt.get_leaves() if l.name==sp)}x")
    for i, sib in enumerate(copy_siblings(gt, sp)):
        print(f"  copy {i} siblings: {sib}")

    mul_bb = extract_backbone(gt)
    bn = {l.name: l for l in mul_bb.get_leaves()}.get(sp)
    if bn is not None and bn.up is not None:
        exp = set()
        for ch in bn.up.get_children():
            if ch is not bn:
                exp |= set(ch.get_leaf_names())
        print(f"\nMUL-derived backbone: {sp} sibling = {sorted(exp - {sp})}")

    bbn = {l.name: l for l in bb.get_leaves()}.get(sp)
    if bbn is not None and bbn.up is not None:
        exp = set()
        for ch in bbn.up.get_children():
            if ch is not bbn:
                exp |= set(ch.get_leaf_names())
        print(f"BUILD backbone ({args.backbone}): {sp} sibling = {sorted(exp - {sp})}")

    events = decompose_mul_tree(gt)
    print(f"\ndecomposed events touching {sp}:")
    for e in events:
        if sp in e.wgd_edge_clade or sp in e.partner_edge_clade:
            print(f"  wgd_clade={sorted(e.wgd_edge_clade)}  partner_clade={sorted(e.partner_edge_clade)}")

    built = build_mul_tree(bb, events)
    print(f"\nREBUILT tree: {sp} appears {sum(1 for l in built.get_leaves() if l.name==sp)}x")
    for i, sib in enumerate(copy_siblings(built, sp)):
        print(f"  copy {i} siblings: {sib}")

    print("\nCompare GT copy-siblings vs REBUILT copy-siblings above = the ret_sisters mismatch.")


if __name__ == "__main__":
    main()
