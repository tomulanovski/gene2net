"""Check whether training samples align with metadata_{idx} / mul_tree_{idx}.

Read-only. For each index it compares three independent views of the polyploid
species:
  - metadata_{idx}.json  "polyploid_species"   (ground truth from packaging)
  - mul_tree_{idx}.nex   leaves seen >1 time   (ground truth from the MUL-tree)
  - sample copy counts   mode copies >= 2      (from the sample's own gene trees)

If meta == mul and the sample's copy-count polyploids are a subset of meta, the
index is aligned and any earlier guard failure was from the (noisy) label-derived
heuristic. If the sample's copy-count polyploids are NOT a subset of meta, the
sample/metadata indices genuinely do not correspond.

Usage:
    python scripts/check_metadata_alignment.py --subdir training/ils_low \
        --indices 0000 0011 0085 0500 1000
"""
import argparse
import json
import os
import pickle
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree
from gene2net_gnn.data.mul_tree_generator import get_polyploid_species


def load_tree_any(path):
    txt = open(path).read()
    if "#nexus" in txt.lower() or "begin trees" in txt.lower():
        for line in txt.splitlines():
            s = line.strip()
            if s.lower().startswith("tree") and "=" in s:
                return Tree(s.split("=", 1)[1].strip(), format=1)
    return Tree(txt.strip(), format=1)


def copy_count_polyploids(sample, min_mode=2):
    names = sample["species_tree_node_names"]
    is_leaf = sample["species_tree_is_leaf"]
    feats = sample["species_tree_node_features"]
    species = set(sample["species_list"])
    out = set()
    for i, nm in enumerate(names):
        if bool(is_leaf[i]) and nm in species and float(feats[i, 2]) >= min_mode:  # col 2 = mode_copies
            out.add(nm)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--subdir", default="training/ils_low")
    ap.add_argument("--indices", nargs="+", default=["0000", "0011", "0085", "0500", "1000"])
    args = ap.parse_args()

    for idx in args.indices:
        md = json.load(open(os.path.join(args.data_root, f"metadata_{idx}.json")))
        meta_poly = set(md.get("polyploid_species", {}))
        mul = load_tree_any(os.path.join(args.data_root, f"mul_tree_{idx}.nex"))
        mul_poly = set(get_polyploid_species(mul))
        with open(os.path.join(args.data_root, args.subdir, f"sample_{idx}", "sample.pkl"), "rb") as f:
            sample = pickle.load(f)
        cc_poly = copy_count_polyploids(sample)

        print(f"[{idx}] n_species(meta)={md.get('n_species')} n_species(sample)={len(sample['species_list'])}")
        print(f"      meta     = {sorted(meta_poly)}")
        print(f"      mul_tree = {sorted(mul_poly)}   meta==mul: {meta_poly == mul_poly}")
        print(f"      sample(copy>=2) = {sorted(cc_poly)}   sample subset of meta: {cc_poly.issubset(meta_poly)}")
        if not cc_poly.issubset(meta_poly):
            print(f"      !! sample-only species: {sorted(cc_poly - meta_poly)}")


if __name__ == "__main__":
    main()
