"""Repair taxon labels in existing reconstructions, so they can be scored again.

Two defects made the mu-distance undefined on some networks and silently inflated the
reticulation Jaccards on the same ones, because a mismatched label is simply a different
taxon to a Jaccard.

1. output.tre kept the substring-fixed names. taxa_map.txt records SUBSTRING replacements
   (fix_substrings_grampa.py rewrites a label that is a substring of another), but
   benchmark_networks.rename_leaves matched whole names only, so a leaf like "filixX-mas"
   was never renamed back to "filix-mas". Fixed in benchmark_networks.py for future runs.

2. ground_truth.nex is a copy of simulations/networks/<net>.tre taken at inference time.
   Reconstructions produced before that file was corrected still hold the old spelling,
   for instance "jalapaensis" where the ground truth now reads "jalapaënsis".

This repairs both in place so the existing reconstructions can be re-scored without
re-running inference. Reports by default; pass --apply to write.

Run from ml_method/:
    python scripts/repair_labels.py
    python scripts/repair_labels.py --apply
"""
import argparse
import os
import sys

from ete3 import Tree

BASE = "output/reconstruct_final/final"
SIM = os.path.join("..", "simulations", "simulations")
NETWORKS = os.path.join("..", "simulations", "networks")
STRATEGIES = ("bound_driven", "detect_only")


def inverse_taxa_map(path):
    """REPLACEMENT -> ORIGINAL from taxa_map.txt.

    prep_grampa writes this file only when a label collided as a substring of
    another, so its absence means no renaming is needed for that network and an
    empty map is correct. A malformed line still raises.
    """
    if not os.path.exists(path):
        return {}
    inv = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"malformed line in {path}: {line!r}")
            inv[parts[1].strip()] = parts[0].strip()
    return inv


def repair_output(tre_path, inv, apply_changes):
    tree = Tree(tre_path, format=9)
    keys = sorted(inv, key=len, reverse=True)
    changed = []
    for leaf in tree.get_leaves():
        if leaf.name in inv:
            changed.append((leaf.name, inv[leaf.name]))
            leaf.name = inv[leaf.name]
            continue
        for k in keys:
            if k in leaf.name:
                new = leaf.name.replace(k, inv[k], 1)
                changed.append((leaf.name, new))
                leaf.name = new
                break
    if changed and apply_changes:
        tree.write(outfile=tre_path, format=9)
    return changed


def repair_ground_truth(gt_copy, net, apply_changes):
    src = os.path.join(NETWORKS, f"{net}.tre")
    if not os.path.exists(src):
        raise FileNotFoundError(f"reference network not found: {src}")
    want = open(src, encoding="utf-8").read()
    have = open(gt_copy, encoding="utf-8").read() if os.path.exists(gt_copy) else None
    if have == want:
        return False
    if apply_changes:
        with open(gt_copy, "w", encoding="utf-8") as f:
            f.write(want)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="write the repairs (default: report only)")
    args = ap.parse_args()

    if not os.path.isdir(BASE):
        raise SystemExit(f"not found: {BASE} (run from ml_method/)")

    n_out, n_gt, seen = 0, 0, set()
    for config in sorted(os.listdir(BASE)):
        for rep in range(1, 6):
            for strat in STRATEGIES:
                d = os.path.join(BASE, config, f"rep{rep}", strat)
                if not os.path.isdir(d):
                    continue
                for net in sorted(os.listdir(d)):
                    net_dir = os.path.join(d, net)
                    tre = os.path.join(net_dir, "output.tre")
                    if not os.path.isdir(net_dir) or not os.path.exists(tre):
                        continue
                    tmap = os.path.join(SIM, net, "processed", config,
                                        "grampa_input", f"replicate_{rep}", "taxa_map.txt")
                    inv = inverse_taxa_map(tmap)
                    for old, new in repair_output(tre, inv, args.apply):
                        n_out += 1
                        if (net, old, new) not in seen:
                            seen.add((net, old, new))
                            print(f"  label  {net:22} {old!r} -> {new!r}")
                    if repair_ground_truth(os.path.join(net_dir, "ground_truth.nex"),
                                           net, args.apply):
                        n_gt += 1
                        if (net, "gt") not in seen:
                            seen.add((net, "gt"))
                            print(f"  truth  {net:22} stale copy refreshed from networks/")

    verb = "repaired" if args.apply else "would repair"
    print(f"\n{verb}: {n_out} leaf labels, {n_gt} ground-truth copies")
    if not args.apply and (n_out or n_gt):
        print("re-run with --apply to write, then re-run jobs/finalize_score.sh")


if __name__ == "__main__":
    main()
