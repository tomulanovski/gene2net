"""Diagnose whether the ASTRAL backbone is the problem: compare each ASTRAL
species tree to the TRUE species tree (the backbone the MUL-tree was built on),
and report the normalized Robinson-Foulds distance.

If RF is high, ASTRAL is genuinely off on these (polyploid, duplicate-label)
gene trees -> that's the edit_distance limiter. If RF is ~0, the oracle's 0.857
is a build/clade-matching artifact, not ASTRAL.

Run in any env with ete3.
"""
import argparse
import os
import sys

from ete3 import Tree


def load_nexus_tree(path):
    for line in open(path).read().split("\n"):
        line = line.strip()
        if line.lower().startswith("tree") and "=" in line:
            return Tree(line.split("=", 1)[1].strip(), format=1)
    raise ValueError(f"No tree in {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mul-trees-dir", required=True)
    parser.add_argument("--config", required=True, help="e.g. ils_low")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--replicate", type=int, default=1)
    args = parser.parse_args()

    norm_rfs = []
    for idx in range(args.start, args.start + args.n):
        idx_str = f"{idx:04d}"
        true_path = os.path.join(args.mul_trees_dir, f"species_tree_{idx_str}.nex")
        astral_path = os.path.join(args.mul_trees_dir, "simphy", args.config, idx_str,
                                   f"replicate_{args.replicate}", "astral_species.tre")
        if not (os.path.exists(true_path) and os.path.exists(astral_path)):
            continue
        true_t = load_nexus_tree(true_path)
        astral_t = Tree(open(astral_path).read().strip(), format=1)

        # restrict to shared leaves, just in case
        shared = set(true_t.get_leaf_names()) & set(astral_t.get_leaf_names())
        if len(shared) < 4:
            continue
        try:
            res = true_t.compare(astral_t, unrooted=True)
            norm_rfs.append(res["norm_rf"])
        except Exception as e:
            print(f"  [{idx_str}] compare failed: {e}")

    if not norm_rfs:
        print("No comparable samples found.")
        return
    mean = sum(norm_rfs) / len(norm_rfs)
    norm_rfs.sort()
    median = norm_rfs[len(norm_rfs) // 2]
    print(f"Compared {len(norm_rfs)} samples")
    print(f"normalized RF (ASTRAL vs true species tree): mean={mean:.3f}  median={median:.3f}  "
          f"min={norm_rfs[0]:.3f}  max={norm_rfs[-1]:.3f}")
    print("\nHigh (>0.3) -> ASTRAL topology is genuinely off (the limiter).")
    print("Low  (<0.1) -> ASTRAL is fine; the oracle 0.857 is a build/clade-match artifact.")


if __name__ == "__main__":
    main()
