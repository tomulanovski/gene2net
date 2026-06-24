"""Test whether the ASTRAL-vs-true backbone edit-distance gap is a ROOTING issue.

backbone_error_analysis showed the gap persists even when ASTRAL has the exact
true UNROOTED topology (rf=0). get_edit_distance_multree is rooting-sensitive (it
builds a directed graph), and ASTRAL output is unrooted while the true backbone is
rooted at the real root. So the suspect is rooting, not topology or branch length.

For each sample we compute the normalized RF between ASTRAL and the true backbone
two ways: unrooted and rooted. If many samples are unrooted-identical (rf=0) but
rooted-different (rf>0), and rooted RF tracks the edit gap, the gap is a rooting
artifact, and the cheap fix is to root ASTRAL correctly before building.

Run in the gene2net env (ete3 works there).

Usage:
    python scripts/backbone_rooting_check.py \
        --astral-dir output/reconstruct_allo/backbone_exp/ils_low_astral/t0.9 \
        --true-dir   output/reconstruct_allo/backbone_exp/ils_low_true/t0.9 \
        --mul-trees-dir data/mul_trees_2k --config ils_low
"""
import argparse
import os

import pandas as pd
from ete3 import Tree


def load_nexus_tree(path):
    for line in open(path).read().split("\n"):
        line = line.strip()
        if line.lower().startswith("tree") and "=" in line:
            return Tree(line.split("=", 1)[1].strip(), format=1)
    raise ValueError(f"No tree in {path}")


def rf_both(ta, tt):
    """Return (unrooted_norm_rf, rooted_norm_rf)."""
    try:
        u = ta.compare(tt, unrooted=True)["norm_rf"]
    except Exception:
        u = None
    try:
        r = ta.compare(tt, unrooted=False)["norm_rf"]
    except Exception:
        r = None
    return u, r


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--astral-dir", required=True)
    parser.add_argument("--true-dir", required=True)
    parser.add_argument("--mul-trees-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--metric", default="edit_distance_multree")
    args = parser.parse_args()

    a = pd.read_csv(os.path.join(args.astral_dir, "reconstruction_scores.csv"))
    t = pd.read_csv(os.path.join(args.true_dir, "reconstruction_scores.csv"))
    m = args.metric
    df = a[["sample", m]].merge(t[["sample", m]], on="sample", suffixes=("_astral", "_true"))
    df["gap"] = df[f"{m}_astral"] - df[f"{m}_true"]

    rows = []
    for _, r in df.iterrows():
        idx = r["sample"].replace("sample_", "")
        ap = os.path.join(args.mul_trees_dir, "simphy", args.config, idx,
                          f"replicate_{args.replicate}", "astral_species.tre")
        tp = os.path.join(args.mul_trees_dir, f"species_tree_{idx}.nex")
        if not (os.path.exists(ap) and os.path.exists(tp)):
            continue
        try:
            ta = Tree(open(ap).read().strip(), format=1)
            tt = load_nexus_tree(tp)
            u, rr = rf_both(ta, tt)
        except Exception:
            continue
        if u is None or rr is None:
            continue
        rows.append({"sample": r["sample"], "unrooted_rf": u, "rooted_rf": rr,
                     "gap": r["gap"]})

    res = pd.DataFrame(rows)
    if res.empty:
        print("No samples — check paths.")
        return

    print(f"\nRooting check on {len(res)} samples, metric={m}\n")
    print(f"mean unrooted RF: {res['unrooted_rf'].mean():.3f}")
    print(f"mean rooted   RF: {res['rooted_rf'].mean():.3f}")
    print(f"corr(unrooted RF, gap) = {res['unrooted_rf'].corr(res['gap']):.3f}")
    print(f"corr(rooted   RF, gap) = {res['rooted_rf'].corr(res['gap']):.3f}")

    same_unrooted = res[res["unrooted_rf"] == 0]
    print(f"\nUnrooted-identical (unrooted_rf=0): {len(same_unrooted)}/{len(res)}")
    if len(same_unrooted):
        n_rooted_diff = (same_unrooted["rooted_rf"] > 0).sum()
        print(f"  of those, rooted-DIFFERENT (rooted_rf>0): {n_rooted_diff}/{len(same_unrooted)} "
              f"({100*n_rooted_diff/len(same_unrooted):.0f}%)")
        print(f"  their mean rooted RF: {same_unrooted['rooted_rf'].mean():.3f}")
        print(f"  their mean edit gap:  {same_unrooted['gap'].mean():.3f}")
        print("\n  -> if most are rooted-different and rooted RF tracks the gap,")
        print("     the gap is a ROOTING artifact and the fix is to re-root ASTRAL.")


if __name__ == "__main__":
    main()
