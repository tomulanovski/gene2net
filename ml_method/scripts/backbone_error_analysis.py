"""Decide HOW to attack the backbone: is the ASTRAL-vs-true edit-distance gap
explained by ASTRAL's topology error, and is that error small/localized?

We have two per-sample score files from the backbone 2x2 (same predicted events,
different backbone):
    <astral-dir>/reconstruction_scores.csv   (built on ASTRAL)
    <true-dir>/reconstruction_scores.csv     (built on the true backbone)

For each sample we add the normalized Robinson-Foulds distance between the ASTRAL
tree and the true diploid backbone, then ask:

  - Among samples where ASTRAL == true backbone (RF = 0), is the edit gap ~0?
    If yes, ALL the gap is the ASTRAL topology, so a better backbone fully closes
    it. If the gap is large even at RF = 0, the build itself adds error.
  - Does the edit gap scale with RF? Strong correlation + usually-small RF ->
    a few wrong edges drive it -> locally refining ASTRAL (route 1) is high-yield.
    Weak/diffuse -> the star-tree / joint rebuild (route 2) is needed.

Run in the final_project env (needs ete3 + pandas).

Usage:
    python scripts/backbone_error_analysis.py \
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


def norm_rf(t_astral, t_true):
    """Normalized RF between two trees on the same leaf set; None if incomparable."""
    try:
        res = t_astral.compare(t_true, unrooted=True)
        return res["norm_rf"], int(res["rf"]), int(res["max_rf"])
    except Exception:
        return None, None, None


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
        astral_path = os.path.join(args.mul_trees_dir, "simphy", args.config, idx,
                                   f"replicate_{args.replicate}", "astral_species.tre")
        true_path = os.path.join(args.mul_trees_dir, f"species_tree_{idx}.nex")
        if not (os.path.exists(astral_path) and os.path.exists(true_path)):
            continue
        try:
            ta = Tree(open(astral_path).read().strip(), format=1)
            tt = load_nexus_tree(true_path)
            nrf, rf, maxrf = norm_rf(ta, tt)
        except Exception:
            continue
        if nrf is None:
            continue
        rows.append({"sample": r["sample"], "norm_rf": nrf, "rf": rf,
                     f"{m}_astral": r[f"{m}_astral"], f"{m}_true": r[f"{m}_true"],
                     "gap": r["gap"]})

    res = pd.DataFrame(rows)
    if res.empty:
        print("No samples with both trees found — check paths.")
        return

    print(f"\nBackbone error analysis on {len(res)} samples, metric={m}\n")
    print(f"ASTRAL edit (mean): {res[f'{m}_astral'].mean():.3f}")
    print(f"TRUE   edit (mean): {res[f'{m}_true'].mean():.3f}")
    print(f"gap   (mean):       {res['gap'].mean():.3f}")
    print(f"ASTRAL norm_rf vs true backbone (mean): {res['norm_rf'].mean():.3f} "
          f"(median {res['norm_rf'].median():.3f})\n")

    corr = res["norm_rf"].corr(res["gap"])
    print(f"correlation(norm_rf, edit gap) = {corr:.3f}")

    # samples where ASTRAL is topologically correct
    exact = res[res["rf"] == 0]
    print(f"\nASTRAL exactly correct (rf=0): {len(exact)}/{len(res)} samples")
    if len(exact):
        print(f"  their mean ASTRAL edit: {exact[f'{m}_astral'].mean():.3f}")
        print(f"  their mean TRUE   edit: {exact[f'{m}_true'].mean():.3f}")
        print(f"  their mean gap:         {exact['gap'].mean():.3f}")
        print("  -> if this gap ~0, the whole gap is ASTRAL topology (backbone fully "
              "explains it). If >0, the build itself adds error even on a correct tree.")

    # gap by RF bucket
    print("\nedit gap by ASTRAL error level:")
    res["rf_bucket"] = pd.cut(res["norm_rf"], [-0.001, 0.0, 0.05, 0.1, 0.2, 1.0],
                              labels=["0", "(0,0.05]", "(0.05,0.1]", "(0.1,0.2]", ">0.2"])
    g = res.groupby("rf_bucket", observed=True).agg(
        n=("gap", "size"), mean_gap=("gap", "mean"),
        mean_astral=(f"{m}_astral", "mean"))
    print(g.round(3).to_string())

    print("\nDecision guide:")
    print("  high corr + small typical RF + gap~0 at rf=0  -> route 1 (refine ASTRAL).")
    print("  gap large even at low RF, or weak corr        -> route 2 (joint/star rebuild).")


if __name__ == "__main__":
    main()
