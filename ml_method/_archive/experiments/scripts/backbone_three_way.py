"""Three-way backbone comparison split by ASTRAL topology correctness.

Merges the per-sample scores from the three backbone builds (astral, rerooted,
true) and splits by whether ASTRAL had the exact true UNROOTED topology (RF=0).

The decisive question: on the RF=0 samples (ASTRAL topology is exactly right),
does the re-rooted build reach the true-backbone score?
  - rerooted ~ true on RF=0  -> rooting IS the fix for correct-topology cases;
    the residual gap is just ASTRAL's topology errors on the other samples.
  - rerooted still >> true on RF=0 -> the true-backbone advantage is structural
    (the GT MUL-tree was built from that exact backbone), i.e. partly artifactual,
    and improving ASTRAL will not reach it.

Run in the gene2net env (ete3).

Usage:
    python scripts/backbone_three_way.py \
        --astral-dir   output/reconstruct_allo/backbone_exp/ils_low_astral/t0.9 \
        --rerooted-dir output/reconstruct_allo/backbone_exp/ils_low_rerooted/t0.9 \
        --true-dir     output/reconstruct_allo/backbone_exp/ils_low_true/t0.9 \
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--astral-dir", required=True)
    parser.add_argument("--rerooted-dir", required=True)
    parser.add_argument("--true-dir", required=True)
    parser.add_argument("--mul-trees-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--metric", default="edit_distance_multree")
    args = parser.parse_args()

    m = args.metric

    def load(d):
        return pd.read_csv(os.path.join(d, "reconstruction_scores.csv"))[["sample", m]]

    a = load(args.astral_dir).rename(columns={m: "astral"})
    rr = load(args.rerooted_dir).rename(columns={m: "rerooted"})
    t = load(args.true_dir).rename(columns={m: "true"})
    df = a.merge(rr, on="sample").merge(t, on="sample")

    # per-sample unrooted RF (ASTRAL vs true backbone)
    rfs = []
    for _, r in df.iterrows():
        idx = r["sample"].replace("sample_", "")
        ap = os.path.join(args.mul_trees_dir, "simphy", args.config, idx,
                          f"replicate_{args.replicate}", "astral_species.tre")
        tp = os.path.join(args.mul_trees_dir, f"species_tree_{idx}.nex")
        rf = None
        if os.path.exists(ap) and os.path.exists(tp):
            try:
                ta = Tree(open(ap).read().strip(), format=1)
                tt = load_nexus_tree(tp)
                rf = ta.compare(tt, unrooted=True)["norm_rf"]
            except Exception:
                rf = None
        rfs.append(rf)
    df["unrooted_rf"] = rfs
    df = df.dropna(subset=["unrooted_rf"])

    def summary(sub, label):
        print(f"\n{label} (n={len(sub)}):")
        print(f"  astral   : {sub['astral'].mean():.3f}")
        print(f"  rerooted : {sub['rerooted'].mean():.3f}")
        print(f"  true     : {sub['true'].mean():.3f}")
        if len(sub):
            closed = (sub['astral'].mean() - sub['rerooted'].mean())
            total = (sub['astral'].mean() - sub['true'].mean())
            if total > 1e-9:
                print(f"  rerooting closes {100*closed/total:.0f}% of the astral->true gap")

    print(f"\nThree-way backbone comparison, metric={m}")
    summary(df, "ALL samples")
    summary(df[df["unrooted_rf"] == 0], "ASTRAL topology EXACT (unrooted RF=0)")
    summary(df[df["unrooted_rf"] > 0], "ASTRAL topology WRONG (unrooted RF>0)")

    print("\nVerdict:")
    print("  rerooted ~ true on the RF=0 rows -> rooting is the fix for correct-topology")
    print("  cases. rerooted still >> true there -> the true-backbone edge is structural")
    print("  (GT built from that backbone) and not reachable by improving ASTRAL.")


if __name__ == "__main__":
    main()
