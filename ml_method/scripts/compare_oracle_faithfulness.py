"""Expose build/scoring convention mismatches on an ORACLE sample.

For an oracle sample (true backbone + true events), the reconstruction SHOULD equal
the ground truth. Any nonzero ret_sisters/edit is therefore a build+fold convention
artifact, not a prediction error. This script folds the built output.tre and the
ground_truth.nex the same way the scorer does (ReticulateTree), then prints, for each
reticulation, its leaf-set (descendants) and its two sister clades side by side, and
flags where the sisters differ.

Reading the output: if the reticulation LEAVES match (same descendants) but the
SISTERS differ, the build is attaching the WGD copy with different local structure
than the ground truth writes it -> a convention bug that penalizes every config.

Run in the gene2net env (needs the comparison core).

Usage:
  # diff one sample:
  python scripts/compare_oracle_faithfulness.py --sample-dir output/oracle/true_ils_low/sample_0002
  # or auto-pick the simplest sample with a nonzero ret_sisters:
  python scripts/compare_oracle_faithfulness.py --recon-dir output/oracle/true_ils_low --auto
"""
import argparse
import os
import sys


def newick_from_file(path):
    txt = open(path).read().strip()
    low = txt.lower()
    if low.startswith("#nexus") or "begin trees" in low:
        for line in txt.splitlines():
            s = line.strip()
            if s.lower().startswith("tree") and "=" in s:
                return s.split("=", 1)[1].strip()
    return txt


def build_rt(tree_str, ReticulateTree):
    rt = ReticulateTree(tree_str)
    if rt.check_duplicated():
        return ReticulateTree(tree_str, is_multree=True)
    return rt


def fold_sample(sdir, ReticulateTree):
    gt = newick_from_file(os.path.join(sdir, "ground_truth.nex"))
    inf = newick_from_file(os.path.join(sdir, "output.tre"))
    return build_rt(gt, ReticulateTree), build_rt(inf, ReticulateTree)


def summarize(rt):
    """{frozenset(leaves): [set,set] sisters} keyed by reticulation leaf-set."""
    leaves = rt.get_reticulation_leaves()
    sisters = rt.get_reticulation_sisters()
    out = {}
    for rid, leaf_list in leaves.items():
        out[frozenset(leaf_list)] = sisters.get(rid, [])
    return out


def show(sdir, ReticulateTree):
    rt_gt, rt_inf = fold_sample(sdir, ReticulateTree)
    gt = summarize(rt_gt)
    inf = summarize(rt_inf)
    print(f"\n=== {os.path.basename(sdir)} ===")
    print(f"reticulations: ground-truth={len(gt)}  built={len(inf)}")

    matched = 0
    for i, (leafset, gt_sis) in enumerate(sorted(gt.items(), key=lambda kv: len(kv[0]))):
        short = sorted(leafset)
        short = short if len(short) <= 6 else short[:6] + [f"...(+{len(short)-6})"]
        print(f"\n[ret {i}] descendants ({len(leafset)}): {short}")
        inf_sis = inf.get(leafset)
        if inf_sis is None:
            print("  built: NO reticulation with these exact descendants (leaf-set differs)")
            continue
        matched += 1
        gt_sides = sorted([tuple(sorted(s)) for s in gt_sis])
        inf_sides = sorted([tuple(sorted(s)) for s in inf_sis])
        same = gt_sides == inf_sides
        print(f"  ground-truth sisters: {[list(s) for s in gt_sides]}")
        print(f"  built        sisters: {[list(s) for s in inf_sides]}")
        print(f"  -> sisters {'MATCH' if same else 'MISMATCH  <-- convention artifact'}")
    print(f"\nsummary: {matched}/{len(gt)} reticulations leaf-matched; "
          f"look for MISMATCH lines above (identical descendants, different sisters).")


def auto_pick(recon_dir, ReticulateTree, max_try=40):
    import pandas as pd
    csv = os.path.join(recon_dir, "reconstruction_scores.csv")
    df = pd.read_csv(csv)
    cand = df[(df["ret_sisters_jaccard"] > 0.05) & (df["ret_sisters_jaccard"] < 0.6)]
    cand = cand.sort_values("ret_sisters_jaccard")
    best = None
    for _, r in cand.head(max_try).iterrows():
        sdir = os.path.join(recon_dir, str(r["sample"]))
        if not os.path.isdir(sdir):
            continue
        try:
            rt_gt, _ = fold_sample(sdir, ReticulateTree)
            n = rt_gt.get_reticulation_count()
        except Exception:
            continue
        if 1 <= n <= 2:
            return sdir
        if best is None:
            best = sdir
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir")
    ap.add_argument("--recon-dir")
    ap.add_argument("--auto", action="store_true")
    ap.add_argument("--sim-scripts",
                    default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts")
    args = ap.parse_args()

    if args.sim_scripts not in sys.path:
        sys.path.insert(0, args.sim_scripts)
    from reticulate_tree import ReticulateTree

    if args.sample_dir:
        show(args.sample_dir, ReticulateTree)
    elif args.recon_dir and args.auto:
        sdir = auto_pick(args.recon_dir, ReticulateTree)
        if not sdir:
            print("No suitable sample found (nonzero ret_sisters).")
            return
        print(f"auto-picked simplest mismatching sample: {sdir}")
        show(sdir, ReticulateTree)
    else:
        ap.error("give --sample-dir, or --recon-dir with --auto")


if __name__ == "__main__":
    main()
