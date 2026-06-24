"""Find out WHY the astral-built and true-built reconstructions score so
differently against the GT, even on samples where ASTRAL has the exact topology.

Internal nodes match freely in get_edit_distance_multree (label None == None), so
the 'scaffold artifact' theory is suspect. This inspects the actual built trees
(already saved, no model needed): for RF=0 samples it reports, per sample,
  edit(astral_recon, GT), edit(true_recon, GT), edit(astral_recon, true_recon),
plus node/leaf counts.

If edit(astral_recon, true_recon) is large -> the two reconstructions are
genuinely different graphs (real structural difference to chase). If ~0 -> they
are the same tree and the GT-comparison difference is a scoring inconsistency.

Run in the gene2net env.

Usage:
    python scripts/backbone_inspect_sample.py \
        --astral-dir   output/reconstruct_allo/backbone_exp/ils_low_astral/t0.9 \
        --true-dir     output/reconstruct_allo/backbone_exp/ils_low_true/t0.9 \
        --mul-trees-dir data/mul_trees_2k --config ils_low --n-show 8
"""
import argparse
import os
import sys

from ete3 import Tree

SIM = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts"
if SIM not in sys.path:
    sys.path.insert(0, SIM)
from reticulate_tree import ReticulateTree
from compare_reticulations import pairwise_compare


def newick_from_file(path):
    txt = open(path).read().strip()
    if txt.lower().startswith("#nexus") or "begin trees" in txt.lower():
        for line in txt.splitlines():
            s = line.strip()
            if s.lower().startswith("tree") and "=" in s:
                return s.split("=", 1)[1].strip()
    return txt


def load_nexus_tree(path):
    for line in open(path).read().split("\n"):
        line = line.strip()
        if line.lower().startswith("tree") and "=" in line:
            return Tree(line.split("=", 1)[1].strip(), format=1)
    raise ValueError(f"No tree in {path}")


def build_rt(tree_str):
    rt = ReticulateTree(tree_str)
    if rt.check_duplicated():
        return ReticulateTree(tree_str, is_multree=True)
    return rt


def counts(newick):
    t = Tree(newick, format=9)
    leaves = t.get_leaves()
    total = len(list(t.traverse()))
    return len(leaves), total


def edit(a_str, b_str):
    try:
        comp = pairwise_compare(build_rt(a_str), build_rt(b_str))
        return comp.get("edit_distance_multree")
    except Exception as e:
        return f"ERR:{type(e).__name__}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--astral-dir", required=True)
    parser.add_argument("--true-dir", required=True)
    parser.add_argument("--mul-trees-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--n-show", type=int, default=8)
    args = parser.parse_args()

    samples = sorted(d for d in os.listdir(args.astral_dir)
                     if os.path.isdir(os.path.join(args.astral_dir, d)))

    shown = 0
    for s in samples:
        idx = s.replace("sample_", "")
        ap = os.path.join(args.mul_trees_dir, "simphy", args.config, idx,
                          f"replicate_{args.replicate}", "astral_species.tre")
        tp = os.path.join(args.mul_trees_dir, f"species_tree_{idx}.nex")
        if not (os.path.exists(ap) and os.path.exists(tp)):
            continue
        try:
            rf = Tree(open(ap).read().strip(), format=1).compare(
                load_nexus_tree(tp), unrooted=True)["norm_rf"]
        except Exception:
            continue
        if rf != 0:
            continue  # only RF=0 samples

        a_out = os.path.join(args.astral_dir, s, "output.tre")
        t_out = os.path.join(args.true_dir, s, "output.tre")
        gt = os.path.join(args.astral_dir, s, "ground_truth.nex")
        if not (os.path.exists(a_out) and os.path.exists(t_out) and os.path.exists(gt)):
            continue

        a_str = newick_from_file(a_out)
        t_str = newick_from_file(t_out)
        g_str = newick_from_file(gt)

        al, an = counts(a_str)
        tl, tn = counts(t_str)
        gl, gn = counts(g_str)

        print(f"\n=== {s} (unrooted RF=0) ===")
        print(f"  leaves/nodes  astral_recon={al}/{an}  true_recon={tl}/{tn}  GT={gl}/{gn}")
        print(f"  edit(astral_recon, GT)        = {edit(a_str, g_str)}")
        print(f"  edit(true_recon,   GT)        = {edit(t_str, g_str)}")
        print(f"  edit(astral_recon, true_recon)= {edit(a_str, t_str)}")

        shown += 1
        if shown >= args.n_show:
            break

    print("\nIf edit(astral_recon, true_recon) is large -> genuinely different trees.")
    print("If ~0 -> same tree, and the GT-comparison gap is a scoring inconsistency.")


if __name__ == "__main__":
    main()
