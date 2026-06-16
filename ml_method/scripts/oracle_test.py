"""Oracle reconstruction test: feed the TRUE WGD events through our build, to
separate 'is the build/folding faithful' from 'are the model predictions bad'.

For each simulated sample:
  - read the ground-truth MUL-tree (mul_tree_NNNN.nex),
  - decompose it into true WGD events (decompose_mul_tree),
  - rebuild a MUL-tree with build_mul_tree on a chosen backbone,
  - write output.tre + ground_truth for scoring.

Backbones:
  --backbone true   : the true species tree (species_tree_NNNN.nex). Isolates the
                      build/fold only. Should reproduce the ground truth (edit ~0)
                      if build_mul_tree is a faithful inverse of decompose.
  --backbone astral : the ASTRAL tree (what our real pipeline uses). The gap vs
                      'true' is the ASTRAL-backbone error; this is the ceiling our
                      real predictions could ever reach.

Then score with score_reconstructions.py. Interpretation:
  - oracle(true) ~ 0          -> build is faithful.
  - oracle(astral) low        -> our gap is the MODEL predictions (detection/partner).
  - oracle(astral) ~ 0.9      -> build/backbone is the limit, not predictions.

Run in the gene2net env (no torch needed).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree

from gene2net_gnn.data.label_extractor import decompose_mul_tree
from gene2net_gnn.inference.mul_tree_builder import build_mul_tree


def load_nexus_tree(path):
    content = open(path).read()
    for line in content.split("\n"):
        line = line.strip()
        if line.lower().startswith("tree") and "=" in line:
            return Tree(line.split("=", 1)[1].strip(), format=1)
    raise ValueError(f"No tree found in {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mul-trees-dir", required=True)
    parser.add_argument("--config", required=True, help="e.g. ils_low")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--backbone", choices=["true", "astral"], default="true")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    done = skipped = 0
    for idx in range(args.start, args.start + args.n):
        idx_str = f"{idx:04d}"
        mul_path = os.path.join(args.mul_trees_dir, f"mul_tree_{idx_str}.nex")
        if not os.path.exists(mul_path):
            skipped += 1
            continue

        if args.backbone == "true":
            bb_path = os.path.join(args.mul_trees_dir, f"species_tree_{idx_str}.nex")
            if not os.path.exists(bb_path):
                skipped += 1
                continue
            backbone = load_nexus_tree(bb_path)
        else:
            bb_path = os.path.join(args.mul_trees_dir, "simphy", args.config, idx_str,
                                   f"replicate_{args.replicate}", "astral_species.tre")
            if not os.path.exists(bb_path):
                skipped += 1
                continue
            backbone = Tree(open(bb_path).read().strip(), format=1)

        gt_mul = load_nexus_tree(mul_path)
        try:
            events = decompose_mul_tree(gt_mul)
            oracle = build_mul_tree(backbone, events)
        except Exception as e:
            print(f"  SKIP [{idx_str}]: {type(e).__name__}: {e}")
            skipped += 1
            continue

        case_dir = os.path.join(args.out_dir, f"sample_{idx_str}")
        os.makedirs(case_dir, exist_ok=True)
        oracle.write(outfile=os.path.join(case_dir, "output.tre"), format=9)
        with open(mul_path) as f:
            gt = f.read()
        with open(os.path.join(case_dir, "ground_truth.nex"), "w") as f:
            f.write(gt)

        print(f"[{idx_str}] events={len(events)} -> {case_dir}/output.tre")
        done += 1

    print(f"\nDone: {done} oracle trees ({args.backbone} backbone), {skipped} skipped.")
    print(f"Score: python scripts/score_reconstructions.py --recon-dir {args.out_dir} --workers 8")


if __name__ == "__main__":
    main()
