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
from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.inference.mul_tree_builder import build_mul_tree, WGDEvent


def preorder_clades(tree):
    return [frozenset(n.get_leaf_names()) for n in tree.traverse("preorder") if not n.is_root()]


def events_from_labels(labels, clades):
    """True events expressed as ASTRAL edges (the labels) -> WGDEvents using the
    ASTRAL clades. These always build cleanly (no exact-match dropping)."""
    events = []
    mask = labels.mask or []
    for k in range(len(labels.wgd_edges)):
        if k < len(mask) and not mask[k]:
            continue
        w = labels.wgd_edges[k]
        p = labels.partner_edges[k] if k < len(labels.partner_edges) else w
        if 0 <= w < len(clades) and 0 <= p < len(clades):
            events.append(WGDEvent(wgd_edge_clade=clades[w], partner_edge_clade=clades[p], confidence=1.0))
    return events


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
    parser.add_argument("--events", choices=["decompose", "labels"], default="decompose",
                        help="decompose: true-tree clades (exact match). labels: true events "
                             "mapped to ASTRAL edges (the fair ceiling; implies --backbone astral).")
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

        gt_mul = load_nexus_tree(mul_path)
        try:
            if args.events == "labels":
                # Fair ceiling: true events mapped to ASTRAL edges, built on ASTRAL.
                astral_path = os.path.join(args.mul_trees_dir, "simphy", args.config, idx_str,
                                           f"replicate_{args.replicate}", "astral_species.tre")
                sample_dir = os.path.join(args.mul_trees_dir, "training", args.config, f"sample_{idx_str}")
                if not (os.path.exists(astral_path) and os.path.isdir(sample_dir)):
                    skipped += 1
                    continue
                backbone = Tree(open(astral_path).read().strip(), format=1)
                sample = Gene2NetSample.load(sample_dir)
                if sample.labels is None:
                    skipped += 1
                    continue
                events = events_from_labels(sample.labels, preorder_clades(backbone))
            else:
                if args.backbone == "true":
                    bb_path = os.path.join(args.mul_trees_dir, f"species_tree_{idx_str}.nex")
                else:
                    bb_path = os.path.join(args.mul_trees_dir, "simphy", args.config, idx_str,
                                           f"replicate_{args.replicate}", "astral_species.tre")
                if not os.path.exists(bb_path):
                    skipped += 1
                    continue
                backbone = (load_nexus_tree(bb_path) if args.backbone == "true"
                            else Tree(open(bb_path).read().strip(), format=1))
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
