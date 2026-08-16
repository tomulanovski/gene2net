"""Detection-vs-multiset diagnostic: does the detection head fire on true events
that the inferred multiset MISSED (copies deleted by fractionation)?

This is the measurement that decides whether a detection-driven decode is worth
building. It runs the trained model on each benchmark network under a (retention)
config, and for every species compares its detection signal across three groups.

For each network:
  - run the model -> wgd_prob per ASTRAL edge
  - inferred copy bound per species (the multiset: infer_copy_bound, or Polyphest's
    multi_set.txt with --copy-bound multiset)
  - true copy number per species from the ground-truth MUL-tree

Classify each species by (true copies, inferred bound):
  - CAUGHT polyploid:  true > 1 AND bound >= 2   (multiset sees the duplication)
  - MISSED polyploid:  true > 1 AND bound == 1   (fractionation hid it from copies)
  - DIPLOID control:   true == 1                 (no duplication)

Detection signal for a species = max wgd_prob over the ASTRAL edges whose clade
contains that species (the strongest WGD call anywhere on its lineage).

Reading:
  - MISSED detection >> DIPLOID detection  -> the detection head has signal beyond
    copy count; a detection-driven decode can recover fractionated events the
    multiset lost, and beat a ploidy-only method (Polyphest) on fractionated data.
  - MISSED ~ DIPLOID                        -> detection is copy-bound too; nothing
    recovers them (the honest information-theoretic limit of fractionation).

Run in the final_project env (needs torch).

Usage:
    python scripts/detection_signal_diagnostic.py \
        --model-dir output/reconstruct_final \
        --config conf_dup_loss_medium_10M_ne1M_fix050 \
        --model-config configs/reconstruct_final.yaml
"""
import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import yaml
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.inference.build_strategies import infer_copy_bound
from scripts.reconstruct_mul_tree import load_model, model_inputs_for, preorder_edge_clades
from scripts.benchmark_networks import (
    NETWORKS, load_gene_trees, load_inverse_taxa_map, copy_bound_from_multiset,
)


def true_copies_from_gt(gt_path, forward):
    """Per-species true copy count from the ground-truth MUL-tree, renamed from
    original to the substring-fixed (replacement) names via forward."""
    t = Tree(open(gt_path).read().strip(), format=1)
    c = Counter(t.get_leaf_names())
    return {forward.get(name, name): cnt for name, cnt in c.items()}


def summarize(xs):
    if not xs:
        return "n=0"
    xs2 = sorted(xs)
    n = len(xs2)
    return (f"n={n:>4}  mean={sum(xs)/n:.3f}  median={xs2[n//2]:.3f}  "
            f"p25={xs2[n//4]:.3f}  p75={xs2[(3*n)//4]:.3f}  "
            f"frac>=0.5={sum(1 for x in xs if x >= 0.5)/n:.2f}")


def main():
    ap = argparse.ArgumentParser(description="Detection-vs-multiset diagnostic.")
    ap.add_argument("--model-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--model-config", default=None)
    ap.add_argument("--copy-bound", choices=["infer", "multiset"], default="infer",
                    help="Which multiset to classify against: our infer_copy_bound, "
                         "or Polyphest's exact multi_set.txt.")
    ap.add_argument("--sim-base",
                    default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations")
    ap.add_argument("--networks-dir",
                    default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks")
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--max-gene-trees", type=int, default=500)
    args = ap.parse_args()

    base = os.path.join(os.path.dirname(__file__), "..")
    cfg_path = args.model_config or os.path.join(base, "configs", "reconstruct_final.yaml")
    model_config = yaml.safe_load(open(cfg_path)).get("model", {})
    # Feature ablations must match how the model was trained (as in benchmark_networks).
    from gene2net_gnn.training.trainer_reconstruct import set_feature_opts
    set_feature_opts(
        coclust_condition_on_dup=model_config.get("coclust_condition_on_dup", False),
        use_n_eff=model_config.get("use_n_eff", False),
    )
    device = "cpu"
    model = load_model(args.model_dir, model_config, device)

    caught, missed, diploid = [], [], []
    n_net = 0
    for net in NETWORKS:
        rep_dir = os.path.join(args.sim_base, net, "processed", args.config,
                               "grampa_input", f"replicate_{args.replicate}")
        gtp = os.path.join(rep_dir, "clean_trees.tre")
        asp = os.path.join(rep_dir, "species.tre")
        gt_path = os.path.join(args.networks_dir, f"{net}.tre")
        if not (os.path.exists(gtp) and os.path.exists(asp) and os.path.exists(gt_path)):
            continue
        try:
            gene_trees = load_gene_trees(gtp, args.max_gene_trees)
            astral = Tree(open(asp).read().strip(), format=1)
        except Exception as e:
            print(f"  skip {net}: parse error ({type(e).__name__}: {e})")
            continue
        species_list = sorted(set(astral.get_leaf_names()))
        inv_map = load_inverse_taxa_map(os.path.join(rep_dir, "taxa_map.txt"))
        forward = {orig: repl for repl, orig in inv_map.items()}

        from gene2net_gnn.data.rooting import hybrid_root
        astral = hybrid_root(astral, gene_trees)

        sample = Gene2NetSample.from_trees(
            species_tree=astral, gene_trees=gene_trees, species_list=species_list)
        clades = preorder_edge_clades(astral)
        with torch.no_grad():
            inputs = model_inputs_for(sample, device)
            wgd_logits, _ = model(**inputs)
            wgd_prob = torch.softmax(wgd_logits, dim=-1)[:, 1].tolist()
        if len(clades) != len(wgd_prob):
            print(f"  skip {net}: clade/edge mismatch {len(clades)} vs {len(wgd_prob)}")
            continue

        if args.copy_bound == "multiset":
            msp = os.path.join(rep_dir.replace("grampa_input", "polyphest_input"),
                               "multi_set.txt")
            if not os.path.exists(msp):
                print(f"  skip {net}: multiset missing {msp}")
                continue
            copy_bound = copy_bound_from_multiset(msp, inv_map)
        else:
            copy_bound = infer_copy_bound(gene_trees)

        true_copies = true_copies_from_gt(gt_path, forward)

        # Precompute, per species, the max wgd_prob over edges whose clade contains it.
        for s in species_list:
            tc = int(true_copies.get(s, 1))
            b = int(copy_bound.get(s, 1))
            det = max((wgd_prob[i] for i in range(len(clades)) if s in clades[i]),
                      default=0.0)
            if tc > 1 and b >= 2:
                caught.append(det)
            elif tc > 1 and b == 1:
                missed.append(det)
            elif tc == 1:
                diploid.append(det)
        n_net += 1

    print(f"\n=== Detection-vs-multiset diagnostic: {args.config} "
          f"({n_net} networks, copy-bound={args.copy_bound}) ===\n")
    print(f"CAUGHT polyploids (true>1, bound>=2): {summarize(caught)}")
    print(f"MISSED polyploids (true>1, bound==1): {summarize(missed)}")
    print(f"DIPLOID control   (true==1)         : {summarize(diploid)}")
    print("\nReading: if MISSED >> DIPLOID, detection fires on fractionated events the")
    print("multiset lost -> a detection-driven decode can recover them and beat a")
    print("ploidy-only method. If MISSED ~ DIPLOID, detection is copy-bound too and the")
    print("fractionation limit is real for every method.")


if __name__ == "__main__":
    main()
