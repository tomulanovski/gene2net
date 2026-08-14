"""Dump the model's flagged WGD edges for one benchmark network so we can see
the actual structure of the predictions (clade sizes, whether flagged edges are
tips / siblings / nested, and their predicted partners). This is what we need
to design the right MUL-tree build.

Usage (final_project env):
    python scripts/debug_predictions.py --model-dir output/reconstruct_allo \
        --config conf_ils_low_10M --network Zhao_2021 --threshold 0.9
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import yaml
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.inference.build_strategies import build_parent_edge_map
from scripts.reconstruct_mul_tree import (
    load_model, model_inputs_for, preorder_edge_clades, build_pairwise_feat,
)
from scripts.benchmark_networks import load_gene_trees


def fmt(clade):
    s = sorted(clade)
    return f"{{{', '.join(s[:6])}{', ...' if len(s) > 6 else ''}}}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--network", required=True)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--threshold", type=float, default=0.9)
    parser.add_argument("--model-config", default=None)
    parser.add_argument("--sim-base",
                        default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    cfg = args.model_config or os.path.join(base_dir, "configs", "reconstruct.yaml")
    mc = yaml.safe_load(open(cfg)).get("model", {})

    rep = os.path.join(args.sim_base, args.network, "processed", args.config,
                       "grampa_input", f"replicate_{args.replicate}")
    gts = load_gene_trees(os.path.join(rep, "clean_trees.tre"), 500)
    astral = Tree(open(os.path.join(rep, "species.tre")).read().strip(), format=1)
    sl = sorted(set(astral.get_leaf_names()))
    sample = Gene2NetSample.from_trees(species_tree=astral, gene_trees=gts, species_list=sl)
    clades = preorder_edge_clades(astral)
    parent_edge = build_parent_edge_map(astral)

    model = load_model(args.model_dir, mc, "cpu")
    with torch.no_grad():
        inp = model_inputs_for(sample, "cpu")
        logits, emb = model(**inp)
        prob = torch.softmax(logits, dim=-1)[:, 1]
        pw = build_pairwise_feat(sample)

    flagged = [i for i in range(len(clades)) if float(prob[i]) >= args.threshold]
    print(f"\n{args.network}: {len(clades)} edges, {len(flagged)} flagged at threshold {args.threshold}")
    print(f"(species tree has {len(sl)} taxa)\n")

    if not flagged:
        print("No flagged edges.")
        return

    flagged_set = set(flagged)
    q = torch.tensor(flagged, dtype=torch.long)
    rows = model.compute_partner_scores_rows(emb, q, pw)

    n_tip = sum(1 for i in flagged if len(clades[i]) == 1)
    n_nested = sum(1 for i in flagged if parent_edge[i] in flagged_set)
    print(f"Flagged edges that are single species (tips): {n_tip}/{len(flagged)}")
    print(f"Flagged edges whose PARENT edge is also flagged (nested): {n_nested}/{len(flagged)}")
    print(f"-> collapse (parent-chain) can only merge {n_nested} edges\n")

    print("edge | prob | clade_size | clade | partner")
    for r, i in enumerate(flagged):
        j = int(rows[r, :len(clades)].argmax())
        ptype = "AUTO(self)" if j == i else "allo"
        print(f"  {i:>3} | {float(prob[i]):.2f} | {len(clades[i]):>2} | {fmt(clades[i])} | {ptype} {fmt(clades[j]) if j != i else ''}")


if __name__ == "__main__":
    main()
