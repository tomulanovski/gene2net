"""Characterize the model's FALSE POSITIVE WGD edges on the simulated val set.

The FP/FN table shows the model is false-positive-leaning (recall ~1.0, the
errors are over-flagging). This asks WHERE those false positives sit relative to
the true WGD edges, which decides the fix:
  - mostly PARENT/ANCESTOR of a true edge  -> a localization problem (the whole
    ancestor chain lights up; build/collapse logic or a sharper feature helps).
  - mostly SCATTERED (unrelated to any true edge) -> a feature/signal problem.

For each val sample we run the model, threshold at --threshold, and for every
false-positive edge classify its relationship to the true WGD edges using the
species-tree topology (reconstructed from the stored edge_index). Edge i of the
model output is the i-th preorder (parent,child) pair, matching wgd_counts[i].

Run in the final_project env.

Usage:
    python scripts/fp_analysis.py --data-dir data/mul_trees_2k/training/ils_low \
        --model-dir output/phase1_feat9_full --threshold 0.90 --max-samples 400
"""
import argparse
import os
import random
import sys
from collections import Counter

import torch
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.data.tree_io import reorder_edge_index_preorder
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2
from gene2net_gnn.training.trainer_phase1 import prepare_sample


def edge_topology(sample):
    """From the stored edge_index, return per preorder-edge:
        child_node[i], parent_edge[i] (or None), clade[i] (frozenset of leaves).
    Aligned with the model's edge order (preorder (parent,child) pairs)."""
    ei = reorder_edge_index_preorder(sample.species_tree_edge_index)
    names = sample.species_tree_node_names
    is_leaf = sample.species_tree_is_leaf.tolist()

    # Even-indexed pairs are (parent, child) in preorder.
    parents, childs = [], []
    children_map = {}
    for k in range(0, ei.shape[1], 2):
        p = int(ei[0, k]); c = int(ei[1, k])
        parents.append(p); childs.append(c)
        children_map.setdefault(p, []).append(c)

    # Descendant leaves of each node (memoized).
    leaf_cache = {}

    def leaves(node):
        if node in leaf_cache:
            return leaf_cache[node]
        if is_leaf[node] or node not in children_map:
            res = frozenset([names[node]])
        else:
            res = frozenset().union(*(leaves(c) for c in children_map[node]))
        leaf_cache[node] = res
        return res

    child_to_edge = {c: i for i, c in enumerate(childs)}
    n = len(childs)
    clade = [leaves(childs[i]) for i in range(n)]
    parent_edge = [child_to_edge.get(parents[i]) for i in range(n)]
    return child_to_edge, parent_edge, clade, n


def classify_fp(i, true_edges, parent_edge, clade):
    """Relationship of false-positive edge i to the set of true WGD edges."""
    if not true_edges:
        return "no_true_events"
    # direct parent of a true edge (the classic ancestor over-flag)
    for t in true_edges:
        if parent_edge[t] == i:
            return "parent_of_true"
    # direct child of a true edge
    if parent_edge[i] in true_edges:
        return "child_of_true"
    # shares a parent edge with a true edge
    for t in true_edges:
        if parent_edge[i] is not None and parent_edge[i] == parent_edge[t]:
            return "sibling_of_true"
    # higher ancestor (clade strictly contains a true clade)
    for t in true_edges:
        if clade[i] > clade[t]:
            return "ancestor_of_true"
    # deeper descendant (clade strictly inside a true clade)
    for t in true_edges:
        if clade[i] < clade[t]:
            return "descendant_of_true"
    return "scattered"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, nargs="+")
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--config", default=None)
    parser.add_argument("--threshold", type=float, default=0.90)
    parser.add_argument("--max-samples", type=int, default=400)
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    cfg = args.config or os.path.join(base_dir, "configs", "phase1.yaml")
    mc = yaml.safe_load(open(cfg)).get("model", {})
    expected_edge_dim = int(mc.get("edge_feat_dim", 9))
    device = "cpu"

    model = SpeciesTreeGNNv2(
        node_feat_dim=int(mc.get("node_feat_dim", 13)),
        edge_feat_dim=expected_edge_dim,
        hidden_dim=int(mc.get("hidden_dim", 64)),
        n_gat_layers=int(mc.get("n_gat_layers", 3)),
        n_gat_heads=int(mc.get("n_gat_heads", 4)),
        dropout=float(mc.get("dropout", 0.2)),
    )
    for name in ["best_f1_model.pt", "best_model.pt"]:
        p = os.path.join(args.model_dir, name)
        if os.path.exists(p):
            model.load_state_dict(
                torch.load(p, map_location=device, weights_only=True), strict=False)
            print(f"Loaded {name}")
            break
    model = model.to(device).eval()

    datasets = [Gene2NetDataset(d) for d in args.data_dir]
    all_pairs = [(ds, i) for ds in datasets for i in range(len(ds))]
    indices = list(range(len(all_pairs)))
    random.seed(42)
    random.shuffle(indices)
    n_val = int(len(all_pairs) * 0.2)
    val_indices = indices[:min(n_val, args.max_samples)]

    rel_counts = Counter()
    fp_clade_sizes = []
    total_fp = total_tp = total_fn = 0

    for idx in val_indices:
        ds, li = all_pairs[idx]
        try:
            s = ds[li]
        except Exception:
            continue
        if s.labels is None or s.species_tree_edge_features is None:
            continue
        if s.species_tree_edge_features.shape[1] != expected_edge_dim:
            continue

        prepared = prepare_sample(s, torch.device(device))
        if prepared is None:
            continue
        inputs, targets, mask = prepared
        with torch.no_grad():
            logits, _ = model(**inputs)
            prob = torch.softmax(logits, dim=-1)[:, 1]

        _, parent_edge, clade, n = edge_topology(s)
        n = min(n, prob.shape[0], targets.shape[0])
        preds = (prob[:n] >= args.threshold).tolist()
        true = [int(targets[i].item()) == 1 for i in range(n)]
        true_edges = [i for i in range(n) if true[i] and bool(mask[i])]

        for i in range(n):
            if not bool(mask[i]):
                continue
            if preds[i] and true[i]:
                total_tp += 1
            elif preds[i] and not true[i]:
                total_fp += 1
                rel_counts[classify_fp(i, true_edges, parent_edge, clade)] += 1
                fp_clade_sizes.append(len(clade[i]))
            elif (not preds[i]) and true[i]:
                total_fn += 1

    print(f"\nThreshold {args.threshold} on {len(val_indices)} val samples")
    print(f"TP={total_tp}  FP={total_fp}  FN={total_fn}\n")
    if total_fp == 0:
        print("No false positives.")
        return
    print("False-positive edge -> relationship to nearest true WGD edge:")
    for rel, c in rel_counts.most_common():
        print(f"  {rel:>20} : {c:>5}  ({100*c/total_fp:.1f}%)")
    avg_sz = sum(fp_clade_sizes) / len(fp_clade_sizes)
    n_tip = sum(1 for z in fp_clade_sizes if z == 1)
    print(f"\nFP clade size: mean {avg_sz:.1f}, tips (size 1) {n_tip}/{total_fp} "
          f"({100*n_tip/total_fp:.1f}%)")
    print("\nIf parent/ancestor dominate -> localization fix (collapse/sharper feature).")
    print("If scattered/tips dominate  -> feature/signal fix.")


if __name__ == "__main__":
    main()
