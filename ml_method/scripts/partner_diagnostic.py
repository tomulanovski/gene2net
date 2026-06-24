"""Diagnose the partner head's failure modes on the simulated val set.

The oracle showed reticulation-leaf accuracy has real headroom that lives in the
PREDICTIONS, and the partner head is the prediction that drives it. This script
isolates partner quality from detection: it conditions on the TRUE WGD edges
(from the labels) and asks, for each, whether the predicted partner is right and
if not, what kind of error it is:

  - auto vs allo confusion (predicted self when truth is a partner, or vice versa)
  - allo partner wrong but NEARBY the true partner (clade overlaps) vs FAR (disjoint)
  - predicted partner is a tip (single species)
  - reciprocal pairs: two WGD edges that each predict the other as partner

This tells us whether the fix is a tie-breaking / consistency rule (reciprocal,
near-misses) or a feature/signal problem (far misses).

Run in the final_project env.

Usage:
    python scripts/partner_diagnostic.py \
        --data-dir data/mul_trees_2k/training/ils_low [...all 9...] \
        --model-dir output/reconstruct_allo --max-samples 1000
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
from gene2net_gnn.training.trainer_reconstruct import build_pairwise_feat, build_partner_targets
from scripts.reconstruct_mul_tree import load_model, model_inputs_for


def edge_clades(sample):
    """Clade (frozenset of leaf names) per preorder edge, aligned with model edges."""
    ei = reorder_edge_index_preorder(sample.species_tree_edge_index)
    names = sample.species_tree_node_names
    is_leaf = sample.species_tree_is_leaf.tolist()
    childs, children_map = [], {}
    for k in range(0, ei.shape[1], 2):
        p, c = int(ei[0, k]), int(ei[1, k])
        childs.append(c)
        children_map.setdefault(p, []).append(c)
    cache = {}

    def leaves(node):
        if node in cache:
            return cache[node]
        if is_leaf[node] or node not in children_map:
            res = frozenset([names[node]])
        else:
            res = frozenset().union(*(leaves(c) for c in children_map[node]))
        cache[node] = res
        return res

    return [leaves(c) for c in childs]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, nargs="+")
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--model-config", default=None)
    parser.add_argument("--max-samples", type=int, default=1000)
    args = parser.parse_args()

    base = os.path.join(os.path.dirname(__file__), "..")
    cfg = args.model_config or os.path.join(base, "configs", "reconstruct.yaml")
    mc = yaml.safe_load(open(cfg)).get("model", {})
    edge_dim = int(mc.get("edge_feat_dim", 9))
    device = "cpu"
    model = load_model(args.model_dir, mc, device)

    datasets = [Gene2NetDataset(d) for d in args.data_dir]
    pairs = [(ds, i) for ds in datasets for i in range(len(ds))]
    idxs = list(range(len(pairs)))
    random.seed(42)
    random.shuffle(idxs)
    val = idxs[:min(int(len(pairs) * 0.2), args.max_samples)]

    auto_total = auto_correct = 0
    allo_total = allo_correct = 0
    err = Counter()       # error-type histogram (allo errors)
    reciprocal = 0
    n_partner_preds = 0
    pred_tip = 0

    for k in val:
        ds, li = pairs[k]
        try:
            s = ds[li]
        except Exception:
            continue
        if s.labels is None or s.species_tree_edge_features is None:
            continue
        if s.species_tree_edge_features.shape[1] != edge_dim:
            continue

        n_edges = s.species_tree_edge_index.shape[1] // 2
        targets = build_partner_targets(s, n_edges, torch.device(device))
        wgd_edges = [i for i in range(n_edges) if int(targets[i]) >= 0]
        if not wgd_edges:
            continue

        clades = edge_clades(s)
        with torch.no_grad():
            inputs = model_inputs_for(s, device)
            _, emb = model(**inputs)
            pw = build_pairwise_feat(s)
            q = torch.tensor(wgd_edges, dtype=torch.long)
            rows = model.compute_partner_scores_rows(emb, q, pw)  # [Q, E]

        pred = {}
        for r, w in enumerate(wgd_edges):
            j = int(rows[r, :n_edges].argmax())
            pred[w] = j
            n_partner_preds += 1
            if j < len(clades) and len(clades[j]) == 1:
                pred_tip += 1

        for r, w in enumerate(wgd_edges):
            p = int(targets[w])
            j = pred[w]
            true_auto = (p == w)
            if true_auto:
                auto_total += 1
                if j == w:
                    auto_correct += 1
            else:
                allo_total += 1
                if j == p:
                    allo_correct += 1
                    continue
                # classify the allo error
                if j == w:
                    err["predicted_auto_(self)"] += 1
                elif pred.get(j) == w:
                    err["reciprocal"] += 1
                elif j < len(clades) and p < len(clades) and (clades[j] & clades[p]):
                    err["near_(clade_overlaps_true)"] += 1
                else:
                    err["far_(disjoint_from_true)"] += 1

        # count reciprocal pairs once
        seen = set()
        for w, j in pred.items():
            if pred.get(j) == w and (j, w) not in seen:
                reciprocal += 1
                seen.add((w, j))

    print(f"\nPartner diagnostic on {len(val)} val samples (conditioned on TRUE WGD edges)\n")
    print(f"AUTO: {auto_correct}/{auto_total} = {auto_correct/max(auto_total,1):.3f}")
    print(f"ALLO: {allo_correct}/{allo_total} = {allo_correct/max(allo_total,1):.3f}")
    print(f"\nReciprocal predicted pairs: {reciprocal}")
    print(f"Predicted partner is a tip: {pred_tip}/{n_partner_preds} "
          f"({100*pred_tip/max(n_partner_preds,1):.1f}%)")
    print(f"\nAllo error breakdown (of {allo_total - allo_correct} allo errors):")
    for kind, c in err.most_common():
        print(f"  {kind:>28} : {c:>5} ({100*c/max(allo_total-allo_correct,1):.1f}%)")
    print("\nReading: 'near' / 'reciprocal' / 'predicted_auto' -> a consistency or "
          "tie-breaking fix. 'far' dominating -> a feature/signal problem.")


if __name__ == "__main__":
    main()
