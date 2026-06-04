"""Permutation importance for Phase 1 edge features.

For each edge-feature column, shuffle its values across all validation edges
(breaking its association with the labels while keeping its distribution), rerun
the model, and measure the F1 drop at a fixed threshold. A large drop means the
feature is important to the trained model.

This is the in-GNN way to rank features and decide what to prune — no separate
model needed.

Usage:
    python scripts/permutation_importance.py \
        --data-dir data/mul_trees_2k/training/ils_low [...all 9...] \
        --model-dir output/phase1_feat9 --threshold 0.88 --max-samples 1000
"""
import argparse
import os
import random
import sys

import numpy as np
import torch
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2
from gene2net_gnn.training.trainer_phase1 import prepare_sample

EDGE_COLS = [
    "concordance", "branch_length", "clade_size", "depth",
    "dup_synchrony", "mirrored_sister_frac",
    "copy_pair_div_mean", "copy_pair_div_cv", "frac_clade_dup",
]


def f1_at(probs, targets, thresh):
    preds = (probs >= thresh).astype(int)
    tp = int(((preds == 1) & (targets == 1)).sum())
    fp = int(((preds == 1) & (targets == 0)).sum())
    fn = int(((preds == 0) & (targets == 1)).sum())
    prec = tp / max(tp + fp, 1)
    rec = tp / max(tp + fn, 1)
    return 2 * prec * rec / max(prec + rec, 1e-8)


def run_model(model, samples, device, override_col=None, perm=None, pool=None, slots=None):
    """Run the model over all samples, optionally overriding one edge-feature
    column with globally-permuted values. Returns (probs, targets) arrays."""
    probs_all, targets_all = [], []
    row_cursor = 0
    with torch.no_grad():
        for s in samples:
            prepared = prepare_sample(s, torch.device(device))
            if prepared is None:
                continue
            inputs, targets, mask = prepared
            n_rows = inputs["edge_features"].shape[0]

            if override_col is not None:
                ef = inputs["edge_features"].clone()
                # rows for this sample within the global pool
                idx = slots[row_cursor:row_cursor + n_rows]
                ef[:, override_col] = torch.tensor(
                    pool[perm][idx, override_col], dtype=ef.dtype, device=ef.device
                )
                inputs = {**inputs, "edge_features": ef}
            row_cursor += n_rows

            logits, _ = model(**inputs)
            p = torch.softmax(logits, dim=-1)[:, 1]
            probs_all.append(p[mask].cpu().numpy())
            targets_all.append(targets[mask].cpu().numpy())

    return np.concatenate(probs_all), np.concatenate(targets_all)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, nargs="+")
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--config", default=None)
    parser.add_argument("--threshold", type=float, default=0.88)
    parser.add_argument("--max-samples", type=int, default=1000)
    parser.add_argument("--n-repeats", type=int, default=3, help="Permutations averaged per feature")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    config_path = args.config or os.path.join(base_dir, "configs", "phase1.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)
    mc = config.get("model", {})
    expected_edge_dim = int(mc.get("edge_feat_dim", 4))
    device = "cuda" if torch.cuda.is_available() else "cpu"

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
            model.load_state_dict(torch.load(p, map_location=device, weights_only=True))
            print(f"Loaded {name}")
            break
    model = model.to(device).eval()

    # Collect val samples (same split logic as tune_threshold), guarding dim.
    datasets = [Gene2NetDataset(d) for d in args.data_dir]
    all_pairs = [(ds, i) for ds in datasets for i in range(len(ds))]
    indices = list(range(len(all_pairs)))
    random.seed(42)
    random.shuffle(indices)
    n_val = int(len(all_pairs) * 0.2)
    val_indices = indices[:min(n_val, args.max_samples)]

    samples = []
    for idx in val_indices:
        ds, li = all_pairs[idx]
        try:
            s = ds[li]
        except Exception:
            continue
        if s.labels is None:
            continue
        ef = s.species_tree_edge_features
        if ef is None or ef.shape[1] != expected_edge_dim:
            continue
        # Phase 1 ignores gene trees — drop them to keep memory small.
        s.gene_tree_edge_indices = []
        s.gene_tree_species_ids = []
        s.gene_tree_branch_lengths = []
        s.gene_tree_leaf_masks = []
        samples.append(s)
    print(f"Val samples used: {len(samples)}")

    # Build global edge-feature pool (rows in sample order) for permutation.
    pool_rows, slots = [], []
    for s in samples:
        ef = s.species_tree_edge_features
        start = len(pool_rows)
        for r in range(ef.shape[0]):
            pool_rows.append(ef[r].numpy())
        slots.extend(range(start, start + ef.shape[0]))
    pool = np.array(pool_rows)
    slots = np.array(slots)

    # Baseline
    probs, targets = run_model(model, samples, device)
    base_f1 = f1_at(probs, targets, args.threshold)
    print(f"\nBaseline F1 @ {args.threshold}: {base_f1:.4f}\n")

    print(f"{'feature':>20} | {'F1 permuted':>11} | {'drop':>7}")
    print("-" * 46)
    results = []
    for c in range(expected_edge_dim):
        drops = []
        for rep in range(args.n_repeats):
            rng = np.random.default_rng(1000 + rep)
            perm = rng.permutation(len(pool))
            p, t = run_model(model, samples, device,
                             override_col=c, perm=perm, pool=pool, slots=slots)
            drops.append(base_f1 - f1_at(p, t, args.threshold))
        mean_drop = float(np.mean(drops))
        results.append((EDGE_COLS[c], base_f1 - mean_drop, mean_drop))

    for name, f1p, drop in sorted(results, key=lambda x: -x[2]):
        print(f"{name:>20} | {f1p:>11.4f} | {drop:>7.4f}")

    print("\nLarger drop = more important. Near-zero or negative drop = prunable.")


if __name__ == "__main__":
    main()
