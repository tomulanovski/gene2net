"""Permutation importance for Phase 1 features (edge AND node).

For each feature column, shuffle its values across all validation rows (breaking
its association with the labels while keeping its marginal distribution), rerun
the model, and measure the F1 drop at a fixed threshold. A large drop means the
feature is important to the trained model. Near-zero or negative drop = prunable.

This is the in-GNN way to rank features and decide what to prune. Crucially it
also lets us justify keeping node features in academia: if a node feature (the
co-clustering summaries especially) shows little detection importance here but
the partner head still needs it, that is the evidence for the design choice.

Node features are permuted over LEAF rows only, then the internal-node
propagation is recomputed, because internal nodes are filled by propagation from
the leaves (their raw input is overwritten).

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
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2, propagate_to_internal
from gene2net_gnn.training.trainer_phase1 import prepare_sample

EDGE_COLS = [
    "concordance", "branch_length", "clade_size", "depth",
    "dup_synchrony", "mirrored_sister_frac",
    "copy_pair_div_mean", "copy_pair_div_cv", "frac_clade_dup",
]
NODE_COLS = [
    "mean_copies", "var_copies", "mode_copies", "p_absent",
    "p_1_copy", "p_2_copies", "p_3plus_copies", "max_copies",
    "coclust_mean", "coclust_std", "coclust_max", "coclust_min", "coclust_median",
]


def f1_at(probs, targets, thresh):
    preds = (probs >= thresh).astype(int)
    tp = int(((preds == 1) & (targets == 1)).sum())
    fp = int(((preds == 1) & (targets == 0)).sum())
    fn = int(((preds == 0) & (targets == 1)).sum())
    prec = tp / max(tp + fp, 1)
    rec = tp / max(tp + fn, 1)
    return 2 * prec * rec / max(prec + rec, 1e-8)


def run_model(model, samples, device, override=None):
    """Run the model over all samples. `override` optionally replaces one feature
    column with globally-permuted values:
        override = dict(kind="edge"|"node", col=int, perm=ndarray,
                        pool=ndarray, slots=ndarray)
    For node overrides only leaf rows are replaced and the propagation is
    recomputed (internal nodes are derived from the leaves).
    """
    probs_all, targets_all = [], []
    cursor = 0
    with torch.no_grad():
        for s in samples:
            prepared = prepare_sample(s, torch.device(device))
            if prepared is None:
                continue
            inputs, targets, mask = prepared

            if override is not None and override["kind"] == "edge":
                ef = inputs["edge_features"].clone()
                n_rows = ef.shape[0]
                idx = override["slots"][cursor:cursor + n_rows]
                ef[:, override["col"]] = torch.tensor(
                    override["pool"][override["perm"]][idx, override["col"]],
                    dtype=ef.dtype, device=ef.device,
                )
                inputs = {**inputs, "edge_features": ef}
                cursor += n_rows

            elif override is not None and override["kind"] == "node":
                nf = inputs["node_features"].clone()
                is_leaf = inputs["is_leaf"]
                leaf_rows = torch.nonzero(is_leaf, as_tuple=False).flatten().tolist()
                n_rows = len(leaf_rows)
                idx = override["slots"][cursor:cursor + n_rows]
                vals = override["pool"][override["perm"]][idx, override["col"]]
                for k, r in enumerate(leaf_rows):
                    nf[r, override["col"]] = float(vals[k])
                prop = propagate_to_internal(
                    nf, inputs["edge_index"], is_leaf, nf.shape[1]
                ).to(device)
                inputs = {**inputs, "node_features": nf,
                          "node_features_propagated": prop}
                cursor += n_rows

            logits, _ = model(**inputs)
            p = torch.softmax(logits, dim=-1)[:, 1]
            probs_all.append(p[mask].cpu().numpy())
            targets_all.append(targets[mask].cpu().numpy())

    return np.concatenate(probs_all), np.concatenate(targets_all)


def build_pool_and_slots(samples, kind):
    """Build a (rows x features) pool and a sample-ordered slots index.
    kind="edge": all edge rows. kind="node": leaf node rows only."""
    pool_rows, slots = [], []
    for s in samples:
        if kind == "edge":
            mat = s.species_tree_edge_features.numpy()
            idx = list(range(mat.shape[0]))
        else:
            mat = s.species_tree_node_features.numpy()
            is_leaf = s.species_tree_is_leaf.numpy().astype(bool)
            idx = list(np.nonzero(is_leaf)[0])
        start = len(pool_rows)
        for r in idx:
            pool_rows.append(mat[r])
        slots.extend(range(start, start + len(idx)))
    return np.array(pool_rows), np.array(slots)


def importance_for(model, samples, device, kind, cols, base_f1, threshold, n_repeats):
    pool, slots = build_pool_and_slots(samples, kind)
    results = []
    for c in range(len(cols)):
        drops = []
        for rep in range(n_repeats):
            rng = np.random.default_rng(1000 + rep)
            perm = rng.permutation(len(pool))
            p, t = run_model(model, samples, device,
                             override=dict(kind=kind, col=c, perm=perm,
                                           pool=pool, slots=slots))
            drops.append(base_f1 - f1_at(p, t, threshold))
        mean_drop = float(np.mean(drops))
        results.append((cols[c], base_f1 - mean_drop, mean_drop))
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, nargs="+")
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--config", default=None)
    parser.add_argument("--threshold", type=float, default=0.88)
    parser.add_argument("--max-samples", type=int, default=1000)
    parser.add_argument("--n-repeats", type=int, default=3, help="Permutations averaged per feature")
    parser.add_argument("--features", choices=["edge", "node", "both"], default="both")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    config_path = args.config or os.path.join(base_dir, "configs", "phase1.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)
    mc = config.get("model", {})
    expected_edge_dim = int(mc.get("edge_feat_dim", 9))
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
            # strict=False: detection-only checkpoints lack the partner_head keys
            # the model class now builds; that head is unused here.
            model.load_state_dict(
                torch.load(p, map_location=device, weights_only=True), strict=False)
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

    # Baseline
    probs, targets = run_model(model, samples, device)
    base_f1 = f1_at(probs, targets, args.threshold)
    print(f"\nBaseline F1 @ {args.threshold}: {base_f1:.4f}\n")

    results = []
    if args.features in ("edge", "both"):
        results += importance_for(model, samples, device, "edge", EDGE_COLS,
                                  base_f1, args.threshold, args.n_repeats)
    if args.features in ("node", "both"):
        results += importance_for(model, samples, device, "node", NODE_COLS,
                                  base_f1, args.threshold, args.n_repeats)

    print(f"{'feature':>20} | {'F1 permuted':>11} | {'drop':>7}")
    print("-" * 46)
    for name, f1p, drop in sorted(results, key=lambda x: -x[2]):
        print(f"{name:>20} | {f1p:>11.4f} | {drop:>7.4f}")

    print("\nLarger drop = more important. Near-zero or negative drop = prunable.")
    print("Note: this ranks features by their effect on DETECTION F1 only. A node")
    print("feature can be unimportant here yet still drive partner/allo accuracy.")


if __name__ == "__main__":
    main()
