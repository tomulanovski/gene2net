"""Permutation importance for PARTNER accuracy (complements permutation_importance.py,
which measures detection F1). Answers the design-justification question: do the features
that are unimportant for detection -- especially the co-clustering features -- drive the
allopolyploid partner prediction, which is the reason we keep them.

For each node feature, each edge feature, and each pairwise partner-feature channel, we
shuffle it on the validation split, rerun the partner head, and measure the drop in
allopolyploid partner accuracy (argmax predicted parent == the away-parent target). Node
and edge features act through the shared trunk and edge embeddings; the pairwise channels
feed the partner head directly. Uses the away labels the model was trained on.

Run in final_project (torch_geometric). Example:
  python scripts/partner_permutation_importance.py --data-dir data/mul_trees_2k/training_rooted/ils_low ... \
      --model-dir output/reconstruct_op_away_clade --config configs/reconstruct_oneparter.yaml
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
from gene2net_gnn.training.trainer_reconstruct import (
    build_pairwise_feat, build_partner_targets, set_feature_opts,
)
from permutation_importance import NODE_COLS, EDGE_COLS, build_pool_and_slots

PAIR_CHANNELS = ["pair_coclust_mean", "pair_coclust_max",
                 "pair_clustersupport_intensity", "pair_clustersupport_peak"]


def _apply_override(inputs, override, cursor, device):
    """Return (inputs, new_cursor) with one node/edge column globally permuted.
    Mirrors permutation_importance.run_model."""
    if override is None:
        return inputs, cursor
    if override["kind"] == "edge":
        ef = inputs["edge_features"].clone()
        n = ef.shape[0]
        idx = override["slots"][cursor:cursor + n]
        ef[:, override["col"]] = torch.tensor(
            override["pool"][override["perm"]][idx, override["col"]],
            dtype=ef.dtype, device=ef.device)
        return {**inputs, "edge_features": ef}, cursor + n
    # node
    nf = inputs["node_features"].clone()
    is_leaf = inputs["is_leaf"]
    leaf_rows = torch.nonzero(is_leaf, as_tuple=False).flatten().tolist()
    n = len(leaf_rows)
    idx = override["slots"][cursor:cursor + n]
    vals = override["pool"][override["perm"]][idx, override["col"]]
    for k, r in enumerate(leaf_rows):
        nf[r, override["col"]] = float(vals[k])
    prop = propagate_to_internal(nf, inputs["edge_index"], is_leaf, nf.shape[1]).to(device)
    return {**inputs, "node_features": nf, "node_features_propagated": prop}, cursor + n


def partner_acc(model, samples, device, override=None, pair_col=None, rng=None):
    """Allopolyploid partner accuracy over the samples, optionally with one node/edge
    column (override) or one pairwise channel (pair_col) permuted."""
    correct = total = 0
    cursor = 0
    with torch.no_grad():
        for s in samples:
            prepared = prepare_sample(s, torch.device(device))
            if prepared is None:
                continue
            inputs, wgd_targets, _ = prepared
            n_edges = wgd_targets.shape[0]
            inputs, cursor = _apply_override(inputs, override, cursor, device)
            _, edge_emb = model(**inputs)

            targets = build_partner_targets(s, n_edges, device)   # away partner per wgd edge
            qsel = (targets >= 0).nonzero(as_tuple=True)[0]
            if qsel.numel() == 0:
                continue
            allo = targets[qsel] != qsel                          # allo events only
            qsel = qsel[allo]
            if qsel.numel() == 0:
                continue

            pf = build_pairwise_feat(s).to(device)
            if pair_col is not None:
                pf = pf.clone()
                E = pf.shape[0]
                flat = pf[:, :, pair_col].reshape(-1).cpu().numpy()
                perm = rng.permutation(flat.shape[0])
                pf[:, :, pair_col] = torch.tensor(
                    flat[perm].reshape(E, E), dtype=pf.dtype, device=pf.device)

            scores = model.compute_partner_scores_rows(edge_emb, qsel, pf).squeeze(-1)
            pred = scores.argmax(dim=-1)
            correct += int((pred == targets[qsel]).sum())
            total += int(qsel.numel())
    return correct / max(total, 1), total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True, nargs="+")
    ap.add_argument("--model-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--max-samples", type=int, default=1000)
    ap.add_argument("--n-repeats", type=int, default=3)
    args = ap.parse_args()

    with open(args.config) as f:
        mc = yaml.safe_load(f).get("model", {})
    device = "cuda" if torch.cuda.is_available() else "cpu"
    set_feature_opts(coclust_condition_on_dup=mc.get("coclust_condition_on_dup", False),
                     use_n_eff=mc.get("use_n_eff", False))

    model = SpeciesTreeGNNv2(
        node_feat_dim=int(mc.get("node_feat_dim", 13)),
        edge_feat_dim=int(mc.get("edge_feat_dim", 9)),
        hidden_dim=int(mc.get("hidden_dim", 64)),
        n_gat_layers=int(mc.get("n_gat_layers", 3)),
        n_gat_heads=int(mc.get("n_gat_heads", 4)),
        dropout=float(mc.get("dropout", 0.2)),
        partner_pair_feat_dim=int(mc.get("partner_pair_feat_dim", 4)),
        n_parents=int(mc.get("n_parents", 1)),
    )
    ckpt = os.path.join(args.model_dir, "best_model.pt")
    model.load_state_dict(torch.load(ckpt, map_location=device, weights_only=True))
    model = model.to(device).eval()

    # away labels: the partner targets the model was trained on
    datasets = [Gene2NetDataset(d, away_labels=True) for d in args.data_dir]
    all_pairs = [ds[i] for ds in datasets for i in range(len(ds))]
    random.seed(42)
    random.shuffle(all_pairs)
    n_val = int(len(all_pairs) * 0.2)
    val = all_pairs[:n_val][:args.max_samples]

    base, n_allo = partner_acc(model, val, device)
    print(f"Val samples: {len(val)} | allo events scored: {n_allo}")
    print(f"Baseline allo partner accuracy: {base:.4f}\n")
    print(f"{'feature':<32}{'acc permuted':>14}{'drop':>10}")
    print("-" * 56)

    rows = []
    # node + edge (through the trunk)
    for kind, cols in (("node", NODE_COLS), ("edge", EDGE_COLS)):
        pool, slots = build_pool_and_slots(val, kind)
        for c, name in enumerate(cols):
            drops = []
            for rep in range(args.n_repeats):
                perm = np.random.default_rng(2000 + rep).permutation(len(pool))
                acc, _ = partner_acc(model, val, device,
                                     override=dict(kind=kind, col=c, perm=perm,
                                                   pool=pool, slots=slots))
                drops.append(base - acc)
            rows.append((name, base - float(np.mean(drops)), float(np.mean(drops))))
    # pairwise channels (direct to the partner head)
    for c, name in enumerate(PAIR_CHANNELS):
        drops = []
        for rep in range(args.n_repeats):
            rng = np.random.default_rng(3000 + rep)
            acc, _ = partner_acc(model, val, device, pair_col=c, rng=rng)
            drops.append(base - acc)
        rows.append((name, base - float(np.mean(drops)), float(np.mean(drops))))

    for name, acc, drop in sorted(rows, key=lambda r: -r[2]):
        print(f"{name:<32}{acc:>14.4f}{drop:>10.4f}")
    print("\nLarger drop = more important for PARTNER. Compare to the detection ranking:")
    print("features dead for detection but important here earn their place via the partner head.")


if __name__ == "__main__":
    main()
