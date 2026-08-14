"""Shared inference helpers for the reconstruction pipeline.

This was once a standalone decode script, but its `main()` was a redundant second
decode path that only supported the arbitrary-rooted ASTRAL backbone (a rooting
footgun). Decoding now goes exclusively through `benchmark_networks.py`
(hybrid-rooted, the real path) and `score_val_reconstruction.py`. This module now
exposes only the helpers those scripts import:

  - preorder_edge_clades(tree)             : per-edge clades in the model's edge order
  - model_inputs_for(sample, device)       : build the model forward() kwargs from a sample
  - load_model(model_dir, model_config, dev): load a checkpoint, inferring n_parents /
                                              pair_dim from the checkpoint shape
  - build_pairwise_feat                     : re-exported from trainer_reconstruct

Edge i of the model output is the i-th non-root node in preorder (alignment fixed in
reorder_edge_index_preorder), so clades are enumerated the same way.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch

from gene2net_gnn.data.tree_io import reorder_edge_index_preorder
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2, propagate_to_internal
# Re-exported so callers can `from scripts.reconstruct_mul_tree import build_pairwise_feat`.
from gene2net_gnn.training.trainer_reconstruct import build_pairwise_feat  # noqa: F401


def preorder_edge_clades(tree):
    clades = []
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        clades.append(frozenset(node.get_leaf_names()))
    return clades


def model_inputs_for(sample, device):
    from gene2net_gnn.training.trainer_reconstruct import _FEATURE_OPTS, augment_node_features_neff
    ei = reorder_edge_index_preorder(sample.species_tree_edge_index)
    sample._edge_index_pre = ei  # used by build_pairwise_feat
    if _FEATURE_OPTS["use_n_eff"]:
        augment_node_features_neff(sample)  # widen node features to 14 before prop
    prop = propagate_to_internal(
        sample.species_tree_node_features, ei,
        sample.species_tree_is_leaf, sample.species_tree_node_features.shape[1],
    )
    return {
        "node_features": sample.species_tree_node_features.to(device),
        "edge_index": ei.to(device),
        "edge_features": sample.species_tree_edge_features.to(device),
        "is_leaf": sample.species_tree_is_leaf.to(device),
        "node_features_propagated": prop.to(device),
    }


def load_model(model_dir, model_config, device):
    ckpt_path = None
    for name in ["best_model.pt", "best_partner_model.pt"]:
        p = os.path.join(model_dir, name)
        if os.path.exists(p):
            ckpt_path = p
            break
    if ckpt_path is None:
        raise FileNotFoundError(
            f"No checkpoint (best_model.pt or best_partner_model.pt) found in {model_dir!r}. "
            "Refusing to evaluate a randomly initialized model."
        )
    state = torch.load(ckpt_path, map_location=device, weights_only=True)

    hidden_dim = int(model_config.get("hidden_dim", 64))
    # Infer the partner-head shape FROM THE CHECKPOINT so a trained model loads
    # regardless of config drift: n_parents = final-layer out-features (1=one-partner,
    # 2=two-parent); pair_dim = first-layer in-features minus the two edge embeddings.
    n_parents = int(model_config.get("n_parents", 1))
    pair_dim = int(model_config.get("partner_pair_feat_dim", 2))
    if "partner_head.3.weight" in state:
        n_parents = int(state["partner_head.3.weight"].shape[0])
    if "partner_head.0.weight" in state:
        pair_dim = int(state["partner_head.0.weight"].shape[1]) - 2 * hidden_dim

    model = SpeciesTreeGNNv2(
        node_feat_dim=int(model_config.get("node_feat_dim", 13)),
        edge_feat_dim=int(model_config.get("edge_feat_dim", 9)),
        hidden_dim=hidden_dim,
        n_gat_layers=int(model_config.get("n_gat_layers", 3)),
        n_gat_heads=int(model_config.get("n_gat_heads", 4)),
        dropout=float(model_config.get("dropout", 0.2)),
        partner_pair_feat_dim=pair_dim,
        n_parents=n_parents,
    )
    model.load_state_dict(state)
    return model.to(device).eval()
