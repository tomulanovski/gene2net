"""End-to-end inference: gene trees → MUL-tree."""
import torch
from typing import List, Optional
from ete3 import Tree

from gene2net_gnn.data.tree_io import load_gene_trees_from_file, get_species_set, tree_to_edge_index
from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.model.gene2net_model import Gene2NetModel
from gene2net_gnn.inference.mul_tree_builder import build_mul_tree, WGDEvent


def predict_mul_tree(
    model: Gene2NetModel,
    gene_trees: List[Tree],
    species_tree: Tree,
    device: torch.device = torch.device("cpu"),
    threshold: float = 0.5,
    max_events: Optional[int] = None,
) -> tuple:
    """Run inference to predict a MUL-tree from gene trees.

    Args:
        model: trained Gene2NetModel
        gene_trees: list of ETE3 gene trees
        species_tree: ETE3 species tree (e.g., from ASTRAL-Pro)
        device: torch device
        threshold: confidence threshold for WGD predictions
        max_events: if set, take top-N events by confidence instead of thresholding

    Returns:
        (mul_tree, events, raw_predictions) tuple
    """
    model.eval()

    species_list = sorted(set(species_tree.get_leaf_names()))
    sp_to_idx = {sp: i for i, sp in enumerate(species_list)}

    # Build sample
    sample = Gene2NetSample.from_trees(species_tree, gene_trees, species_list)

    # Prepare model inputs
    gene_tree_data = []
    for i in range(len(sample.gene_tree_edge_indices)):
        gene_tree_data.append({
            "edge_index": sample.gene_tree_edge_indices[i].to(device),
            "species_ids": sample.gene_tree_species_ids[i].to(device),
            "branch_lengths": sample.gene_tree_branch_lengths[i].to(device),
            "leaf_mask": sample.gene_tree_leaf_masks[i].to(device),
        })

    species_ids = torch.tensor(
        [sp_to_idx.get(n, -1) for n in sample.species_tree_node_names],
        dtype=torch.long, device=device
    )

    species_tree_data = {
        "edge_index": sample.species_tree_edge_index.to(device),
        "node_features": sample.species_tree_node_features.to(device),
        "is_leaf": sample.species_tree_is_leaf.to(device),
        "species_ids": species_ids,
    }
    if sample.species_tree_edge_features is not None:
        species_tree_data["edge_features"] = sample.species_tree_edge_features.to(device)

    # Forward pass
    with torch.no_grad():
        wgd_logits, partner_logits, edge_emb = model(gene_tree_data, species_tree_data)

    # Decode predictions
    wgd_probs = torch.softmax(wgd_logits, dim=-1)  # [E, 4]
    wgd_pred = wgd_probs.argmax(dim=-1)  # [E] predicted event count per edge
    wgd_confidence = 1.0 - wgd_probs[:, 0]  # confidence = 1 - P(no event)

    # Get edge→clade mapping from species tree
    edge_index = sample.species_tree_edge_index
    node_names = sample.species_tree_node_names
    edge_clades = _get_edge_clades(species_tree, edge_index, node_names)

    # Select WGD events
    events = []
    if max_events is not None:
        # Take top-N by confidence
        top_indices = wgd_confidence.argsort(descending=True)[:max_events]
        for idx in top_indices:
            idx = idx.item()
            if wgd_pred[idx] > 0:
                partner_idx = partner_logits[idx].argmax().item()
                events.append(WGDEvent(
                    wgd_edge_clade=frozenset(edge_clades[idx]),
                    partner_edge_clade=frozenset(edge_clades[partner_idx]),
                    confidence=wgd_confidence[idx].item(),
                ))
    else:
        # Threshold-based
        for idx in range(len(wgd_pred)):
            if wgd_pred[idx] > 0 and wgd_confidence[idx] > threshold:
                partner_idx = partner_logits[idx].argmax().item()
                events.append(WGDEvent(
                    wgd_edge_clade=frozenset(edge_clades[idx]),
                    partner_edge_clade=frozenset(edge_clades[partner_idx]),
                    confidence=wgd_confidence[idx].item(),
                ))

    # Build MUL-tree
    mul_tree = build_mul_tree(species_tree, events)

    raw_predictions = {
        "wgd_logits": wgd_logits.cpu(),
        "partner_logits": partner_logits.cpu(),
        "wgd_probs": wgd_probs.cpu(),
        "wgd_confidence": wgd_confidence.cpu(),
        "edge_clades": edge_clades,
    }

    return mul_tree, events, raw_predictions


def _get_edge_clades(species_tree, edge_index, node_names):
    """Map each undirected edge index to the leaf set below it.

    Each undirected edge corresponds to a parent→child edge.
    The clade is the set of leaf names below the child node.
    """
    # Build node→leaves mapping using the species tree
    node_to_idx = {}
    for i, node in enumerate(species_tree.traverse("preorder")):
        node_to_idx[id(node)] = i

    idx_to_leaves = {}
    for i, node in enumerate(species_tree.traverse("preorder")):
        idx_to_leaves[i] = set(node.get_leaf_names())

    # Each undirected edge: take the child node (even-indexed edges are parent→child)
    edge_clades = {}
    n_undirected = edge_index.shape[1] // 2
    for e in range(n_undirected):
        child_idx = edge_index[1, e * 2].item()  # child node of parent→child edge
        edge_clades[e] = idx_to_leaves.get(child_idx, set())

    return edge_clades
