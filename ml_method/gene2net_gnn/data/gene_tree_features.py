"""Gene tree summary features for Phase 2.

Computes per-edge features that capture clade-level duplication patterns
in gene trees. These supplement the 4 existing edge features (concordance,
branch length, clade size, depth) with 5 new features derived from gene
tree copy number patterns.
"""
import torch
from typing import Dict, List


def get_edge_clades(
    edge_index: torch.Tensor,
    is_leaf: torch.Tensor,
) -> Dict[int, List[int]]:
    """Identify which leaf nodes belong to the clade defined by each edge.

    For each undirected edge (stored as parent->child at even indices),
    the clade is the set of leaf nodes reachable from the child node.

    Args:
        edge_index: [2, 2E] undirected edges (parent->child at even indices).
        is_leaf: [N] boolean mask indicating leaf nodes.

    Returns:
        Dict mapping edge index (0..E-1) to list of leaf node indices in the clade.
    """
    n_edges = edge_index.shape[1] // 2

    # Build parent->children adjacency from even-indexed (parent->child) edges
    parent_to_children: Dict[int, List[int]] = {}
    for i in range(0, edge_index.shape[1], 2):
        parent = edge_index[0, i].item()
        child = edge_index[1, i].item()
        parent_to_children.setdefault(parent, []).append(child)

    # For each edge, find all descendant leaves of the child node
    def get_descendant_leaves(node: int) -> List[int]:
        if is_leaf[node]:
            return [node]
        leaves = []
        for child in parent_to_children.get(node, []):
            leaves.extend(get_descendant_leaves(child))
        return leaves

    edge_clades = {}
    for edge_idx in range(n_edges):
        child_node = edge_index[1, edge_idx * 2].item()
        edge_clades[edge_idx] = get_descendant_leaves(child_node)

    return edge_clades


def compute_gene_tree_edge_features(
    edge_clades: Dict[int, List[int]],
    species_ids_list: List[torch.Tensor],
    leaf_masks_list: List[torch.Tensor],
    species_tree_is_leaf: torch.Tensor,
    species_tree_node_names: List[str],
    species_list: List[str],
    n_edges: int,
) -> torch.Tensor:
    """Compute 5 gene tree summary features per species tree edge.

    For each edge (defining a clade), analyzes copy number patterns across
    gene trees to detect WGD signatures.

    Features computed per edge:
        1. clade_duplication_freq: fraction of gene trees where >=1 clade
           species has >1 copy.
        2. clade_mean_extra_copies: average extra copies (above 1) for
           clade species, averaged across gene trees.
        3. clade_copy_ratio: mean copies of clade species / mean copies of
           non-clade species, averaged across gene trees.
        4. clade_multi_copy_coherence: when one clade species is multi-copy,
           what fraction of other clade species are also multi-copy?
        5. clade_copy_cv: coefficient of variation of copy counts across
           clade species, averaged across gene trees.

    Args:
        edge_clades: Dict mapping edge index to list of leaf node indices.
        species_ids_list: List of [N_gt] tensors with species ID per gene tree node.
        leaf_masks_list: List of [N_gt] boolean tensors marking gene tree leaves.
        species_tree_is_leaf: [N_sp] boolean mask for species tree.
        species_tree_node_names: List of species tree node names.
        species_list: Sorted list of species names.
        n_edges: Number of undirected edges in species tree.

    Returns:
        Tensor [n_edges, 5] of new edge features.
    """
    sp_to_idx = {sp: i for i, sp in enumerate(species_list)}
    n_species = len(species_list)
    n_gene_trees = len(species_ids_list)

    # Map species tree leaf node indices to species indices
    node_to_sp_idx = {}
    for node_idx, name in enumerate(species_tree_node_names):
        if species_tree_is_leaf[node_idx] and name in sp_to_idx:
            node_to_sp_idx[node_idx] = sp_to_idx[name]

    # Precompute per-gene-tree copy counts: [n_gene_trees, n_species]
    copy_counts = torch.zeros(n_gene_trees, n_species, dtype=torch.float)
    for gt_idx in range(n_gene_trees):
        sp_ids = species_ids_list[gt_idx]
        leaf_mask = leaf_masks_list[gt_idx]
        leaf_species = sp_ids[leaf_mask]
        for sp_id in leaf_species:
            sp_val = sp_id.item()
            if 0 <= sp_val < n_species:
                copy_counts[gt_idx, sp_val] += 1

    # For each edge, compute the 5 features
    features = torch.zeros(n_edges, 5, dtype=torch.float)
    all_sp_indices = set(range(n_species))

    for edge_idx in range(n_edges):
        clade_nodes = edge_clades.get(edge_idx, [])
        clade_sp = set()
        for node_idx in clade_nodes:
            if node_idx in node_to_sp_idx:
                clade_sp.add(node_to_sp_idx[node_idx])

        if len(clade_sp) == 0:
            continue

        non_clade_sp = all_sp_indices - clade_sp
        clade_sp_list = sorted(clade_sp)
        non_clade_sp_list = sorted(non_clade_sp)

        # Per-gene-tree statistics
        dup_count = 0
        extra_copies_sum = 0.0
        ratio_sum = 0.0
        ratio_count = 0
        coherence_sum = 0.0
        coherence_count = 0
        cv_sum = 0.0
        cv_count = 0

        for gt_idx in range(n_gene_trees):
            clade_copies = copy_counts[gt_idx, clade_sp_list]

            # Feature 1: duplication frequency
            if (clade_copies > 1).any():
                dup_count += 1

            # Feature 2: mean extra copies
            extra = (clade_copies - 1).clamp(min=0).mean().item()
            extra_copies_sum += extra

            # Feature 3: copy ratio (clade vs non-clade)
            clade_mean = clade_copies.mean().item()
            if len(non_clade_sp_list) > 0:
                non_clade_mean = copy_counts[gt_idx, non_clade_sp_list].mean().item()
                if non_clade_mean > 0:
                    ratio_sum += clade_mean / non_clade_mean
                    ratio_count += 1

            # Feature 4: multi-copy coherence
            if len(clade_sp_list) > 1:
                multi_copy = (clade_copies > 1)
                n_multi = multi_copy.sum().item()
                if n_multi > 0:
                    coherence_sum += n_multi / len(clade_sp_list)
                    coherence_count += 1

            # Feature 5: copy count CV within clade
            if len(clade_sp_list) > 1:
                mean_c = clade_copies.mean().item()
                if mean_c > 0:
                    std_c = clade_copies.std().item()
                    cv_sum += std_c / mean_c
                    cv_count += 1

        # Aggregate across gene trees
        features[edge_idx, 0] = dup_count / max(n_gene_trees, 1)
        features[edge_idx, 1] = extra_copies_sum / max(n_gene_trees, 1)
        features[edge_idx, 2] = ratio_sum / max(ratio_count, 1)
        features[edge_idx, 3] = coherence_sum / max(coherence_count, 1)
        features[edge_idx, 4] = cv_sum / max(cv_count, 1)

    return features
