"""Gene tree topology features for Phase 3.

Computes per-edge features from gene tree STRUCTURE — specifically,
where paralogous copies are placed relative to each other. This captures
the key WGD signature that copy-count-based features miss:

WGD on edge leading to clade {A,B,C} produces gene trees like:
    ((A1,B1,C1),(A2,B2,C2))  ← copies co-located as sister clades

Gene-level duplication of just A produces:
    ((A1,B,C), ..., A2)  ← A2 scattered elsewhere

Three features per species tree edge:
    1. paralog_colocation: Are extra copies of different clade species
       co-located in the same gene tree subtree?
    2. clade_monophyly: Do clade species form a monophyletic group?
    3. duplication_depth_consistency: Is the duplication point at a
       consistent depth across gene trees?
"""
import torch
from typing import Dict, List, Optional, Set, Tuple


def _build_tree_maps(edge_index: torch.Tensor) -> Tuple[Dict[int, int], Dict[int, List[int]], int]:
    """Build parent and children maps from edge_index tensor.

    Args:
        edge_index: [2, 2E] with parent->child at even indices.

    Returns:
        (parent_map, children_map, root_node)
    """
    parent = {}
    children: Dict[int, List[int]] = {}
    all_nodes = set()

    for i in range(0, edge_index.shape[1], 2):
        p = edge_index[0, i].item()
        c = edge_index[1, i].item()
        parent[c] = p
        children.setdefault(p, []).append(c)
        all_nodes.add(p)
        all_nodes.add(c)

    root = None
    for n in all_nodes:
        if n not in parent:
            root = n
            break

    return parent, children, root


def _get_descendant_leaves(node: int, children: Dict[int, List[int]], is_leaf_set: Set[int]) -> Set[int]:
    """Get all leaf descendants of a node."""
    if node in is_leaf_set:
        return {node}
    result = set()
    for c in children.get(node, []):
        result |= _get_descendant_leaves(c, children, is_leaf_set)
    return result


def _find_lca(nodes: List[int], parent: Dict[int, int]) -> Optional[int]:
    """Find the lowest common ancestor of a set of nodes."""
    if not nodes:
        return None
    if len(nodes) == 1:
        return nodes[0]

    # Build ancestor path for first node
    path = []
    n = nodes[0]
    while True:
        path.append(n)
        if n not in parent:
            break
        n = parent[n]
    path_set = set(path)

    # For each subsequent node, walk up until hitting an ancestor
    for node in nodes[1:]:
        curr = node
        while curr not in path_set:
            curr = parent[curr]
        # curr is the LCA so far; restrict path_set to curr's ancestors
        new_path_set = set()
        c = curr
        while True:
            new_path_set.add(c)
            if c not in parent:
                break
            c = parent[c]
        path_set = new_path_set

    # The LCA is the deepest node still in path_set
    # Walk from any input node upward; first hit in path_set is LCA
    curr = nodes[0]
    while curr not in path_set:
        curr = parent[curr]
    return curr


def _node_depth(node: int, parent: Dict[int, int]) -> int:
    """Compute depth of a node (edges from root)."""
    depth = 0
    curr = node
    while curr in parent:
        depth += 1
        curr = parent[curr]
    return depth


def compute_topology_features(
    edge_clades: Dict[int, List[int]],
    gene_tree_edge_indices: List[torch.Tensor],
    gene_tree_species_ids: List[torch.Tensor],
    gene_tree_leaf_masks: List[torch.Tensor],
    species_tree_is_leaf: torch.Tensor,
    species_tree_node_names: List[str],
    species_list: List[str],
    n_edges: int,
) -> torch.Tensor:
    """Compute 3 topology features per species tree edge from gene trees.

    Args:
        edge_clades: Dict mapping edge index to list of species tree leaf node indices.
        gene_tree_edge_indices: List of [2, 2E_gt] tensors per gene tree.
        gene_tree_species_ids: List of [N_gt] tensors with species ID per node.
        gene_tree_leaf_masks: List of [N_gt] boolean tensors.
        species_tree_is_leaf: [N_sp] boolean.
        species_tree_node_names: List of species tree node names.
        species_list: Sorted species names.
        n_edges: Number of undirected edges in species tree.

    Returns:
        Tensor [n_edges, 3] of topology features.
    """
    sp_to_idx = {sp: i for i, sp in enumerate(species_list)}
    n_gene_trees = len(gene_tree_edge_indices)

    # Map species tree leaf nodes to species indices
    node_to_sp_idx = {}
    for node_idx, name in enumerate(species_tree_node_names):
        if species_tree_is_leaf[node_idx] and name in sp_to_idx:
            node_to_sp_idx[node_idx] = sp_to_idx[name]

    # Convert edge_clades to species index sets
    edge_clade_species: Dict[int, Set[int]] = {}
    for edge_idx, leaf_nodes in edge_clades.items():
        sp_set = set()
        for node_idx in leaf_nodes:
            if node_idx in node_to_sp_idx:
                sp_set.add(node_to_sp_idx[node_idx])
        edge_clade_species[edge_idx] = sp_set

    # Preprocess all gene trees: build tree structures and species->leaf maps
    gt_structures = []
    for gt_idx in range(n_gene_trees):
        ei = gene_tree_edge_indices[gt_idx]
        sp_ids = gene_tree_species_ids[gt_idx]
        leaf_mask = gene_tree_leaf_masks[gt_idx]

        parent, children_map, root = _build_tree_maps(ei)
        leaf_set = set(i for i in range(len(leaf_mask)) if leaf_mask[i])

        # Map species index -> list of gene tree leaf nodes
        species_to_leaves: Dict[int, List[int]] = {}
        for leaf_node in leaf_set:
            sp_val = sp_ids[leaf_node].item()
            if sp_val >= 0:
                species_to_leaves.setdefault(sp_val, []).append(leaf_node)

        gt_structures.append((parent, children_map, root, leaf_set, species_to_leaves))

    features = torch.zeros(n_edges, 3, dtype=torch.float)

    for edge_idx in range(n_edges):
        clade_sp = edge_clade_species.get(edge_idx, set())
        if len(clade_sp) < 2:
            # Single-species clade: no meaningful topology signal
            continue

        colocation_sum = 0.0
        colocation_count = 0
        monophyly_sum = 0.0
        monophyly_count = 0
        depth_values = []

        for gt_idx in range(n_gene_trees):
            parent, children_map, root, leaf_set, sp_to_leaves = gt_structures[gt_idx]

            # Collect all gene tree leaves belonging to clade species
            clade_leaves = []
            multi_copy_species = []
            for sp in clade_sp:
                leaves = sp_to_leaves.get(sp, [])
                clade_leaves.extend(leaves)
                if len(leaves) > 1:
                    multi_copy_species.append(sp)

            if not clade_leaves:
                continue

            # --- Feature 1: Paralog co-location ---
            # When 2+ clade species have duplicates, check if their extra
            # copies end up in the same subtree
            if len(multi_copy_species) >= 2:
                # Find LCA of ALL clade leaves
                lca = _find_lca(clade_leaves, parent)
                if lca is not None and lca in children_map and len(children_map[lca]) >= 2:
                    # Split leaves by which child subtree they fall into
                    child_subtrees = {}
                    for child_node in children_map[lca]:
                        desc = _get_descendant_leaves(child_node, children_map, leaf_set)
                        child_subtrees[child_node] = desc

                    # For each multi-copy species, check if copies land in
                    # different subtrees (= coordinated duplication = WGD)
                    coordinated = 0
                    total_pairs = 0
                    for sp in multi_copy_species:
                        leaves = sp_to_leaves[sp]
                        # Which subtrees contain copies of this species?
                        subtrees_with_sp = set()
                        for leaf in leaves:
                            for child_node, desc in child_subtrees.items():
                                if leaf in desc:
                                    subtrees_with_sp.add(child_node)
                                    break
                        if len(subtrees_with_sp) >= 2:
                            coordinated += 1
                        total_pairs += 1

                    if total_pairs > 0:
                        colocation_sum += coordinated / total_pairs
                        colocation_count += 1

            # --- Feature 2: Clade monophyly ---
            # Take one copy per species, check if they form a monophyletic group
            representative_leaves = []
            for sp in clade_sp:
                leaves = sp_to_leaves.get(sp, [])
                if leaves:
                    representative_leaves.append(leaves[0])

            if len(representative_leaves) >= 2:
                lca_rep = _find_lca(representative_leaves, parent)
                if lca_rep is not None:
                    all_desc = _get_descendant_leaves(lca_rep, children_map, leaf_set)
                    # How many non-clade leaves are mixed in?
                    clade_desc = set(representative_leaves)
                    n_total = len(all_desc)
                    n_clade = len(clade_desc & all_desc)
                    if n_total > 0:
                        monophyly_sum += n_clade / n_total
                        monophyly_count += 1

            # --- Feature 3: Duplication depth consistency ---
            # For species with >1 copy, find the depth of their LCA
            # (where the "duplication" happened in the gene tree)
            if multi_copy_species:
                for sp in multi_copy_species:
                    leaves = sp_to_leaves[sp]
                    dup_lca = _find_lca(leaves, parent)
                    if dup_lca is not None:
                        depth_values.append(_node_depth(dup_lca, parent))

        # Aggregate
        features[edge_idx, 0] = colocation_sum / max(colocation_count, 1)
        features[edge_idx, 1] = monophyly_sum / max(monophyly_count, 1)

        # Depth consistency: use 1 - normalized_std (high = consistent)
        if len(depth_values) >= 2:
            t = torch.tensor(depth_values, dtype=torch.float)
            mean_d = t.mean()
            if mean_d > 0:
                features[edge_idx, 2] = 1.0 - (t.std() / mean_d).clamp(max=1.0).item()
            else:
                features[edge_idx, 2] = 1.0
        elif len(depth_values) == 1:
            features[edge_idx, 2] = 1.0  # single observation = perfectly consistent

    return features
