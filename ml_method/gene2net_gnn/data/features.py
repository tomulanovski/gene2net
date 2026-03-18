"""Hand-crafted features for gene tree discordance and copy number patterns."""
from collections import Counter
from typing import Dict, List, Set

from ete3 import Tree


def compute_copy_count_features(gene_trees: List[Tree], species: str) -> dict:
    """Compute copy count statistics for a species across all gene trees.

    Parameters
    ----------
    gene_trees : list of Tree
        Gene trees to analyse.
    species : str
        Target species name.

    Returns
    -------
    dict with keys:
        mean_copies, var_copies, mode_copies, p_absent, p_1_copy,
        p_2_copies, p_3plus_copies, max_copies
    """
    counts = []
    for tree in gene_trees:
        n = sum(1 for leaf in tree.get_leaves() if leaf.name == species)
        counts.append(n)

    n_trees = len(gene_trees)
    if n_trees == 0:
        return {
            "mean_copies": 0.0,
            "var_copies": 0.0,
            "mode_copies": 0,
            "p_absent": 1.0,
            "p_1_copy": 0.0,
            "p_2_copies": 0.0,
            "p_3plus_copies": 0.0,
            "max_copies": 0,
        }

    mean_copies = sum(counts) / n_trees
    var_copies = sum((c - mean_copies) ** 2 for c in counts) / n_trees

    counter = Counter(counts)
    mode_copies = counter.most_common(1)[0][0]

    p_absent = sum(1 for c in counts if c == 0) / n_trees
    p_1_copy = sum(1 for c in counts if c == 1) / n_trees
    p_2_copies = sum(1 for c in counts if c == 2) / n_trees
    p_3plus_copies = sum(1 for c in counts if c >= 3) / n_trees
    max_copies = max(counts) if counts else 0

    return {
        "mean_copies": mean_copies,
        "var_copies": var_copies,
        "mode_copies": mode_copies,
        "p_absent": p_absent,
        "p_1_copy": p_1_copy,
        "p_2_copies": p_2_copies,
        "p_3plus_copies": p_3plus_copies,
        "max_copies": max_copies,
    }


def compute_clustering_profile(
    gene_trees: List[Tree],
    species: str,
    all_species: Set[str],
) -> Dict[str, float]:
    """Compute co-clustering frequencies between a target species and all others.

    For each other species S, compute the fraction of gene trees where at least
    one leaf of 'species' has a sibling leaf that is species S (i.e., they share
    a direct common parent as leaves).

    Parameters
    ----------
    gene_trees : list of Tree
    species : str
        Target species.
    all_species : set of str
        All species to consider.

    Returns
    -------
    dict mapping each other species name to its co-clustering frequency (float in [0,1]).
    """
    other_species = all_species - {species}
    co_cluster_counts: Dict[str, int] = {s: 0 for s in other_species}

    for tree in gene_trees:
        # For this tree, find which other species co-cluster (are sisters) with any copy of target
        tree_co_clusters: Set[str] = set()
        for leaf in tree.get_leaves():
            if leaf.name != species:
                continue
            parent = leaf.up
            if parent is None:
                continue
            for sibling in parent.get_children():
                if sibling is leaf:
                    continue
                if sibling.is_leaf() and sibling.name in other_species:
                    tree_co_clusters.add(sibling.name)

        for s in tree_co_clusters:
            co_cluster_counts[s] += 1

    n_trees = len(gene_trees)
    if n_trees == 0:
        return {s: 0.0 for s in other_species}

    return {s: co_cluster_counts[s] / n_trees for s in other_species}


def compute_clustering_summary(
    gene_trees: List[Tree],
    species: str,
    all_species: Set[str],
) -> List[float]:
    """Compute fixed-size summary statistics of the clustering profile.

    Returns 5 values: [mean, std, max, min, median] of co-clustering
    frequencies. This is independent of the number of species, unlike
    the full clustering profile.
    """
    profile = compute_clustering_profile(gene_trees, species, all_species)
    values = list(profile.values())
    if not values:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / n
    std = var ** 0.5
    sorted_vals = sorted(values)
    median = sorted_vals[n // 2] if n % 2 == 1 else (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2
    return [mean, std, max(values), min(values), median]


def compute_species_tree_edge_features(
    species_tree: Tree,
    gene_trees: List[Tree],
) -> Dict[int, dict]:
    """Compute edge-level features for a species tree using gene tree concordance.

    For each edge (parent -> child) in the species tree, compute:
        - concordance_factor: fraction of gene trees supporting this edge's bipartition
        - branch_length: the edge's branch length
        - clade_size: number of leaves below this edge
        - depth: number of edges from root to this edge

    Bipartition concordance is computed at the species level: group all copies of
    a species in a gene tree to the same side as the species appears in the
    species tree bipartition.

    Parameters
    ----------
    species_tree : Tree
        The species tree (ete3 Tree with branch lengths).
    gene_trees : list of Tree

    Returns
    -------
    dict mapping edge_index (int) to a dict of features.
    """
    all_species = frozenset(species_tree.get_leaf_names())
    n_trees = len(gene_trees)

    # Build list of (edge_index, node) for all non-root nodes
    edge_nodes = []
    idx = 0
    for node in species_tree.traverse("preorder"):
        if node.is_root():
            continue
        edge_nodes.append((idx, node))
        idx += 1

    # Pre-compute species-level bipartitions for the species tree edges
    # bipartition: (below_set, above_set)
    species_bipartitions = []
    for edge_idx, node in edge_nodes:
        below = frozenset(node.get_leaf_names())
        above = all_species - below
        species_bipartitions.append((edge_idx, below, above))

    # For each gene tree, compute its species-level bipartitions
    def gene_tree_species_bipartitions(gtree: Tree) -> Set[frozenset]:
        """Return the set of bipartitions (as frozensets of below-species) for a gene tree."""
        gt_all_species = frozenset(leaf.name for leaf in gtree.get_leaves())
        biparts: Set[frozenset] = set()
        for node in gtree.traverse("preorder"):
            if node.is_root():
                continue
            below_sp = frozenset(leaf.name for leaf in node.get_leaves())
            # Only include if it's a proper bipartition (not the full set)
            if below_sp and below_sp != gt_all_species:
                biparts.add(below_sp)
        return biparts

    # Count concordance for each edge
    concordance_counts = [0] * len(edge_nodes)
    for gtree in gene_trees:
        gt_biparts = gene_tree_species_bipartitions(gtree)
        gt_species = frozenset(leaf.name for leaf in gtree.get_leaves())
        for i, (edge_idx, below, above) in enumerate(species_bipartitions):
            # Restrict bipartition to the species present in this gene tree
            restricted_below = below & gt_species
            restricted_above = above & gt_species
            if not restricted_below or not restricted_above:
                # Uninformative for this gene tree
                continue
            # Check if this restricted bipartition is supported
            if restricted_below in gt_biparts or restricted_above in gt_biparts:
                concordance_counts[i] += 1

    # Compute depth for each edge (number of edges from root)
    def compute_depth(node: Tree) -> int:
        depth = 0
        current = node
        while not current.is_root():
            depth += 1
            current = current.up
        return depth

    # Build output
    result: Dict[int, dict] = {}
    for i, (edge_idx, node) in enumerate(edge_nodes):
        n_informative = sum(
            1 for gtree in gene_trees
            if (species_bipartitions[i][1] & frozenset(leaf.name for leaf in gtree.get_leaves()))
            and (species_bipartitions[i][2] & frozenset(leaf.name for leaf in gtree.get_leaves()))
        )
        cf = concordance_counts[i] / n_informative if n_informative > 0 else 0.0

        result[edge_idx] = {
            "concordance_factor": cf,
            "branch_length": node.dist,
            "clade_size": len(node.get_leaves()),
            "depth": compute_depth(node),
        }

    return result
