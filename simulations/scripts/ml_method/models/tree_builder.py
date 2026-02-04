"""
MUL-Tree reconstruction from bipartition scores.

Given predicted bipartition scores, reconstruct a valid MUL-tree
by selecting compatible bipartitions and building the tree structure.
"""

from collections import Counter
from typing import List, Dict, Set, Tuple, Optional, FrozenSet
import numpy as np

try:
    from ete3 import Tree
except ImportError:
    raise ImportError("Please install ete3: pip install ete3")

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data.bipartition_extractor import bipartition_to_species_sets, Bipartition


def are_bipartitions_compatible(bp1: Bipartition, bp2: Bipartition) -> bool:
    """
    Check if two bipartitions are compatible (can coexist in the same tree).

    Two bipartitions are compatible if, when considering all four intersections
    of their parts, at least one intersection is empty.
    """
    left1, right1 = bp1
    left2, right2 = bp2

    # Convert to species sets (ignoring copy indices)
    species_left1 = frozenset(sp for sp, _ in left1)
    species_right1 = frozenset(sp for sp, _ in right1)
    species_left2 = frozenset(sp for sp, _ in left2)
    species_right2 = frozenset(sp for sp, _ in right2)

    # Check the four intersections
    i1 = species_left1 & species_left2
    i2 = species_left1 & species_right2
    i3 = species_right1 & species_left2
    i4 = species_right1 & species_right2

    # Compatible if at least one is empty
    return not (i1 and i2 and i3 and i4)


def select_compatible_bipartitions(
    bipartitions: List[Tuple[str, Bipartition, float]],
    max_bipartitions: Optional[int] = None
) -> List[Tuple[str, Bipartition, float]]:
    """
    Select a maximal compatible set of bipartitions.

    Uses a greedy algorithm: sort by score, add if compatible with all selected.

    Args:
        bipartitions: List of (key, bipartition, score) tuples
        max_bipartitions: Maximum number to select (optional)

    Returns:
        List of selected (key, bipartition, score) tuples
    """
    # Sort by score descending
    sorted_bps = sorted(bipartitions, key=lambda x: -x[2])

    selected = []

    for key, bp, score in sorted_bps:
        # Check compatibility with all selected
        compatible = True
        for _, sel_bp, _ in selected:
            if not are_bipartitions_compatible(bp, sel_bp):
                compatible = False
                break

        if compatible:
            selected.append((key, bp, score))

            if max_bipartitions and len(selected) >= max_bipartitions:
                break

    return selected


def build_tree_from_bipartitions(
    bipartitions: List[Tuple[str, Bipartition, float]],
    all_taxa: Set[Tuple[str, int]]
) -> Tree:
    """
    Build a tree from a compatible set of bipartitions.

    Args:
        bipartitions: List of (key, bipartition, score) tuples (must be compatible)
        all_taxa: Set of all (species, copy_index) tuples

    Returns:
        ete3 Tree object
    """
    if not bipartitions:
        # Return star tree
        tree = Tree()
        for species, copy_idx in all_taxa:
            leaf_name = f"{species}" if copy_idx == 0 else f"{species}_copy{copy_idx + 1}"
            tree.add_child(name=leaf_name)
        return tree

    # Sort bipartitions by size of smaller side (larger first = deeper splits)
    sorted_bps = sorted(bipartitions, key=lambda x: -min(len(x[1][0]), len(x[1][1])))

    # Build tree iteratively
    # Start with a tree containing all taxa as leaves under root
    tree = Tree()

    # Create leaf nodes
    taxa_to_node = {}
    for species, copy_idx in all_taxa:
        leaf_name = f"{species}" if copy_idx == 0 else f"{species}_copy{copy_idx + 1}"
        leaf = tree.add_child(name=leaf_name)
        taxa_to_node[(species, copy_idx)] = leaf

    # Process each bipartition
    for key, (left, right), score in sorted_bps:
        # Find nodes for left side taxa
        left_nodes = [taxa_to_node.get(taxon) for taxon in left]
        left_nodes = [n for n in left_nodes if n is not None]

        if not left_nodes:
            continue

        # Find their common ancestor
        if len(left_nodes) == 1:
            mrca = left_nodes[0]
        else:
            mrca = tree.get_common_ancestor(left_nodes)

        # Check if these leaves are already in their own clade
        mrca_leaves = set(mrca.get_leaves())
        target_leaves = set(left_nodes)

        if mrca_leaves == target_leaves:
            # Already in correct clade, no need to modify
            continue

        # Need to create a new internal node
        # Detach left_nodes from their parents and attach to new node
        new_node = Tree()
        new_node.support = score

        for node in left_nodes:
            if node.up is not None:
                parent = node.up
                node.detach()
                new_node.add_child(node)

        # Add new_node to the tree
        tree.add_child(new_node)

        # Update taxa_to_node
        for taxon in left:
            if taxon in taxa_to_node:
                taxa_to_node[taxon] = new_node.search_nodes(name=taxa_to_node[taxon].name)[0] if taxa_to_node[taxon].name else taxa_to_node[taxon]

    return tree


def reconstruct_mul_tree(
    bipartition_scores: Dict[str, float],
    bipartition_info: Dict[str, Dict],
    threshold: float = 0.5,
    max_bipartitions: Optional[int] = None
) -> Tree:
    """
    Reconstruct MUL-tree from predicted bipartition scores.

    Args:
        bipartition_scores: Dict mapping bipartition key to predicted probability
        bipartition_info: Dict mapping bipartition key to info (from extract_bipartitions_from_trees)
        threshold: Probability threshold for including bipartition
        max_bipartitions: Maximum number of bipartitions to include

    Returns:
        Reconstructed MUL-tree as ete3 Tree
    """
    # Filter by threshold
    candidates = [
        (key, bipartition_info[key]['bipartition'], score)
        for key, score in bipartition_scores.items()
        if score >= threshold and key in bipartition_info
    ]

    if not candidates:
        print(f"Warning: No bipartitions above threshold {threshold}")
        # Return empty tree
        return Tree()

    # Select compatible bipartitions
    selected = select_compatible_bipartitions(candidates, max_bipartitions)

    # Get all taxa from bipartitions
    all_taxa = set()
    for _, bp, _ in selected:
        all_taxa.update(bp[0])
        all_taxa.update(bp[1])

    # Build tree
    tree = build_tree_from_bipartitions(selected, all_taxa)

    return tree


def newick_from_bipartitions(
    bipartition_scores: Dict[str, float],
    bipartition_info: Dict[str, Dict],
    threshold: float = 0.5
) -> str:
    """
    Get Newick string from bipartition scores.

    Args:
        bipartition_scores: Predicted scores
        bipartition_info: Bipartition metadata
        threshold: Score threshold

    Returns:
        Newick string
    """
    tree = reconstruct_mul_tree(bipartition_scores, bipartition_info, threshold)
    return tree.write(format=1)


class MULTreeBuilder:
    """
    Class for building MUL-trees from bipartition scores.

    Combines bipartition selection with copy number prediction.
    """

    def __init__(
        self,
        threshold: float = 0.5,
        use_copy_numbers: bool = True
    ):
        """
        Initialize the tree builder.

        Args:
            threshold: Score threshold for bipartition selection
            use_copy_numbers: Whether to use separate copy number predictions
        """
        self.threshold = threshold
        self.use_copy_numbers = use_copy_numbers

    def build(
        self,
        bipartition_scores: Dict[str, float],
        bipartition_info: Dict[str, Dict],
        copy_numbers: Optional[Dict[str, int]] = None
    ) -> Tree:
        """
        Build MUL-tree from predictions.

        Args:
            bipartition_scores: Predicted bipartition probabilities
            bipartition_info: Bipartition metadata
            copy_numbers: Optional dict mapping species to copy count

        Returns:
            Reconstructed MUL-tree
        """
        # Reconstruct base tree from bipartitions
        tree = reconstruct_mul_tree(
            bipartition_scores,
            bipartition_info,
            threshold=self.threshold
        )

        # TODO: Integrate copy number predictions
        # This would expand species with predicted copies

        return tree

    def build_and_compare(
        self,
        bipartition_scores: Dict[str, float],
        bipartition_info: Dict[str, Dict],
        ground_truth: Tree
    ) -> Dict[str, float]:
        """
        Build tree and compare to ground truth.

        Args:
            bipartition_scores: Predicted scores
            bipartition_info: Bipartition metadata
            ground_truth: True MUL-tree

        Returns:
            Dict of comparison metrics
        """
        predicted = self.build(bipartition_scores, bipartition_info)

        # Compute comparison metrics
        metrics = {}

        # Robinson-Foulds distance
        try:
            rf, max_rf, _, _, _, _, _ = predicted.robinson_foulds(
                ground_truth,
                unrooted_trees=True
            )
            metrics['rf_distance'] = rf
            metrics['rf_normalized'] = rf / max_rf if max_rf > 0 else 0
        except Exception as e:
            metrics['rf_distance'] = -1
            metrics['rf_normalized'] = -1

        return metrics


# Command-line testing
if __name__ == "__main__":
    print("Testing MUL-tree reconstruction...")

    # Create simple test case
    # Taxa: A (2 copies), B, C
    taxa = {('A', 0), ('A', 1), ('B', 0), ('C', 0)}

    # Bipartitions
    bp1 = (frozenset({('A', 0), ('A', 1)}), frozenset({('B', 0), ('C', 0)}))
    bp2 = (frozenset({('B', 0)}), frozenset({('A', 0), ('A', 1), ('C', 0)}))

    bipartitions = [
        ('bp1', bp1, 0.9),
        ('bp2', bp2, 0.7),
    ]

    # Test compatibility
    print(f"BP1 and BP2 compatible: {are_bipartitions_compatible(bp1, bp2)}")

    # Test selection
    selected = select_compatible_bipartitions(bipartitions)
    print(f"\nSelected {len(selected)} compatible bipartitions")

    # Test tree building
    tree = build_tree_from_bipartitions(selected, taxa)
    print(f"\nReconstructed tree: {tree.write(format=9)}")

    print("\nTest passed!")
