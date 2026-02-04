"""
Bipartition extraction from gene trees with polyploidy support.

A bipartition (split) divides the taxa into two groups based on an internal edge.
For MUL-trees with polyploidy, taxa can appear multiple times, so we use multisets.
"""

from collections import Counter
from typing import List, Tuple, Set, FrozenSet, Dict, Optional
import os
from pathlib import Path

try:
    from ete3 import Tree
except ImportError:
    raise ImportError("Please install ete3: pip install ete3")


# Type aliases
Multiset = FrozenSet[Tuple[str, int]]  # Species name + copy index
Bipartition = Tuple[Multiset, Multiset]  # Two sides of the split


def parse_tree(newick_string: str) -> Tree:
    """Parse a Newick string into an ete3 Tree object."""
    # Try different formats
    for fmt in [0, 1, 5, 9]:
        try:
            tree = Tree(newick_string, format=fmt)
            return tree
        except:
            continue
    raise ValueError(f"Could not parse tree: {newick_string[:100]}...")


def get_species_from_leaf(leaf_name: str) -> str:
    """
    Extract species name from leaf label.

    Handles common naming conventions:
    - "SpeciesA" -> "SpeciesA"
    - "SpeciesA_copy1" -> "SpeciesA"
    - "SpeciesA_1" -> "SpeciesA"
    """
    # If the leaf name ends with _copy# or _#, strip it
    name = leaf_name

    # Handle _copy1, _copy2, etc.
    if "_copy" in name:
        name = name.split("_copy")[0]

    # Handle trailing _1, _2 etc. (but be careful with species that have numbers)
    # Only strip if it's a single digit or two digits at the end
    parts = name.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit() and len(parts[1]) <= 2:
        name = parts[0]

    return name


def get_leaf_taxa(tree: Tree) -> List[Tuple[str, str]]:
    """
    Get list of (leaf_name, species_name) tuples for all leaves.

    Returns:
        List of (original_leaf_name, species_name) tuples
    """
    taxa = []
    for leaf in tree.get_leaves():
        species = get_species_from_leaf(leaf.name)
        taxa.append((leaf.name, species))
    return taxa


def extract_bipartitions_from_tree(
    tree: Tree,
    all_species: Optional[Set[str]] = None
) -> List[Tuple[Bipartition, float]]:
    """
    Extract all bipartitions from a single tree.

    Args:
        tree: ete3 Tree object
        all_species: Set of all species (for consistent encoding). If None, uses tree's species.

    Returns:
        List of (bipartition, branch_support) tuples.
        Each bipartition is a tuple of two frozensets (multisets represented as sets of (species, copy_index) tuples).
    """
    # Get leaf taxa with species mapping
    leaf_taxa = get_leaf_taxa(tree)

    # Build mapping: leaf_name -> (species, copy_index within this tree)
    species_copy_counter: Dict[str, int] = {}
    leaf_to_species_copy: Dict[str, Tuple[str, int]] = {}

    for leaf_name, species in leaf_taxa:
        copy_idx = species_copy_counter.get(species, 0)
        leaf_to_species_copy[leaf_name] = (species, copy_idx)
        species_copy_counter[species] = copy_idx + 1

    # All taxa in this tree (as multiset)
    all_taxa = frozenset(leaf_to_species_copy.values())

    bipartitions = []

    for node in tree.traverse("postorder"):
        # Skip leaves and root
        if node.is_leaf() or node.is_root():
            continue

        # Get taxa on this side of the split
        descendant_leaves = node.get_leaves()
        left_taxa = frozenset(
            leaf_to_species_copy[leaf.name]
            for leaf in descendant_leaves
        )

        # Other side is complement
        right_taxa = all_taxa - left_taxa

        # Skip trivial splits (one taxon on either side)
        if len(left_taxa) <= 1 or len(right_taxa) <= 1:
            continue

        # Canonical form: smaller set first (by size, then lexicographically)
        if len(left_taxa) < len(right_taxa):
            bipart = (left_taxa, right_taxa)
        elif len(left_taxa) > len(right_taxa):
            bipart = (right_taxa, left_taxa)
        else:
            # Same size: use lexicographic order
            if sorted(left_taxa) < sorted(right_taxa):
                bipart = (left_taxa, right_taxa)
            else:
                bipart = (right_taxa, left_taxa)

        # Get branch support (if available)
        support = node.support if hasattr(node, 'support') and node.support is not None else 1.0

        bipartitions.append((bipart, support))

    return bipartitions


def bipartition_to_species_sets(bipart: Bipartition) -> Tuple[Counter, Counter]:
    """
    Convert bipartition (with copy indices) to species counts (multisets).

    Args:
        bipart: Tuple of two frozensets of (species, copy_index) tuples

    Returns:
        Tuple of two Counters (species -> count)
    """
    left_counts = Counter(species for species, _ in bipart[0])
    right_counts = Counter(species for species, _ in bipart[1])
    return left_counts, right_counts


def bipartition_to_string(bipart: Bipartition) -> str:
    """Convert bipartition to human-readable string."""
    left_counts, right_counts = bipartition_to_species_sets(bipart)

    def format_side(counts: Counter) -> str:
        items = []
        for species, count in sorted(counts.items()):
            if count == 1:
                items.append(species)
            else:
                items.append(f"{species}×{count}")
        return ",".join(items)

    return f"{{{format_side(left_counts)}}} | {{{format_side(right_counts)}}}"


def extract_bipartitions_from_trees(
    trees: List[Tree],
    min_frequency: float = 0.0
) -> Dict[str, Dict]:
    """
    Extract all bipartitions from a list of gene trees.

    Args:
        trees: List of ete3 Tree objects
        min_frequency: Minimum frequency threshold (0-1) to include bipartition

    Returns:
        Dictionary mapping bipartition string to info dict containing:
        - 'bipartition': The Bipartition tuple
        - 'count': Number of trees containing this bipartition
        - 'frequency': count / len(trees)
        - 'supports': List of branch supports when present
        - 'mean_support': Average branch support
    """
    n_trees = len(trees)

    # Track bipartitions across trees
    # Key: canonical bipartition string (for deduplication)
    # Value: info dict
    bipart_info: Dict[str, Dict] = {}

    for tree in trees:
        tree_biparts = extract_bipartitions_from_tree(tree)

        for bipart, support in tree_biparts:
            # Convert to string key (for consistent hashing)
            bipart_str = bipartition_to_string(bipart)

            if bipart_str not in bipart_info:
                bipart_info[bipart_str] = {
                    'bipartition': bipart,
                    'count': 0,
                    'supports': [],
                    'species_left': bipartition_to_species_sets(bipart)[0],
                    'species_right': bipartition_to_species_sets(bipart)[1],
                }

            bipart_info[bipart_str]['count'] += 1
            bipart_info[bipart_str]['supports'].append(support)

    # Compute derived statistics
    for bipart_str, info in bipart_info.items():
        info['frequency'] = info['count'] / n_trees
        info['mean_support'] = sum(info['supports']) / len(info['supports']) if info['supports'] else 0.0

    # Filter by minimum frequency
    if min_frequency > 0:
        bipart_info = {
            k: v for k, v in bipart_info.items()
            if v['frequency'] >= min_frequency
        }

    return bipart_info


def load_trees_from_file(filepath: str) -> List[Tree]:
    """
    Load trees from a file (one Newick string per line).

    Args:
        filepath: Path to file containing Newick trees

    Returns:
        List of ete3 Tree objects
    """
    trees = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                try:
                    tree = parse_tree(line)
                    trees.append(tree)
                except Exception as e:
                    print(f"Warning: Could not parse tree: {e}")
    return trees


def extract_ground_truth_bipartitions(mul_tree: Tree) -> Set[str]:
    """
    Extract bipartitions from a ground truth MUL-tree.

    Args:
        mul_tree: The ground truth MUL-tree

    Returns:
        Set of bipartition strings that are in the true tree
    """
    biparts = extract_bipartitions_from_tree(mul_tree)
    return {bipartition_to_string(bp) for bp, _ in biparts}


# Command-line interface for testing
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python bipartition_extractor.py <trees_file> [ground_truth_tree]")
        sys.exit(1)

    trees_file = sys.argv[1]
    print(f"Loading trees from {trees_file}...")
    trees = load_trees_from_file(trees_file)
    print(f"Loaded {len(trees)} trees")

    print("\nExtracting bipartitions...")
    bipart_info = extract_bipartitions_from_trees(trees, min_frequency=0.01)

    print(f"\nFound {len(bipart_info)} unique bipartitions (frequency >= 1%)")

    # Show top 20 by frequency
    sorted_biparts = sorted(bipart_info.items(), key=lambda x: -x[1]['frequency'])
    print("\nTop 20 bipartitions by frequency:")
    for bipart_str, info in sorted_biparts[:20]:
        print(f"  {info['frequency']:.1%} (support={info['mean_support']:.2f}): {bipart_str}")

    # If ground truth provided, compute overlap
    if len(sys.argv) >= 3:
        gt_file = sys.argv[2]
        print(f"\nLoading ground truth from {gt_file}...")
        with open(gt_file, 'r') as f:
            gt_newick = f.read().strip()
        gt_tree = parse_tree(gt_newick)
        gt_biparts = extract_ground_truth_bipartitions(gt_tree)

        print(f"Ground truth has {len(gt_biparts)} bipartitions")

        # Check overlap
        predicted = set(bipart_info.keys())
        true_positive = predicted & gt_biparts
        false_positive = predicted - gt_biparts
        false_negative = gt_biparts - predicted

        print(f"\nOverlap analysis:")
        print(f"  True positives (in both): {len(true_positive)}")
        print(f"  False positives (in predicted, not in GT): {len(false_positive)}")
        print(f"  False negatives (in GT, not in predicted): {len(false_negative)}")

        if len(gt_biparts) > 0:
            recall = len(true_positive) / len(gt_biparts)
            print(f"  Recall: {recall:.1%}")
        if len(predicted) > 0:
            precision = len(true_positive) / len(predicted)
            print(f"  Precision: {precision:.1%}")
