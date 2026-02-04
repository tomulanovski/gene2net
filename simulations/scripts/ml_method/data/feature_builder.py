"""
Feature computation for bipartitions.

Computes features for each candidate bipartition extracted from gene trees.
Features are designed to predict whether a bipartition is in the true MUL-tree.
"""

from collections import Counter
from typing import List, Dict, Set, Tuple, Optional
import numpy as np
from pathlib import Path

try:
    from ete3 import Tree
except ImportError:
    raise ImportError("Please install ete3: pip install ete3")

from .bipartition_extractor import (
    extract_bipartitions_from_tree,
    extract_bipartitions_from_trees,
    bipartition_to_string,
    bipartition_to_species_sets,
    Bipartition,
    get_leaf_taxa,
)


def are_bipartitions_compatible(bp1: Bipartition, bp2: Bipartition) -> bool:
    """
    Check if two bipartitions are compatible (can coexist in the same tree).

    Two bipartitions are compatible if, for each pair of sets from bp1 and bp2,
    at least one of the four intersections is empty.
    """
    left1, right1 = bp1
    left2, right2 = bp2

    # Convert to species sets (ignoring copy indices for compatibility check)
    species_left1 = frozenset(sp for sp, _ in left1)
    species_right1 = frozenset(sp for sp, _ in right1)
    species_left2 = frozenset(sp for sp, _ in left2)
    species_right2 = frozenset(sp for sp, _ in right2)

    # Check the four intersections
    # Compatible if at least one is empty
    i1 = species_left1 & species_left2
    i2 = species_left1 & species_right2
    i3 = species_right1 & species_left2
    i4 = species_right1 & species_right2

    # All four non-empty means incompatible
    return not (i1 and i2 and i3 and i4)


def compute_conflict_rate(
    target_bipart: str,
    all_biparts: Dict[str, Dict],
    target_info: Dict
) -> float:
    """
    Compute the fraction of bipartitions that conflict with the target.

    Args:
        target_bipart: String key of target bipartition
        all_biparts: Dict of all bipartitions
        target_info: Info dict for target bipartition

    Returns:
        Fraction of bipartitions that are incompatible (0-1)
    """
    target_bp = target_info['bipartition']
    n_conflicts = 0
    n_total = 0

    for other_key, other_info in all_biparts.items():
        if other_key == target_bipart:
            continue

        other_bp = other_info['bipartition']
        n_total += 1

        if not are_bipartitions_compatible(target_bp, other_bp):
            # Weight by frequency of conflicting bipartition
            n_conflicts += other_info['frequency']

    return n_conflicts / n_total if n_total > 0 else 0.0


def compute_bipartition_features(
    trees: List[Tree],
    bipart_info: Dict[str, Dict],
    all_species: Optional[Set[str]] = None
) -> Dict[str, np.ndarray]:
    """
    Compute features for all bipartitions.

    Args:
        trees: List of gene trees
        bipart_info: Dict from extract_bipartitions_from_trees()
        all_species: Set of all species (optional)

    Returns:
        Dict mapping bipartition string to feature vector
    """
    n_trees = len(trees)

    # Pre-compute species info if not provided
    if all_species is None:
        all_species = set()
        for tree in trees:
            taxa = get_leaf_taxa(tree)
            all_species.update(sp for _, sp in taxa)

    n_species = len(all_species)

    # Pre-compute per-tree bipartitions for consistency checks
    tree_biparts = []
    for tree in trees:
        tree_bps = extract_bipartitions_from_tree(tree)
        tree_bipart_strs = {bipartition_to_string(bp) for bp, _ in tree_bps}
        tree_biparts.append(tree_bipart_strs)

    features = {}

    for bipart_str, info in bipart_info.items():
        feat = compute_single_bipartition_features(
            bipart_str, info, bipart_info,
            trees, tree_biparts, n_species
        )
        features[bipart_str] = feat

    return features


def compute_single_bipartition_features(
    bipart_str: str,
    info: Dict,
    all_biparts: Dict[str, Dict],
    trees: List[Tree],
    tree_biparts: List[Set[str]],
    n_species: int
) -> np.ndarray:
    """
    Compute feature vector for a single bipartition.

    Returns:
        Feature vector as numpy array
    """
    n_trees = len(trees)

    # ========== FREQUENCY-BASED FEATURES ==========

    # Basic frequency
    frequency = info['frequency']

    # Weighted frequency (by support) - already computed
    supports = info['supports']
    if supports:
        weighted_freq = sum(s * (1/n_trees) for s in supports)
        weighted_freq = weighted_freq / frequency if frequency > 0 else 0
    else:
        weighted_freq = frequency

    # Presence variance (variance in presence across trees)
    presence = [1 if bipart_str in tb else 0 for tb in tree_biparts]
    presence_var = np.var(presence)

    # ========== SUPPORT-BASED FEATURES ==========

    if supports:
        mean_support = info['mean_support']
        min_support = min(supports)
        max_support = max(supports)
        support_std = np.std(supports)
    else:
        mean_support = 0.0
        min_support = 0.0
        max_support = 0.0
        support_std = 0.0

    # ========== CONSISTENCY FEATURES ==========

    # Conflict rate with other bipartitions
    conflict_rate = compute_conflict_rate(bipart_str, all_biparts, info)

    # Consistency: how often this bipartition co-occurs with other high-freq bipartitions
    high_freq_biparts = [k for k, v in all_biparts.items()
                        if v['frequency'] > 0.5 and k != bipart_str]
    cooccurrence = 0
    for other_str in high_freq_biparts:
        # Count trees where both appear
        co_trees = sum(1 for tb in tree_biparts
                      if bipart_str in tb and other_str in tb)
        cooccurrence += co_trees / n_trees
    consistency_score = cooccurrence / len(high_freq_biparts) if high_freq_biparts else 0.0

    # ========== POLYPLOIDY-SPECIFIC FEATURES ==========

    left_counts, right_counts = info['species_left'], info['species_right']

    # Does this bipartition involve polyploid species (species with >1 copy)?
    left_polyploid = sum(1 for c in left_counts.values() if c > 1)
    right_polyploid = sum(1 for c in right_counts.values() if c > 1)
    involves_polyploid = 1 if (left_polyploid > 0 or right_polyploid > 0) else 0

    # Total copy count on each side
    left_total = sum(left_counts.values())
    right_total = sum(right_counts.values())

    # Copy ratio (balance of copies between sides)
    copy_ratio = min(left_total, right_total) / max(left_total, right_total) if max(left_total, right_total) > 0 else 0

    # Species count on each side
    left_species = len(left_counts)
    right_species = len(right_counts)

    # Monophyly indicator: if polyploid, are all copies on the same side?
    # (This suggests autopolyploidy if true)
    all_species_in_bipart = set(left_counts.keys()) | set(right_counts.keys())
    monophyly_score = 0
    n_polyploids_checked = 0
    for sp in all_species_in_bipart:
        left_c = left_counts.get(sp, 0)
        right_c = right_counts.get(sp, 0)
        if left_c + right_c > 1:  # Polyploid
            n_polyploids_checked += 1
            if left_c == 0 or right_c == 0:
                monophyly_score += 1  # All copies on one side
    monophyly_score = monophyly_score / n_polyploids_checked if n_polyploids_checked > 0 else 0

    # ========== STRUCTURAL FEATURES ==========

    # Split size (smaller side)
    split_size = min(left_total, right_total)
    split_size_ratio = split_size / (left_total + right_total) if (left_total + right_total) > 0 else 0

    # Depth indicator (shallow splits have larger sides)
    # Normalized by total species count
    depth_indicator = 1 - (split_size / n_species) if n_species > 0 else 0

    # ========== ASSEMBLE FEATURE VECTOR ==========

    feature_vector = np.array([
        # Frequency features (3)
        frequency,
        weighted_freq,
        presence_var,

        # Support features (4)
        mean_support,
        min_support,
        max_support,
        support_std,

        # Consistency features (2)
        conflict_rate,
        consistency_score,

        # Polyploidy features (5)
        involves_polyploid,
        copy_ratio,
        left_polyploid,
        right_polyploid,
        monophyly_score,

        # Structural features (4)
        split_size,
        split_size_ratio,
        depth_indicator,
        left_species + right_species,  # Total unique species
    ], dtype=np.float32)

    return feature_vector


def get_feature_names() -> List[str]:
    """Get names of all features in the feature vector."""
    return [
        # Frequency features
        "frequency",
        "weighted_frequency",
        "presence_variance",

        # Support features
        "mean_support",
        "min_support",
        "max_support",
        "support_std",

        # Consistency features
        "conflict_rate",
        "consistency_score",

        # Polyploidy features
        "involves_polyploid",
        "copy_ratio",
        "left_polyploid_count",
        "right_polyploid_count",
        "monophyly_score",

        # Structural features
        "split_size",
        "split_size_ratio",
        "depth_indicator",
        "total_species",
    ]


def build_dataset(
    trees: List[Tree],
    ground_truth_biparts: Optional[Set[str]] = None,
    min_frequency: float = 0.0
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Build feature matrix and labels for training.

    Args:
        trees: List of gene trees
        ground_truth_biparts: Set of bipartition strings that are in the true tree
        min_frequency: Minimum frequency to include bipartition

    Returns:
        X: Feature matrix [n_bipartitions, n_features]
        y: Labels [n_bipartitions] (1 if in ground truth, 0 otherwise)
        bipart_keys: List of bipartition strings (for reference)
    """
    # Extract bipartitions
    bipart_info = extract_bipartitions_from_trees(trees, min_frequency=min_frequency)

    # Compute features
    features = compute_bipartition_features(trees, bipart_info)

    # Build arrays
    bipart_keys = list(features.keys())
    X = np.stack([features[k] for k in bipart_keys])

    if ground_truth_biparts is not None:
        y = np.array([1 if k in ground_truth_biparts else 0 for k in bipart_keys], dtype=np.float32)
    else:
        y = np.zeros(len(bipart_keys), dtype=np.float32)

    return X, y, bipart_keys


# Command-line interface for testing
if __name__ == "__main__":
    import sys
    from .bipartition_extractor import load_trees_from_file, parse_tree, extract_ground_truth_bipartitions

    if len(sys.argv) < 2:
        print("Usage: python feature_builder.py <trees_file> [ground_truth_tree]")
        sys.exit(1)

    trees_file = sys.argv[1]
    print(f"Loading trees from {trees_file}...")
    trees = load_trees_from_file(trees_file)
    print(f"Loaded {len(trees)} trees")

    # Load ground truth if provided
    ground_truth = None
    if len(sys.argv) >= 3:
        gt_file = sys.argv[2]
        print(f"Loading ground truth from {gt_file}...")
        with open(gt_file, 'r') as f:
            gt_newick = f.read().strip()
        gt_tree = parse_tree(gt_newick)
        ground_truth = extract_ground_truth_bipartitions(gt_tree)
        print(f"Ground truth has {len(ground_truth)} bipartitions")

    print("\nBuilding dataset...")
    X, y, bipart_keys = build_dataset(trees, ground_truth, min_frequency=0.01)

    print(f"\nDataset shape: {X.shape}")
    print(f"Positive examples: {y.sum():.0f} / {len(y)}")

    # Show feature statistics
    feature_names = get_feature_names()
    print("\nFeature statistics:")
    for i, name in enumerate(feature_names):
        print(f"  {name}: mean={X[:, i].mean():.3f}, std={X[:, i].std():.3f}")

    # Show class balance
    if ground_truth is not None:
        print(f"\nClass balance: {y.sum():.0f} positive, {len(y) - y.sum():.0f} negative")
        print(f"Positive rate: {y.mean():.1%}")
