"""Data loading and feature extraction utilities."""

from .bipartition_extractor import (
    parse_tree,
    extract_bipartitions_from_tree,
    extract_bipartitions_from_trees,
    extract_ground_truth_bipartitions,
    load_trees_from_file,
    bipartition_to_string,
    bipartition_to_species_sets,
    get_leaf_taxa,
    get_species_from_leaf,
)

from .feature_builder import (
    compute_bipartition_features,
    compute_single_bipartition_features,
    build_dataset,
    get_feature_names,
    are_bipartitions_compatible,
)

__all__ = [
    'parse_tree',
    'extract_bipartitions_from_tree',
    'extract_bipartitions_from_trees',
    'extract_ground_truth_bipartitions',
    'load_trees_from_file',
    'bipartition_to_string',
    'bipartition_to_species_sets',
    'get_leaf_taxa',
    'get_species_from_leaf',
    'compute_bipartition_features',
    'compute_single_bipartition_features',
    'build_dataset',
    'get_feature_names',
    'are_bipartitions_compatible',
]
