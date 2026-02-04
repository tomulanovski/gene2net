"""
ML Method for MUL-Tree Inference from Gene Trees.

This package implements a machine learning approach to infer MUL-trees
(multi-labeled trees representing polyploid species relationships)
from collections of gene trees.

Key modules:
- data: Bipartition extraction and feature computation
- models: XGBoost scorer and tree reconstruction
- training: Training and evaluation utilities
- analysis: Interpretability and visualization

Main approach:
1. Extract bipartitions from gene trees
2. Compute features for each bipartition
3. Train a classifier to predict P(bipartition ∈ true MUL-tree)
4. Reconstruct tree from high-scoring compatible bipartitions

Usage:
    from ml_method.data import load_trees_from_file, build_dataset
    from ml_method.models import BipartitionScorer, reconstruct_mul_tree
    from ml_method.training import train_on_networks
"""

__version__ = "0.1.0"
