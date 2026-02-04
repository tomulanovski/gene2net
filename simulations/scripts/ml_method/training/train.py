"""
Main training script for the MUL-tree inference model.

Trains a bipartition scorer on simulated data and evaluates on held-out networks.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from data.bipartition_extractor import (
    load_trees_from_file,
    parse_tree,
    extract_bipartitions_from_trees,
    extract_ground_truth_bipartitions,
)
from data.feature_builder import build_dataset, get_feature_names
from models.xgb_scorer import BipartitionScorer, train_scorer
from models.tree_builder import MULTreeBuilder, reconstruct_mul_tree


def load_network_data(
    network_dir: str,
    gene_trees_file: str = "gene_trees.tre",
    ground_truth_file: str = "ground_truth.tre"
) -> Tuple[List, set, str]:
    """
    Load gene trees and ground truth for a single network.

    Args:
        network_dir: Directory containing network data
        gene_trees_file: Filename of gene trees
        ground_truth_file: Filename of ground truth MUL-tree

    Returns:
        (gene_trees, ground_truth_bipartitions, ground_truth_newick)
    """
    network_path = Path(network_dir)

    # Load gene trees
    gt_path = network_path / gene_trees_file
    if not gt_path.exists():
        raise FileNotFoundError(f"Gene trees not found: {gt_path}")
    gene_trees = load_trees_from_file(str(gt_path))

    # Load ground truth
    truth_path = network_path / ground_truth_file
    if not truth_path.exists():
        raise FileNotFoundError(f"Ground truth not found: {truth_path}")
    with open(truth_path, 'r') as f:
        gt_newick = f.read().strip()
    gt_tree = parse_tree(gt_newick)
    gt_biparts = extract_ground_truth_bipartitions(gt_tree)

    return gene_trees, gt_biparts, gt_newick


def train_on_networks(
    network_dirs: List[str],
    gene_trees_file: str = "gene_trees.tre",
    ground_truth_file: str = "ground_truth.tre",
    min_frequency: float = 0.01,
    **model_kwargs
) -> Tuple[BipartitionScorer, Dict]:
    """
    Train bipartition scorer on multiple networks.

    Args:
        network_dirs: List of directories containing network data
        gene_trees_file: Filename of gene trees in each directory
        ground_truth_file: Filename of ground truth in each directory
        min_frequency: Minimum bipartition frequency to include
        **model_kwargs: Arguments for BipartitionScorer

    Returns:
        (trained_scorer, training_info)
    """
    all_X = []
    all_y = []
    feature_names = get_feature_names()

    print(f"Loading data from {len(network_dirs)} networks...")

    for i, network_dir in enumerate(network_dirs):
        print(f"  [{i+1}/{len(network_dirs)}] {Path(network_dir).name}")

        try:
            gene_trees, gt_biparts, _ = load_network_data(
                network_dir, gene_trees_file, ground_truth_file
            )

            X, y, bipart_keys = build_dataset(
                gene_trees, gt_biparts, min_frequency=min_frequency
            )

            all_X.append(X)
            all_y.append(y)

            print(f"    {len(gene_trees)} trees, {len(X)} bipartitions, {y.sum():.0f} positive")

        except Exception as e:
            print(f"    Error: {e}")
            continue

    if not all_X:
        raise ValueError("No data loaded successfully")

    # Concatenate all data
    X_train = np.vstack(all_X)
    y_train = np.concatenate(all_y)

    print(f"\nTotal training data: {X_train.shape[0]} bipartitions, {y_train.sum():.0f} positive ({y_train.mean():.1%})")

    # Train model
    print("\nTraining model...")
    scorer = BipartitionScorer(**model_kwargs)
    scorer.fit(X_train, y_train, feature_names=feature_names)

    # Get feature importance
    importance = scorer.get_feature_importance()
    sorted_importance = sorted(importance.items(), key=lambda x: -x[1])

    print("\nTop 10 features by importance:")
    for name, value in sorted_importance[:10]:
        print(f"  {name}: {value:.2f}")

    training_info = {
        'n_networks': len(network_dirs),
        'n_bipartitions': len(X_train),
        'n_positive': int(y_train.sum()),
        'positive_rate': float(y_train.mean()),
        'feature_importance': dict(sorted_importance),
    }

    return scorer, training_info


def evaluate_on_network(
    scorer: BipartitionScorer,
    network_dir: str,
    gene_trees_file: str = "gene_trees.tre",
    ground_truth_file: str = "ground_truth.tre",
    min_frequency: float = 0.01,
    threshold: float = 0.5
) -> Dict:
    """
    Evaluate trained scorer on a single network.

    Args:
        scorer: Trained BipartitionScorer
        network_dir: Directory containing network data
        gene_trees_file: Filename of gene trees
        ground_truth_file: Filename of ground truth
        min_frequency: Minimum bipartition frequency
        threshold: Classification threshold

    Returns:
        Dict of evaluation metrics
    """
    # Load data
    gene_trees, gt_biparts, gt_newick = load_network_data(
        network_dir, gene_trees_file, ground_truth_file
    )

    # Extract bipartitions and features
    bipart_info = extract_bipartitions_from_trees(gene_trees, min_frequency=min_frequency)
    X, y, bipart_keys = build_dataset(gene_trees, gt_biparts, min_frequency=min_frequency)

    # Predict
    probs = scorer.predict_proba(X)
    preds = (probs >= threshold).astype(int)

    # Compute metrics
    tp = ((preds == 1) & (y == 1)).sum()
    fp = ((preds == 1) & (y == 0)).sum()
    fn = ((preds == 0) & (y == 1)).sum()

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    # Build predicted tree and compare
    bipart_scores = {k: float(p) for k, p in zip(bipart_keys, probs)}
    builder = MULTreeBuilder(threshold=threshold)

    metrics = {
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'n_true_positive': int(tp),
        'n_false_positive': int(fp),
        'n_false_negative': int(fn),
        'n_ground_truth_biparts': int(y.sum()),
        'n_predicted_biparts': int(preds.sum()),
    }

    return metrics


def cross_validate_networks(
    network_dirs: List[str],
    gene_trees_file: str = "gene_trees.tre",
    ground_truth_file: str = "ground_truth.tre",
    min_frequency: float = 0.01,
    **model_kwargs
) -> Dict:
    """
    Leave-one-network-out cross-validation.

    Args:
        network_dirs: List of network directories
        gene_trees_file: Gene trees filename
        ground_truth_file: Ground truth filename
        min_frequency: Minimum bipartition frequency
        **model_kwargs: Arguments for BipartitionScorer

    Returns:
        Dict of CV results
    """
    results = []

    for i, test_dir in enumerate(network_dirs):
        print(f"\n=== Fold {i+1}/{len(network_dirs)}: Testing on {Path(test_dir).name} ===")

        # Train on all except test
        train_dirs = [d for j, d in enumerate(network_dirs) if j != i]

        try:
            scorer, _ = train_on_networks(
                train_dirs,
                gene_trees_file=gene_trees_file,
                ground_truth_file=ground_truth_file,
                min_frequency=min_frequency,
                **model_kwargs
            )

            # Evaluate on test
            metrics = evaluate_on_network(
                scorer, test_dir,
                gene_trees_file=gene_trees_file,
                ground_truth_file=ground_truth_file,
                min_frequency=min_frequency
            )

            metrics['network'] = Path(test_dir).name
            results.append(metrics)

            print(f"  Precision: {metrics['precision']:.3f}")
            print(f"  Recall: {metrics['recall']:.3f}")
            print(f"  F1: {metrics['f1']:.3f}")

        except Exception as e:
            print(f"  Error: {e}")
            continue

    # Aggregate results
    if results:
        avg_metrics = {
            'mean_precision': np.mean([r['precision'] for r in results]),
            'mean_recall': np.mean([r['recall'] for r in results]),
            'mean_f1': np.mean([r['f1'] for r in results]),
            'std_f1': np.std([r['f1'] for r in results]),
        }
    else:
        avg_metrics = {}

    return {
        'per_network': results,
        'aggregate': avg_metrics
    }


def main():
    parser = argparse.ArgumentParser(description='Train bipartition scorer for MUL-tree inference')

    parser.add_argument('--data-dir', type=str, required=True,
                       help='Directory containing network subdirectories')
    parser.add_argument('--output-dir', type=str, default='output',
                       help='Directory for output files')
    parser.add_argument('--gene-trees-file', type=str, default='gene_trees.tre',
                       help='Filename of gene trees in each network directory')
    parser.add_argument('--ground-truth-file', type=str, default='ground_truth.tre',
                       help='Filename of ground truth in each network directory')
    parser.add_argument('--min-frequency', type=float, default=0.01,
                       help='Minimum bipartition frequency')
    parser.add_argument('--cv', action='store_true',
                       help='Run leave-one-out cross-validation')
    parser.add_argument('--n-estimators', type=int, default=100,
                       help='Number of XGBoost estimators')
    parser.add_argument('--max-depth', type=int, default=6,
                       help='Maximum tree depth')

    args = parser.parse_args()

    # Find network directories
    data_dir = Path(args.data_dir)
    network_dirs = [str(d) for d in data_dir.iterdir() if d.is_dir()]

    if not network_dirs:
        print(f"No network directories found in {data_dir}")
        sys.exit(1)

    print(f"Found {len(network_dirs)} networks")

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model_kwargs = {
        'n_estimators': args.n_estimators,
        'max_depth': args.max_depth,
    }

    if args.cv:
        # Cross-validation
        results = cross_validate_networks(
            network_dirs,
            gene_trees_file=args.gene_trees_file,
            ground_truth_file=args.ground_truth_file,
            min_frequency=args.min_frequency,
            **model_kwargs
        )

        print("\n=== Cross-Validation Results ===")
        print(f"Mean F1: {results['aggregate'].get('mean_f1', 0):.3f} ± {results['aggregate'].get('std_f1', 0):.3f}")

        # Save results
        with open(output_dir / 'cv_results.json', 'w') as f:
            json.dump(results, f, indent=2)

    else:
        # Train on all data
        scorer, training_info = train_on_networks(
            network_dirs,
            gene_trees_file=args.gene_trees_file,
            ground_truth_file=args.ground_truth_file,
            min_frequency=args.min_frequency,
            **model_kwargs
        )

        # Save model
        scorer.save(str(output_dir / 'bipartition_scorer.pkl'))

        # Save training info
        with open(output_dir / 'training_info.json', 'w') as f:
            json.dump(training_info, f, indent=2)

        print(f"\nModel saved to {output_dir}")


if __name__ == "__main__":
    main()
