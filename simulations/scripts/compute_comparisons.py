#!/usr/bin/env python3
"""
compute_comparisons.py - Comparison Engine

Computes all metrics between ground truth and inferred networks, with caching.

Usage:
    python compute_comparisons.py INVENTORY_CSV CACHE_DIR [--force-recompute] [--export CSV]
    python compute_comparisons.py inventory.csv cache/ --export comparisons.csv
"""

import os
import sys
import argparse
import pickle
import hashlib
from pathlib import Path
from typing import Dict, Optional, Tuple
from datetime import datetime
import pandas as pd
import numpy as np

# Import existing comparison tools
from reticulate_tree import ReticulateTree
from compare_reticulations import pairwise_compare


class ComparisonEngine:
    """Wrapper around existing comparison tools with caching"""

    def __init__(self, cache_dir: str, force_recompute: bool = False):
        """
        Initialize comparison engine

        Args:
            cache_dir: Directory for caching comparison results
            force_recompute: If True, ignore cache and recompute all
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.force_recompute = force_recompute

        # Statistics
        self.stats = {
            'total': 0,
            'success': 0,
            'from_cache': 0,
            'computed': 0,
            'failed': 0,
            'errors': []
        }

    def _file_hash(self, filepath: str) -> str:
        """Compute SHA256 hash of file for cache invalidation"""
        try:
            with open(filepath, 'rb') as f:
                return hashlib.sha256(f.read()).hexdigest()[:16]
        except Exception:
            return "missing"

    def _cache_key(self, network: str, config: str, method: str, replicate: int) -> str:
        """Generate cache filename"""
        return f"{network}_{config}_{method}_rep{replicate}.pkl"

    def _is_cache_valid(self, cache_path: Path, gt_path: str, inf_path: str) -> bool:
        """
        Check if cache is still valid

        Cache is valid if:
        1. Cache file exists
        2. Ground truth and inferred files haven't changed (based on file hashes)
        """
        if not cache_path.exists():
            return False

        try:
            with open(cache_path, 'rb') as f:
                cache_data = pickle.load(f)

            # Check file hashes
            current_gt_hash = self._file_hash(gt_path)
            current_inf_hash = self._file_hash(inf_path)

            return (cache_data['gt_hash'] == current_gt_hash and
                   cache_data['inf_hash'] == current_inf_hash)
        except Exception:
            return False

    def load_network(self, path: str, is_multree: bool = False) -> Optional[ReticulateTree]:
        """
        Load network from file using ReticulateTree

        Args:
            path: Path to network file
            is_multree: If True, treat as MUL-tree and convert to network

        Returns:
            ReticulateTree object or None if failed
        """
        try:
            with open(path, 'r') as f:
                newick_str = f.read().strip()

            # ReticulateTree handles format detection and MUL-tree conversion
            tree = ReticulateTree(newick_str)
            return tree
        except Exception as e:
            return None

    def compare_pair(self, gt_path: str, inf_path: str, network: str,
                    method: str) -> Dict:
        """
        Compare ground truth vs inferred network

        Args:
            gt_path: Path to ground truth network
            inf_path: Path to inferred network
            network: Network name (for error reporting)
            method: Method name (for error reporting)

        Returns:
            Dictionary with comparison results or error info
        """
        try:
            # Load networks
            gt_tree = self.load_network(gt_path)
            if gt_tree is None:
                return {
                    'status': 'ERROR',
                    'error': f'Failed to load ground truth: {gt_path}',
                    'metrics': None
                }

            inf_tree = self.load_network(inf_path, is_multree=True)
            if inf_tree is None:
                return {
                    'status': 'ERROR',
                    'error': f'Failed to load inferred network: {inf_path}',
                    'metrics': None
                }

            # Run comparison using existing tool
            metrics = pairwise_compare(gt_tree, inf_tree)

            return {
                'status': 'SUCCESS',
                'error': None,
                'metrics': metrics
            }

        except Exception as e:
            return {
                'status': 'ERROR',
                'error': str(e),
                'metrics': None
            }

    def flatten_metrics(self, metrics: Dict) -> Dict:
        """
        Flatten nested metric dictionaries

        Input:
            {
                'edit_distance': 0.234,
                'num_rets_diff': 2,
                'ploidy_diff': {'dist': 0.1, 'FP': 0.2, 'FN': 0.1, 'TP': 0.7},
                ...
            }

        Output:
            {
                'edit_distance': 0.234,
                'num_rets_diff': 2,
                'ploidy_diff.dist': 0.1,
                'ploidy_diff.FP': 0.2,
                ...
            }
        """
        flattened = {}

        for key, value in metrics.items():
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    flattened[f"{key}.{sub_key}"] = sub_value
            else:
                flattened[key] = value

        return flattened

    def compute_all_comparisons(self, inventory: pd.DataFrame) -> pd.DataFrame:
        """
        Compute comparisons for all valid entries in inventory

        Args:
            inventory: DataFrame from collect_results.py

        Returns:
            DataFrame with columns: network, config, method, replicate, metric, value, status
        """
        results = []
        self.stats['total'] = len(inventory)

        # Filter to only valid combinations (both files exist)
        valid_inventory = inventory[inventory['gt_exists'] & inventory['inferred_exists']].copy()

        print(f"\nProcessing {len(valid_inventory)}/{len(inventory)} valid combinations...")

        for idx, row in valid_inventory.iterrows():
            network = row['network']
            config = row['config']
            method = row['method']
            replicate = row['replicate']
            gt_path = row['gt_path']
            inf_path = row['inferred_path']

            # Check cache
            cache_filename = self._cache_key(network, config, method, replicate)
            cache_path = self.cache_dir / cache_filename

            if not self.force_recompute and self._is_cache_valid(cache_path, gt_path, inf_path):
                # Load from cache
                try:
                    with open(cache_path, 'rb') as f:
                        cache_data = pickle.load(f)
                    comparison = cache_data['comparison']
                    self.stats['from_cache'] += 1
                    source = 'cache'
                except Exception as e:
                    # Cache read failed, recompute
                    comparison = self.compare_pair(gt_path, inf_path, network, method)
                    self.stats['computed'] += 1
                    source = 'computed'
            else:
                # Compute comparison
                comparison = self.compare_pair(gt_path, inf_path, network, method)
                self.stats['computed'] += 1
                source = 'computed'

                # Save to cache if successful
                if comparison['status'] == 'SUCCESS':
                    try:
                        cache_data = {
                            'network': network,
                            'config': config,
                            'method': method,
                            'replicate': replicate,
                            'timestamp': datetime.now().isoformat(),
                            'gt_hash': self._file_hash(gt_path),
                            'inf_hash': self._file_hash(inf_path),
                            'comparison': comparison
                        }
                        with open(cache_path, 'wb') as f:
                            pickle.dump(cache_data, f)
                    except Exception as e:
                        pass  # Cache write failed, but computation succeeded

            # Record result
            if comparison['status'] == 'SUCCESS':
                self.stats['success'] += 1

                # Flatten metrics and add to results
                flat_metrics = self.flatten_metrics(comparison['metrics'])

                for metric_name, value in flat_metrics.items():
                    results.append({
                        'network': network,
                        'config': config,
                        'method': method,
                        'replicate': replicate,
                        'metric': metric_name,
                        'value': value,
                        'status': 'SUCCESS'
                    })
            else:
                self.stats['failed'] += 1
                self.stats['errors'].append({
                    'network': network,
                    'config': config,
                    'method': method,
                    'replicate': replicate,
                    'error': comparison['error']
                })

                # Add failed entry with NaN value
                results.append({
                    'network': network,
                    'config': config,
                    'method': method,
                    'replicate': replicate,
                    'metric': 'FAILED',
                    'value': np.nan,
                    'status': 'FAILED'
                })

            # Progress update
            if (idx + 1) % 20 == 0 or (idx + 1) == len(valid_inventory):
                print(f"  Processed {idx + 1}/{len(valid_inventory)} "
                      f"(Success: {self.stats['success']}, Failed: {self.stats['failed']}, "
                      f"Cached: {self.stats['from_cache']})")

        return pd.DataFrame(results)

    def print_statistics(self):
        """Print comparison statistics"""
        print(f"\n{'='*80}")
        print(f"Comparison Engine Statistics")
        print(f"{'='*80}")
        print(f"Total entries:      {self.stats['total']}")
        print(f"Successful:         {self.stats['success']} ({self.stats['success']/self.stats['total']*100:.1f}%)")
        print(f"Failed:             {self.stats['failed']} ({self.stats['failed']/self.stats['total']*100:.1f}%)")
        print(f"From cache:         {self.stats['from_cache']}")
        print(f"Newly computed:     {self.stats['computed']}")

        if self.stats['errors']:
            print(f"\nErrors ({len(self.stats['errors'])}):")
            for i, error in enumerate(self.stats['errors'][:10], 1):  # Show first 10
                print(f"  {i}. {error['network']} / {error['method']} / rep_{error['replicate']}")
                print(f"     Error: {error['error']}")

            if len(self.stats['errors']) > 10:
                print(f"  ... and {len(self.stats['errors']) - 10} more")

        print(f"{'='*80}\n")

    def write_comparison_report(self, output_file: str):
        """Write detailed comparison report to file"""
        with open(output_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("Comparison Report\n")
            f.write("="*80 + "\n\n")

            f.write(f"Total comparisons: {self.stats['total']}\n")
            f.write(f"Successful:        {self.stats['success']} ({self.stats['success']/self.stats['total']*100:.1f}%)\n")
            f.write(f"Failed:            {self.stats['failed']} ({self.stats['failed']/self.stats['total']*100:.1f}%)\n")
            f.write(f"From cache:        {self.stats['from_cache']}\n")
            f.write(f"Newly computed:    {self.stats['computed']}\n\n")

            if self.stats['errors']:
                f.write(f"\nFailed Comparisons ({len(self.stats['errors'])}):\n")
                f.write("-"*80 + "\n")
                for error in self.stats['errors']:
                    f.write(f"\n{error['network']} / {error['config']} / {error['method']} / rep_{error['replicate']}\n")
                    f.write(f"Error: {error['error']}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Compute network comparison metrics with caching',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compute comparisons from inventory
  %(prog)s inventory.csv cache/

  # Force recompute (ignore cache)
  %(prog)s inventory.csv cache/ --force-recompute

  # Export results and report
  %(prog)s inventory.csv cache/ --export comparisons.csv --report report.txt
        """
    )

    parser.add_argument('inventory', help='Inventory CSV file from collect_results.py')
    parser.add_argument('cache_dir', help='Directory for caching comparison results')
    parser.add_argument('--force-recompute', action='store_true',
                       help='Force recompute all comparisons (ignore cache)')
    parser.add_argument('--export', metavar='FILE',
                       help='Export comparison results to CSV file')
    parser.add_argument('--report', metavar='FILE',
                       help='Write detailed comparison report to file')

    args = parser.parse_args()

    # Load inventory
    try:
        inventory = pd.read_csv(args.inventory)
        print(f"Loaded inventory: {len(inventory)} combinations")
    except Exception as e:
        print(f"Error loading inventory: {e}", file=sys.stderr)
        sys.exit(1)

    # Create comparison engine
    engine = ComparisonEngine(args.cache_dir, force_recompute=args.force_recompute)

    # Compute comparisons
    comparisons_df = engine.compute_all_comparisons(inventory)

    # Print statistics
    engine.print_statistics()

    # Export if requested
    if args.export:
        comparisons_df.to_csv(args.export, index=False)
        print(f"Comparisons exported to: {args.export}")
        print(f"Total rows: {len(comparisons_df)}")

    # Write report if requested
    if args.report:
        engine.write_comparison_report(args.report)
        print(f"Comparison report written to: {args.report}")

    return comparisons_df


if __name__ == '__main__':
    comparisons_df = main()
