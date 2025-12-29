#!/usr/bin/env python3
"""
compute_comparisons.py - Pairwise Comparison Engine for Real Data

Computes all pairwise metrics between different methods (and published networks), with caching.

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

# Import existing comparison tools from simulations/scripts
sys.path.insert(0, str(Path(__file__).parent.parent / 'simulations' / 'scripts'))
from reticulate_tree import ReticulateTree
from compare_reticulations import pairwise_compare


class ComparisonEngine:
    """Wrapper around existing comparison tools with caching for pairwise method comparisons"""

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

    def _cache_key(self, network: str, method1: str, method2: str) -> str:
        """Generate cache filename for pairwise comparison"""
        # Sort methods alphabetically to ensure consistent cache keys
        methods = tuple(sorted([method1, method2]))
        return f"{network}_{methods[0]}_vs_{methods[1]}.pkl"

    def _is_cache_valid(self, cache_path: Path, path1: str, path2: str) -> bool:
        """
        Check if cache is still valid

        Cache is valid if:
        1. Cache file exists
        2. Both network files haven't changed (based on file hashes)
        """
        if not cache_path.exists():
            return False

        try:
            with open(cache_path, 'rb') as f:
                cache_data = pickle.load(f)

            # Check file hashes
            current_hash1 = self._file_hash(path1)
            current_hash2 = self._file_hash(path2)

            return (cache_data['hash1'] == current_hash1 and
                   cache_data['hash2'] == current_hash2)
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

    def compare_pair(self, path1: str, path2: str, network: str,
                    method1: str, method2: str) -> Dict:
        """
        Compare two method outputs (pairwise comparison)

        Args:
            path1: Path to first method's network
            path2: Path to second method's network
            network: Network name (for error reporting)
            method1: First method name (for error reporting)
            method2: Second method name (for error reporting)

        Returns:
            Dictionary with comparison results or error info
        """
        try:
            # Load networks
            tree1 = self.load_network(path1)
            if tree1 is None:
                return {
                    'status': 'ERROR',
                    'error': f'Failed to load {method1} network: {path1}',
                    'metrics': None
                }

            tree2 = self.load_network(path2)
            if tree2 is None:
                return {
                    'status': 'ERROR',
                    'error': f'Failed to load {method2} network: {path2}',
                    'metrics': None
                }

            # Run comparison using existing tool
            metrics = pairwise_compare(tree1, tree2)

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
        Compute pairwise comparisons for all valid method pairs per network

        Args:
            inventory: DataFrame from collect_results.py with columns:
                      network, method, network_path, exists

        Returns:
            DataFrame with columns: network, method1, method2, metric, value, status
        """
        results = []
        self.stats['total'] = 0

        # Group by network
        for network in inventory['network'].unique():
            network_data = inventory[inventory['network'] == network]
            
            # Get available methods for this network
            available_methods = network_data[network_data['exists']].copy()
            
            if len(available_methods) < 2:
                # Need at least 2 methods to compare
                continue

            # Create all pairwise comparisons
            method_list = available_methods['method'].tolist()
            method_paths = dict(zip(available_methods['method'], available_methods['network_path']))

            for i, method1 in enumerate(method_list):
                for method2 in method_list[i+1:]:
                    path1 = method_paths[method1]
                    path2 = method_paths[method2]

                    self.stats['total'] += 1

                    # Check cache
                    cache_filename = self._cache_key(network, method1, method2)
                    cache_path = self.cache_dir / cache_filename

                    if not self.force_recompute and self._is_cache_valid(cache_path, path1, path2):
                        # Load from cache
                        try:
                            with open(cache_path, 'rb') as f:
                                cache_data = pickle.load(f)
                            comparison = cache_data['comparison']
                            self.stats['from_cache'] += 1
                            source = 'cache'
                        except Exception as e:
                            # Cache read failed, recompute
                            comparison = self.compare_pair(path1, path2, network, method1, method2)
                            self.stats['computed'] += 1
                            source = 'computed'
                    else:
                        # Compute comparison
                        comparison = self.compare_pair(path1, path2, network, method1, method2)
                        self.stats['computed'] += 1
                        source = 'computed'

                        # Save to cache if successful
                        if comparison['status'] == 'SUCCESS':
                            try:
                                cache_data = {
                                    'network': network,
                                    'method1': method1,
                                    'method2': method2,
                                    'timestamp': datetime.now().isoformat(),
                                    'hash1': self._file_hash(path1),
                                    'hash2': self._file_hash(path2),
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
                                'method1': method1,
                                'method2': method2,
                                'metric': metric_name,
                                'value': value,
                                'status': 'SUCCESS'
                            })
                    else:
                        self.stats['failed'] += 1
                        self.stats['errors'].append({
                            'network': network,
                            'method1': method1,
                            'method2': method2,
                            'error': comparison['error']
                        })

                        # Add failed entry with NaN value
                        results.append({
                            'network': network,
                            'method1': method1,
                            'method2': method2,
                            'metric': 'FAILED',
                            'value': np.nan,
                            'status': 'FAILED'
                        })

        # Progress update
        print(f"\nProcessed {self.stats['total']} pairwise comparisons")
        print(f"  Success: {self.stats['success']}, Failed: {self.stats['failed']}, "
              f"Cached: {self.stats['from_cache']}")

        return pd.DataFrame(results)

    def print_statistics(self):
        """Print comparison statistics"""
        print(f"\n{'='*80}")
        print(f"Comparison Engine Statistics")
        print(f"{'='*80}")
        print(f"Total comparisons:   {self.stats['total']}")
        if self.stats['total'] > 0:
            print(f"Successful:          {self.stats['success']} ({self.stats['success']/self.stats['total']*100:.1f}%)")
            print(f"Failed:              {self.stats['failed']} ({self.stats['failed']/self.stats['total']*100:.1f}%)")
        print(f"From cache:         {self.stats['from_cache']}")
        print(f"Newly computed:     {self.stats['computed']}")

        if self.stats['errors']:
            print(f"\nErrors ({len(self.stats['errors'])}):")
            for i, error in enumerate(self.stats['errors'][:10], 1):  # Show first 10
                print(f"  {i}. {error['network']} / {error['method1']} vs {error['method2']}")
                print(f"     Error: {error['error']}")

            if len(self.stats['errors']) > 10:
                print(f"  ... and {len(self.stats['errors']) - 10} more")

        print(f"{'='*80}\n")

    def write_comparison_report(self, output_file: str):
        """Write detailed comparison report to file"""
        with open(output_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("Pairwise Comparison Report\n")
            f.write("="*80 + "\n\n")

            f.write(f"Total comparisons: {self.stats['total']}\n")
            if self.stats['total'] > 0:
                f.write(f"Successful:        {self.stats['success']} ({self.stats['success']/self.stats['total']*100:.1f}%)\n")
                f.write(f"Failed:            {self.stats['failed']} ({self.stats['failed']/self.stats['total']*100:.1f}%)\n")
            f.write(f"From cache:        {self.stats['from_cache']}\n")
            f.write(f"Newly computed:    {self.stats['computed']}\n\n")

            if self.stats['errors']:
                f.write(f"\nFailed Comparisons ({len(self.stats['errors'])}):\n")
                f.write("-"*80 + "\n")
                for error in self.stats['errors']:
                    f.write(f"\n{error['network']} / {error['method1']} vs {error['method2']}\n")
                    f.write(f"Error: {error['error']}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Compute pairwise network comparison metrics with caching',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compute comparisons from inventory
  %(prog)s inventory.csv cache/

  # Export to CSV
  %(prog)s inventory.csv cache/ --export comparisons.csv

  # Force recompute (ignore cache)
  %(prog)s inventory.csv cache/ --force-recompute
        """
    )

    parser.add_argument('inventory', help='Path to inventory CSV file')
    parser.add_argument('cache_dir', help='Directory for caching comparison results')
    parser.add_argument('--force-recompute', action='store_true',
                       help='Force recompute all comparisons (ignore cache)')
    parser.add_argument('--export', metavar='FILE',
                       help='Export comparisons to CSV file')
    parser.add_argument('--report', metavar='FILE',
                       help='Write detailed comparison report to file')

    args = parser.parse_args()

    # Load inventory
    try:
        inventory = pd.read_csv(args.inventory)
    except FileNotFoundError:
        print(f"Error: Inventory file not found: {args.inventory}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading inventory: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate inventory columns
    required_cols = ['network', 'method', 'network_path', 'exists']
    missing_cols = [col for col in required_cols if col not in inventory.columns]
    if missing_cols:
        print(f"Error: Inventory missing required columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)

    # Create comparison engine
    engine = ComparisonEngine(args.cache_dir, force_recompute=args.force_recompute)

    # Compute comparisons
    print(f"Computing pairwise comparisons from: {args.inventory}")
    comparisons = engine.compute_all_comparisons(inventory)

    # Print statistics
    engine.print_statistics()

    # Export if requested
    if args.export:
        comparisons.to_csv(args.export, index=False)
        print(f"Comparisons exported to: {args.export}")
        print(f"Total rows: {len(comparisons)}")

    # Write report if requested
    if args.report:
        engine.write_comparison_report(args.report)
        print(f"Comparison report written to: {args.report}")

    return comparisons


if __name__ == '__main__':
    comparisons = main()

