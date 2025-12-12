#!/usr/bin/env python3
"""
collect_results.py - Data Collection Module

Scans directory structure and creates inventory of all available simulation results
for phylogenetic network inference evaluation.

Usage:
    python collect_results.py CONFIG [--config YAML] [--export CSV]
    python collect_results.py conf_ils_low_10M --export inventory.csv
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import yaml


class ResultInventory:
    """Scans directory structure and catalogs all available results"""

    def __init__(self, config: str, config_dict: Dict):
        """
        Initialize inventory scanner

        Args:
            config: Configuration name (e.g., 'conf_ils_low_10M')
            config_dict: Configuration dictionary from YAML
        """
        self.config = config
        self.base_dir = Path(config_dict['base_dir'])
        self.networks_dir = Path(config_dict['networks_dir'])
        self.networks = config_dict['networks']
        self.methods = config_dict['methods']
        self.num_replicates = config_dict['num_replicates']

    def scan_ground_truth_networks(self) -> pd.DataFrame:
        """
        Scan for ground truth network files

        Returns:
            DataFrame with columns: network, gt_path, gt_exists
        """
        results = []

        for network in self.networks:
            gt_path = self.networks_dir / f"{network}.tre"

            results.append({
                'network': network,
                'gt_path': str(gt_path),
                'gt_exists': gt_path.exists()
            })

        return pd.DataFrame(results)

    def scan_method_results(self, method: str) -> pd.DataFrame:
        """
        Scan for result files for a specific method

        Args:
            method: Method name (e.g., 'grampa', 'polyphest_p50')

        Returns:
            DataFrame with columns: network, replicate, inferred_path, inferred_exists, file_size
        """
        results = []
        method_config = self.methods[method]
        method_dir = method_config['directory']
        output_file = method_config['output_file']

        for network in self.networks:
            for replicate in range(1, self.num_replicates + 1):
                # Construct path: {base}/{network}/results/{config}/{method_dir}/replicate_{N}/{output_file}
                result_path = (self.base_dir / network / "results" / self.config /
                              method_dir / f"replicate_{replicate}" / output_file)

                file_size = result_path.stat().st_size if result_path.exists() else 0

                results.append({
                    'network': network,
                    'method': method,
                    'replicate': replicate,
                    'inferred_path': str(result_path),
                    'inferred_exists': result_path.exists() and file_size > 0,
                    'file_size': file_size
                })

        return pd.DataFrame(results)

    def create_full_inventory(self) -> pd.DataFrame:
        """
        Create master inventory of all combinations

        Returns:
            DataFrame with columns: network, config, method, replicate,
                                   ground_truth_path, inferred_path,
                                   gt_exists, inferred_exists, file_size
        """
        # Get ground truth info
        gt_df = self.scan_ground_truth_networks()

        # Collect all method results
        method_dfs = []
        for method in self.methods.keys():
            method_df = self.scan_method_results(method)
            method_dfs.append(method_df)

        # Combine all method results
        all_methods_df = pd.concat(method_dfs, ignore_index=True)

        # Merge with ground truth
        inventory = all_methods_df.merge(gt_df, on='network', how='left')

        # Add config column
        inventory.insert(1, 'config', self.config)

        # Reorder columns
        inventory = inventory[['network', 'config', 'method', 'replicate',
                              'gt_path', 'inferred_path',
                              'gt_exists', 'inferred_exists', 'file_size']]

        return inventory

    def get_completion_report(self, inventory: pd.DataFrame) -> Dict:
        """
        Generate summary of data availability

        Args:
            inventory: Master inventory DataFrame

        Returns:
            Dictionary with completion statistics
        """
        total_combinations = len(inventory)

        # Overall stats
        both_exist = (inventory['gt_exists'] & inventory['inferred_exists']).sum()
        gt_missing = (~inventory['gt_exists']).sum()
        inferred_missing = (~inventory['inferred_exists']).sum()

        # Per-method stats
        method_stats = {}
        for method in inventory['method'].unique():
            method_data = inventory[inventory['method'] == method]
            expected = len(method_data)
            available = (method_data['gt_exists'] & method_data['inferred_exists']).sum()
            missing = expected - available

            method_stats[method] = {
                'expected': expected,
                'available': available,
                'missing': missing,
                'completion_rate': (available / expected * 100) if expected > 0 else 0
            }

        # Per-network stats
        network_stats = {}
        for network in inventory['network'].unique():
            network_data = inventory[inventory['network'] == network]
            expected = len(network_data)
            available = (network_data['gt_exists'] & network_data['inferred_exists']).sum()
            missing = expected - available

            network_stats[network] = {
                'expected': expected,
                'available': available,
                'missing': missing,
                'completion_rate': (available / expected * 100) if expected > 0 else 0
            }

        return {
            'total_combinations': total_combinations,
            'both_exist': both_exist,
            'gt_missing': gt_missing,
            'inferred_missing': inferred_missing,
            'completion_rate': (both_exist / total_combinations * 100) if total_combinations > 0 else 0,
            'method_stats': method_stats,
            'network_stats': network_stats
        }

    def print_completion_report(self, report: Dict):
        """Print formatted completion report"""
        print(f"\n{'='*80}")
        print(f"Data Collection Report: {self.config}")
        print(f"{'='*80}")

        print(f"\nOverall Summary:")
        print(f"  Total combinations: {report['total_combinations']}")
        print(f"  Both files exist:   {report['both_exist']} ({report['completion_rate']:.1f}%)")
        print(f"  GT missing:         {report['gt_missing']}")
        print(f"  Inferred missing:   {report['inferred_missing']}")

        print(f"\nPer-Method Completion:")
        print(f"  {'Method':<20} {'Expected':<10} {'Available':<10} {'Missing':<10} {'Rate':<10}")
        print(f"  {'-'*70}")
        for method, stats in sorted(report['method_stats'].items()):
            print(f"  {method:<20} {stats['expected']:<10} {stats['available']:<10} "
                  f"{stats['missing']:<10} {stats['completion_rate']:>6.1f}%")

        # Show networks with issues
        incomplete_networks = {net: stats for net, stats in report['network_stats'].items()
                              if stats['completion_rate'] < 100}

        if incomplete_networks:
            print(f"\nNetworks with Missing Data ({len(incomplete_networks)}):")
            for network, stats in sorted(incomplete_networks.items(),
                                        key=lambda x: x[1]['completion_rate']):
                print(f"  {network:<30} {stats['available']}/{stats['expected']} "
                      f"({stats['completion_rate']:.1f}%)")

        print(f"\n{'='*80}\n")


def load_config(config_path: str) -> Dict:
    """Load configuration from YAML file"""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(
        description='Collect and inventory phylogenetic network inference results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Scan results for a configuration
  %(prog)s conf_ils_low_10M

  # Export inventory to CSV
  %(prog)s conf_ils_low_10M --export inventory.csv

  # Use custom config file
  %(prog)s conf_ils_low_10M --config custom_config.yaml
        """
    )

    parser.add_argument('configuration', help='Configuration name (e.g., conf_ils_low_10M)')
    parser.add_argument('--config', default='simulations/summary_config.yaml',
                       help='Path to configuration YAML file (default: simulations/summary_config.yaml)')
    parser.add_argument('--export', metavar='FILE',
                       help='Export inventory to CSV file')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Show detailed file paths')

    args = parser.parse_args()

    # Load configuration
    try:
        config_dict = load_config(args.config)
    except FileNotFoundError:
        print(f"Error: Configuration file not found: {args.config}", file=sys.stderr)
        print(f"Current directory: {os.getcwd()}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading configuration: {e}", file=sys.stderr)
        sys.exit(1)

    # Create inventory
    print(f"Scanning results for configuration: {args.configuration}")
    inventory_scanner = ResultInventory(args.configuration, config_dict)

    # Scan and create inventory
    inventory = inventory_scanner.create_full_inventory()

    # Generate and print report
    report = inventory_scanner.get_completion_report(inventory)
    inventory_scanner.print_completion_report(report)

    # Show sample of inventory if verbose
    if args.verbose:
        print("Sample inventory (first 10 rows):")
        print(inventory.head(10).to_string())
        print()

    # Export if requested
    if args.export:
        inventory.to_csv(args.export, index=False)
        print(f"Inventory exported to: {args.export}")
        print(f"Total rows: {len(inventory)}")

    # Return inventory for programmatic use
    return inventory, report


if __name__ == '__main__':
    inventory, report = main()
