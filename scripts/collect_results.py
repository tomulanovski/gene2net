#!/usr/bin/env python3
"""
collect_results.py - Data Collection Module for Real Data

Scans papers directory structure and creates inventory of all available method outputs
for phylogenetic network inference comparison.

Usage:
    python collect_results.py [--config YAML] [--export CSV]
    python collect_results.py --export inventory.csv
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List
import pandas as pd
import yaml


class ResultInventory:
    """Scans papers directory structure and catalogs all available method outputs"""

    def __init__(self, config_dict: Dict):
        """
        Initialize inventory scanner

        Args:
            config_dict: Configuration dictionary from YAML
        """
        self.papers_dir = Path(config_dict['papers_dir'])
        self.networks = config_dict['networks']
        self.methods = config_dict['methods']

    def scan_method_results(self, network: str, method: str) -> Dict:
        """
        Scan for result file for a specific method and network

        Args:
            network: Network name (e.g., "Bendiksby_2011")
            method: Method name (e.g., "grampa", "polyphest", "paper")

        Returns:
            Dictionary with network, method, path, exists, file_size
        """
        method_config = self.methods[method]
        method_dir = method_config['directory']
        output_file = method_config['output_file']

        # Construct path: papers/{network}/networks/{method_dir}/{output_file}
        result_path = self.papers_dir / network / "networks" / method_dir / output_file

        file_size = result_path.stat().st_size if result_path.exists() else 0

        return {
            'network': network,
            'method': method,
            'network_path': str(result_path),
            'exists': result_path.exists() and file_size > 0,
            'file_size': file_size
        }

    def create_full_inventory(self) -> pd.DataFrame:
        """
        Create master inventory of all network-method combinations

        Returns:
            DataFrame with columns: network, method, network_path, exists, file_size
        """
        results = []

        for network in self.networks:
            for method in self.methods.keys():
                result = self.scan_method_results(network, method)
                results.append(result)

        return pd.DataFrame(results)

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
        available = inventory['exists'].sum()
        missing = (~inventory['exists']).sum()

        # Per-method stats
        method_stats = {}
        for method in inventory['method'].unique():
            method_data = inventory[inventory['method'] == method]
            expected = len(method_data)
            method_available = method_data['exists'].sum()
            method_missing = expected - method_available

            method_stats[method] = {
                'expected': expected,
                'available': method_available,
                'missing': method_missing,
                'completion_rate': (method_available / expected * 100) if expected > 0 else 0
            }

        # Per-network stats
        network_stats = {}
        for network in inventory['network'].unique():
            network_data = inventory[inventory['network'] == network]
            expected = len(network_data)
            network_available = network_data['exists'].sum()
            network_missing = expected - network_available

            network_stats[network] = {
                'expected': expected,
                'available': network_available,
                'missing': network_missing,
                'completion_rate': (network_available / expected * 100) if expected > 0 else 0
            }

        return {
            'total_combinations': total_combinations,
            'available': available,
            'missing': missing,
            'completion_rate': (available / total_combinations * 100) if total_combinations > 0 else 0,
            'method_stats': method_stats,
            'network_stats': network_stats
        }

    def print_completion_report(self, report: Dict):
        """Print formatted completion report"""
        print(f"\n{'='*80}")
        print(f"Data Collection Report: Papers Analysis")
        print(f"{'='*80}")

        print(f"\nOverall Summary:")
        print(f"  Total combinations: {report['total_combinations']}")
        print(f"  Available:         {report['available']} ({report['completion_rate']:.1f}%)")
        print(f"  Missing:            {report['missing']}")

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
        description='Collect and inventory phylogenetic network inference results from papers directory',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Scan results using default config
  %(prog)s

  # Export inventory to CSV
  %(prog)s --export inventory.csv

  # Use custom config file
  %(prog)s --config custom_config.yaml
        """
    )

    parser.add_argument('--config', default='scripts/papers_config.yaml',
                       help='Path to configuration YAML file (default: scripts/papers_config.yaml)')
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
    print(f"Scanning papers directory: {config_dict['papers_dir']}")
    inventory_scanner = ResultInventory(config_dict)

    # Scan and create inventory
    inventory = inventory_scanner.create_full_inventory()

    # Generate and print report
    report = inventory_scanner.get_completion_report(inventory)
    inventory_scanner.print_completion_report(report)

    # Show sample of inventory if verbose
    if args.verbose:
        print("Sample inventory (first 20 rows):")
        print(inventory.head(20).to_string())
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

