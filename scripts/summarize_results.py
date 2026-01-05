#!/usr/bin/env python3
"""
summarize_results.py - Summary Tables Generator for Real Data

Generates summary tables from pairwise comparisons:
- Method availability
- Pairwise comparison summaries (aggregated across networks)
- Per-network comparisons
- Method rankings

Usage:
    python summarize_results.py COMPARISONS_CSV INVENTORY_CSV OUTPUT_DIR
    python summarize_results.py comparisons.csv inventory.csv analysis/
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Dict
import pandas as pd
import numpy as np


class ResultSummarizer:
    """Generate summary tables from pairwise comparisons"""

    def __init__(self, comparisons_df: pd.DataFrame, inventory_df: pd.DataFrame,
                 comparable_networks: list = None):
        """
        Initialize summarizer

        Args:
            comparisons_df: DataFrame from compute_comparisons.py
            inventory_df: DataFrame from collect_results.py
            comparable_networks: List of network names to use for completion rate calculation.
                               If None, uses all networks. If provided, completion rates are
                               calculated as percentage of these networks that have results.
        """
        self.comparisons_df = comparisons_df
        self.inventory_df = inventory_df
        self.comparable_networks = comparable_networks

        # Filter out failed comparisons
        if not comparisons_df.empty and 'status' in comparisons_df.columns:
            self.valid_df = comparisons_df[comparisons_df['status'] == 'SUCCESS'].copy()
        else:
            self.valid_df = pd.DataFrame()

    def generate_method_availability(self) -> pd.DataFrame:
        """
        Generate method availability table

        Returns:
            DataFrame with columns: method, total_networks, available, missing, completion_rate,
                                   comparable_networks, available_comparable
        """
        availability = []
        
        # Determine which networks to use for completion rate calculation
        if self.comparable_networks:
            # Filter inventory to only comparable networks
            comparable_inventory = self.inventory_df[
                self.inventory_df['network'].isin(self.comparable_networks)
            ]
            total_for_completion = len(self.comparable_networks)
        else:
            comparable_inventory = self.inventory_df
            total_for_completion = len(self.inventory_df['network'].unique())

        for method in self.inventory_df['method'].unique():
            method_data = self.inventory_df[self.inventory_df['method'] == method]
            total = len(method_data)
            available = method_data['exists'].sum()
            missing = total - available
            
            # Calculate completion rate based on comparable networks only
            if self.comparable_networks:
                method_comparable = comparable_inventory[
                    comparable_inventory['method'] == method
                ]
                available_comparable = method_comparable['exists'].sum()
                completion_rate = (available_comparable / total_for_completion * 100) if total_for_completion > 0 else 0
            else:
                available_comparable = available
                completion_rate = (available / total * 100) if total > 0 else 0

            availability.append({
                'method': method,
                'total_networks': total,
                'available': available,
                'missing': missing,
                'completion_rate': completion_rate,
                'comparable_networks': total_for_completion,
                'available_comparable': available_comparable
            })

        return pd.DataFrame(availability).sort_values('completion_rate', ascending=False)

    def generate_pairwise_summary(self) -> pd.DataFrame:
        """
        Generate pairwise comparison summary aggregated across all networks

        Returns:
            DataFrame with columns: method1, method2, metric, mean, std, min, max, n_networks
        """
        if self.valid_df.empty:
            return pd.DataFrame()

        # Group by (method1, method2, metric) and aggregate across networks
        grouped = self.valid_df.groupby(['method1', 'method2', 'metric'])['value']

        summary = grouped.agg([
            ('mean', 'mean'),
            ('std', 'std'),
            ('min', 'min'),
            ('max', 'max'),
            ('n_networks', 'count')
        ]).reset_index()

        # Fill NaN std (when only 1 network) with 0
        summary['std'] = summary['std'].fillna(0)

        # Sort by method pair
        summary = summary.sort_values(['method1', 'method2', 'metric'])

        return summary

    def generate_per_network_comparisons(self) -> pd.DataFrame:
        """
        Generate detailed per-network comparison table

        Returns:
            Wide-format DataFrame with one row per network and columns for each method pair × metric
        """
        if self.valid_df.empty:
            return pd.DataFrame()

        # Create method pair identifier
        self.valid_df['method_pair'] = self.valid_df.apply(
            lambda row: f"{row['method1']}_vs_{row['method2']}", axis=1
        )

        # Pivot: network × (method_pair_metric)
        pivot_data = self.valid_df.pivot_table(
            index='network',
            columns=['method_pair', 'metric'],
            values='value',
            aggfunc='first'  # Should only be one value per combination
        )

        # Flatten column names
        pivot_data.columns = [f"{pair}_{metric}" for pair, metric in pivot_data.columns]

        return pivot_data.reset_index()

    def generate_method_rankings(self) -> pd.DataFrame:
        """
        Generate method rankings based on pairwise comparisons

        For each method, compute average distance to all other methods.
        Lower distance = better (methods are more similar to others).

        Returns:
            DataFrame with columns: method, avg_edit_distance_multree, avg_rf_distance, 
                                   avg_num_rets_diff, rank_edit_distance, rank_rf_distance, etc.
        """
        if self.valid_df.empty:
            return pd.DataFrame()

        # Focus on key metrics for ranking
        key_metrics = ['edit_distance_multree', 'rf_distance', 'num_rets_diff']

        rankings = []

        # Get all unique methods
        all_methods = set(self.valid_df['method1'].unique()) | set(self.valid_df['method2'].unique())

        for method in all_methods:
            # Get all comparisons involving this method
            method_comparisons = self.valid_df[
                (self.valid_df['method1'] == method) | (self.valid_df['method2'] == method)
            ]

            method_stats = {'method': method}

            # Compute average for each key metric
            for metric in key_metrics:
                metric_data = method_comparisons[method_comparisons['metric'] == metric]['value']
                if len(metric_data) > 0:
                    method_stats[f'avg_{metric}'] = metric_data.mean()
                    method_stats[f'std_{metric}'] = metric_data.std()
                    method_stats[f'n_{metric}'] = len(metric_data)
                else:
                    method_stats[f'avg_{metric}'] = np.nan
                    method_stats[f'std_{metric}'] = np.nan
                    method_stats[f'n_{metric}'] = 0

            rankings.append(method_stats)

        rankings_df = pd.DataFrame(rankings)

        # Add rankings (lower is better for distances)
        for metric in key_metrics:
            col = f'avg_{metric}'
            if col in rankings_df.columns:
                rankings_df[f'rank_{metric}'] = rankings_df[col].rank(ascending=True, method='min')

        return rankings_df.sort_values('rank_edit_distance_multree' if 'rank_edit_distance_multree' in rankings_df.columns else 'method')

    def generate_all_summaries(self, output_dir: Path):
        """
        Generate all summary tables and save to output directory

        Args:
            output_dir: Directory to save summary files
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*80}")
        print(f"Generating Summary Tables")
        print(f"{'='*80}\n")

        # Method availability
        print("Generating method availability table...")
        if self.comparable_networks:
            print(f"  Using {len(self.comparable_networks)} comparable networks for completion rate calculation")
        availability = self.generate_method_availability()
        availability_file = output_dir / "method_availability.csv"
        availability.to_csv(availability_file, index=False)
        print(f"  Saved: {availability_file}")
        print(f"  Methods: {len(availability)}")

        # Pairwise summary
        print("\nGenerating pairwise comparison summary...")
        pairwise_summary = self.generate_pairwise_summary()
        if not pairwise_summary.empty:
            pairwise_file = output_dir / "pairwise_summary.csv"
            pairwise_summary.to_csv(pairwise_file, index=False)
            print(f"  Saved: {pairwise_file}")
            print(f"  Method pairs: {pairwise_summary[['method1', 'method2']].drop_duplicates().shape[0]}")
            print(f"  Metrics: {pairwise_summary['metric'].nunique()}")
        else:
            print("  WARNING: No valid comparisons to summarize")

        # Per-network comparisons
        print("\nGenerating per-network comparison table...")
        per_network = self.generate_per_network_comparisons()
        if not per_network.empty:
            per_network_file = output_dir / "per_network_comparisons.csv"
            per_network.to_csv(per_network_file, index=False)
            print(f"  Saved: {per_network_file}")
            print(f"  Networks: {len(per_network)}")
            print(f"  Comparison columns: {len(per_network.columns) - 1}")
        else:
            print("  WARNING: No valid comparisons for per-network table")

        # Method rankings
        print("\nGenerating method rankings...")
        rankings = self.generate_method_rankings()
        if not rankings.empty:
            rankings_file = output_dir / "method_rankings.csv"
            rankings.to_csv(rankings_file, index=False)
            print(f"  Saved: {rankings_file}")
            print(f"  Methods ranked: {len(rankings)}")
        else:
            print("  WARNING: No valid comparisons for rankings")

        print(f"\n{'='*80}")
        print(f"Summary generation complete!")
        print(f"Output directory: {output_dir}")
        print(f"{'='*80}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate summary tables from pairwise comparisons',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate all summaries
  %(prog)s comparisons.csv inventory.csv analysis/

  # Custom output directory
  %(prog)s comparisons.csv inventory.csv analysis/latest/
        """
    )

    parser.add_argument('comparisons', help='Path to comparisons CSV file')
    parser.add_argument('inventory', help='Path to inventory CSV file')
    parser.add_argument('output_dir', help='Output directory for summary tables')

    args = parser.parse_args()

    # Load data
    try:
        comparisons_df = pd.read_csv(args.comparisons)
        print(f"Loaded comparisons: {len(comparisons_df)} rows")
    except FileNotFoundError:
        print(f"Error: Comparisons file not found: {args.comparisons}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading comparisons: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        inventory_df = pd.read_csv(args.inventory)
        print(f"Loaded inventory: {len(inventory_df)} rows")
    except FileNotFoundError:
        print(f"Error: Inventory file not found: {args.inventory}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading inventory: {e}", file=sys.stderr)
        sys.exit(1)

    # Create summarizer
    summarizer = ResultSummarizer(comparisons_df, inventory_df)

    # Generate all summaries
    summarizer.generate_all_summaries(args.output_dir)


if __name__ == '__main__':
    main()

