#!/usr/bin/env python3
"""
aggregate_and_summarize.py - Aggregation and Multi-Level Summary Module

Aggregates replicate results and generates multiple analytical views.

Usage:
    python aggregate_and_summarize.py COMPARISONS_CSV OUTPUT_DIR [--network-stats CSV]
    python aggregate_and_summarize.py comparisons.csv summary/ --network-stats networks/mul_tree_final_stats.csv
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Dict, Optional
import pandas as pd
import numpy as np
from scipy.stats import pearsonr


class MetricAggregator:
    """Aggregate metrics across replicates"""

    def __init__(self, comparisons_df: pd.DataFrame):
        """
        Initialize aggregator

        Args:
            comparisons_df: DataFrame from compute_comparisons.py
        """
        self.comparisons_df = comparisons_df

        # Filter out failed comparisons
        if not comparisons_df.empty and 'status' in comparisons_df.columns:
            self.valid_df = comparisons_df[comparisons_df['status'] == 'SUCCESS'].copy()
        else:
            # No comparisons or missing status column - create empty DataFrame with expected columns
            self.valid_df = pd.DataFrame(columns=['network', 'config', 'method', 'metric', 'value'])

    def aggregate_replicates(self) -> pd.DataFrame:
        """
        Aggregate across 5 replicates

        Returns:
            DataFrame with columns: network, config, method, metric, mean, std, min, max, n_valid
        """
        # Handle empty dataframe
        if self.valid_df.empty:
            print("\nWARNING: No valid comparisons to aggregate (all methods missing or failed)")
            return pd.DataFrame(columns=['network', 'config', 'method', 'metric', 'mean', 'std', 'min', 'max', 'n_valid'])

        # Group by (network, config, method, metric)
        grouped = self.valid_df.groupby(['network', 'config', 'method', 'metric'])['value']

        # Compute statistics
        aggregated = grouped.agg([
            ('mean', 'mean'),
            ('std', 'std'),
            ('min', 'min'),
            ('max', 'max'),
            ('n_valid', 'count')
        ]).reset_index()

        # Fill NaN std (when only 1 replicate) with 0
        aggregated['std'] = aggregated['std'].fillna(0)

        return aggregated


class MultiLevelSummary:
    """Generate different summary views"""

    def __init__(self, aggregated_df: pd.DataFrame, config: str):
        """
        Initialize summary generator

        Args:
            aggregated_df: Aggregated metrics DataFrame
            config: Configuration name
        """
        self.aggregated_df = aggregated_df
        self.config = config

    def level1_detailed_per_network(self, output_dir: Path):
        """
        Generate detailed per-network pivot tables for each metric

        Creates files:
            - edit_distance.csv
            - num_rets_diff.csv
            - ploidy_diff.csv (for ploidy_diff.dist)
            - ret_leaf_jaccard.csv (for ret_leaf_jaccard.dist)
            - ret_sisters_jaccard.csv (for ret_sisters_jaccard.dist)
        """
        level1_dir = output_dir / "level1_detailed_per_network"
        level1_dir.mkdir(parents=True, exist_ok=True)

        # Define metrics to generate pivot tables for
        # PRIMARY METRICS (MUL-tree based)
        metrics_to_pivot = {
            'edit_distance_multree': 'edit_distance_multree',  # PRIMARY: MUL-tree edit distance
            'rf_distance': 'rf_distance',  # PRIMARY: RF distance on MUL-trees
            'num_rets_diff': 'num_rets_diff',
            'num_rets_bias': 'num_rets_bias',  # Signed difference (bias)
            'ploidy_diff.dist': 'ploidy_diff',
            'ret_leaf_jaccard.dist': 'ret_leaf_jaccard',
            'ret_sisters_jaccard.dist': 'ret_sisters_jaccard',
            'edit_distance': 'edit_distance_network'  # LEGACY: Network edit distance (for comparison)
        }

        for metric_name, output_name in metrics_to_pivot.items():
            # Filter to this metric
            metric_data = self.aggregated_df[self.aggregated_df['metric'] == metric_name].copy()

            if len(metric_data) == 0:
                print(f"Warning: No data for metric '{metric_name}'")
                continue

            # Create pivot table: rows=networks, columns=methods, values=mean±std
            pivot_mean = metric_data.pivot(index='network', columns='method', values='mean')
            pivot_std = metric_data.pivot(index='network', columns='method', values='std')
            pivot_n = metric_data.pivot(index='network', columns='method', values='n_valid')

            # Combine mean ± std
            pivot_combined = pivot_mean.copy()
            for col in pivot_combined.columns:
                # Format as "mean ± std (n)"
                pivot_combined[col] = pivot_mean[col].apply(lambda x: f"{x:.4f}" if not np.isnan(x) else "N/A") + \
                                     " ± " + \
                                     pivot_std[col].apply(lambda x: f"{x:.4f}" if not np.isnan(x) else "0") + \
                                     " (" + \
                                     pivot_n[col].apply(lambda x: f"{int(x)}" if not np.isnan(x) else "0") + \
                                     ")"

            # Save to CSV
            output_file = level1_dir / f"{output_name}.csv"
            pivot_combined.to_csv(output_file)

            # Also save separate mean/std tables
            pivot_mean.to_csv(level1_dir / f"{output_name}_mean.csv")
            pivot_std.to_csv(level1_dir / f"{output_name}_std.csv")
            pivot_n.to_csv(level1_dir / f"{output_name}_n.csv")

            print(f"  Created Level 1: {output_name}.csv ({len(pivot_combined)} networks × {len(pivot_combined.columns)} methods)")

    def level2_method_rankings(self, output_file: Path):
        """
        Generate overall method rankings across all networks

        Output columns: method, metric, avg_value, avg_rank, num_best, num_worst
        """
        # Metrics to rank (lower is better for distances)
        # PRIMARY METRICS first, then secondary
        # Note: num_rets_bias is NOT ranked (it's signed, not an error magnitude)
        distance_metrics = ['edit_distance_multree', 'rf_distance',  # PRIMARY MUL-tree metrics
                           'num_rets_diff', 'ploidy_diff.dist',
                           'ret_leaf_jaccard.dist', 'ret_sisters_jaccard.dist',
                           'edit_distance']  # LEGACY network edit distance

        rankings = []

        for metric in distance_metrics:
            # Filter to this metric
            metric_data = self.aggregated_df[self.aggregated_df['metric'] == metric].copy()

            if len(metric_data) == 0:
                continue

            # For each network, rank methods (1 = best = lowest distance)
            network_ranks = []

            for network in metric_data['network'].unique():
                network_data = metric_data[metric_data['network'] == network].copy()
                network_data = network_data.sort_values('mean')

                # Assign ranks
                network_data['rank'] = range(1, len(network_data) + 1)

                # Identify best method
                best_method = network_data.iloc[0]['method']
                worst_method = network_data.iloc[-1]['method']

                # Record ranks
                for _, row in network_data.iterrows():
                    network_ranks.append({
                        'network': network,
                        'method': row['method'],
                        'metric': metric,
                        'value': row['mean'],
                        'rank': row['rank'],
                        'is_best': row['method'] == best_method,
                        'is_worst': row['method'] == worst_method
                    })

            rank_df = pd.DataFrame(network_ranks)

            # Aggregate across networks
            for method in rank_df['method'].unique():
                method_data = rank_df[rank_df['method'] == method]

                rankings.append({
                    'method': method,
                    'metric': metric,
                    'avg_value': method_data['value'].mean(),
                    'avg_rank': method_data['rank'].mean(),
                    'num_best': method_data['is_best'].sum(),
                    'num_worst': method_data['is_worst'].sum(),
                    'num_networks': len(method_data)
                })

        rankings_df = pd.DataFrame(rankings)

        # Sort by metric and avg_rank
        rankings_df = rankings_df.sort_values(['metric', 'avg_rank'])

        # Save
        rankings_df.to_csv(output_file, index=False)

        print(f"  Created Level 2: method_rankings.csv ({len(rankings_df)} rows)")

        return rankings_df

    def level3_network_correlations(self, output_file: Path, network_stats_csv: Optional[str] = None):
        """
        Correlate network properties with reconstruction difficulty

        Args:
            output_file: Output CSV file path
            network_stats_csv: Path to mul_tree_final_stats.csv

        Output columns: network_property, method, correlation, p_value
        """
        if network_stats_csv is None or not Path(network_stats_csv).exists():
            print(f"  Skipping Level 3: Network stats file not provided or not found")
            return None

        # Load network characteristics
        try:
            network_stats = pd.read_csv(network_stats_csv)

            # Clean up network names (remove .tre extension if present)
            network_stats['network'] = network_stats['Filename'].str.replace('.tre', '', regex=False)

            # Select relevant columns
            properties = {
                'Num_Species': 'num_species',
                'Num_Polyploids': 'num_polyploids',
                'Max_Copies': 'max_copies',
                'H_Strict': 'num_reticulations'
            }

            network_stats = network_stats[['network'] + list(properties.keys())]

        except Exception as e:
            print(f"  Error loading network stats: {e}")
            return None

        # Merge with aggregated metrics
        merged = self.aggregated_df.merge(network_stats, on='network', how='left')

        # Compute correlations for each (method, metric, property) combination
        correlations = []

        # Focus on main distance metrics (bias is signed, not a distance metric for correlation)
        # PRIMARY METRICS first
        main_metrics = ['edit_distance_multree', 'rf_distance',  # PRIMARY MUL-tree metrics
                       'num_rets_diff', 'num_rets_bias']

        for method in merged['method'].unique():
            for metric in main_metrics:
                method_metric_data = merged[(merged['method'] == method) &
                                           (merged['metric'] == metric)].copy()

                if len(method_metric_data) < 3:  # Need at least 3 points for correlation
                    continue

                for prop_col, prop_name in properties.items():
                    # Remove NaN values
                    clean_data = method_metric_data[['mean', prop_col]].dropna()

                    if len(clean_data) < 3:
                        continue

                    # Compute Pearson correlation
                    try:
                        corr, p_value = pearsonr(clean_data[prop_col], clean_data['mean'])

                        correlations.append({
                            'method': method,
                            'metric': metric,
                            'network_property': prop_name,
                            'correlation': corr,
                            'p_value': p_value,
                            'significant': p_value < 0.05,
                            'n_networks': len(clean_data)
                        })
                    except Exception:
                        continue

        if len(correlations) == 0:
            print(f"  Skipping Level 3: No valid correlations computed")
            return None

        corr_df = pd.DataFrame(correlations)

        # Sort by significance and absolute correlation
        corr_df = corr_df.sort_values(['significant', 'correlation'], ascending=[False, False])

        # Save
        corr_df.to_csv(output_file, index=False)

        print(f"  Created Level 3: network_correlations.csv ({len(corr_df)} rows)")

        return corr_df

    def generate_all_summaries(self, output_dir: Path, network_stats_csv: Optional[str] = None):
        """Generate all summary levels"""
        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\nGenerating multi-level summaries...")

        # Level 1: Detailed per-network
        self.level1_detailed_per_network(output_dir)

        # Level 2: Method rankings
        self.level2_method_rankings(output_dir / "level2_method_rankings.csv")

        # Level 3: Network correlations
        self.level3_network_correlations(output_dir / "level3_network_correlations.csv",
                                        network_stats_csv)

        print(f"\nAll summaries generated in: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate replicate results and generate multi-level summaries',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate all summaries
  %(prog)s comparisons.csv summary/

  # Include network characteristics correlation
  %(prog)s comparisons.csv summary/ --network-stats networks/mul_tree_final_stats.csv
        """
    )

    parser.add_argument('comparisons', help='Comparisons CSV file from compute_comparisons.py')
    parser.add_argument('output_dir', help='Output directory for summary files')
    parser.add_argument('--network-stats', metavar='CSV',
                       help='Network characteristics CSV (mul_tree_final_stats.csv)')
    parser.add_argument('--export-aggregated', metavar='FILE',
                       help='Export aggregated metrics to CSV file')

    args = parser.parse_args()

    # Load comparisons
    try:
        comparisons_df = pd.read_csv(args.comparisons)
        print(f"Loaded comparisons: {len(comparisons_df)} rows")
    except Exception as e:
        print(f"Error loading comparisons: {e}", file=sys.stderr)
        sys.exit(1)

    # Aggregate replicates
    print("\nAggregating replicates...")
    aggregator = MetricAggregator(comparisons_df)
    aggregated_df = aggregator.aggregate_replicates()

    print(f"Aggregated to {len(aggregated_df)} unique (network, config, method, metric) combinations")

    # Export aggregated if requested
    if args.export_aggregated:
        aggregated_df.to_csv(args.export_aggregated, index=False)
        print(f"Aggregated metrics exported to: {args.export_aggregated}")

    # Generate multi-level summaries
    # Get config from first row (all should be the same)
    config = aggregated_df['config'].iloc[0] if len(aggregated_df) > 0 else "unknown"

    summary_gen = MultiLevelSummary(aggregated_df, config)
    summary_gen.generate_all_summaries(Path(args.output_dir), args.network_stats)

    return aggregated_df


if __name__ == '__main__':
    aggregated_df = main()
