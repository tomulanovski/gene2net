#!/usr/bin/env python3
"""
create_analysis_figures.py - Complete Analysis Figure Generation

Generates comprehensive analysis for phylogenetic network inference methods.
Organizes outputs by configuration with proper directory structure.

Usage:
    # Single configuration
    python create_analysis_figures.py --config conf_ils_low_10M

    # Multiple configurations
    python create_analysis_figures.py --config conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M

    # Custom network stats path (optional)
    python create_analysis_figures.py --config conf_ils_low_10M --network-stats /custom/path/stats.csv
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr
from typing import List, Dict
import warnings
warnings.filterwarnings('ignore')

# Publication settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['font.family'] = 'sans-serif'

# Color schemes
METHOD_COLORS = {
    'grampa': '#1f77b4',
    'polyphest_p50': '#ff7f0e',
    'polyphest_p70': '#2ca02c',
    'polyphest_p90': '#d62728',
    'mpsugar': '#9467bd',
    'padre': '#8c564b'
}


class ConfigurationAnalyzer:
    """Analyze a single configuration and generate all figures"""

    def __init__(self, config: str, network_stats_file: str):
        self.config = config
        self.config_name = config  # e.g., "conf_ils_low_10M"

        # Set up directories
        self.base_dir = Path(f"simulations/analysis/summary/{config}")
        self.plots_dir = self.base_dir / "plots"
        self.tables_dir = self.base_dir / "tables"

        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.tables_dir.mkdir(parents=True, exist_ok=True)

        # Load network statistics
        self.network_stats = pd.read_csv(network_stats_file)
        self.network_stats['network'] = self.network_stats['Filename'].str.replace('.tre', '')

        # Load configuration data
        self.inventory = None
        self.aggregated = None
        self.comparisons = None
        self._load_data()

        # Extract ILS level
        self.ils_level = self._extract_ils_level()

    def _extract_ils_level(self) -> str:
        """Extract ILS level from config name"""
        if 'low' in self.config.lower():
            return 'low'
        elif 'medium' in self.config.lower():
            return 'medium'
        elif 'high' in self.config.lower():
            return 'high'
        return 'unknown'

    def _load_data(self):
        """Load all data for this configuration"""
        if (self.base_dir / "inventory.csv").exists():
            self.inventory = pd.read_csv(self.base_dir / "inventory.csv")

        if (self.base_dir / "aggregated_metrics.csv").exists():
            self.aggregated = pd.read_csv(self.base_dir / "aggregated_metrics.csv")

        if (self.base_dir / "comparisons_raw.csv").exists():
            self.comparisons = pd.read_csv(self.base_dir / "comparisons_raw.csv")

    # ========================================================================
    # CORE PLOTS: Success vs Network Characteristics
    # ========================================================================

    def plot_success_vs_reticulations(self):
        """Plot number of successful runs as a function of number of reticulations"""
        if self.inventory is None:
            print("  WARNING: No inventory data, skipping success vs reticulations")
            return

        fig, ax = plt.subplots(figsize=(10, 6))

        # Merge with network stats
        inv = self.inventory.merge(
            self.network_stats[['network', 'H_Strict']],
            on='network', how='left'
        )

        # Calculate success rate per (H_Strict, method)
        for method in sorted(inv['method'].unique()):
            method_inv = inv[inv['method'] == method]

            # Group by H_Strict
            grouped = method_inv.groupby('H_Strict').apply(
                lambda x: pd.Series({
                    'success_rate': (x['status'] == 'exists').sum() / len(x) * 100,
                    'num_networks': len(x['network'].unique())
                })
            ).reset_index()

            ax.plot(grouped['H_Strict'], grouped['success_rate'],
                   'o-', label=method, color=METHOD_COLORS.get(method, '#000000'),
                   markersize=8, linewidth=2.5, alpha=0.8)

        ax.set_xlabel('Number of Reticulations (H_Strict)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Success Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title(f'Success Rate vs Number of Reticulations\n{self.config_name} (ILS {self.ils_level})',
                    fontsize=13, fontweight='bold', pad=15)
        ax.legend(frameon=True, loc='best', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_ylim(-5, 105)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "01_success_vs_reticulations.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "01_success_vs_reticulations.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_success_vs_polyploids(self):
        """Plot number of successful runs as a function of number of polyploids"""
        if self.inventory is None:
            print("  WARNING: No inventory data, skipping success vs polyploids")
            return

        fig, ax = plt.subplots(figsize=(10, 6))

        # Merge with network stats
        inv = self.inventory.merge(
            self.network_stats[['network', 'Num_Polyploids']],
            on='network', how='left'
        )

        # Calculate success rate per (Num_Polyploids, method)
        for method in sorted(inv['method'].unique()):
            method_inv = inv[inv['method'] == method]

            grouped = method_inv.groupby('Num_Polyploids').apply(
                lambda x: pd.Series({
                    'success_rate': (x['status'] == 'exists').sum() / len(x) * 100,
                    'num_networks': len(x['network'].unique())
                })
            ).reset_index()

            ax.plot(grouped['Num_Polyploids'], grouped['success_rate'],
                   'o-', label=method, color=METHOD_COLORS.get(method, '#000000'),
                   markersize=8, linewidth=2.5, alpha=0.8)

        ax.set_xlabel('Number of Polyploid Species', fontsize=12, fontweight='bold')
        ax.set_ylabel('Success Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title(f'Success Rate vs Number of Polyploids\n{self.config_name} (ILS {self.ils_level})',
                    fontsize=13, fontweight='bold', pad=15)
        ax.legend(frameon=True, loc='best', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_ylim(-5, 105)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "02_success_vs_polyploids.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "02_success_vs_polyploids.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_success_vs_wgd(self):
        """Plot number of successful runs as a function of total WGD events"""
        if self.inventory is None:
            print("  WARNING: No inventory data, skipping success vs WGD")
            return

        fig, ax = plt.subplots(figsize=(10, 6))

        # Merge with network stats
        inv = self.inventory.merge(
            self.network_stats[['network', 'Total_WGD']],
            on='network', how='left'
        )

        # Calculate success rate per (Total_WGD, method)
        for method in sorted(inv['method'].unique()):
            method_inv = inv[inv['method'] == method]

            grouped = method_inv.groupby('Total_WGD').apply(
                lambda x: pd.Series({
                    'success_rate': (x['status'] == 'exists').sum() / len(x) * 100,
                    'num_networks': len(x['network'].unique())
                })
            ).reset_index()

            ax.plot(grouped['Total_WGD'], grouped['success_rate'],
                   'o-', label=method, color=METHOD_COLORS.get(method, '#000000'),
                   markersize=8, linewidth=2.5, alpha=0.8)

        ax.set_xlabel('Total WGD Events (Autopolyploidization + Reticulations)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Success Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title(f'Success Rate vs Total WGD Events\n{self.config_name} (ILS {self.ils_level})',
                    fontsize=13, fontweight='bold', pad=15)
        ax.legend(frameon=True, loc='best', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_ylim(-5, 105)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "03_success_vs_wgd.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "03_success_vs_wgd.png", bbox_inches='tight', dpi=300)
        plt.close()

    # ========================================================================
    # ADDITIONAL ANALYSIS PLOTS
    # ========================================================================

    def plot_method_performance_overview(self):
        """Overall method performance: edit distance and success rate"""
        if self.aggregated is None or self.inventory is None:
            print("  WARNING: Missing data, skipping performance overview")
            return

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Panel 1: Edit distance boxplot
        edit_data = self.aggregated[self.aggregated['metric'] == 'edit_distance']

        if not edit_data.empty:
            methods = sorted(edit_data['method'].unique())
            data_for_box = []
            colors = []

            for method in methods:
                method_data = edit_data[edit_data['method'] == method]['mean'].values
                data_for_box.append(method_data)
                colors.append(METHOD_COLORS.get(method, '#cccccc'))

            bp = ax1.boxplot(data_for_box, labels=methods, patch_artist=True,
                            showfliers=True, widths=0.6)

            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            # Overlay points
            for i, (method, data) in enumerate(zip(methods, data_for_box), 1):
                x = np.random.normal(i, 0.04, size=len(data))
                ax1.scatter(x, data, alpha=0.4, s=30, color=colors[i-1],
                           edgecolors='black', linewidths=0.5)

        ax1.set_ylabel('Edit Distance (lower = better)', fontsize=11, fontweight='bold')
        ax1.set_title('Edit Distance Distribution', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')
        ax1.tick_params(axis='x', rotation=45)

        # Panel 2: Success rate
        success_rates = []
        for method in self.inventory['method'].unique():
            method_inv = self.inventory[self.inventory['method'] == method]
            total = len(method_inv)
            success = (method_inv['status'] == 'exists').sum()
            success_rate = success / total * 100 if total > 0 else 0
            success_rates.append({
                'method': method,
                'success_rate': success_rate,
                'total': total,
                'success': success
            })

        sr_df = pd.DataFrame(success_rates)
        methods = sr_df['method'].values
        rates = sr_df['success_rate'].values
        colors = [METHOD_COLORS.get(m, '#cccccc') for m in methods]

        bars = ax2.bar(range(len(methods)), rates, color=colors, alpha=0.7, edgecolor='black', linewidth=1)
        ax2.set_xticks(range(len(methods)))
        ax2.set_xticklabels(methods, rotation=45, ha='right')
        ax2.set_ylabel('Success Rate (%)', fontsize=11, fontweight='bold')
        ax2.set_title('Success Rate', fontsize=12, fontweight='bold')
        ax2.set_ylim(0, 105)
        ax2.grid(True, alpha=0.3, axis='y')

        # Add counts on bars
        for i, (bar, row) in enumerate(zip(bars, sr_df.itertuples())):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 2,
                    f'{row.success}/{row.total}',
                    ha='center', va='bottom', fontsize=8)

        fig.suptitle(f'Method Performance Overview - {self.config_name} (ILS {self.ils_level})',
                    fontsize=14, fontweight='bold', y=1.02)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "04_method_performance_overview.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "04_method_performance_overview.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_method_network_heatmap(self):
        """Heatmap showing edit distance for each (network, method) combination"""
        if self.aggregated is None:
            print("  WARNING: No aggregated data, skipping heatmap")
            return

        # Filter to edit_distance metric
        edit_data = self.aggregated[self.aggregated['metric'] == 'edit_distance']

        if edit_data.empty:
            print("  WARNING: No edit distance data, skipping heatmap")
            return

        # Pivot table
        pivot = edit_data.pivot_table(
            index='network',
            columns='method',
            values='mean',
            aggfunc='first'
        )

        # Sort by mean difficulty
        pivot['mean_difficulty'] = pivot.mean(axis=1)
        pivot = pivot.sort_values('mean_difficulty')
        pivot = pivot.drop('mean_difficulty', axis=1)

        # Plot
        fig, ax = plt.subplots(figsize=(10, 14))

        sns.heatmap(pivot, cmap='RdYlGn_r', center=0.5, vmin=0, vmax=1.5,
                   annot=False, fmt='.2f', linewidths=0.5,
                   cbar_kws={'label': 'Edit Distance'},
                   ax=ax)

        ax.set_title(f'Method × Network Performance Heatmap\n{self.config_name} (ILS {self.ils_level})',
                    fontsize=13, fontweight='bold', pad=15)
        ax.set_xlabel('Method', fontsize=11, fontweight='bold')
        ax.set_ylabel('Network (sorted by difficulty)', fontsize=11, fontweight='bold')

        plt.tight_layout()
        fig.savefig(self.plots_dir / "05_method_network_heatmap.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "05_method_network_heatmap.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_reticulation_accuracy(self):
        """Accuracy of reticulation count inference"""
        if self.aggregated is None:
            print("  WARNING: No aggregated data, skipping reticulation accuracy")
            return

        ret_data = self.aggregated[self.aggregated['metric'] == 'num_rets_diff']

        if ret_data.empty:
            print("  WARNING: No reticulation difference data, skipping")
            return

        fig, ax = plt.subplots(figsize=(10, 6))

        methods = sorted(ret_data['method'].unique())
        data_for_box = []
        colors = []

        for method in methods:
            method_data = ret_data[ret_data['method'] == method]['mean'].values
            data_for_box.append(method_data)
            colors.append(METHOD_COLORS.get(method, '#cccccc'))

        bp = ax.boxplot(data_for_box, labels=methods, patch_artist=True,
                       showfliers=True, widths=0.6)

        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.axhline(y=0, color='green', linestyle='--', linewidth=2, alpha=0.7,
                  label='Perfect inference')

        ax.set_ylabel('|Inferred H - True H|', fontsize=11, fontweight='bold')
        ax.set_title(f'Reticulation Inference Accuracy\n{self.config_name} (ILS {self.ils_level})',
                    fontsize=13, fontweight='bold', pad=15)
        ax.grid(True, alpha=0.3, axis='y')
        ax.tick_params(axis='x', rotation=45)
        ax.legend(fontsize=10)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "06_reticulation_accuracy.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "06_reticulation_accuracy.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_difficulty_correlations(self):
        """Correlation between network characteristics and edit distance"""
        if self.aggregated is None:
            print("  WARNING: No aggregated data, skipping correlations")
            return

        # Get edit distance data
        edit_data = self.aggregated[self.aggregated['metric'] == 'edit_distance']

        if edit_data.empty:
            return

        # Merge with network stats
        edit_with_stats = edit_data.merge(
            self.network_stats[['network', 'Num_Species', 'Num_Polyploids', 'Max_Copies',
                               'H_Strict', 'Num_Autopolyploidization_Events', 'Total_WGD']],
            on='network', how='left'
        )

        network_props = ['Num_Species', 'Num_Polyploids', 'Max_Copies', 'H_Strict',
                        'Num_Autopolyploidization_Events', 'Total_WGD']

        methods = sorted(edit_with_stats['method'].unique())

        # Build correlation matrix
        corr_matrix = pd.DataFrame(index=methods, columns=network_props)

        for method in methods:
            method_data = edit_with_stats[edit_with_stats['method'] == method]
            for prop in network_props:
                valid = method_data[[prop, 'mean']].dropna()
                if len(valid) > 3:
                    corr, _ = pearsonr(valid[prop], valid['mean'])
                    corr_matrix.loc[method, prop] = corr
                else:
                    corr_matrix.loc[method, prop] = np.nan

        corr_matrix = corr_matrix.astype(float)

        # Plot
        fig, ax = plt.subplots(figsize=(10, 7))

        sns.heatmap(corr_matrix, annot=True, cmap='RdYlGn_r', center=0,
                   vmin=-1, vmax=1, fmt='.2f', linewidths=0.5,
                   cbar_kws={'label': 'Pearson Correlation'},
                   ax=ax, annot_kws={'size': 9})

        ax.set_title(f'Network Characteristics vs Edit Distance\n{self.config_name} (ILS {self.ils_level})\n(Red = harder with more, Green = easier with more)',
                    fontsize=13, fontweight='bold', pad=15)
        ax.set_xlabel('Network Characteristics', fontsize=11, fontweight='bold')
        ax.set_ylabel('Method', fontsize=11, fontweight='bold')

        plt.tight_layout()
        fig.savefig(self.plots_dir / "07_difficulty_correlations.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "07_difficulty_correlations.png", bbox_inches='tight', dpi=300)
        plt.close()

    # ========================================================================
    # SUMMARY TABLE
    # ========================================================================

    def create_summary_table(self):
        """Create comprehensive summary table for this configuration"""
        if self.inventory is None or self.aggregated is None:
            print("  WARNING: Missing data, skipping summary table")
            return

        rows = []

        for method in self.inventory['method'].unique():
            # Success rate
            method_inv = self.inventory[self.inventory['method'] == method]
            total = len(method_inv)
            success = (method_inv['status'] == 'exists').sum()
            success_rate = success / total * 100 if total > 0 else 0

            # Edit distance
            method_edit = self.aggregated[
                (self.aggregated['method'] == method) &
                (self.aggregated['metric'] == 'edit_distance')
            ]

            if not method_edit.empty:
                mean_edit = method_edit['mean'].mean()
                std_edit = method_edit['std'].mean()
                min_edit = method_edit['mean'].min()
                max_edit = method_edit['mean'].max()
            else:
                mean_edit = std_edit = min_edit = max_edit = np.nan

            # Reticulation error
            method_ret = self.aggregated[
                (self.aggregated['method'] == method) &
                (self.aggregated['metric'] == 'num_rets_diff')
            ]

            if not method_ret.empty:
                mean_ret_error = method_ret['mean'].mean()
            else:
                mean_ret_error = np.nan

            rows.append({
                'Configuration': self.config_name,
                'ILS_Level': self.ils_level,
                'Method': method,
                'Total_Runs': total,
                'Successful_Runs': success,
                'Success_Rate_%': success_rate,
                'Mean_Edit_Distance': mean_edit,
                'Std_Edit_Distance': std_edit,
                'Min_Edit_Distance': min_edit,
                'Max_Edit_Distance': max_edit,
                'Mean_Reticulation_Error': mean_ret_error
            })

        summary_df = pd.DataFrame(rows)

        # Save
        summary_df.to_csv(self.tables_dir / "summary_table.csv",
                         index=False, float_format='%.4f')

        # Also create a simplified version for quick viewing
        simple_df = summary_df[['Method', 'Success_Rate_%', 'Mean_Edit_Distance',
                                'Mean_Reticulation_Error']].copy()
        simple_df.to_csv(self.tables_dir / "summary_simple.csv",
                        index=False, float_format='%.2f')

        return summary_df

    # ========================================================================
    # MAIN GENERATION
    # ========================================================================

    def generate_all_figures(self):
        """Generate all figures for this configuration"""
        print(f"\n{'='*80}")
        print(f"Generating Analysis for: {self.config_name}")
        print(f"ILS Level: {self.ils_level}")
        print(f"Output: {self.base_dir}")
        print(f"{'='*80}\n")

        print("[1/8] Success vs Reticulations...")
        self.plot_success_vs_reticulations()

        print("[2/8] Success vs Polyploids...")
        self.plot_success_vs_polyploids()

        print("[3/8] Success vs WGD Events...")
        self.plot_success_vs_wgd()

        print("[4/8] Method Performance Overview...")
        self.plot_method_performance_overview()

        print("[5/8] Method × Network Heatmap...")
        self.plot_method_network_heatmap()

        print("[6/8] Reticulation Accuracy...")
        self.plot_reticulation_accuracy()

        print("[7/8] Difficulty Correlations...")
        self.plot_difficulty_correlations()

        print("[8/8] Summary Tables...")
        self.create_summary_table()

        print(f"\n{'='*80}")
        print(f"Complete! All outputs saved to:")
        print(f"  Plots:  {self.plots_dir}")
        print(f"  Tables: {self.tables_dir}")
        print(f"{'='*80}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate comprehensive analysis figures for phylogenetic network inference',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single configuration
  %(prog)s --config conf_ils_low_10M

  # Multiple configurations (processes each separately)
  %(prog)s --config conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M

  # Custom network stats path (if needed)
  %(prog)s --config conf_ils_low_10M --network-stats /custom/path/stats.csv

Output structure (per configuration):
  simulations/analysis/summary/{config}/
  ├── plots/
  │   ├── 01_success_vs_reticulations.pdf/png
  │   ├── 02_success_vs_polyploids.pdf/png
  │   ├── 03_success_vs_wgd.pdf/png
  │   ├── 04_method_performance_overview.pdf/png
  │   ├── 05_method_network_heatmap.pdf/png
  │   ├── 06_reticulation_accuracy.pdf/png
  │   └── 07_difficulty_correlations.pdf/png
  └── tables/
      ├── summary_table.csv
      └── summary_simple.csv
        """
    )

    parser.add_argument('--config', nargs='+', required=True,
                       help='Configuration name(s) to analyze')
    parser.add_argument('--network-stats',
                       default='simulations/networks/mul_tree_final_stats.csv',
                       help='Path to network characteristics CSV (default: simulations/networks/mul_tree_final_stats.csv)')

    args = parser.parse_args()

    # Process each configuration
    for config in args.config:
        analyzer = ConfigurationAnalyzer(
            config=config,
            network_stats_file=args.network_stats
        )
        analyzer.generate_all_figures()


if __name__ == '__main__':
    main()
