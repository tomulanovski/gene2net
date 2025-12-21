#!/usr/bin/env python3
"""
create_analysis_figures.py - Publication-Quality Analysis Figures

Generates comprehensive, clean visualizations for phylogenetic network inference methods.
Creates both single and faceted plots for maximum clarity.

Usage:
    # Single configuration
    python create_analysis_figures.py --config conf_ils_low_10M

    # Multiple configurations
    python create_analysis_figures.py --config conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - no X11 required
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Dict
import warnings
warnings.filterwarnings('ignore')

# Publication-quality settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['legend.fontsize'] = 11
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['grid.linewidth'] = 0.8
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['lines.markersize'] = 7

# Enhanced color palette - distinct and colorblind-friendly
METHOD_COLORS = {
    'grampa': '#0173B2',        # Blue
    'polyphest': '#DE8F05',     # Orange
    'polyphest_p50': '#DE8F05', # Orange
    'polyphest_p70': '#029E73', # Teal
    'polyphest_p90': '#CC78BC', # Pink
    'mpsugar': '#CA9161',       # Tan
    'padre': '#ECE133'          # Yellow
}

# Markers for better distinction
METHOD_MARKERS = {
    'grampa': 'o',
    'polyphest': 's',
    'polyphest_p50': 's',
    'polyphest_p70': '^',
    'polyphest_p90': 'D',
    'mpsugar': 'v',
    'padre': 'P'
}


class ConfigurationAnalyzer:
    """Analyze and visualize results for a single configuration"""

    def __init__(self, config: str, network_stats_file: str):
        self.config = config

        # Extract ILS level from config name
        if 'low' in config.lower():
            self.ils_level = 'Low'
        elif 'medium' in config.lower():
            self.ils_level = 'Medium'
        elif 'high' in config.lower():
            self.ils_level = 'High'
        else:
            self.ils_level = 'Unknown'

        self.config_name = config.replace('conf_', '').replace('_10M', '')

        # Setup directories
        self.base_dir = Path(f"simulations/analysis/summary/{config}")
        self.plots_dir = self.base_dir / "plots"
        self.tables_dir = self.base_dir / "tables"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.tables_dir.mkdir(parents=True, exist_ok=True)

        # Load data
        self.network_stats = pd.read_csv(network_stats_file)
        # Remove .tre extension from network names if present
        self.network_stats['network'] = self.network_stats['Filename'].str.replace('.tre', '')

        # Load inventory
        inventory_file = self.base_dir / "inventory.csv"
        self.inventory = pd.read_csv(inventory_file) if inventory_file.exists() else None

        # Load aggregated metrics
        metrics_file = self.base_dir / "aggregated_metrics.csv"
        self.metrics = pd.read_csv(metrics_file) if metrics_file.exists() else None

        # Load comparisons
        comparisons_file = self.base_dir / "comparisons_raw.csv"
        self.comparisons = pd.read_csv(comparisons_file) if comparisons_file.exists() else None

        print(f"\nLoaded data for {config}:")
        print(f"  Networks: {len(self.network_stats)}")
        print(f"  Inventory: {len(self.inventory) if self.inventory is not None else 0}")
        print(f"  Metrics: {len(self.metrics) if self.metrics is not None else 0}")

    def generate_all_figures(self):
        """Generate all analysis figures"""
        print(f"\n{'='*80}")
        print(f"Generating Analysis for: {self.config}")
        print(f"ILS Level: {self.ils_level}")
        print(f"Output: {self.base_dir}")
        print(f"{'='*80}\n")

        # Core completion rate plots (both single and faceted)
        print("[1/12] Completion Rate vs Holm Fold (single)...")
        self.plot_completion_vs_characteristic('H_Strict', 'Holm Fold', '01a')

        print("[2/12] Completion Rate vs Holm Fold (faceted)...")
        self.plot_completion_vs_characteristic_faceted('H_Strict', 'Holm Fold', '01b')

        print("[3/12] Completion Rate vs Polyphest Fold (single)...")
        self.plot_completion_vs_characteristic('H_Relaxed', 'Polyphest Fold (threshold=0.2)', '02a')

        print("[4/12] Completion Rate vs Polyphest Fold (faceted)...")
        self.plot_completion_vs_characteristic_faceted('H_Relaxed', 'Polyphest Fold (threshold=0.2)', '02b')

        print("[5/12] Folding Method Comparison...")
        self.plot_folding_comparison()

        print("[6/12] Completion Rate vs Polyploids (single)...")
        self.plot_completion_vs_characteristic('Num_Polyploids', 'Number of Polyploid Species', '03a')

        print("[7/12] Completion Rate vs Polyploids (faceted)...")
        self.plot_completion_vs_characteristic_faceted('Num_Polyploids', 'Number of Polyploid Species', '03b')

        print("[8/12] Completion Rate vs Total WGD (single)...")
        self.plot_completion_vs_characteristic('Total_WGD', 'Total WGD Events', '04a')

        print("[9/12] Completion Rate vs Total WGD (faceted)...")
        self.plot_completion_vs_characteristic_faceted('Total_WGD', 'Total WGD Events', '04b')

        # Performance and accuracy plots
        print("[10/12] Edit Distance Distribution...")
        self.plot_edit_distance_distribution()

        print("[11/12] Accuracy Metrics Heatmap...")
        self.plot_accuracy_heatmap()

        print("[12/12] Method Performance Summary...")
        self.plot_method_summary()

        # Generate summary table
        print("[Tables] Generating summary tables...")
        self.generate_summary_table()

        print(f"\n{'='*80}")
        print(f"✓ Analysis complete! Outputs in:")
        print(f"  Plots: {self.plots_dir}")
        print(f"  Tables: {self.tables_dir}")
        print(f"{'='*80}\n")

    def plot_completion_vs_characteristic(self, char_col: str, char_label: str, fig_num: str):
        """Plot completion rate vs network characteristic (single plot, all methods)"""
        if self.inventory is None:
            print(f"  WARNING: No inventory data, skipping {char_label}")
            return

        # Check if column exists
        if char_col not in self.network_stats.columns:
            print(f"  WARNING: Column {char_col} not found, skipping")
            return

        fig, ax = plt.subplots(figsize=(12, 7))

        # Merge with network stats
        inv = self.inventory.merge(
            self.network_stats[['network', char_col]],
            on='network', how='left'
        )

        # Remove NaN values
        inv = inv.dropna(subset=[char_col])

        # Calculate completion rate per characteristic value per method
        for method in sorted(inv['method'].unique()):
            method_inv = inv[inv['method'] == method]

            grouped = method_inv.groupby(char_col).apply(
                lambda x: pd.Series({
                    'completion_rate': x['inferred_exists'].sum() / len(x) * 100,
                    'n_runs': len(x)
                })
            ).reset_index()

            # Only plot if we have data
            if len(grouped) > 0:
                ax.plot(grouped[char_col], grouped['completion_rate'],
                       marker=METHOD_MARKERS.get(method, 'o'),
                       label=method,
                       color=METHOD_COLORS.get(method, '#000000'),
                       linewidth=2.5,
                       markersize=8,
                       alpha=0.85,
                       markeredgewidth=1.5,
                       markeredgecolor='white')

        ax.set_xlabel(f'{char_label}', fontsize=14, fontweight='bold')
        ax.set_ylabel('Completion Rate (%)', fontsize=14, fontweight='bold')
        ax.set_title(f'Completion Rate vs {char_label}\nILS {self.ils_level}',
                    fontsize=15, fontweight='bold', pad=20)
        ax.legend(frameon=True, loc='best', fontsize=12, framealpha=0.9)
        ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.8)
        ax.set_ylim(-5, 105)

        # Integer x-axis if characteristic is integer
        if inv[char_col].dtype in ['int64', 'int32']:
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        plt.tight_layout()
        fig.savefig(self.plots_dir / f"{fig_num}_completion_vs_{char_col.lower()}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"{fig_num}_completion_vs_{char_col.lower()}.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_completion_vs_characteristic_faceted(self, char_col: str, char_label: str, fig_num: str):
        """Plot completion rate vs characteristic (faceted by method)"""
        if self.inventory is None:
            return

        if char_col not in self.network_stats.columns:
            return

        # Merge with network stats
        inv = self.inventory.merge(
            self.network_stats[['network', char_col]],
            on='network', how='left'
        )
        inv = inv.dropna(subset=[char_col])

        methods = sorted(inv['method'].unique())
        n_methods = len(methods)

        # Create subplots
        fig, axes = plt.subplots(1, n_methods, figsize=(5*n_methods, 5), sharey=True)
        if n_methods == 1:
            axes = [axes]

        for ax, method in zip(axes, methods):
            method_inv = inv[inv['method'] == method]

            grouped = method_inv.groupby(char_col).apply(
                lambda x: pd.Series({
                    'completion_rate': x['inferred_exists'].sum() / len(x) * 100,
                    'n_runs': len(x)
                })
            ).reset_index()

            if len(grouped) > 0:
                ax.plot(grouped[char_col], grouped['completion_rate'],
                       marker=METHOD_MARKERS.get(method, 'o'),
                       color=METHOD_COLORS.get(method, '#000000'),
                       linewidth=3,
                       markersize=10,
                       markeredgewidth=1.5,
                       markeredgecolor='white')

                ax.set_xlabel(f'{char_label}', fontsize=13, fontweight='bold')
                ax.set_title(method, fontsize=14, fontweight='bold', pad=10)
                ax.grid(True, alpha=0.25, linestyle='--')
                ax.set_ylim(-5, 105)

                if inv[char_col].dtype in ['int64', 'int32']:
                    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        axes[0].set_ylabel('Completion Rate (%)', fontsize=13, fontweight='bold')
        fig.suptitle(f'Completion Rate vs {char_label} (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=1.02)

        plt.tight_layout()
        fig.savefig(self.plots_dir / f"{fig_num}_completion_vs_{char_col.lower()}_faceted.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"{fig_num}_completion_vs_{char_col.lower()}_faceted.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_folding_comparison(self):
        """Compare Holm Fold vs Polyphest Fold"""
        if self.inventory is None:
            return

        # Check both columns exist
        if 'H_Strict' not in self.network_stats.columns or 'H_Relaxed' not in self.network_stats.columns:
            print("  WARNING: Missing H_Strict or H_Relaxed, skipping folding comparison")
            return

        inv = self.inventory.merge(
            self.network_stats[['network', 'H_Strict', 'H_Relaxed']],
            on='network', how='left'
        )

        methods = sorted(inv['method'].unique())
        n_methods = len(methods)

        fig, axes = plt.subplots(1, n_methods, figsize=(5*n_methods, 5), sharey=True)
        if n_methods == 1:
            axes = [axes]

        for ax, method in zip(axes, methods):
            method_inv = inv[inv['method'] == method]

            # Holm Fold
            grouped_strict = method_inv.groupby('H_Strict').apply(
                lambda x: pd.Series({'completion_rate': x['inferred_exists'].sum() / len(x) * 100})
            ).reset_index()

            # Polyphest Fold
            grouped_relaxed = method_inv.groupby('H_Relaxed').apply(
                lambda x: pd.Series({'completion_rate': x['inferred_exists'].sum() / len(x) * 100})
            ).reset_index()

            if len(grouped_strict) > 0:
                ax.plot(grouped_strict['H_Strict'], grouped_strict['completion_rate'],
                       'o-', label='Holm Fold', color='#0173B2', linewidth=2.5, markersize=8,
                       markeredgewidth=1.5, markeredgecolor='white')

            if len(grouped_relaxed) > 0:
                ax.plot(grouped_relaxed['H_Relaxed'], grouped_relaxed['completion_rate'],
                       's--', label='Polyphest Fold (θ=0.2)', color='#DE8F05', linewidth=2.5, markersize=8,
                       markeredgewidth=1.5, markeredgecolor='white')

            ax.set_xlabel('Number of Reticulations', fontsize=13, fontweight='bold')
            ax.set_title(method, fontsize=14, fontweight='bold', pad=10)
            ax.grid(True, alpha=0.25, linestyle='--')
            ax.set_ylim(-5, 105)
            ax.legend(fontsize=11, framealpha=0.9)
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        axes[0].set_ylabel('Completion Rate (%)', fontsize=13, fontweight='bold')
        fig.suptitle(f'Folding Method Comparison (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=1.02)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "05_folding_comparison.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "05_folding_comparison.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_edit_distance_distribution(self):
        """Plot edit distance distribution for each method"""
        if self.metrics is None:
            print("  WARNING: No metrics data, skipping edit distance")
            return

        # Filter to edit_distance metric
        edit_data = self.metrics[self.metrics['metric'] == 'edit_distance'].copy()

        if len(edit_data) == 0:
            print("  WARNING: No edit_distance data found")
            return

        methods = sorted(edit_data['method'].unique())

        fig, ax = plt.subplots(figsize=(12, 7))

        # Prepare data for boxplot
        data_by_method = []
        labels = []
        colors = []

        for method in methods:
            method_data = edit_data[edit_data['method'] == method]['mean'].dropna()
            if len(method_data) > 0:
                data_by_method.append(method_data)
                labels.append(method)
                colors.append(METHOD_COLORS.get(method, '#000000'))

        # Create boxplot
        bp = ax.boxplot(data_by_method, labels=labels, patch_artist=True,
                       widths=0.6, showfliers=True,
                       boxprops=dict(linewidth=1.5),
                       whiskerprops=dict(linewidth=1.5),
                       capprops=dict(linewidth=1.5),
                       medianprops=dict(linewidth=2.5, color='red'))

        # Color boxes
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_ylabel('Edit Distance', fontsize=14, fontweight='bold')
        ax.set_xlabel('Method', fontsize=14, fontweight='bold')
        ax.set_title(f'Edit Distance Distribution (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')
        plt.xticks(rotation=0, fontsize=12)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "06_edit_distance_distribution.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "06_edit_distance_distribution.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_accuracy_heatmap(self):
        """Plot heatmap of edit distance for each method × network"""
        if self.metrics is None:
            return

        edit_data = self.metrics[self.metrics['metric'] == 'edit_distance'].copy()

        if len(edit_data) == 0:
            return

        # Pivot to method × network matrix
        pivot = edit_data.pivot(index='method', columns='network', values='mean')

        # Sort networks by H_Strict for better visualization
        network_order = self.network_stats.sort_values('H_Strict')['network'].tolist()
        pivot = pivot[[col for col in network_order if col in pivot.columns]]

        fig, ax = plt.subplots(figsize=(16, 6))

        # Create heatmap
        sns.heatmap(pivot, annot=False, cmap='RdYlGn_r', center=0.5,
                   cbar_kws={'label': 'Edit Distance'},
                   linewidths=0.5, linecolor='gray',
                   ax=ax, vmin=0, vmax=1)

        ax.set_xlabel('Network', fontsize=13, fontweight='bold')
        ax.set_ylabel('Method', fontsize=13, fontweight='bold')
        ax.set_title(f'Edit Distance Heatmap (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', pad=20)

        plt.xticks(rotation=90, fontsize=9)
        plt.yticks(rotation=0, fontsize=11)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "07_accuracy_heatmap.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "07_accuracy_heatmap.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_method_summary(self):
        """Summary bar plot: completion rate and mean edit distance"""
        if self.inventory is None or self.metrics is None:
            return

        methods = sorted(self.inventory['method'].unique())

        # Calculate completion rates
        completion_rates = []
        edit_distances = []

        for method in methods:
            method_inv = self.inventory[self.inventory['method'] == method]
            comp_rate = method_inv['inferred_exists'].sum() / len(method_inv) * 100
            completion_rates.append(comp_rate)

            # Get mean edit distance
            method_edit = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'edit_distance')
            ]
            if len(method_edit) > 0:
                edit_distances.append(method_edit['mean'].mean())
            else:
                edit_distances.append(np.nan)

        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Completion rate bars
        colors1 = [METHOD_COLORS.get(m, '#000000') for m in methods]
        bars1 = ax1.bar(methods, completion_rates, color=colors1, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax1.set_ylabel('Completion Rate (%)', fontsize=13, fontweight='bold')
        ax1.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax1.set_title('Completion Rate by Method', fontsize=14, fontweight='bold', pad=15)
        ax1.set_ylim(0, 105)
        ax1.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax1.tick_params(axis='x', rotation=0)

        # Add value labels on bars
        for bar, val in zip(bars1, completion_rates):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Edit distance bars
        colors2 = [METHOD_COLORS.get(m, '#000000') for m in methods]
        bars2 = ax2.bar(methods, edit_distances, color=colors2, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax2.set_ylabel('Mean Edit Distance', fontsize=13, fontweight='bold')
        ax2.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax2.set_title('Accuracy by Method (lower is better)', fontsize=14, fontweight='bold', pad=15)
        ax2.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax2.tick_params(axis='x', rotation=0)

        # Add value labels on bars
        for bar, val in zip(bars2, edit_distances):
            if not np.isnan(val):
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        fig.suptitle(f'Method Performance Summary (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=1.00)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "08_method_summary.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "08_method_summary.png", bbox_inches='tight', dpi=300)
        plt.close()

    def generate_summary_table(self):
        """Generate summary statistics table"""
        if self.inventory is None or self.metrics is None:
            return

        methods = sorted(self.inventory['method'].unique())

        summary_data = []
        for method in methods:
            method_inv = self.inventory[self.inventory['method'] == method]
            total = len(method_inv)
            successful = method_inv['inferred_exists'].sum()
            comp_rate = successful / total * 100 if total > 0 else 0

            # Edit distance stats
            method_edit = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'edit_distance')
            ]

            if len(method_edit) > 0:
                mean_ed = method_edit['mean'].mean()
                std_ed = method_edit['mean'].std()
                min_ed = method_edit['mean'].min()
                max_ed = method_edit['mean'].max()
            else:
                mean_ed = std_ed = min_ed = max_ed = np.nan

            summary_data.append({
                'method': method,
                'total_runs': total,
                'completed_runs': successful,
                'completion_rate_%': comp_rate,
                'mean_edit_distance': mean_ed,
                'std_edit_distance': std_ed,
                'min_edit_distance': min_ed,
                'max_edit_distance': max_ed
            })

        df = pd.DataFrame(summary_data)
        df.to_csv(self.tables_dir / "summary_table.csv", index=False, float_format='%.4f')

        print(f"\n  Summary Table:")
        print(df.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(
        description='Generate publication-quality analysis figures',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output structure (per configuration):
  simulations/analysis/summary/{config}/
  ├── plots/
  │   ├── 01a_completion_vs_h_strict.pdf/png       # Single plot
  │   ├── 01b_completion_vs_h_strict_faceted.pdf/png
  │   ├── 02a_completion_vs_h_relaxed.pdf/png
  │   ├── 02b_completion_vs_h_relaxed_faceted.pdf/png
  │   ├── 03a_completion_vs_num_polyploids.pdf/png
  │   ├── 03b_completion_vs_num_polyploids_faceted.pdf/png
  │   ├── 04a_completion_vs_total_wgd.pdf/png
  │   ├── 04b_completion_vs_total_wgd_faceted.pdf/png
  │   ├── 05_folding_comparison.pdf/png
  │   ├── 06_edit_distance_distribution.pdf/png
  │   ├── 07_accuracy_heatmap.pdf/png
  │   └── 08_method_summary.pdf/png
  └── tables/
      └── summary_table.csv
        """
    )

    parser.add_argument('--config', nargs='+', required=True,
                       help='Configuration name(s) to analyze')

    # Default path relative to this script's location
    default_stats_path = Path(__file__).parent.parent / 'networks' / 'mul_tree_final_stats.csv'
    parser.add_argument('--network-stats',
                       default=str(default_stats_path),
                       help=f'Path to network characteristics CSV (default: {default_stats_path})')

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
