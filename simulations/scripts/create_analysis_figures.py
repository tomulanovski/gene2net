#!/usr/bin/env python3
"""
create_analysis_figures.py - Publication-Quality Analysis Figures

Generates comprehensive, clean visualizations for phylogenetic network inference methods.
Creates combined plots, per-method plots, and folding method comparisons.

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
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['lines.markersize'] = 8

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
        self.plots_individual_dir = self.plots_dir / "individual_methods"
        self.tables_dir = self.base_dir / "tables"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.plots_individual_dir.mkdir(parents=True, exist_ok=True)
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

        # Core completion rate plots - combined and per-method
        print("[1/15] Completion Rate vs Holm Fold (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'H_Strict', 'Number of Reticulations (Holm Fold)', '01_combined')

        print("[2/15] Completion Rate vs Holm Fold (per method)...")
        self.plot_completion_vs_characteristic_individual(
            'H_Strict', 'Number of Reticulations (Holm Fold)', '01_individual')

        print("[3/15] Completion Rate vs Polyphest Fold (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'H_Relaxed', 'Number of Reticulations (Polyphest Fold)', '02_combined')

        print("[4/15] Completion Rate vs Polyphest Fold (per method)...")
        self.plot_completion_vs_characteristic_individual(
            'H_Relaxed', 'Number of Reticulations (Polyphest Fold)', '02_individual')

        print("[5/15] Folding Method Comparison (completion rates)...")
        self.plot_folding_comparison()

        print("[6/15] Folding Method Accuracy Comparison...")
        self.plot_folding_accuracy_comparison()

        print("[7/15] Completion Rate vs Polyploids (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'Num_Polyploids', 'Number of Polyploid Species', '03_combined')

        print("[8/15] Completion Rate vs Polyploids (per method)...")
        self.plot_completion_vs_characteristic_individual(
            'Num_Polyploids', 'Number of Polyploid Species', '03_individual')

        print("[9/15] Completion Rate vs Total WGD (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'Total_WGD', 'Total WGD Events', '04_combined')

        print("[10/15] Completion Rate vs Total WGD (per method)...")
        self.plot_completion_vs_characteristic_individual(
            'Total_WGD', 'Total WGD Events', '04_individual')

        # Performance and accuracy plots
        print("[11/15] Reticulation Error Distribution...")
        self.plot_reticulation_error_distribution()

        print("[12/15] Edit Distance Distribution...")
        self.plot_edit_distance_distribution()

        print("[13/15] Accuracy Heatmap...")
        self.plot_accuracy_heatmap()

        print("[14/15] Per-Network Completion Breakdown...")
        self.plot_per_network_breakdown()

        print("[15/15] Method Performance Summary...")
        self.plot_method_summary()

        # Generate summary table
        print("[Tables] Generating summary tables...")
        self.generate_summary_table()

        print(f"\n{'='*80}")
        print(f"✓ Analysis complete! Outputs in:")
        print(f"  Combined plots: {self.plots_dir}")
        print(f"  Individual plots: {self.plots_individual_dir}")
        print(f"  Tables: {self.tables_dir}")
        print(f"{'='*80}\n")

    def plot_completion_vs_characteristic_combined(self, char_col: str, char_label: str, fig_prefix: str):
        """Plot completion rate vs network characteristic (all methods combined with error bars)"""
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
        inv = inv.dropna(subset=[char_col])

        # Plot each method
        for method in sorted(inv['method'].unique()):
            method_inv = inv[inv['method'] == method]

            # Calculate completion rate and variability per characteristic value
            grouped = []
            for char_val in sorted(method_inv[char_col].unique()):
                char_data = method_inv[method_inv[char_col] == char_val]

                # Calculate per-network completion rates
                network_rates = []
                for network in char_data['network'].unique():
                    net_data = char_data[char_data['network'] == network]
                    rate = net_data['inferred_exists'].sum() / len(net_data) * 100
                    network_rates.append(rate)

                grouped.append({
                    char_col: char_val,
                    'completion_rate': np.mean(network_rates),
                    'std_err': np.std(network_rates) / np.sqrt(len(network_rates)) if len(network_rates) > 1 else 0,
                    'n_networks': len(network_rates),
                    'n_runs': len(char_data)
                })

            grouped_df = pd.DataFrame(grouped)

            if len(grouped_df) > 0:
                # Plot with error bars
                ax.errorbar(grouped_df[char_col], grouped_df['completion_rate'],
                           yerr=grouped_df['std_err'],
                           marker=METHOD_MARKERS.get(method, 'o'),
                           label=method,
                           color=METHOD_COLORS.get(method, '#000000'),
                           linewidth=2.5,
                           markersize=9,
                           capsize=5,
                           capthick=2,
                           alpha=0.85,
                           markeredgewidth=1.5,
                           markeredgecolor='white')

        ax.set_xlabel(char_label, fontsize=14, fontweight='bold')
        ax.set_ylabel('Completion Rate (%)', fontsize=14, fontweight='bold')
        ax.set_title(f'Completion Rate vs {char_label}\nILS {self.ils_level}',
                    fontsize=15, fontweight='bold', pad=20)
        ax.legend(frameon=True, loc='best', fontsize=12, framealpha=0.9)
        ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.8)
        ax.set_ylim(-5, 105)

        if inv[char_col].dtype in ['int64', 'int32']:
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        plt.tight_layout()
        fig.savefig(self.plots_dir / f"{fig_prefix}_{char_col.lower()}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"{fig_prefix}_{char_col.lower()}.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_completion_vs_characteristic_individual(self, char_col: str, char_label: str, fig_prefix: str):
        """Plot completion rate vs characteristic - one file per method"""
        if self.inventory is None:
            return

        if char_col not in self.network_stats.columns:
            return

        inv = self.inventory.merge(
            self.network_stats[['network', char_col]],
            on='network', how='left'
        )
        inv = inv.dropna(subset=[char_col])

        methods = sorted(inv['method'].unique())

        for method in methods:
            fig, ax = plt.subplots(figsize=(10, 6))

            method_inv = inv[inv['method'] == method]

            # Calculate stats per characteristic value
            grouped = []
            for char_val in sorted(method_inv[char_col].unique()):
                char_data = method_inv[method_inv[char_col] == char_val]

                network_rates = []
                for network in char_data['network'].unique():
                    net_data = char_data[char_data['network'] == network]
                    rate = net_data['inferred_exists'].sum() / len(net_data) * 100
                    network_rates.append(rate)

                grouped.append({
                    char_col: char_val,
                    'completion_rate': np.mean(network_rates),
                    'std_err': np.std(network_rates) / np.sqrt(len(network_rates)) if len(network_rates) > 1 else 0,
                    'n_networks': len(network_rates)
                })

            grouped_df = pd.DataFrame(grouped)

            if len(grouped_df) > 0:
                ax.errorbar(grouped_df[char_col], grouped_df['completion_rate'],
                           yerr=grouped_df['std_err'],
                           marker=METHOD_MARKERS.get(method, 'o'),
                           color=METHOD_COLORS.get(method, '#000000'),
                           linewidth=3,
                           markersize=10,
                           capsize=6,
                           capthick=2.5,
                           markeredgewidth=2,
                           markeredgecolor='white')

                # Add network count annotations
                for _, row in grouped_df.iterrows():
                    ax.annotate(f"n={int(row['n_networks'])}",
                               xy=(row[char_col], row['completion_rate']),
                               xytext=(0, 10), textcoords='offset points',
                               ha='center', fontsize=9, color='gray')

            ax.set_xlabel(char_label, fontsize=13, fontweight='bold')
            ax.set_ylabel('Completion Rate (%)', fontsize=13, fontweight='bold')
            ax.set_title(f'{method}\nCompletion Rate vs {char_label} (ILS {self.ils_level})',
                        fontsize=14, fontweight='bold', pad=15)
            ax.grid(True, alpha=0.25, linestyle='--')
            ax.set_ylim(-5, 105)

            if inv[char_col].dtype in ['int64', 'int32']:
                ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

            plt.tight_layout()
            safe_method = method.replace('_', '-')
            fig.savefig(self.plots_individual_dir / f"{fig_prefix}_{char_col.lower()}_{safe_method}.pdf", bbox_inches='tight')
            fig.savefig(self.plots_individual_dir / f"{fig_prefix}_{char_col.lower()}_{safe_method}.png", bbox_inches='tight', dpi=300)
            plt.close()

    def plot_folding_comparison(self):
        """Compare Holm Fold vs Polyphest Fold - completion rates"""
        if self.inventory is None:
            return

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
            ax.legend(fontsize=10, framealpha=0.9)
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        axes[0].set_ylabel('Completion Rate (%)', fontsize=13, fontweight='bold')
        fig.suptitle(f'Folding Method Comparison: Completion Rates (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=1.02)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "05_folding_completion_comparison.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "05_folding_completion_comparison.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_folding_accuracy_comparison(self):
        """Compare folding methods: which produces more accurate reticulation counts?"""
        if self.metrics is None:
            print("  WARNING: No metrics data, skipping folding accuracy comparison")
            return

        # Get reticulation error data
        ret_error = self.metrics[self.metrics['metric'] == 'num_rets_diff'].copy()

        if len(ret_error) == 0:
            print("  WARNING: No num_rets_diff data found")
            return

        # We need to determine which folding was used for each method's output
        # This requires knowing the folding method used in network comparison
        # For now, we'll show absolute reticulation error by method

        methods = sorted(ret_error['method'].unique())
        n_methods = len(methods)

        fig, axes = plt.subplots(1, n_methods, figsize=(5*n_methods, 5), sharey=True)
        if n_methods == 1:
            axes = [axes]

        for ax, method in zip(axes, methods):
            method_data = ret_error[ret_error['method'] == method]

            # Get absolute mean error
            abs_errors = method_data['mean'].abs()

            ax.hist(abs_errors, bins=20, color=METHOD_COLORS.get(method, '#000000'),
                   alpha=0.7, edgecolor='black')

            mean_error = abs_errors.mean()
            ax.axvline(mean_error, color='red', linestyle='--', linewidth=2.5,
                      label=f'Mean = {mean_error:.2f}')

            ax.set_xlabel('Absolute Reticulation Count Error', fontsize=12, fontweight='bold')
            ax.set_title(method, fontsize=13, fontweight='bold')
            ax.legend(fontsize=10)
            ax.grid(True, alpha=0.25, axis='y')

        axes[0].set_ylabel('Frequency', fontsize=12, fontweight='bold')
        fig.suptitle(f'Reticulation Count Accuracy by Method (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', y=1.02)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "06_reticulation_accuracy.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "06_reticulation_accuracy.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_reticulation_error_distribution(self):
        """Boxplot of reticulation count errors"""
        if self.metrics is None:
            return

        ret_error = self.metrics[self.metrics['metric'] == 'num_rets_diff'].copy()

        if len(ret_error) == 0:
            return

        methods = sorted(ret_error['method'].unique())

        fig, ax = plt.subplots(figsize=(12, 7))

        data_by_method = []
        labels = []
        colors = []

        for method in methods:
            method_data = ret_error[ret_error['method'] == method]['mean'].dropna()
            if len(method_data) > 0:
                data_by_method.append(method_data)
                labels.append(method)
                colors.append(METHOD_COLORS.get(method, '#000000'))

        bp = ax.boxplot(data_by_method, labels=labels, patch_artist=True,
                       widths=0.6, showfliers=True,
                       boxprops=dict(linewidth=1.5),
                       whiskerprops=dict(linewidth=1.5),
                       capprops=dict(linewidth=1.5),
                       medianprops=dict(linewidth=2.5, color='red'))

        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.5, label='Perfect accuracy')
        ax.set_ylabel('Reticulation Count Difference\n(Inferred - True)', fontsize=14, fontweight='bold')
        ax.set_xlabel('Method', fontsize=14, fontweight='bold')
        ax.set_title(f'Reticulation Count Error Distribution (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax.legend(fontsize=11)
        plt.xticks(rotation=0, fontsize=12)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "07_reticulation_error_boxplot.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "07_reticulation_error_boxplot.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_edit_distance_distribution(self):
        """Plot edit distance distribution for each method"""
        if self.metrics is None:
            return

        edit_data = self.metrics[self.metrics['metric'] == 'edit_distance'].copy()

        if len(edit_data) == 0:
            return

        methods = sorted(edit_data['method'].unique())

        fig, ax = plt.subplots(figsize=(12, 7))

        data_by_method = []
        labels = []
        colors = []

        for method in methods:
            method_data = edit_data[edit_data['method'] == method]['mean'].dropna()
            if len(method_data) > 0:
                data_by_method.append(method_data)
                labels.append(method)
                colors.append(METHOD_COLORS.get(method, '#000000'))

        bp = ax.boxplot(data_by_method, labels=labels, patch_artist=True,
                       widths=0.6, showfliers=True,
                       boxprops=dict(linewidth=1.5),
                       whiskerprops=dict(linewidth=1.5),
                       capprops=dict(linewidth=1.5),
                       medianprops=dict(linewidth=2.5, color='red'))

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
        fig.savefig(self.plots_dir / "08_edit_distance_boxplot.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "08_edit_distance_boxplot.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_accuracy_heatmap(self):
        """Plot heatmap of edit distance for each method × network"""
        if self.metrics is None:
            return

        edit_data = self.metrics[self.metrics['metric'] == 'edit_distance'].copy()

        if len(edit_data) == 0:
            return

        pivot = edit_data.pivot(index='method', columns='network', values='mean')

        network_order = self.network_stats.sort_values('H_Strict')['network'].tolist()
        pivot = pivot[[col for col in network_order if col in pivot.columns]]

        fig, ax = plt.subplots(figsize=(18, 6))

        sns.heatmap(pivot, annot=False, cmap='RdYlGn_r', center=0.5,
                   cbar_kws={'label': 'Edit Distance'},
                   linewidths=0.5, linecolor='gray',
                   ax=ax, vmin=0, vmax=1)

        ax.set_xlabel('Network (sorted by H_Strict)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Method', fontsize=13, fontweight='bold')
        ax.set_title(f'Edit Distance Heatmap (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', pad=20)

        plt.xticks(rotation=90, fontsize=8)
        plt.yticks(rotation=0, fontsize=11)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "09_accuracy_heatmap.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "09_accuracy_heatmap.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_per_network_breakdown(self):
        """Show per-network completion rates to visualize aggregation"""
        if self.inventory is None:
            return

        # Merge with H_Strict for sorting
        inv = self.inventory.merge(
            self.network_stats[['network', 'H_Strict']],
            on='network', how='left'
        )

        methods = sorted(inv['method'].unique())
        networks_sorted = self.network_stats.sort_values('H_Strict')['network'].tolist()

        # Calculate completion rate per method per network
        data = []
        for method in methods:
            for network in networks_sorted:
                net_data = inv[(inv['method'] == method) & (inv['network'] == network)]
                if len(net_data) > 0:
                    comp_rate = net_data['inferred_exists'].sum() / len(net_data) * 100
                    h_strict = self.network_stats[self.network_stats['network'] == network]['H_Strict'].values[0]
                    data.append({
                        'method': method,
                        'network': network,
                        'completion_rate': comp_rate,
                        'H_Strict': h_strict
                    })

        df = pd.DataFrame(data)

        fig, ax = plt.subplots(figsize=(18, 6))

        # Plot grouped bars
        x = np.arange(len(networks_sorted))
        width = 0.8 / len(methods)

        for i, method in enumerate(methods):
            method_data = df[df['method'] == method]
            method_data = method_data.set_index('network').reindex(networks_sorted).reset_index()

            ax.bar(x + i*width, method_data['completion_rate'],
                  width, label=method,
                  color=METHOD_COLORS.get(method, '#000000'),
                  alpha=0.8, edgecolor='black', linewidth=0.5)

        ax.set_xlabel('Network (sorted by H_Strict)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Completion Rate (%)', fontsize=13, fontweight='bold')
        ax.set_title(f'Per-Network Completion Rates (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', pad=20)
        ax.set_xticks(x + width * len(methods) / 2)
        ax.set_xticklabels(networks_sorted, rotation=90, fontsize=8)
        ax.legend(fontsize=10, ncol=len(methods))
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax.set_ylim(0, 105)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "10_per_network_breakdown.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "10_per_network_breakdown.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_method_summary(self):
        """Summary bar plot: completion rate and mean edit distance"""
        if self.inventory is None or self.metrics is None:
            return

        methods = sorted(self.inventory['method'].unique())

        completion_rates = []
        edit_distances = []
        ret_errors = []

        for method in methods:
            method_inv = self.inventory[self.inventory['method'] == method]
            comp_rate = method_inv['inferred_exists'].sum() / len(method_inv) * 100
            completion_rates.append(comp_rate)

            method_edit = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'edit_distance')
            ]
            if len(method_edit) > 0:
                edit_distances.append(method_edit['mean'].mean())
            else:
                edit_distances.append(np.nan)

            method_ret = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'num_rets_diff')
            ]
            if len(method_ret) > 0:
                ret_errors.append(method_ret['mean'].abs().mean())
            else:
                ret_errors.append(np.nan)

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

        colors = [METHOD_COLORS.get(m, '#000000') for m in methods]

        # Completion rate
        bars1 = ax1.bar(methods, completion_rates, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax1.set_ylabel('Completion Rate (%)', fontsize=13, fontweight='bold')
        ax1.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax1.set_title('Completion Rate', fontsize=14, fontweight='bold', pad=15)
        ax1.set_ylim(0, 105)
        ax1.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax1.tick_params(axis='x', rotation=0)

        for bar, val in zip(bars1, completion_rates):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Edit distance
        bars2 = ax2.bar(methods, edit_distances, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax2.set_ylabel('Mean Edit Distance', fontsize=13, fontweight='bold')
        ax2.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax2.set_title('Edit Distance (lower = better)', fontsize=14, fontweight='bold', pad=15)
        ax2.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax2.tick_params(axis='x', rotation=0)

        for bar, val in zip(bars2, edit_distances):
            if not np.isnan(val):
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Reticulation error
        bars3 = ax3.bar(methods, ret_errors, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax3.set_ylabel('Mean Absolute Reticulation Error', fontsize=13, fontweight='bold')
        ax3.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax3.set_title('Reticulation Accuracy (lower = better)', fontsize=14, fontweight='bold', pad=15)
        ax3.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax3.tick_params(axis='x', rotation=0)

        for bar, val in zip(bars3, ret_errors):
            if not np.isnan(val):
                height = bar.get_height()
                ax3.text(bar.get_x() + bar.get_width()/2., height,
                        f'{val:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        fig.suptitle(f'Method Performance Summary (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=1.00)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "11_method_summary.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "11_method_summary.png", bbox_inches='tight', dpi=300)
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

            method_edit = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'edit_distance')
            ]

            if len(method_edit) > 0:
                mean_ed = method_edit['mean'].mean()
                std_ed = method_edit['mean'].std()
            else:
                mean_ed = std_ed = np.nan

            method_ret = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'num_rets_diff')
            ]

            if len(method_ret) > 0:
                mean_ret_err = method_ret['mean'].abs().mean()
            else:
                mean_ret_err = np.nan

            summary_data.append({
                'method': method,
                'total_runs': total,
                'completed_runs': successful,
                'completion_rate_%': comp_rate,
                'mean_edit_distance': mean_ed,
                'std_edit_distance': std_ed,
                'mean_reticulation_error': mean_ret_err
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
  │   ├── 01_combined_*.pdf/png           # All methods together
  │   ├── 02_combined_*.pdf/png
  │   ├── 03_combined_*.pdf/png
  │   ├── 04_combined_*.pdf/png
  │   ├── 05_folding_completion_comparison.pdf/png
  │   ├── 06_reticulation_accuracy.pdf/png
  │   ├── 07_reticulation_error_boxplot.pdf/png
  │   ├── 08_edit_distance_boxplot.pdf/png
  │   ├── 09_accuracy_heatmap.pdf/png
  │   ├── 10_per_network_breakdown.pdf/png
  │   ├── 11_method_summary.pdf/png
  │   └── individual_methods/             # Separate file per method
  │       ├── 01_individual_*_grampa.pdf/png
  │       ├── 01_individual_*_polyphest.pdf/png
  │       └── ...
  └── tables/
      └── summary_table.csv
        """
    )

    parser.add_argument('--config', nargs='+', required=True,
                       help='Configuration name(s) to analyze')

    default_stats_path = Path(__file__).parent.parent / 'networks' / 'mul_tree_final_stats.csv'
    parser.add_argument('--network-stats',
                       default=str(default_stats_path),
                       help=f'Path to network characteristics CSV (default: {default_stats_path})')

    args = parser.parse_args()

    for config in args.config:
        analyzer = ConfigurationAnalyzer(
            config=config,
            network_stats_file=args.network_stats
        )
        analyzer.generate_all_figures()


if __name__ == '__main__':
    main()
