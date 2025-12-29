#!/usr/bin/env python3
"""
create_analysis_figures.py - Publication-Quality Analysis Figures

Generates comprehensive, clean visualizations for phylogenetic network inference methods.
Creates combined plots and faceted per-method views.

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
    'mpsugar': '#8B4513',       # Saddle Brown (darker, distinct from orange)
    'padre': '#ECE133',         # Yellow
    'alloppnet': '#DC143C'      # Crimson/Red
}

# Markers for better distinction
METHOD_MARKERS = {
    'grampa': 'o',
    'polyphest': 's',
    'polyphest_p50': 's',
    'polyphest_p70': '^',
    'polyphest_p90': 'D',
    'mpsugar': 'v',
    'padre': 'P',
    'alloppnet': 'X'
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

        # Enrich network stats with derived metrics
        self._prepare_enriched_stats()

        print(f"\nLoaded data for {config}:")
        print(f"  Networks: {len(self.network_stats)}")
        print(f"  Inventory: {len(self.inventory) if self.inventory is not None else 0}")
        print(f"  Metrics: {len(self.metrics) if self.metrics is not None else 0}")

    def _clean_output_directories(self):
        """Clean plots and tables directories before generating new figures"""
        print("Cleaning output directories...")
        
        # Clean plots directory
        if self.plots_dir.exists():
            plot_files = list(self.plots_dir.glob("*"))
            plot_files = [f for f in plot_files if f.is_file()]  # Only files, not subdirectories
            for file in plot_files:
                file.unlink()
            print(f"  Cleaned {len(plot_files)} files from plots/")
        
        # Clean individual_methods subdirectory
        if self.plots_individual_dir.exists():
            individual_files = list(self.plots_individual_dir.glob("*"))
            individual_files = [f for f in individual_files if f.is_file()]
            for file in individual_files:
                file.unlink()
            print(f"  Cleaned {len(individual_files)} files from plots/individual_methods/")
        
        # Clean tables directory
        if self.tables_dir.exists():
            table_files = list(self.tables_dir.glob("*"))
            table_files = [f for f in table_files if f.is_file()]
            for file in table_files:
                file.unlink()
            print(f"  Cleaned {len(table_files)} files from tables/")
        
        print("  ✓ Output directories cleaned (preserved run_full_summary files)\n")

    def _prepare_enriched_stats(self):
        """Add derived columns to network_stats for additional analyses"""
        # Polyploid ratio: proportion of species that are polyploid
        self.network_stats['Polyploid_Ratio'] = (
            self.network_stats['Num_Polyploids'] / self.network_stats['Num_Species']
        ).fillna(0)

        # Reticulation density: reticulations per species
        self.network_stats['Ret_Density'] = (
            self.network_stats['H_Strict'] / self.network_stats['Num_Species']
        ).fillna(0)

    def generate_all_figures(self):
        """Generate all analysis figures - comprehensive suite"""
        print(f"\n{'='*80}")
        print(f"Generating COMPREHENSIVE Analysis for: {self.config}")
        print(f"ILS Level: {self.ils_level}")
        print(f"Output: {self.base_dir}")
        print(f"{'='*80}\n")

        # Clean old plots and tables before generating new ones
        self._clean_output_directories()

        total_plots = 48  # Updated: added Jaccard distribution boxplots, per-method correlations, and per-network bias
        plot_num = 0

        # ========================================================================
        # CATEGORY 1: COMPLETION RATE vs NETWORK CHARACTERISTICS
        # ========================================================================
        print("\n" + "="*80)
        print("CATEGORY 1: Completion Rate vs Network Characteristics")
        print("="*80)

        # Original characteristics
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Holm Fold (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'H_Strict', 'Number of Reticulations (Holm Fold)', '01_combined_completion_vs_h_strict')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Holm Fold (faceted)...")
        self.plot_completion_vs_characteristic_faceted(
            'H_Strict', 'Number of Reticulations (Holm Fold)', '01_faceted_completion_vs_h_strict')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Polyphest Fold (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'H_Relaxed', 'Number of Reticulations (Polyphest Fold)', '02_combined_completion_vs_h_relaxed')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Polyphest Fold (faceted)...")
        self.plot_completion_vs_characteristic_faceted(
            'H_Relaxed', 'Number of Reticulations (Polyphest Fold)', '02_faceted_completion_vs_h_relaxed')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Polyploids (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'Num_Polyploids', 'Number of Polyploid Species', '03_combined_completion_vs_polyploids')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Polyploids (faceted)...")
        self.plot_completion_vs_characteristic_faceted(
            'Num_Polyploids', 'Number of Polyploid Species', '03_faceted_completion_vs_polyploids')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Total WGD (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'Total_WGD', 'Total WGD Events', '04_combined_completion_vs_total_wgd')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Total WGD (faceted)...")
        self.plot_completion_vs_characteristic_faceted(
            'Total_WGD', 'Total WGD Events', '04_faceted_completion_vs_total_wgd')

        # NEW: Additional network characteristics
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Num Species (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'Num_Species', 'Number of Species', '05_combined_completion_vs_num_species')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Num Species (faceted)...")
        self.plot_completion_vs_characteristic_faceted(
            'Num_Species', 'Number of Species', '05_faceted_completion_vs_num_species')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Max Copies (combined)...")
        self.plot_completion_vs_characteristic_combined(
            'Max_Copies', 'Maximum Copies per Species', '06_combined_completion_vs_max_copies')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion Rate vs Max Copies (faceted)...")
        self.plot_completion_vs_characteristic_faceted(
            'Max_Copies', 'Maximum Copies per Species', '06_faceted_completion_vs_max_copies')

        # ========================================================================
        # CATEGORY 2: EDIT DISTANCE (ACCURACY) vs NETWORK CHARACTERISTICS
        # ========================================================================
        print("\n" + "="*80)
        print("CATEGORY 2: MUL-tree Accuracy Metrics vs Network Characteristics")
        print("="*80)

        # MUL-tree Edit Distance plots (PRIMARY METRIC)
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] MUL-tree Edit Distance vs Num Species (combined)...")
        self.plot_accuracy_vs_characteristic_combined(
            'Num_Species', 'Number of Species', 'edit_distance_multree',
            'Edit Distance (MUL-tree)', '11_combined_editdist_multree_vs_num_species')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] MUL-tree Edit Distance vs Num Species (faceted)...")
        self.plot_accuracy_vs_characteristic_faceted(
            'Num_Species', 'Number of Species', 'edit_distance_multree',
            'Edit Distance (MUL-tree)', '11_faceted_editdist_multree_vs_num_species')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] MUL-tree Edit Distance vs H_Strict (combined)...")
        self.plot_accuracy_vs_characteristic_combined(
            'H_Strict', 'Number of Reticulations (Holm Fold)', 'edit_distance_multree',
            'Edit Distance (MUL-tree)', '12_combined_editdist_multree_vs_h_strict')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] MUL-tree Edit Distance vs H_Strict (faceted)...")
        self.plot_accuracy_vs_characteristic_faceted(
            'H_Strict', 'Number of Reticulations (Holm Fold)', 'edit_distance_multree',
            'Edit Distance (MUL-tree)', '12_faceted_editdist_multree_vs_h_strict')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] MUL-tree Edit Distance vs Num Polyploids (combined)...")
        self.plot_accuracy_vs_characteristic_combined(
            'Num_Polyploids', 'Number of Polyploid Species', 'edit_distance_multree',
            'Edit Distance (MUL-tree)', '13_combined_editdist_multree_vs_polyploids')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] MUL-tree Edit Distance vs Num Polyploids (faceted)...")
        self.plot_accuracy_vs_characteristic_faceted(
            'Num_Polyploids', 'Number of Polyploid Species', 'edit_distance_multree',
            'Edit Distance (MUL-tree)', '13_faceted_editdist_multree_vs_polyploids')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] MUL-tree Edit Distance vs Max Copies (combined)...")
        self.plot_accuracy_vs_characteristic_combined(
            'Max_Copies', 'Maximum Copies per Species', 'edit_distance_multree',
            'Edit Distance (MUL-tree)', '14_combined_editdist_multree_vs_max_copies')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] MUL-tree Edit Distance vs Max Copies (faceted)...")
        self.plot_accuracy_vs_characteristic_faceted(
            'Max_Copies', 'Maximum Copies per Species', 'edit_distance_multree',
            'Edit Distance (MUL-tree)', '14_faceted_editdist_multree_vs_max_copies')

        # RF Distance plots (PRIMARY METRIC)
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] RF Distance vs Num Species (combined)...")
        self.plot_accuracy_vs_characteristic_combined(
            'Num_Species', 'Number of Species', 'rf_distance',
            'Robinson-Foulds Distance (MUL-tree)', '15_combined_rf_vs_num_species')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] RF Distance vs Num Species (faceted)...")
        self.plot_accuracy_vs_characteristic_faceted(
            'Num_Species', 'Number of Species', 'rf_distance',
            'Robinson-Foulds Distance (MUL-tree)', '15_faceted_rf_vs_num_species')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] RF Distance vs H_Strict (combined)...")
        self.plot_accuracy_vs_characteristic_combined(
            'H_Strict', 'Number of Reticulations (Holm Fold)', 'rf_distance',
            'Robinson-Foulds Distance (MUL-tree)', '16_combined_rf_vs_h_strict')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] RF Distance vs H_Strict (faceted)...")
        self.plot_accuracy_vs_characteristic_faceted(
            'H_Strict', 'Number of Reticulations (Holm Fold)', 'rf_distance',
            'Robinson-Foulds Distance (MUL-tree)', '16_faceted_rf_vs_h_strict')

        # ========================================================================
        # CATEGORY 3: ADVANCED METRICS (Jaccard, Polyploid F1)
        # ========================================================================
        print("\n" + "="*80)
        print("CATEGORY 3: Advanced Performance Metrics")
        print("="*80)

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Reticulation Leaf Jaccard vs H_Strict (combined)...")
        self.plot_jaccard_vs_characteristic_combined(
            'H_Strict', 'Number of Reticulations (Holm Fold)',
            'ret_leaf_jaccard', 'Reticulation Leaf Set Jaccard',
            '21_combined_ret_leaf_jaccard_vs_h_strict')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Reticulation Leaf Jaccard vs H_Strict (faceted)...")
        self.plot_jaccard_vs_characteristic_faceted(
            'H_Strict', 'Number of Reticulations (Holm Fold)',
            'ret_leaf_jaccard', 'Reticulation Leaf Set Jaccard',
            '21_faceted_ret_leaf_jaccard_vs_h_strict')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Sister Relationship Jaccard vs H_Strict (combined)...")
        self.plot_jaccard_vs_characteristic_combined(
            'H_Strict', 'Number of Reticulations (Holm Fold)',
            'ret_sisters_jaccard', 'Sister Relationship Jaccard',
            '22_combined_ret_sisters_jaccard_vs_h_strict')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Sister Relationship Jaccard vs H_Strict (faceted)...")
        self.plot_jaccard_vs_characteristic_faceted(
            'H_Strict', 'Number of Reticulations (Holm Fold)',
            'ret_sisters_jaccard', 'Sister Relationship Jaccard',
            '22_faceted_ret_sisters_jaccard_vs_h_strict')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Polyploid Identification F1 Score...")
        self.plot_polyploid_f1_performance()

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Reticulation Leaf Jaccard Distribution...")
        self.plot_metric_distribution('ret_leaf_jaccard.dist',
                                      'Reticulation Leaf Set Distance',
                                      '08d_ret_leaf_jaccard_distribution')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Sister Relationship Jaccard Distribution...")
        self.plot_metric_distribution('ret_sisters_jaccard.dist',
                                      'Sister Relationship Distance',
                                      '08e_ret_sisters_jaccard_distribution')

        # ========================================================================
        # CATEGORY 4: DISTRIBUTIONS, COMPARISONS, AND SUMMARIES
        # ========================================================================
        print("\n" + "="*80)
        print("CATEGORY 4: Distributions, Comparisons, and Summary Plots")
        print("="*80)

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Folding Method Comparison (completion rates)...")
        self.plot_folding_comparison()

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Folding Method Accuracy Comparison...")
        self.plot_folding_accuracy_comparison()

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Reticulation Error Distribution...")
        self.plot_reticulation_error_distribution()

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] MUL-tree Edit Distance Distribution...")
        self.plot_edit_distance_distribution()
        
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] 3-Way Distance Metric Comparison...")
        self.plot_distance_metrics_comparison()

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Edit Distance MUL-tree Distribution...")
        self.plot_metric_distribution('edit_distance_multree',
                                      'Edit Distance (MUL-tree)',
                                      '08b_edit_distance_multree_distribution')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] RF Distance Distribution...")
        self.plot_metric_distribution('rf_distance',
                                      'Robinson-Foulds Distance (MUL-tree)',
                                      '08c_rf_distance_distribution')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Per-Network Completion Breakdown...")
        self.plot_per_network_breakdown()

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Per-Network Reticulation Bias...")
        self.plot_reticulation_bias_per_network()

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Method Performance Summary...")
        self.plot_method_summary()

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Comprehensive Correlation Heatmap (Aggregated)...")
        self.plot_comprehensive_correlation_heatmap()

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Per-Method Correlation Heatmaps...")
        self.plot_correlation_heatmap_per_method()

        # ========================================================================
        # TABLES
        # ========================================================================
        print("\n" + "="*80)
        print("Generating Summary Tables")
        print("="*80)
        self.generate_summary_tables()

        print(f"\n{'='*80}")
        print(f"✓ COMPREHENSIVE ANALYSIS COMPLETE!")
        print(f"  Generated {plot_num} plots across 4 categories")
        print(f"  Plots: {self.plots_dir}")
        print(f"  Faceted plots: {self.plots_individual_dir}")
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
                # Plot with error bars (scatter plot, no connecting lines - data is discrete)
                ax.errorbar(grouped_df[char_col], grouped_df['completion_rate'],
                           yerr=grouped_df['std_err'],
                           marker=METHOD_MARKERS.get(method, 'o'),
                           label=method,
                           color=METHOD_COLORS.get(method, '#000000'),
                           linestyle='None',  # No connecting lines - data is discrete
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

    def plot_completion_vs_characteristic_faceted(self, char_col: str, char_label: str, fig_prefix: str):
        """Plot completion rate vs characteristic - faceted subplots, one per method"""
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
        n_methods = len(methods)

        # Create faceted plot
        ncols = min(3, n_methods)
        nrows = (n_methods + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows), squeeze=False)
        axes = axes.flatten()

        for idx, method in enumerate(methods):
            ax = axes[idx]
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
                           linestyle='None',  # No connecting lines - data is discrete
                           markersize=10,
                           capsize=6,
                           capthick=2.5,
                           markeredgewidth=2,
                           markeredgecolor='white')

                # Only show network count if multiple networks contribute to same point
                for _, row in grouped_df.iterrows():
                    if row['n_networks'] > 1:
                        ax.annotate(f"n={int(row['n_networks'])}",
                                   xy=(row[char_col], row['completion_rate']),
                                   xytext=(0, 10), textcoords='offset points',
                                   ha='center', fontsize=9, color='dimgray',
                                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.7))

            ax.set_xlabel(char_label, fontsize=12, fontweight='bold')
            ax.set_ylabel('Completion Rate (%)', fontsize=12, fontweight='bold')
            ax.set_title(f'{method}',
                        fontsize=13, fontweight='bold', pad=10)
            ax.grid(True, alpha=0.25, linestyle='--')
            ax.set_ylim(-5, 105)

            if inv[char_col].dtype in ['int64', 'int32']:
                ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        # Hide unused subplots
        for idx in range(n_methods, len(axes)):
            axes[idx].set_visible(False)

        fig.suptitle(f'Completion Rate vs {char_label} (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=0.995)

        plt.tight_layout()
        fig.savefig(self.plots_individual_dir / f"{fig_prefix}_{char_col.lower()}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_individual_dir / f"{fig_prefix}_{char_col.lower()}.png", bbox_inches='tight', dpi=300)
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

        ncols = min(3, n_methods)
        nrows = (n_methods + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows), squeeze=False, sharey=True)
        axes = axes.flatten()

        for idx, method in enumerate(methods):
            ax = axes[idx]
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

            ax.set_xlabel('Number of Reticulations', fontsize=12, fontweight='bold')
            ax.set_title(method, fontsize=13, fontweight='bold', pad=10)
            ax.grid(True, alpha=0.25, linestyle='--')
            ax.set_ylim(-5, 105)
            ax.legend(fontsize=10, framealpha=0.9)
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        # Hide unused subplots
        for idx in range(n_methods, len(axes)):
            axes[idx].set_visible(False)

        axes[0].set_ylabel('Completion Rate (%)', fontsize=12, fontweight='bold')
        fig.suptitle(f'Folding Method Comparison: Completion Rates (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=0.995)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "05_folding_completion_comparison.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "05_folding_completion_comparison.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_folding_accuracy_comparison(self):
        """Compare folding methods: which produces more accurate reticulation counts? Shows bias."""
        if self.metrics is None:
            print("  WARNING: No metrics data, skipping folding accuracy comparison")
            return

        # Get reticulation bias data (signed errors)
        ret_bias = self.metrics[self.metrics['metric'] == 'num_rets_bias'].copy()

        if len(ret_bias) == 0:
            # Fall back to num_rets_diff for backward compatibility
            print("  WARNING: No num_rets_bias data found, falling back to num_rets_diff")
            ret_bias = self.metrics[self.metrics['metric'] == 'num_rets_diff'].copy()
            if len(ret_bias) == 0:
                print("  WARNING: No num_rets_diff data found")
                return

        methods = sorted(ret_bias['method'].unique())
        n_methods = len(methods)

        ncols = min(3, n_methods)
        nrows = (n_methods + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows), squeeze=False, sharey=True)
        axes = axes.flatten()

        for idx, method in enumerate(methods):
            ax = axes[idx]
            method_data = ret_bias[ret_bias['method'] == method]

            # Get signed errors (bias)
            biases = method_data['mean']

            # Create histogram
            ax.hist(biases, bins=20, color=METHOD_COLORS.get(method, '#000000'),
                   alpha=0.7, edgecolor='black')

            mean_bias = biases.mean()
            mae = biases.abs().mean()
            
            # Show both mean bias and MAE
            ax.axvline(mean_bias, color='red', linestyle='--', linewidth=2.5,
                      label=f'Mean bias = {mean_bias:+.2f}')
            ax.axvline(0, color='black', linestyle='-', linewidth=1.5, alpha=0.5)

            ax.set_xlabel('Reticulation Count Bias\n(+ = Over-estimation)', fontsize=11, fontweight='bold')
            ax.set_title(f'{method}\nMAE = {mae:.2f}', fontsize=12, fontweight='bold')
            ax.legend(fontsize=10)
            ax.grid(True, alpha=0.25, axis='y')

        # Hide unused subplots
        for idx in range(n_methods, len(axes)):
            axes[idx].set_visible(False)

        axes[0].set_ylabel('Frequency', fontsize=11, fontweight='bold')
        fig.suptitle(f'Reticulation Count Bias by Method (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', y=0.995)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "06_reticulation_bias_histogram.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "06_reticulation_bias_histogram.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_reticulation_error_distribution(self):
        """Boxplot of reticulation count errors - shows percentage bias (signed)"""
        if self.metrics is None:
            return

        # Use num_rets_bias (signed) instead of num_rets_diff (absolute)
        ret_bias = self.metrics[self.metrics['metric'] == 'num_rets_bias'].copy()

        if len(ret_bias) == 0:
            # Fall back to num_rets_diff for backward compatibility
            ret_bias = self.metrics[self.metrics['metric'] == 'num_rets_diff'].copy()
            if len(ret_bias) == 0:
                return
            metric_name = 'num_rets_diff'
            use_percentage = False
        else:
            metric_name = 'num_rets_bias'
            use_percentage = True

        # Merge with network stats to get true H_Strict for percentage calculation
        ret_bias = ret_bias.merge(
            self.network_stats[['network', 'H_Strict']],
            on='network',
            how='left'
        )

        methods = sorted(ret_bias['method'].unique())

        fig, ax = plt.subplots(figsize=(12, 7))

        data_by_method = []
        labels = []
        colors = []
        mean_biases = []

        for method in methods:
            method_data = ret_bias[ret_bias['method'] == method].copy()
            
            if use_percentage:
                # Calculate percentage bias: (bias / H_Strict) * 100
                # Handle division by zero: if H_Strict=0, use absolute bias
                method_data['bias_pct'] = (
                    method_data['mean'] / method_data['H_Strict'] * 100
                ).replace([np.inf, -np.inf], np.nan)
                
                # For networks with H_Strict=0, use absolute bias (can't calculate percentage)
                zero_h_mask = method_data['H_Strict'] == 0
                if zero_h_mask.any():
                    method_data.loc[zero_h_mask, 'bias_pct'] = method_data.loc[zero_h_mask, 'mean']
                
                method_values = method_data['bias_pct'].dropna()
            else:
                method_values = method_data['mean'].dropna()
            
            if len(method_values) > 0:
                data_by_method.append(method_values)
                labels.append(method)
                colors.append(METHOD_COLORS.get(method, '#000000'))
                mean_biases.append(method_values.mean())

        if len(data_by_method) == 0:
            return

        bp = ax.boxplot(data_by_method, labels=labels, patch_artist=True,
                       widths=0.6, showfliers=True,
                       boxprops=dict(linewidth=1.5),
                       whiskerprops=dict(linewidth=1.5),
                       capprops=dict(linewidth=1.5),
                       medianprops=dict(linewidth=2.5, color='red'))

        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Add mean bias annotations
        for i, (method, mean_bias) in enumerate(zip(labels, mean_biases), 1):
            sign = '+' if mean_bias >= 0 else ''
            ax.text(i, ax.get_ylim()[1] * 0.95, f'{sign}{mean_bias:.1f}%',
                   ha='center', va='top', fontsize=10, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8))

        ax.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.5, label='Perfect accuracy (0%)')
        
        if use_percentage:
            ylabel = 'Reticulation Error (%)\n(Inferred - True) / True × 100\n[Positive = Over-estimation]'
        else:
            ylabel = 'Reticulation Count Bias\n(Inferred - True)\n[Positive = Over-estimation]'
        
        ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
        ax.set_xlabel('Method', fontsize=14, fontweight='bold')
        
        if use_percentage:
            title = f'Reticulation Error Distribution - Percentage Bias (ILS {self.ils_level})\nMean percentage bias shown above each box'
        else:
            title = f'Reticulation Count Bias Distribution (ILS {self.ils_level})\nMean bias shown above each box'
        
        ax.set_title(title, fontsize=15, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax.legend(fontsize=11)
        plt.xticks(rotation=0, fontsize=12)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "07_reticulation_bias_boxplot.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "07_reticulation_bias_boxplot.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_edit_distance_distribution(self):
        """Plot MUL-tree edit distance distribution for each method"""
        if self.metrics is None:
            return

        # Use MUL-tree edit distance (PRIMARY METRIC)
        edit_data = self.metrics[self.metrics['metric'] == 'edit_distance_multree'].copy()
        
        # Fallback to network edit distance if MUL-tree not available
        if len(edit_data) == 0:
            edit_data = self.metrics[self.metrics['metric'] == 'edit_distance'].copy()
            metric_type = 'Network'
            if len(edit_data) == 0:
                return
        else:
            metric_type = 'MUL-tree'

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

        ax.set_ylabel(f'Edit Distance ({metric_type})\n(0 = identical, 1 = very different)', 
                     fontsize=14, fontweight='bold')
        ax.set_xlabel('Method', fontsize=14, fontweight='bold')
        ax.set_title(f'Edit Distance Distribution - {metric_type} (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')
        plt.xticks(rotation=0, fontsize=12)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "08_edit_distance_multree_boxplot.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "08_edit_distance_multree_boxplot.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_distance_metrics_comparison(self):
        """Compare all three distance metrics side-by-side: Network ED, MUL-tree ED, and RF"""
        if self.metrics is None:
            return

        methods = sorted(self.inventory['method'].unique()) if self.inventory is not None else []
        if len(methods) == 0:
            return

        # Collect data for all three metrics
        metrics_to_compare = {
            'edit_distance': 'Network Edit Distance',
            'edit_distance_multree': 'MUL-tree Edit Distance',
            'rf_distance': 'RF Distance (MUL-tree)'
        }
        
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        
        for idx, (metric_name, metric_label) in enumerate(metrics_to_compare.items()):
            ax = axes[idx]
            
            metric_data = self.metrics[self.metrics['metric'] == metric_name]
            
            if len(metric_data) == 0:
                ax.text(0.5, 0.5, f'No data for\n{metric_label}', 
                       ha='center', va='center', fontsize=14, color='gray')
                ax.set_title(metric_label, fontsize=13, fontweight='bold')
                ax.axis('off')
                continue
            
            # Prepare data for boxplot
            data_by_method = []
            labels = []
            colors = []
            means = []
            
            for method in methods:
                method_data = metric_data[metric_data['method'] == method]['mean'].dropna()
                if len(method_data) > 0:
                    data_by_method.append(method_data)
                    labels.append(method)
                    colors.append(METHOD_COLORS.get(method, '#000000'))
                    means.append(method_data.mean())
            
            if len(data_by_method) == 0:
                ax.text(0.5, 0.5, f'No data for\n{metric_label}', 
                       ha='center', va='center', fontsize=14, color='gray')
                ax.set_title(metric_label, fontsize=13, fontweight='bold')
                ax.axis('off')
                continue
            
            # Create boxplot
            bp = ax.boxplot(data_by_method, labels=labels, patch_artist=True,
                           widths=0.5, showfliers=True,
                           boxprops=dict(linewidth=1.5),
                           whiskerprops=dict(linewidth=1.5),
                           capprops=dict(linewidth=1.5),
                           medianprops=dict(linewidth=2.5, color='red'))
            
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            # Add mean values as text
            for i, (label, mean_val) in enumerate(zip(labels, means), 1):
                ax.text(i, ax.get_ylim()[1] * 0.95, f'{mean_val:.3f}',
                       ha='center', va='top', fontsize=9, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                edgecolor='gray', alpha=0.8))
            
            ax.set_ylabel('Distance\n(0 = identical, 1 = very different)', 
                         fontsize=12, fontweight='bold')
            ax.set_xlabel('Method', fontsize=12, fontweight='bold')
            ax.set_title(metric_label, fontsize=13, fontweight='bold', pad=15)
            ax.grid(True, alpha=0.25, axis='y', linestyle='--')
            ax.tick_params(axis='x', rotation=0, labelsize=10)
            
            # Highlight if this is the primary metric
            if metric_name in ['edit_distance_multree', 'rf_distance']:
                ax.patch.set_edgecolor('#2E8B57')
                ax.patch.set_linewidth(3)
        
        fig.suptitle(f'Distance Metrics Comparison (ILS {self.ils_level})\n' + 
                    'Green border = Primary metrics (MUL-tree based)',
                    fontsize=16, fontweight='bold', y=1.02)
        
        plt.tight_layout()
        fig.savefig(self.plots_dir / "08a_distance_metrics_comparison.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "08a_distance_metrics_comparison.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_metric_distribution(self, metric_name: str, metric_label: str, filename_prefix: str):
        """Generic method to plot distribution of any metric as box plots"""
        if self.metrics is None:
            print(f"  WARNING: No metrics data, skipping {metric_label}")
            return

        # Filter for this metric
        metric_data = self.metrics[self.metrics['metric'] == metric_name]

        if len(metric_data) == 0:
            print(f"  WARNING: No data for metric '{metric_name}', skipping")
            return

        methods = sorted(metric_data['method'].unique())

        # Prepare data for box plots
        plot_data = []
        for method in methods:
            method_data = metric_data[metric_data['method'] == method]
            plot_data.append(method_data['mean'].values)

        fig, ax = plt.subplots(figsize=(12, 7))

        # Create box plots
        bp = ax.boxplot(plot_data, labels=methods, patch_artist=True,
                       widths=0.6,
                       boxprops=dict(linewidth=1.5, edgecolor='black'),
                       whiskerprops=dict(linewidth=1.5),
                       capprops=dict(linewidth=1.5),
                       medianprops=dict(linewidth=2, color='red'))

        # Color boxes
        for patch, method in zip(bp['boxes'], methods):
            patch.set_facecolor(METHOD_COLORS.get(method, '#CCCCCC'))
            patch.set_alpha(0.7)

        ax.set_ylabel(metric_label, fontsize=14, fontweight='bold')
        ax.set_xlabel('Method', fontsize=14, fontweight='bold')
        ax.set_title(f'{metric_label} Distribution (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')
        plt.xticks(rotation=0, fontsize=12)

        plt.tight_layout()
        fig.savefig(self.plots_dir / f"{filename_prefix}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"{filename_prefix}.png", bbox_inches='tight', dpi=300)
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
        fig.savefig(self.plots_dir / "09_per_network_breakdown.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "09_per_network_breakdown.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_reticulation_bias_per_network(self):
        """Grouped bar chart showing reticulation bias (percentage) per network for all methods"""
        if self.metrics is None or self.network_stats is None:
            print("  WARNING: Missing data for per-network reticulation bias plot")
            return

        # Get num_rets_bias (signed error)
        ret_bias = self.metrics[self.metrics['metric'] == 'num_rets_bias'].copy()
        
        if len(ret_bias) == 0:
            print("  WARNING: No num_rets_bias data found, skipping per-network bias plot")
            return

        # Merge with network stats to get H_Strict for percentage calculation
        ret_bias = ret_bias.merge(
            self.network_stats[['network', 'H_Strict']],
            on='network',
            how='left'
        )

        # Calculate percentage bias
        ret_bias['bias_pct'] = (
            ret_bias['mean'] / ret_bias['H_Strict'] * 100
        ).replace([np.inf, -np.inf], np.nan)
        
        # For networks with H_Strict=0, use absolute bias
        zero_h_mask = ret_bias['H_Strict'] == 0
        if zero_h_mask.any():
            ret_bias.loc[zero_h_mask, 'bias_pct'] = ret_bias.loc[zero_h_mask, 'mean']

        methods = sorted(ret_bias['method'].unique())
        networks_sorted = self.network_stats.sort_values('H_Strict')['network'].tolist()

        # Prepare data for plotting
        data = []
        for method in methods:
            for network in networks_sorted:
                method_net_data = ret_bias[
                    (ret_bias['method'] == method) & 
                    (ret_bias['network'] == network)
                ]
                if len(method_net_data) > 0:
                    bias_pct = method_net_data['bias_pct'].values[0]
                    h_strict = method_net_data['H_Strict'].values[0]
                    data.append({
                        'method': method,
                        'network': network,
                        'bias_pct': bias_pct,
                        'H_Strict': h_strict
                    })
                else:
                    # No data for this method-network combination
                    data.append({
                        'method': method,
                        'network': network,
                        'bias_pct': np.nan,
                        'H_Strict': self.network_stats[
                            self.network_stats['network'] == network
                        ]['H_Strict'].values[0] if len(self.network_stats[
                            self.network_stats['network'] == network
                        ]) > 0 else 0
                    })

        df = pd.DataFrame(data)

        fig, ax = plt.subplots(figsize=(18, 7))

        # Plot grouped bars
        x = np.arange(len(networks_sorted))
        width = 0.8 / len(methods)

        for i, method in enumerate(methods):
            method_data = df[df['method'] == method]
            method_data = method_data.set_index('network').reindex(networks_sorted).reset_index()
            
            # Get method-specific color
            method_color = METHOD_COLORS.get(method, '#000000')
            
            # Plot each bar individually with method-specific color
            bias_values = method_data['bias_pct'].values
            bars = []
            for j, (network, bias) in enumerate(zip(networks_sorted, bias_values)):
                if np.isnan(bias):
                    # Gray bar for missing data (height 0, just a marker)
                    bar = ax.bar(x[j] + i*width, 0, width,
                               color='#CCCCCC', alpha=0.3,
                               edgecolor='black', linewidth=0.5)
                else:
                    # Use method-specific color for all bars
                    bar = ax.bar(x[j] + i*width, bias, width,
                               color=method_color, alpha=0.8,
                               edgecolor='black', linewidth=0.5)
                bars.append(bar[0])
            
            # Add label only for first bar of each method
            if len(bars) > 0:
                bars[0].set_label(method)

        ax.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.5, label='Perfect accuracy (0%)')
        ax.set_xlabel('Network (sorted by H_Strict)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Reticulation Bias (%)\n(Inferred - True) / True × 100', 
                     fontsize=13, fontweight='bold')
        ax.set_title(f'Per-Network Reticulation Bias (ILS {self.ils_level})',
                    fontsize=15, fontweight='bold', pad=20)
        ax.set_xticks(x + width * len(methods) / 2)
        ax.set_xticklabels(networks_sorted, rotation=90, fontsize=8)
        ax.legend(fontsize=10, ncol=len(methods) + 1)
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')

        plt.tight_layout()
        fig.savefig(self.plots_dir / "09b_per_network_reticulation_bias.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "09b_per_network_reticulation_bias.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_method_summary(self):
        """Summary bar plot: completion rate, edit distance, and reticulation error with bias"""
        if self.inventory is None or self.metrics is None:
            return

        methods = sorted(self.inventory['method'].unique())

        completion_rates = []
        edit_distances = []
        ret_errors = []
        ret_biases = []  # NEW: track bias

        for method in methods:
            method_inv = self.inventory[self.inventory['method'] == method]
            comp_rate = method_inv['inferred_exists'].sum() / len(method_inv) * 100
            completion_rates.append(comp_rate)

            # Use MUL-tree edit distance (PRIMARY METRIC)
            method_edit = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'edit_distance_multree')
            ]
            
            # Fallback to network edit distance if needed
            if len(method_edit) == 0:
                method_edit = self.metrics[
                    (self.metrics['method'] == method) &
                    (self.metrics['metric'] == 'edit_distance')
                ]
            
            if len(method_edit) > 0:
                edit_distances.append(method_edit['mean'].mean())
            else:
                edit_distances.append(np.nan)

            # Absolute error (MAE)
            method_ret = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'num_rets_diff')
            ]
            if len(method_ret) > 0:
                ret_errors.append(method_ret['mean'].mean())  # Already absolute
            else:
                ret_errors.append(np.nan)
            
            # Bias (signed error) - calculate as percentage
            method_bias = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'num_rets_bias')
            ]
            if len(method_bias) > 0:
                # Merge with network stats to calculate percentage bias
                method_bias_with_stats = method_bias.merge(
                    self.network_stats[['network', 'H_Strict']],
                    on='network',
                    how='left'
                )
                # Calculate percentage: (bias / H_Strict) * 100
                method_bias_with_stats['bias_pct'] = (
                    method_bias_with_stats['mean'] / method_bias_with_stats['H_Strict'] * 100
                ).replace([np.inf, -np.inf], np.nan)
                # For H_Strict=0, use absolute bias
                zero_h_mask = method_bias_with_stats['H_Strict'] == 0
                if zero_h_mask.any():
                    method_bias_with_stats.loc[zero_h_mask, 'bias_pct'] = method_bias_with_stats.loc[zero_h_mask, 'mean']
                
                ret_biases.append(method_bias_with_stats['bias_pct'].mean())
            else:
                ret_biases.append(np.nan)

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

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

        # Edit distance (MUL-tree)
        bars2 = ax2.bar(methods, edit_distances, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax2.set_ylabel('Mean Edit Distance (MUL-tree)', fontsize=13, fontweight='bold')
        ax2.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax2.set_title('MUL-tree Accuracy (lower = better)', fontsize=14, fontweight='bold', pad=15)
        ax2.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax2.tick_params(axis='x', rotation=0)

        for bar, val in zip(bars2, edit_distances):
            if not np.isnan(val):
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Reticulation absolute error (MAE)
        bars3 = ax3.bar(methods, ret_errors, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax3.set_ylabel('Mean Absolute Error (MAE)', fontsize=13, fontweight='bold')
        ax3.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax3.set_title('Reticulation Count: Absolute Error', fontsize=14, fontweight='bold', pad=15)
        ax3.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax3.tick_params(axis='x', rotation=0)

        for bar, val in zip(bars3, ret_errors):
            if not np.isnan(val):
                height = bar.get_height()
                ax3.text(bar.get_x() + bar.get_width()/2., height,
                        f'{val:.2f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

        # Reticulation bias (signed error)
        # Color bars based on bias direction: red for over-estimation, blue for under-estimation
        bias_colors = []
        for bias in ret_biases:
            if np.isnan(bias):
                bias_colors.append('#CCCCCC')
            elif bias > 0:
                bias_colors.append('#D62728')  # Red for over-estimation
            else:
                bias_colors.append('#1F77B4')  # Blue for under-estimation
        
        bars4 = ax4.bar(methods, ret_biases, color=bias_colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax4.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.5, label='No bias (0%)')
        ax4.set_ylabel('Mean Bias (%)\n(Signed Error / True × 100)', fontsize=13, fontweight='bold')
        ax4.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax4.set_title('Reticulation Count: Percentage Bias\n[+ = Over-estimation, - = Under-estimation]', 
                     fontsize=14, fontweight='bold', pad=15)
        ax4.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax4.tick_params(axis='x', rotation=0)
        ax4.legend(fontsize=11)

        for bar, val in zip(bars4, ret_biases):
            if not np.isnan(val):
                height = bar.get_height()
                sign = '+' if val >= 0 else ''
                va = 'bottom' if val >= 0 else 'top'
                offset = 0.02 if val >= 0 else -0.02
                ax4.text(bar.get_x() + bar.get_width()/2., height + offset * (ax4.get_ylim()[1] - ax4.get_ylim()[0]),
                        f'{sign}{val:.1f}%', ha='center', va=va, fontsize=9, fontweight='bold')

        fig.suptitle(f'Method Performance Summary (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=0.995)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "10_method_summary.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "10_method_summary.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_accuracy_vs_characteristic_combined(self, char_col: str, char_label: str,
                                                  metric_name: str, metric_label: str, fig_prefix: str):
        """Plot accuracy metric vs network characteristic (all methods combined)"""
        if self.metrics is None:
            print(f"  WARNING: No metrics data, skipping {metric_label}")
            return

        if char_col not in self.network_stats.columns:
            print(f"  WARNING: Column {char_col} not found, skipping")
            return

        fig, ax = plt.subplots(figsize=(12, 7))

        # Merge metrics with network stats
        metrics_with_stats = self.metrics[self.metrics['metric'] == metric_name].merge(
            self.network_stats[['network', char_col]],
            on='network', how='left'
        )
        metrics_with_stats = metrics_with_stats.dropna(subset=[char_col, 'mean'])

        # Plot each method
        for method in sorted(metrics_with_stats['method'].unique()):
            method_data = metrics_with_stats[metrics_with_stats['method'] == method]

            # Calculate mean and std error per characteristic value
            grouped = method_data.groupby(char_col).agg({
                'mean': ['mean', 'std', 'count']
            }).reset_index()
            grouped.columns = [char_col, 'metric_mean', 'metric_std', 'n']
            grouped['std_err'] = grouped['metric_std'] / np.sqrt(grouped['n'])

            if len(grouped) > 0:
                ax.errorbar(grouped[char_col], grouped['metric_mean'],
                           yerr=grouped['std_err'],
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
        ax.set_ylabel(metric_label, fontsize=14, fontweight='bold')
        ax.set_title(f'{metric_label} vs {char_label}\nILS {self.ils_level}',
                    fontsize=15, fontweight='bold', pad=20)
        ax.legend(frameon=True, loc='best', fontsize=12, framealpha=0.9)
        ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.8)

        if metrics_with_stats[char_col].dtype in ['int64', 'int32']:
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        plt.tight_layout()
        fig.savefig(self.plots_dir / f"{fig_prefix}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"{fig_prefix}.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_accuracy_vs_characteristic_faceted(self, char_col: str, char_label: str,
                                                 metric_name: str, metric_label: str, fig_prefix: str):
        """Plot accuracy metric vs characteristic - faceted subplots, one per method"""
        if self.metrics is None:
            return

        if char_col not in self.network_stats.columns:
            return

        metrics_with_stats = self.metrics[self.metrics['metric'] == metric_name].merge(
            self.network_stats[['network', char_col]],
            on='network', how='left'
        )
        metrics_with_stats = metrics_with_stats.dropna(subset=[char_col, 'mean'])

        methods = sorted(metrics_with_stats['method'].unique())
        n_methods = len(methods)

        ncols = min(3, n_methods)
        nrows = (n_methods + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows), squeeze=False)
        axes = axes.flatten()

        for idx, method in enumerate(methods):
            ax = axes[idx]
            method_data = metrics_with_stats[metrics_with_stats['method'] == method]

            grouped = method_data.groupby(char_col).agg({
                'mean': ['mean', 'std', 'count']
            }).reset_index()
            grouped.columns = [char_col, 'metric_mean', 'metric_std', 'n']
            grouped['std_err'] = grouped['metric_std'] / np.sqrt(grouped['n'])

            if len(grouped) > 0:
                ax.errorbar(grouped[char_col], grouped['metric_mean'],
                           yerr=grouped['std_err'],
                           marker=METHOD_MARKERS.get(method, 'o'),
                           color=METHOD_COLORS.get(method, '#000000'),
                           linewidth=2.5,
                           markersize=9,
                           capsize=5,
                           capthick=2,
                           alpha=0.85,
                           markeredgewidth=1.5,
                           markeredgecolor='white')

            ax.set_xlabel(char_label, fontsize=12, fontweight='bold')
            ax.set_ylabel(metric_label, fontsize=12, fontweight='bold')
            ax.set_title(f'{method}', fontsize=13, fontweight='bold', pad=10)
            ax.grid(True, alpha=0.25, linestyle='--')

            if metrics_with_stats[char_col].dtype in ['int64', 'int32']:
                ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        # Hide unused subplots
        for idx in range(n_methods, len(axes)):
            axes[idx].axis('off')

        fig.suptitle(f'{metric_label} vs {char_label} (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()
        fig.savefig(self.plots_individual_dir / f"{fig_prefix}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_individual_dir / f"{fig_prefix}.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_jaccard_vs_characteristic_combined(self, char_col: str, char_label: str,
                                                 jaccard_metric: str, jaccard_label: str, fig_prefix: str):
        """Plot Jaccard similarity (distance metric) vs network characteristic"""
        if self.metrics is None:
            print(f"  WARNING: No metrics data, skipping {jaccard_label}")
            return

        if char_col not in self.network_stats.columns:
            print(f"  WARNING: Column {char_col} not found, skipping")
            return

        # Use the .dist variant which is 1 - Jaccard similarity
        metric_name = f"{jaccard_metric}.dist"

        fig, ax = plt.subplots(figsize=(12, 7))

        metrics_with_stats = self.metrics[self.metrics['metric'] == metric_name].merge(
            self.network_stats[['network', char_col]],
            on='network', how='left'
        )
        metrics_with_stats = metrics_with_stats.dropna(subset=[char_col, 'mean'])

        if len(metrics_with_stats) == 0:
            print(f"  WARNING: No data for {metric_name}, skipping")
            plt.close()
            return

        for method in sorted(metrics_with_stats['method'].unique()):
            method_data = metrics_with_stats[metrics_with_stats['method'] == method]

            grouped = method_data.groupby(char_col).agg({
                'mean': ['mean', 'std', 'count']
            }).reset_index()
            grouped.columns = [char_col, 'metric_mean', 'metric_std', 'n']
            grouped['std_err'] = grouped['metric_std'] / np.sqrt(grouped['n'])

            if len(grouped) > 0:
                # Plot distance directly (don't convert to similarity)
                ax.errorbar(grouped[char_col], grouped['metric_mean'],
                           yerr=grouped['std_err'],
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
        ax.set_ylabel(f'{jaccard_label}\n(0 = perfect match, 1 = no match)', fontsize=14, fontweight='bold')
        ax.set_title(f'{jaccard_label} vs {char_label}\nILS {self.ils_level}',
                    fontsize=15, fontweight='bold', pad=20)
        ax.legend(frameon=True, loc='best', fontsize=12, framealpha=0.9)
        ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.8)
        ax.set_ylim(-0.05, 1.05)

        if metrics_with_stats[char_col].dtype in ['int64', 'int32']:
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        plt.tight_layout()
        fig.savefig(self.plots_dir / f"{fig_prefix}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"{fig_prefix}.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_jaccard_vs_characteristic_faceted(self, char_col: str, char_label: str,
                                                jaccard_metric: str, jaccard_label: str, fig_prefix: str):
        """Plot Jaccard similarity vs characteristic - faceted subplots"""
        if self.metrics is None:
            return

        if char_col not in self.network_stats.columns:
            return

        metric_name = f"{jaccard_metric}.dist"

        metrics_with_stats = self.metrics[self.metrics['metric'] == metric_name].merge(
            self.network_stats[['network', char_col]],
            on='network', how='left'
        )
        metrics_with_stats = metrics_with_stats.dropna(subset=[char_col, 'mean'])

        if len(metrics_with_stats) == 0:
            return

        methods = sorted(metrics_with_stats['method'].unique())
        n_methods = len(methods)

        ncols = min(3, n_methods)
        nrows = (n_methods + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows), squeeze=False)
        axes = axes.flatten()

        for idx, method in enumerate(methods):
            ax = axes[idx]
            method_data = metrics_with_stats[metrics_with_stats['method'] == method]

            grouped = method_data.groupby(char_col).agg({
                'mean': ['mean', 'std', 'count']
            }).reset_index()
            grouped.columns = [char_col, 'metric_mean', 'metric_std', 'n']
            grouped['std_err'] = grouped['metric_std'] / np.sqrt(grouped['n'])

            if len(grouped) > 0:
                # Plot distance directly (don't convert to similarity)
                ax.errorbar(grouped[char_col], grouped['metric_mean'],
                           yerr=grouped['std_err'],
                           marker=METHOD_MARKERS.get(method, 'o'),
                           color=METHOD_COLORS.get(method, '#000000'),
                           linewidth=2.5,
                           markersize=9,
                           capsize=5,
                           capthick=2,
                           alpha=0.85,
                           markeredgewidth=1.5,
                           markeredgecolor='white')

            ax.set_xlabel(char_label, fontsize=12, fontweight='bold')
            ax.set_ylabel(f'{jaccard_label}\n(0 = perfect, 1 = no match)', fontsize=12, fontweight='bold')
            ax.set_title(f'{method}', fontsize=13, fontweight='bold', pad=10)
            ax.grid(True, alpha=0.25, linestyle='--')
            ax.set_ylim(-0.05, 1.05)

            if metrics_with_stats[char_col].dtype in ['int64', 'int32']:
                ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        for idx in range(n_methods, len(axes)):
            axes[idx].axis('off')

        fig.suptitle(f'{jaccard_label} vs {char_label} (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()
        fig.savefig(self.plots_individual_dir / f"{fig_prefix}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_individual_dir / f"{fig_prefix}.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_polyploid_f1_performance(self):
        """Plot F1 score for polyploid identification per method"""
        if self.metrics is None:
            print("  WARNING: No metrics data, skipping Polyploid F1")
            return

        # Calculate F1 scores from TP, FP, FN
        ploidy_metrics = self.metrics[self.metrics['metric'].str.startswith('ploidy_diff.')]

        if len(ploidy_metrics) == 0:
            print("  WARNING: No ploidy metrics found, skipping")
            return

        methods = sorted(ploidy_metrics['method'].unique())
        f1_scores = []
        precisions = []
        recalls = []

        for method in methods:
            method_data = ploidy_metrics[ploidy_metrics['method'] == method]

            tp = method_data[method_data['metric'] == 'ploidy_diff.TP']['mean'].sum()
            fp = method_data[method_data['metric'] == 'ploidy_diff.FP']['mean'].sum()
            fn = method_data[method_data['metric'] == 'ploidy_diff.FN']['mean'].sum()

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

            f1_scores.append(f1)
            precisions.append(precision)
            recalls.append(recall)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        colors = [METHOD_COLORS.get(m, '#000000') for m in methods]

        # F1 scores
        bars1 = ax1.bar(methods, f1_scores, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax1.set_ylabel('F1 Score', fontsize=13, fontweight='bold')
        ax1.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax1.set_title('Polyploid Identification F1 Score', fontsize=14, fontweight='bold', pad=15)
        ax1.set_ylim(0, 1.05)
        ax1.grid(True, alpha=0.25, axis='y', linestyle='--')
        # Add reference line for perfect score only
        ax1.axhline(y=1.0, color='green', linestyle='--', linewidth=1.5, alpha=0.5, label='Perfect (F1=1.0)')
        ax1.legend(fontsize=9)

        for bar, val in zip(bars1, f1_scores):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Precision and Recall
        x = np.arange(len(methods))
        width = 0.35

        bars2 = ax2.bar(x - width/2, precisions, width, label='Precision', color='#0173B2', alpha=0.8, edgecolor='black')
        bars3 = ax2.bar(x + width/2, recalls, width, label='Recall', color='#DE8F05', alpha=0.8, edgecolor='black')

        ax2.set_ylabel('Score', fontsize=13, fontweight='bold')
        ax2.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax2.set_title('Polyploid Identification: Precision vs Recall', fontsize=14, fontweight='bold', pad=15)
        ax2.set_xticks(x)
        ax2.set_xticklabels(methods)
        ax2.set_ylim(0, 1.05)
        ax2.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax2.legend(fontsize=9)
        
        # Add value labels on bars
        for bars, values in [(bars2, precisions), (bars3, recalls)]:
            for bar, val in zip(bars, values):
                if not np.isnan(val):
                    height = bar.get_height()
                    ax2.text(bar.get_x() + bar.get_width()/2., height,
                            f'{val:.2f}', ha='center', va='bottom', fontsize=8, fontweight='bold')

        fig.suptitle(f'Polyploid Identification Performance (ILS {self.ils_level})',
                    fontsize=16, fontweight='bold', y=1.00)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "23_polyploid_f1_performance.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "23_polyploid_f1_performance.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_comprehensive_correlation_heatmap(self):
        """Comprehensive correlation heatmap: all network properties vs all performance metrics"""
        if self.inventory is None or self.metrics is None:
            print("  WARNING: Missing data for correlation heatmap")
            return

        # Prepare data: for each method and network, get all properties and metrics
        correlation_data = []

        for method in self.inventory['method'].unique():
            method_inv = self.inventory[self.inventory['method'] == method]

            for network in method_inv['network'].unique():
                net_inv = method_inv[method_inv['network'] == network]

                # Get network properties
                net_stats = self.network_stats[self.network_stats['network'] == network]
                if len(net_stats) == 0:
                    continue

                row = {
                    'method': method,
                    'network': network,
                    'completion_rate': net_inv['inferred_exists'].sum() / len(net_inv) * 100,
                    'Num_Species': net_stats['Num_Species'].values[0],
                    'H_Strict': net_stats['H_Strict'].values[0],
                    'H_Relaxed': net_stats['H_Relaxed'].values[0],
                    'Num_Polyploids': net_stats['Num_Polyploids'].values[0],
                    'Max_Copies': net_stats['Max_Copies'].values[0],
                    'Total_WGD': net_stats['Total_WGD'].values[0],
                    'Polyploid_Ratio': net_stats['Polyploid_Ratio'].values[0],
                }

                # Get performance metrics
                net_metrics = self.metrics[
                    (self.metrics['method'] == method) &
                    (self.metrics['network'] == network)
                ]

                for _, metric_row in net_metrics.iterrows():
                    metric_name = metric_row['metric']
                    if metric_name in ['edit_distance', 'num_rets_diff']:
                        row[metric_name] = metric_row['mean']

                correlation_data.append(row)

        if len(correlation_data) == 0:
            print("  WARNING: No correlation data available")
            return

        df = pd.DataFrame(correlation_data)

        # Select columns for correlation
        property_cols = ['Num_Species', 'H_Strict', 'H_Relaxed', 'Num_Polyploids',
                        'Max_Copies', 'Total_WGD', 'Polyploid_Ratio']
        metric_cols = ['completion_rate', 'edit_distance', 'num_rets_diff']

        # Keep only available columns
        property_cols = [c for c in property_cols if c in df.columns]
        metric_cols = [c for c in metric_cols if c in df.columns]

        if len(property_cols) == 0 or len(metric_cols) == 0:
            print("  WARNING: Insufficient data for correlation analysis")
            return

        # Calculate correlation matrix
        corr_matrix = df[property_cols + metric_cols].corr()

        # Extract the subset: properties vs metrics
        corr_subset = corr_matrix.loc[property_cols, metric_cols]

        fig, ax = plt.subplots(figsize=(10, 8))

        sns.heatmap(corr_subset, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                   vmin=-1, vmax=1, square=True, linewidths=1, cbar_kws={'label': 'Correlation'},
                   ax=ax, annot_kws={'fontsize': 10, 'fontweight': 'bold'})

        ax.set_xlabel('Performance Metrics', fontsize=13, fontweight='bold')
        ax.set_ylabel('Network Properties', fontsize=13, fontweight='bold')
        ax.set_title(f'Network Properties vs Performance Metrics Correlation (Aggregated Across All Methods)\nILS {self.ils_level}',
                    fontsize=15, fontweight='bold', pad=20)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "31_comprehensive_correlation_heatmap.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "31_comprehensive_correlation_heatmap.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_correlation_heatmap_per_method(self):
        """Create per-method correlation heatmaps showing which network properties affect each method"""
        if self.inventory is None or self.metrics is None:
            print("  WARNING: Missing data for per-method correlation heatmap")
            return

        # Prepare data: for each method and network, get all properties and metrics
        correlation_data = []

        for method in self.inventory['method'].unique():
            method_inv = self.inventory[self.inventory['method'] == method]

            for network in method_inv['network'].unique():
                net_inv = method_inv[method_inv['network'] == network]

                # Get network properties
                net_stats = self.network_stats[self.network_stats['network'] == network]
                if len(net_stats) == 0:
                    continue

                row = {
                    'method': method,
                    'network': network,
                    'completion_rate': net_inv['inferred_exists'].sum() / len(net_inv) * 100,
                    'Num_Species': net_stats['Num_Species'].values[0],
                    'H_Strict': net_stats['H_Strict'].values[0],
                    'H_Relaxed': net_stats['H_Relaxed'].values[0],
                    'Num_Polyploids': net_stats['Num_Polyploids'].values[0],
                    'Max_Copies': net_stats['Max_Copies'].values[0],
                    'Total_WGD': net_stats['Total_WGD'].values[0],
                    'Polyploid_Ratio': net_stats['Polyploid_Ratio'].values[0],
                }

                # Get performance metrics
                net_metrics = self.metrics[
                    (self.metrics['method'] == method) &
                    (self.metrics['network'] == network)
                ]

                for _, metric_row in net_metrics.iterrows():
                    metric_name = metric_row['metric']
                    if metric_name in ['edit_distance_multree', 'edit_distance', 'num_rets_diff']:
                        row[metric_name] = metric_row['mean']

                correlation_data.append(row)

        if len(correlation_data) == 0:
            print("  WARNING: No correlation data available")
            return

        df = pd.DataFrame(correlation_data)

        # Select columns for correlation
        property_cols = ['Num_Species', 'H_Strict', 'H_Relaxed', 'Num_Polyploids',
                        'Max_Copies', 'Total_WGD', 'Polyploid_Ratio']
        metric_cols = ['completion_rate', 'edit_distance_multree', 'edit_distance', 'num_rets_diff']

        # Keep only available columns
        property_cols = [c for c in property_cols if c in df.columns]
        metric_cols = [c for c in metric_cols if c in df.columns]

        if len(property_cols) == 0 or len(metric_cols) == 0:
            print("  WARNING: Insufficient data for correlation analysis")
            return

        methods = sorted(df['method'].unique())
        n_methods = len(methods)

        # Create faceted subplots - one per method
        ncols = min(3, n_methods)
        nrows = (n_methods + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows), squeeze=False)
        axes = axes.flatten()

        for idx, method in enumerate(methods):
            ax = axes[idx]
            method_df = df[df['method'] == method]

            # Calculate correlation matrix for this method
            corr_matrix = method_df[property_cols + metric_cols].corr()
            corr_subset = corr_matrix.loc[property_cols, metric_cols]

            if len(corr_subset) == 0 or corr_subset.isna().all().all():
                ax.text(0.5, 0.5, f'Insufficient data\nfor {method}',
                       ha='center', va='center', fontsize=12, color='gray')
                ax.set_title(method, fontsize=13, fontweight='bold')
                ax.axis('off')
                continue

            # Plot heatmap
            sns.heatmap(corr_subset, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                       vmin=-1, vmax=1, square=True, linewidths=0.5,
                       cbar_kws={'label': 'Correlation'},
                       ax=ax, annot_kws={'fontsize': 8, 'fontweight': 'bold'})

            ax.set_title(method, fontsize=13, fontweight='bold', pad=10)
            if idx % ncols == 0:  # Leftmost column
                ax.set_ylabel('Network Properties', fontsize=11, fontweight='bold')
            if idx >= (nrows - 1) * ncols:  # Bottom row
                ax.set_xlabel('Performance Metrics', fontsize=11, fontweight='bold')
            ax.tick_params(axis='x', labelsize=8, rotation=45)
            ax.tick_params(axis='y', labelsize=8)

        # Hide unused subplots
        for idx in range(n_methods, len(axes)):
            axes[idx].axis('off')

        fig.suptitle(f'Network Properties vs Performance Metrics Correlation (Per Method)\nILS {self.ils_level}',
                    fontsize=16, fontweight='bold', y=0.995)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "32_per_method_correlation_heatmap.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "32_per_method_correlation_heatmap.png", bbox_inches='tight', dpi=300)
        plt.close()

    def generate_summary_tables(self):
        """Generate comprehensive summary tables for publication"""
        if self.inventory is None or self.metrics is None:
            return

        methods = sorted(self.inventory['method'].unique())

        # Table 1: Overall performance summary
        summary_data = []
        for method in methods:
            method_inv = self.inventory[self.inventory['method'] == method]
            total = len(method_inv)
            successful = method_inv['inferred_exists'].sum()
            comp_rate = successful / total * 100 if total > 0 else 0

            # Use MUL-tree edit distance (PRIMARY METRIC)
            method_edit = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'edit_distance_multree')
            ]
            
            # Fallback to network edit distance if needed
            if len(method_edit) == 0:
                method_edit = self.metrics[
                    (self.metrics['method'] == method) &
                    (self.metrics['metric'] == 'edit_distance')
                ]

            if len(method_edit) > 0:
                mean_ed = method_edit['mean'].mean()
                std_ed = method_edit['mean'].std()
                median_ed = method_edit['mean'].median()
            else:
                mean_ed = std_ed = median_ed = np.nan
            
            # Also get RF distance
            method_rf = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'rf_distance')
            ]
            
            if len(method_rf) > 0:
                mean_rf = method_rf['mean'].mean()
                median_rf = method_rf['mean'].median()
            else:
                mean_rf = median_rf = np.nan

            # Reticulation absolute error (MAE)
            method_ret = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'num_rets_diff')
            ]

            if len(method_ret) > 0:
                mean_ret_err = method_ret['mean'].mean()  # Already absolute
                median_ret_err = method_ret['mean'].median()
            else:
                mean_ret_err = median_ret_err = np.nan
            
            # Reticulation bias (signed error)
            method_bias = self.metrics[
                (self.metrics['method'] == method) &
                (self.metrics['metric'] == 'num_rets_bias')
            ]

            if len(method_bias) > 0:
                mean_ret_bias = method_bias['mean'].mean()
                median_ret_bias = method_bias['mean'].median()
            else:
                mean_ret_bias = median_ret_bias = np.nan

            summary_data.append({
                'Method': method,
                'Total_Runs': total,
                'Completed_Runs': successful,
                'Completion_Rate_%': comp_rate,
                'Mean_Edit_Distance_MULtree': mean_ed,
                'Median_Edit_Distance_MULtree': median_ed,
                'Std_Edit_Distance_MULtree': std_ed,
                'Mean_RF_Distance': mean_rf,
                'Median_RF_Distance': median_rf,
                'Mean_Reticulation_MAE': mean_ret_err,
                'Median_Reticulation_MAE': median_ret_err,
                'Mean_Reticulation_Bias': mean_ret_bias,
                'Median_Reticulation_Bias': median_ret_bias
            })

        df = pd.DataFrame(summary_data)
        df.to_csv(self.tables_dir / "01_method_performance_summary.csv", index=False, float_format='%.4f')

        # Table 2: Per-network performance (for supplementary)
        network_data = []
        for network in sorted(self.network_stats['network'].unique()):
            net_stats = self.network_stats[self.network_stats['network'] == network].iloc[0]
            row = {
                'Network': network,
                'H_Strict': net_stats['H_Strict'],
                'H_Relaxed': net_stats['H_Relaxed'],
                'Num_Polyploids': net_stats['Num_Polyploids'],
                'Total_WGD': net_stats['Total_WGD']
            }

            for method in methods:
                # Completion rate
                method_inv = self.inventory[
                    (self.inventory['method'] == method) &
                    (self.inventory['network'] == network)
                ]
                if len(method_inv) > 0:
                    comp_rate = method_inv['inferred_exists'].sum() / len(method_inv) * 100
                    row[f'{method}_CompRate_%'] = comp_rate
                else:
                    row[f'{method}_CompRate_%'] = np.nan

                # Edit distance (MUL-tree, PRIMARY)
                method_edit = self.metrics[
                    (self.metrics['method'] == method) &
                    (self.metrics['network'] == network) &
                    (self.metrics['metric'] == 'edit_distance_multree')
                ]
                
                # Fallback to network edit distance
                if len(method_edit) == 0:
                    method_edit = self.metrics[
                        (self.metrics['method'] == method) &
                        (self.metrics['network'] == network) &
                        (self.metrics['metric'] == 'edit_distance')
                    ]
                
                if len(method_edit) > 0:
                    row[f'{method}_EditDist_MULtree'] = method_edit['mean'].values[0]
                else:
                    row[f'{method}_EditDist_MULtree'] = np.nan
                
                # RF distance
                method_rf = self.metrics[
                    (self.metrics['method'] == method) &
                    (self.metrics['network'] == network) &
                    (self.metrics['metric'] == 'rf_distance')
                ]
                if len(method_rf) > 0:
                    row[f'{method}_RF'] = method_rf['mean'].values[0]
                else:
                    row[f'{method}_RF'] = np.nan

            network_data.append(row)

        df_networks = pd.DataFrame(network_data)
        df_networks.to_csv(self.tables_dir / "02_per_network_performance.csv", index=False, float_format='%.4f')

        print(f"\n  Summary Table (Method Performance):")
        print(df.to_string(index=False))
        print(f"\n  Per-Network Table saved to: {self.tables_dir / '02_per_network_performance.csv'}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate publication-quality analysis figures',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output structure (per configuration):
  simulations/analysis/summary/{config}/
  ├── plots/
  │   ├── 01_combined_*.pdf/png                      # All methods on one plot
  │   ├── 02_combined_*.pdf/png
  │   ├── 03_combined_*.pdf/png
  │   ├── 04_combined_*.pdf/png
  │   ├── 05_folding_completion_comparison.pdf/png
  │   ├── 06_reticulation_accuracy.pdf/png
  │   ├── 07_reticulation_error_boxplot.pdf/png
  │   ├── 08_edit_distance_boxplot.pdf/png
  │   ├── 09_per_network_breakdown.pdf/png
  │   ├── 10_method_summary.pdf/png
  │   └── individual_methods/                        # Faceted plots (all methods in one file)
  │       ├── 01_faceted_h_strict.pdf/png
  │       ├── 02_faceted_h_relaxed.pdf/png
  │       ├── 03_faceted_num_polyploids.pdf/png
  │       └── 04_faceted_total_wgd.pdf/png
  └── tables/
      ├── 01_method_performance_summary.csv           # Main results table
      └── 02_per_network_performance.csv              # Supplementary table
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
