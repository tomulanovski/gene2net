#!/usr/bin/env python3
"""
create_analysis_figures.py - Publication-Quality Analysis Figures for Real Data

Generates comprehensive visualizations for pairwise method comparisons in real data.

Usage:
    python create_analysis_figures.py --comparisons CSV --inventory CSV --output DIR
    python create_analysis_figures.py --comparisons comparisons.csv --inventory inventory.csv --output analysis/plots
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
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
    'padre': '#ECE133',         # Yellow
    'mpallop': '#8B4513',       # Saddle Brown
    'grandma_split': '#CC78BC',  # Muted Purple
    'alloppnet': '#DC143C',     # Crimson/Red
    'paper': '#029E73'          # Teal
}


class RealDataAnalyzer:
    """Analyze and visualize results for real data pairwise comparisons"""

    def __init__(self, comparisons_file: str, inventory_file: str, output_dir: str,
                 comparable_networks: list = None, method_stats_file: str = None):
        """
        Initialize analyzer

        Args:
            comparisons_file: Path to comparisons CSV
            inventory_file: Path to inventory CSV
            output_dir: Output directory for plots
            comparable_networks: List of network names to use for completion rate calculation
            method_stats_file: Path to per-method network stats CSV (reticulation counts etc.)
        """
        self.output_dir = Path(output_dir)
        self.plots_dir = self.output_dir / "plots"
        self.tables_dir = self.output_dir / "tables"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        self.comparable_networks = comparable_networks
        self.method_stats = pd.read_csv(method_stats_file) if method_stats_file else None

        # Load data
        self.comparisons = pd.read_csv(comparisons_file)
        self.inventory = pd.read_csv(inventory_file)
        
        # Exclude 'paper' method from all analyses
        self.inventory = self.inventory[self.inventory['method'] != 'paper'].copy()
        self.comparisons = self.comparisons[
            (self.comparisons['method1'] != 'paper') & 
            (self.comparisons['method2'] != 'paper')
        ].copy()

        # Filter to successful comparisons
        if 'status' in self.comparisons.columns:
            self.valid_comparisons = self.comparisons[self.comparisons['status'] == 'SUCCESS'].copy()
        else:
            self.valid_comparisons = self.comparisons.copy()

        print(f"\nLoaded data:")
        print(f"  Comparisons: {len(self.comparisons)} total, {len(self.valid_comparisons)} successful")
        print(f"  Inventory: {len(self.inventory)} entries")
        print(f"  Networks: {self.inventory['network'].nunique()}")
        print(f"  Methods: {self.inventory['method'].nunique()}")
        if comparable_networks:
            print(f"  Comparable networks (for completion rate): {len(comparable_networks)}")
    
    def _get_metric_label(self, metric: str) -> str:
        """Get human-readable label for a metric"""
        metric_labels = {
            'edit_distance_multree': 'Edit Distance',
            'edit_distance': 'Network Edit Distance',
            # 'rf_distance': 'RF Distance',  # Disabled: RF not well-defined for MUL-trees
            'num_rets_diff': 'Reticulation Count Difference',
            'polyploid_species_jaccard': 'Polyploid Species Distance',
            'ploidy_diff.dist': 'Ploidy Distance',
            'ret_leaf_jaccard.dist': 'Reticulation Leaf Distance',
            'ret_sisters_jaccard.dist': 'Sister-Taxa Distance'
        }
        return metric_labels.get(metric, metric.replace('_', ' ').title())

    def plot_method_availability(self):
        """Bar chart showing completion rate per method"""
        # Use comparable networks for completion rate if specified
        if self.comparable_networks:
            comparable_inventory = self.inventory[
                self.inventory['network'].isin(self.comparable_networks)
            ]
            total_for_completion = len(self.comparable_networks)
            
            availability = comparable_inventory.groupby('method').agg({
                'exists': 'sum'
            }).reset_index()
            availability.columns = ['method', 'available']
            availability['total'] = total_for_completion
            availability['completion_rate'] = (availability['available'] / total_for_completion * 100)
        else:
            availability = self.inventory.groupby('method').agg({
                'exists': ['sum', 'count']
            }).reset_index()
            availability.columns = ['method', 'available', 'total']
            availability['completion_rate'] = (availability['available'] / availability['total'] * 100)

        fig, ax = plt.subplots(figsize=(10, 6))

        colors = [METHOD_COLORS.get(m, '#CCCCCC') for m in availability['method']]
        bars = ax.bar(availability['method'], availability['completion_rate'],
                     color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)

        ax.set_ylabel('Completion Rate (%)', fontsize=13, fontweight='bold')
        ax.set_xlabel('Method', fontsize=13, fontweight='bold')
        if self.comparable_networks:
            title = f'Method Availability Across Networks\n(Based on {len(self.comparable_networks)} Comparable Networks)'
        else:
            title = 'Method Availability Across Networks'
        ax.set_title(title, fontsize=15, fontweight='bold', pad=20)
        ax.set_ylim(0, 105)
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')

        # Add value labels
        for bar, rate in zip(bars, availability['completion_rate']):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{rate:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        fig.savefig(self.plots_dir / "01_method_availability.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "01_method_availability.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_availability_heatmap(self):
        """Heatmap showing method × network availability matrix"""
        # Filter out 'paper' method (already filtered in __init__, but double-check)
        inventory_filtered = self.inventory[self.inventory['method'] != 'paper'].copy()
        
        pivot = inventory_filtered.pivot_table(
            index='network',
            columns='method',
            values='exists',
            aggfunc='first'
        ).astype(float)

        fig, ax = plt.subplots(figsize=(10, max(8, len(pivot) * 0.3)))

        sns.heatmap(pivot, annot=True, fmt='.0f', cmap='RdYlGn', vmin=0, vmax=1,
                   cbar_kws={'label': 'Available (1) / Missing (0)'},
                   ax=ax, linewidths=0.5, linecolor='gray')

        ax.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax.set_ylabel('Network', fontsize=13, fontweight='bold')
        ax.set_title('Method Availability Matrix', fontsize=15, fontweight='bold', pad=20)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "02_availability_heatmap.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "02_availability_heatmap.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_pairwise_heatmap(self, metric: str = 'edit_distance_multree'):
        """Heatmap showing mean pairwise distances between methods"""
        if self.valid_comparisons.empty:
            print(f"  WARNING: No valid comparisons for {metric}")
            return

        metric_data = self.valid_comparisons[self.valid_comparisons['metric'] == metric]

        if metric_data.empty:
            print(f"  WARNING: No data for metric {metric}")
            return

        # Get all unique methods
        all_methods = sorted(set(metric_data['method1'].unique()) | set(metric_data['method2'].unique()))

        # Create symmetric matrix
        matrix = np.full((len(all_methods), len(all_methods)), np.nan)

        for i, method1 in enumerate(all_methods):
            for j, method2 in enumerate(all_methods):
                if i == j:
                    matrix[i, j] = 0.0  # Distance to self is 0
                else:
                    # Get comparisons for this pair (both directions)
                    pair_data = metric_data[
                        ((metric_data['method1'] == method1) & (metric_data['method2'] == method2)) |
                        ((metric_data['method1'] == method2) & (metric_data['method2'] == method1))
                    ]['value']

                    if len(pair_data) > 0:
                        matrix[i, j] = pair_data.mean()

        # Create DataFrame
        df_matrix = pd.DataFrame(matrix, index=all_methods, columns=all_methods)

        fig, ax = plt.subplots(figsize=(10, 8))

        # Use sequential colormap for distances (darker = more similar, lighter = more different)
        # All metrics are non-negative (distances or absolute differences)
        sns.heatmap(df_matrix, annot=True, fmt='.3f', cmap='YlGnBu_r', vmin=0,
                   cbar_kws={'label': self._get_metric_label(metric)},
                   ax=ax, linewidths=1, linecolor='black', square=True)

        # Get proper metric label
        metric_label = self._get_metric_label(metric)
        
        ax.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax.set_ylabel('Method', fontsize=13, fontweight='bold')
        ax.set_title(f'Method Similarity: {metric_label}\n(Lower = More Similar)',
                    fontsize=15, fontweight='bold', pad=20)

        # Add GRAMPA footnote on Jaccard heatmaps
        if 'jaccard' in metric:
            from compare_reticulations import SINGLE_RETICULATION_METHODS
            if set(all_methods) & SINGLE_RETICULATION_METHODS:
                fig.text(0.01, 0.01, '* GRAMPA: best-match only (1 reticulation)',
                         fontsize=9, fontstyle='italic', color='gray')

        plt.tight_layout()
        metric_safe = metric.replace('.', '_')
        fig.savefig(self.plots_dir / f"03_pairwise_{metric_safe}_heatmap.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"03_pairwise_{metric_safe}_heatmap.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_pairwise_boxplots(self, metric: str = 'edit_distance_multree'):
        """Boxplots showing distribution of pairwise distances for each method pair"""
        if self.valid_comparisons.empty:
            return

        metric_data = self.valid_comparisons[self.valid_comparisons['metric'] == metric]

        if metric_data.empty:
            return

        # Create method pair identifier
        metric_data = metric_data.copy()
        metric_data['pair'] = metric_data.apply(
            lambda row: f"{row['method1']} vs {row['method2']}", axis=1
        )

        pairs = sorted(metric_data['pair'].unique())
        data_by_pair = [metric_data[metric_data['pair'] == pair]['value'].values for pair in pairs]

        fig, ax = plt.subplots(figsize=(max(12, len(pairs) * 0.8), 6))

        bp = ax.boxplot(data_by_pair, labels=pairs, patch_artist=True,
                       widths=0.6, showfliers=True,
                       boxprops=dict(linewidth=1.5),
                       whiskerprops=dict(linewidth=1.5),
                       capprops=dict(linewidth=1.5),
                       medianprops=dict(linewidth=2.5, color='red'))

        for patch in bp['boxes']:
            patch.set_facecolor('#0173B2')
            patch.set_alpha(0.7)

        # Get proper metric label
        metric_label = self._get_metric_label(metric)
        
        ax.set_ylabel(metric_label, fontsize=13, fontweight='bold')
        ax.set_xlabel('Method Pair', fontsize=13, fontweight='bold')
        ax.set_title(f'Method Similarity Distribution: {metric_label}\nby Method Pair',
                    fontsize=15, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')

        # Add GRAMPA footnote on Jaccard boxplots
        if 'jaccard' in metric:
            from compare_reticulations import SINGLE_RETICULATION_METHODS
            all_methods = set(metric_data['method1'].unique()) | set(metric_data['method2'].unique())
            if all_methods & SINGLE_RETICULATION_METHODS:
                fig.text(0.01, 0.01, '* GRAMPA: best-match only (1 reticulation)',
                         fontsize=9, fontstyle='italic', color='gray')

        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        metric_safe = metric.replace('.', '_')
        fig.savefig(self.plots_dir / f"04_pairwise_{metric_safe}_boxplot.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"04_pairwise_{metric_safe}_boxplot.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_per_network_comparisons(self, metric: str = 'edit_distance_multree'):
        """Grouped bar chart showing pairwise distances per network"""
        if self.valid_comparisons.empty:
            return

        metric_data = self.valid_comparisons[self.valid_comparisons['metric'] == metric]

        if metric_data.empty:
            return

        # Create method pair identifier
        metric_data = metric_data.copy()
        metric_data['pair'] = metric_data.apply(
            lambda row: f"{row['method1']} vs {row['method2']}", axis=1
        )

        # Pivot: network × pair
        pivot = metric_data.pivot_table(
            index='network',
            columns='pair',
            values='value',
            aggfunc='first'
        )

        fig, ax = plt.subplots(figsize=(max(14, len(pivot) * 0.5), 8))

        pivot.plot(kind='bar', ax=ax, width=0.8, colormap='Set3')

        # Get proper metric label
        metric_label = self._get_metric_label(metric)
        
        ax.set_ylabel(metric_label, fontsize=13, fontweight='bold')
        ax.set_xlabel('Network', fontsize=13, fontweight='bold')
        ax.set_title(f'Method Similarity by Network: {metric_label}',
                    fontsize=15, fontweight='bold', pad=20)
        ax.legend(title='Method Pair', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')

        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        metric_safe = metric.replace('.', '_')
        fig.savefig(self.plots_dir / f"05_per_network_{metric_safe}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"05_per_network_{metric_safe}.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_method_rankings(self):
        """Bar chart showing method rankings based on average distances"""
        # Compute average distance for each method
        all_methods = set(self.valid_comparisons['method1'].unique()) | set(self.valid_comparisons['method2'].unique())

        rankings = []
        for method in all_methods:
            method_comparisons = self.valid_comparisons[
                (self.valid_comparisons['method1'] == method) | (self.valid_comparisons['method2'] == method)
            ]

            # Focus on edit_distance_multree
            ed_data = method_comparisons[method_comparisons['metric'] == 'edit_distance_multree']['value']
            avg_ed = ed_data.mean() if len(ed_data) > 0 else np.nan

            rankings.append({
                'method': method,
                'avg_edit_distance_multree': avg_ed
            })

        rankings_df = pd.DataFrame(rankings).sort_values('avg_edit_distance_multree', ascending=True)

        fig, ax = plt.subplots(figsize=(10, 6))

        colors = [METHOD_COLORS.get(m, '#CCCCCC') for m in rankings_df['method']]
        bars = ax.barh(rankings_df['method'], rankings_df['avg_edit_distance_multree'],
                      color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)

        ax.set_xlabel('Average Edit Distance', fontsize=13, fontweight='bold')
        ax.set_ylabel('Method', fontsize=13, fontweight='bold')
        ax.set_title('Method Similarity Rankings\nAverage Distance to Other Methods (Lower = More Similar)', 
                    fontsize=15, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.25, axis='x', linestyle='--')

        # Add value labels
        for bar, val in zip(bars, rankings_df['avg_edit_distance_multree']):
            if not np.isnan(val):
                width = bar.get_width()
                ax.text(width, bar.get_y() + bar.get_height()/2.,
                       f'{val:.3f}', ha='left', va='center', fontsize=10, fontweight='bold')

        plt.tight_layout()
        fig.savefig(self.plots_dir / "06_method_rankings.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "06_method_rankings.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_distance_metrics_comparison(self):
        """Compare Edit Distance metrics side-by-side"""
        if self.valid_comparisons.empty:
            return

        metrics_to_compare = {
            'edit_distance_multree': 'Edit Distance',
            # 'rf_distance': 'RF Distance',  # Disabled: RF not well-defined for MUL-trees
        }

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        for idx, (metric_name, metric_label) in enumerate(metrics_to_compare.items()):
            ax = axes[idx]

            metric_data = self.valid_comparisons[self.valid_comparisons['metric'] == metric_name]

            if len(metric_data) == 0:
                ax.text(0.5, 0.5, f'No data for\n{metric_label}',
                       ha='center', va='center', fontsize=14, color='gray')
                ax.set_title(metric_label, fontsize=13, fontweight='bold')
                ax.axis('off')
                continue

            # Create method pair identifier
            metric_data = metric_data.copy()
            metric_data['pair'] = metric_data.apply(
                lambda row: f"{row['method1']} vs {row['method2']}", axis=1
            )

            pairs = sorted(metric_data['pair'].unique())
            data_by_pair = [metric_data[metric_data['pair'] == pair]['value'].values for pair in pairs]

            if len(data_by_pair) == 0:
                ax.text(0.5, 0.5, f'No data for\n{metric_label}',
                       ha='center', va='center', fontsize=14, color='gray')
                ax.set_title(metric_label, fontsize=13, fontweight='bold')
                ax.axis('off')
                continue

            # Create boxplot
            bp = ax.boxplot(data_by_pair, labels=pairs, patch_artist=True,
                           widths=0.5, showfliers=True,
                           boxprops=dict(linewidth=1.5),
                           whiskerprops=dict(linewidth=1.5),
                           capprops=dict(linewidth=1.5),
                           medianprops=dict(linewidth=2.5, color='red'))

            for patch in bp['boxes']:
                patch.set_facecolor(METHOD_COLORS.get('grampa', '#0173B2'))
                patch.set_alpha(0.7)

            # Add mean values as text
            for i, pair in enumerate(pairs):
                pair_data = metric_data[metric_data['pair'] == pair]['value']
                if len(pair_data) > 0:
                    mean_val = pair_data.mean()
                    ax.text(i + 1, ax.get_ylim()[1] * 0.95, f'{mean_val:.3f}',
                           ha='center', va='top', fontsize=9, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                    edgecolor='gray', alpha=0.8))

            ax.set_ylabel('Distance\n(0 = identical, 1 = very different)',
                         fontsize=12, fontweight='bold')
            ax.set_xlabel('Method Pair', fontsize=12, fontweight='bold')
            ax.set_title(metric_label, fontsize=13, fontweight='bold', pad=15)
            ax.grid(True, alpha=0.25, axis='y', linestyle='--')
            ax.tick_params(axis='x', rotation=45, labelsize=9)

        fig.suptitle('Method Similarity: Distance Metrics Comparison',
                    fontsize=16, fontweight='bold', y=1.02)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "07_distance_metrics_comparison.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "07_distance_metrics_comparison.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_polyploid_agreement(self):
        """Heatmap showing mean polyploid species Jaccard distance between methods"""
        if self.valid_comparisons.empty:
            return

        ploidy_data = self.valid_comparisons[self.valid_comparisons['metric'] == 'polyploid_species_jaccard']

        if ploidy_data.empty:
            print("  WARNING: No polyploid identification data available")
            return

        # Build symmetric method×method mean distance matrix
        all_methods = sorted(set(ploidy_data['method1'].unique()) | set(ploidy_data['method2'].unique()))
        matrix = np.full((len(all_methods), len(all_methods)), np.nan)

        for i, m1 in enumerate(all_methods):
            for j, m2 in enumerate(all_methods):
                if i == j:
                    matrix[i, j] = 0.0
                else:
                    pair_vals = ploidy_data[
                        ((ploidy_data['method1'] == m1) & (ploidy_data['method2'] == m2)) |
                        ((ploidy_data['method1'] == m2) & (ploidy_data['method2'] == m1))
                    ]['value']
                    if len(pair_vals) > 0:
                        matrix[i, j] = pair_vals.mean()

        df_matrix = pd.DataFrame(matrix, index=all_methods, columns=all_methods)

        fig, ax = plt.subplots(figsize=(10, 8))

        sns.heatmap(df_matrix, annot=True, fmt='.3f', cmap='YlOrRd', vmin=0, vmax=1,
                   cbar_kws={'label': 'Polyploid Species Distance\n(0 = same species, 1 = different)'},
                   ax=ax, linewidths=1, linecolor='black', square=True,
                   annot_kws={'fontsize': 10, 'fontweight': 'bold'})

        ax.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax.set_ylabel('Method', fontsize=13, fontweight='bold')
        ax.set_title('Polyploid Species Distance\n(Lower = More Similar)',
                    fontsize=15, fontweight='bold', pad=20)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "08_polyploid_identification_agreement.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "08_polyploid_identification_agreement.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_dataset_agreement_ranking(self):
        """Rank datasets by method agreement: one panel per metric + one for average"""
        if self.valid_comparisons.empty:
            return

        agreement_metrics = ['polyploid_species_jaccard', 'ret_leaf_jaccard.dist']
        available_metrics = [m for m in agreement_metrics
                            if m in self.valid_comparisons['metric'].values]

        if not available_metrics:
            print("  WARNING: No agreement metrics available for dataset ranking")
            return

        from matplotlib.patches import Patch

        n_panels = len(available_metrics) + 1  # +1 for average
        fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, max(6, 12 * 0.4)),
                                sharey=True)
        if n_panels == 1:
            axes = [axes]

        # Sort datasets by average across all metrics (used for consistent y-axis order)
        all_metric_data = self.valid_comparisons[self.valid_comparisons['metric'].isin(available_metrics)]
        overall_means = all_metric_data.groupby('network')['value'].mean().sort_values()
        dataset_order = overall_means.index.tolist()

        metric_labels = {
            'polyploid_species_jaccard': 'Polyploid Species\nDist.',
            'ret_leaf_jaccard.dist': 'Reticulation Leaf\nDistance',
        }

        def color_for_val(v):
            if v < 0.3:
                return '#2E8B57'
            elif v < 0.6:
                return '#DE8F05'
            else:
                return '#DC143C'

        # Individual metric panels
        for idx, metric in enumerate(available_metrics):
            ax = axes[idx]
            metric_data = self.valid_comparisons[self.valid_comparisons['metric'] == metric]
            network_means = metric_data.groupby('network')['value'].mean()

            # Reindex to consistent order
            network_means = network_means.reindex(dataset_order)

            colors = [color_for_val(v) if not np.isnan(v) else '#CCCCCC'
                      for v in network_means.values]
            bars = ax.barh(range(len(dataset_order)), network_means.values,
                          color=colors, alpha=0.8, edgecolor='black', linewidth=0.8)

            for i, val in enumerate(network_means.values):
                if not np.isnan(val):
                    ax.text(val + 0.01, i, f'{val:.2f}', ha='left', va='center', fontsize=9)

            ax.set_yticks(range(len(dataset_order)))
            if idx == 0:
                ax.set_yticklabels(dataset_order)
            ax.set_xlabel('Mean Pairwise Distance', fontsize=11)
            ax.set_title(metric_labels.get(metric, metric), fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.25, axis='x', linestyle='--')
            ax.set_xlim(0, 1.15)

        # Average panel (last)
        ax = axes[-1]
        colors = [color_for_val(v) for v in overall_means.values]
        bars = ax.barh(range(len(dataset_order)), overall_means.values,
                      color=colors, alpha=0.8, edgecolor='black', linewidth=0.8)

        for i, val in enumerate(overall_means.values):
            ax.text(val + 0.01, i, f'{val:.2f}', ha='left', va='center', fontsize=9)

        ax.set_xlabel('Mean Pairwise Distance', fontsize=11)
        ax.set_title('Average\n(All Metrics)', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.25, axis='x', linestyle='--')
        ax.set_xlim(0, 1.15)

        # Legend
        legend_elements = [
            Patch(facecolor='#2E8B57', edgecolor='black', label='High agreement (<0.3)'),
            Patch(facecolor='#DE8F05', edgecolor='black', label='Moderate (0.3-0.6)'),
            Patch(facecolor='#DC143C', edgecolor='black', label='Low agreement (>0.6)')
        ]
        axes[-1].legend(handles=legend_elements, loc='lower right', fontsize=9)

        fig.suptitle('Dataset Ranking by Method Agreement\n(Lower = Methods Agree More)',
                    fontsize=15, fontweight='bold', y=1.02)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "09_dataset_agreement_ranking.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "09_dataset_agreement_ranking.png", bbox_inches='tight', dpi=300)
        plt.close()

    def plot_reticulation_counts_per_dataset(self):
        """Show actual inferred reticulation count per method per dataset"""
        if self.method_stats is None or self.method_stats.empty:
            print("  WARNING: No method_network_stats data for reticulation count plot")
            return

        df = self.method_stats.copy()

        # Pivot for grouped bar chart
        pivot = df.pivot_table(index='network', columns='method', values='reticulation_count')

        if pivot.empty:
            return

        fig, ax = plt.subplots(figsize=(max(14, len(pivot) * 0.8), 7))

        # Use method colors in consistent order
        method_order = [m for m in METHOD_COLORS if m in pivot.columns]
        remaining = [m for m in pivot.columns if m not in method_order]
        method_order.extend(remaining)
        pivot = pivot[[m for m in method_order if m in pivot.columns]]

        colors = [METHOD_COLORS.get(m, '#CCCCCC') for m in pivot.columns]
        pivot.plot(kind='bar', ax=ax, color=colors, width=0.8, edgecolor='black', linewidth=0.5)

        ax.set_ylabel('Inferred Reticulations', fontsize=13, fontweight='bold')
        ax.set_xlabel('Dataset', fontsize=13, fontweight='bold')
        ax.set_title('Inferred Reticulation Count by Method and Dataset',
                    fontsize=15, fontweight='bold', pad=20)
        ax.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.xticks(rotation=45, ha='right')

        plt.tight_layout()
        fig.savefig(self.plots_dir / "10_reticulation_counts_per_dataset.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "10_reticulation_counts_per_dataset.png", bbox_inches='tight', dpi=300)
        plt.close()

    def generate_all_figures(self):
        """Generate all analysis figures"""
        print(f"\n{'='*80}")
        print(f"Generating Analysis Figures")
        print(f"Output: {self.plots_dir}")
        print(f"{'='*80}\n")

        # Method availability
        print("Plotting method availability...")
        self.plot_method_availability()
        self.plot_availability_heatmap()

        # Pairwise comparisons
        print("Plotting pairwise comparisons...")
        metrics_to_plot = [
            'edit_distance_multree', 'num_rets_diff',
            'polyploid_species_jaccard',
            'ret_leaf_jaccard.dist', 'ret_sisters_jaccard.dist',
            'ploidy_diff.dist'
            # 'rf_distance' disabled: RF not well-defined for MUL-trees
        ]
        for metric in metrics_to_plot:
            print(f"  {metric}...")
            self.plot_pairwise_heatmap(metric)
            self.plot_pairwise_boxplots(metric)

        # Distance comparison
        print("Plotting distance metrics comparison...")
        self.plot_distance_metrics_comparison()

        # Per-network analysis
        print("Plotting per-network comparisons...")
        self.plot_per_network_comparisons('edit_distance_multree')

        # Method rankings
        print("Plotting method rankings...")
        self.plot_method_rankings()

        # Polyploid identification agreement
        print("Plotting polyploid identification agreement...")
        self.plot_polyploid_agreement()

        # Dataset agreement ranking
        print("Plotting dataset agreement ranking...")
        self.plot_dataset_agreement_ranking()

        # Reticulation counts per dataset
        print("Plotting reticulation counts per dataset...")
        self.plot_reticulation_counts_per_dataset()

        print(f"\n{'='*80}")
        print(f"Figure generation complete!")
        print(f"  Plots saved to: {self.plots_dir}")
        print(f"{'='*80}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate publication-quality analysis figures for real data',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--comparisons', required=True,
                       help='Path to comparisons CSV file')
    parser.add_argument('--inventory', required=True,
                       help='Path to inventory CSV file')
    parser.add_argument('--output', required=True,
                       help='Output directory for plots')

    args = parser.parse_args()

    analyzer = RealDataAnalyzer(
        comparisons_file=args.comparisons,
        inventory_file=args.inventory,
        output_dir=args.output
    )
    analyzer.generate_all_figures()


if __name__ == '__main__':
    main()

