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
    'alloppnet': '#DC143C',     # Crimson/Red
    'paper': '#029E73'          # Teal
}


class RealDataAnalyzer:
    """Analyze and visualize results for real data pairwise comparisons"""

    def __init__(self, comparisons_file: str, inventory_file: str, output_dir: str,
                 comparable_networks: list = None):
        """
        Initialize analyzer

        Args:
            comparisons_file: Path to comparisons CSV
            inventory_file: Path to inventory CSV
            output_dir: Output directory for plots
            comparable_networks: List of network names to use for completion rate calculation
        """
        self.output_dir = Path(output_dir)
        self.plots_dir = self.output_dir / "plots"
        self.tables_dir = self.output_dir / "tables"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        self.comparable_networks = comparable_networks

        # Load data
        self.comparisons = pd.read_csv(comparisons_file)
        self.inventory = pd.read_csv(inventory_file)

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
        pivot = self.inventory.pivot_table(
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

        sns.heatmap(df_matrix, annot=True, fmt='.3f', cmap='YlOrRd', vmin=0,
                   cbar_kws={'label': metric.replace('_', ' ').title()},
                   ax=ax, linewidths=1, linecolor='black', square=True)

        ax.set_xlabel('Method', fontsize=13, fontweight='bold')
        ax.set_ylabel('Method', fontsize=13, fontweight='bold')
        ax.set_title(f'Pairwise {metric.replace("_", " ").title()} Between Methods', 
                    fontsize=15, fontweight='bold', pad=20)

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

        ax.set_ylabel(metric.replace('_', ' ').title(), fontsize=13, fontweight='bold')
        ax.set_xlabel('Method Pair', fontsize=13, fontweight='bold')
        ax.set_title(f'Distribution of {metric.replace("_", " ").title()} by Method Pair',
                    fontsize=15, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.25, axis='y', linestyle='--')

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

        ax.set_ylabel(metric.replace('_', ' ').title(), fontsize=13, fontweight='bold')
        ax.set_xlabel('Network', fontsize=13, fontweight='bold')
        ax.set_title(f'{metric.replace("_", " ").title()} by Network and Method Pair',
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

        ax.set_xlabel('Average Edit Distance (MUL-tree)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Method', fontsize=13, fontweight='bold')
        ax.set_title('Method Rankings: Average Distance to Other Methods\n(Lower = More Similar)', 
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
        metrics_to_plot = ['edit_distance_multree', 'rf_distance', 'num_rets_diff']
        for metric in metrics_to_plot:
            print(f"  {metric}...")
            self.plot_pairwise_heatmap(metric)
            self.plot_pairwise_boxplots(metric)

        # Per-network analysis
        print("Plotting per-network comparisons...")
        self.plot_per_network_comparisons('edit_distance_multree')

        # Method rankings
        print("Plotting method rankings...")
        self.plot_method_rankings()

        print(f"\n{'='*80}")
        print(f"✓ Figure generation complete!")
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

