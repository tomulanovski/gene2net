#!/usr/bin/env python3
"""
create_comprehensive_figures.py - Comprehensive Method Evaluation

Generates complete analysis of phylogenetic network inference methods covering:
1. Overall accuracy (edit distance)
2. Success rates (completion)
3. Reticulation inference accuracy
4. Polyploid identification
5. Network complexity effects
6. Method-specific strengths/weaknesses

Usage:
    python create_comprehensive_figures.py \
      --configs conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
      --network-stats ../networks/mul_tree_final_stats.csv
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from typing import List, Dict, Optional
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
plt.rcParams['figure.titlesize'] = 14
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

ILS_COLORS = {
    'low': '#27ae60',
    'medium': '#f39c12',
    'high': '#e74c3c'
}

ILS_ORDER = ['low', 'medium', 'high']


class ComprehensiveAnalyzer:
    """Comprehensive analysis of network inference methods"""

    def __init__(self, config_names: List[str], network_stats_file: str,
                 output_dir: str = "simulations/analysis/comprehensive_figures"):
        self.config_names = config_names
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load network statistics
        self.network_stats = pd.read_csv(network_stats_file)
        self.network_stats['network'] = self.network_stats['Filename'].str.replace('.tre', '')

        # Load data for all configs
        self.data = {}
        for config in config_names:
            self.data[config] = self._load_config_data(config)

    def _load_config_data(self, config: str) -> Dict[str, pd.DataFrame]:
        """Load all data for a configuration"""
        summary_dir = Path(f"simulations/analysis/summary/{config}")
        data = {}

        if (summary_dir / "inventory.csv").exists():
            data['inventory'] = pd.read_csv(summary_dir / "inventory.csv")
        if (summary_dir / "aggregated_metrics.csv").exists():
            data['aggregated'] = pd.read_csv(summary_dir / "aggregated_metrics.csv")
        if (summary_dir / "comparisons_raw.csv").exists():
            data['comparisons'] = pd.read_csv(summary_dir / "comparisons_raw.csv")

        return data

    def _extract_ils_level(self, config: str) -> str:
        """Extract ILS level from config name"""
        if 'low' in config.lower():
            return 'low'
        elif 'medium' in config.lower():
            return 'medium'
        elif 'high' in config.lower():
            return 'high'
        return 'unknown'

    def create_master_dataframe(self):
        """
        Create unified dataframe with all metrics and network characteristics

        Columns: config, ils_level, network, method, replicate,
                 edit_distance, num_rets_diff, ploidy_diff, success,
                 true_H_Strict, true_Num_Auto, true_Total_WGD, true_Num_Polyploids, etc.
        """
        all_rows = []

        for config in self.config_names:
            if 'comparisons' not in self.data[config]:
                continue

            comp = self.data[config]['comparisons'].copy()
            ils_level = self._extract_ils_level(config)

            # Merge with network stats
            comp = comp.merge(
                self.network_stats[['network', 'Num_Species', 'Num_Polyploids', 'Max_Copies',
                                   'H_Strict', 'Num_Autopolyploidization_Events', 'Total_WGD']],
                on='network', how='left', suffixes=('', '_true')
            )

            # Add ILS level
            comp['ils_level'] = ils_level
            comp['config'] = config

            # Mark success
            comp['success'] = comp['status'] == 'SUCCESS'

            all_rows.append(comp)

        if not all_rows:
            return pd.DataFrame()

        master_df = pd.concat(all_rows, ignore_index=True)
        return master_df

    # ========================================================================
    # FIGURE 1: Overall Method Performance (The Big Picture)
    # ========================================================================

    def plot_overall_method_performance(self, master_df):
        """
        3-panel figure showing overall method performance:
        A) Edit distance by method and ILS level
        B) Success rate by method and ILS level
        C) Combined: Success rate vs Mean edit distance (scatter)
        """
        fig = plt.figure(figsize=(18, 5))
        gs = fig.add_gridspec(1, 3, hspace=0.3, wspace=0.3)

        # Panel A: Edit Distance (Accuracy)
        ax1 = fig.add_subplot(gs[0, 0])

        # Prepare data for grouped boxplot
        success_df = master_df[master_df['success']].copy()
        edit_data = success_df[success_df['metric'] == 'edit_distance']

        if not edit_data.empty:
            # Create grouped data
            methods = sorted(edit_data['method'].unique())
            ils_levels = sorted(edit_data['ils_level'].unique(), key=lambda x: ILS_ORDER.index(x) if x in ILS_ORDER else 999)

            x_pos = []
            labels = []
            colors_list = []
            data_list = []

            for i, method in enumerate(methods):
                for j, ils in enumerate(ils_levels):
                    data = edit_data[(edit_data['method'] == method) &
                                    (edit_data['ils_level'] == ils)]['value'].values
                    if len(data) > 0:
                        x_pos.append(i * (len(ils_levels) + 0.5) + j)
                        data_list.append(data)
                        colors_list.append(ILS_COLORS.get(ils, '#888888'))
                        if j == 0:
                            labels.append(method)

            bp = ax1.boxplot(data_list, positions=x_pos, widths=0.6, patch_artist=True,
                            showfliers=False, medianprops=dict(color='black', linewidth=1.5))

            for patch, color in zip(bp['boxes'], colors_list):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            # Set x-axis labels
            ax1.set_xticks([i * (len(ils_levels) + 0.5) + (len(ils_levels)-1)/2 for i in range(len(methods))])
            ax1.set_xticklabels(labels, rotation=45, ha='right')

        ax1.set_ylabel('Edit Distance (lower = better)', fontsize=11, fontweight='bold')
        ax1.set_title('A) Accuracy: Edit Distance by Method', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')

        # Legend for ILS levels
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=ILS_COLORS[ils], alpha=0.7, label=f'ILS {ils}')
                          for ils in ils_levels if ils in ILS_COLORS]
        ax1.legend(handles=legend_elements, loc='upper right', frameon=True)

        # Panel B: Success Rate
        ax2 = fig.add_subplot(gs[0, 1])

        if 'inventory' in self.data[self.config_names[0]]:
            success_rates = []

            for config in self.config_names:
                if 'inventory' not in self.data[config]:
                    continue

                inv = self.data[config]['inventory']
                ils = self._extract_ils_level(config)

                for method in inv['method'].unique():
                    method_inv = inv[inv['method'] == method]
                    success_rate = (method_inv['status'] == 'exists').sum() / len(method_inv) * 100
                    success_rates.append({
                        'method': method,
                        'ils_level': ils,
                        'success_rate': success_rate
                    })

            sr_df = pd.DataFrame(success_rates)

            if not sr_df.empty:
                # Grouped bar chart
                methods = sorted(sr_df['method'].unique())
                ils_levels = sorted(sr_df['ils_level'].unique(),
                                   key=lambda x: ILS_ORDER.index(x) if x in ILS_ORDER else 999)

                x = np.arange(len(methods))
                width = 0.25

                for i, ils in enumerate(ils_levels):
                    ils_data = sr_df[sr_df['ils_level'] == ils]
                    heights = [ils_data[ils_data['method'] == m]['success_rate'].values[0]
                              if len(ils_data[ils_data['method'] == m]) > 0 else 0
                              for m in methods]

                    ax2.bar(x + i*width, heights, width, label=f'ILS {ils}',
                           color=ILS_COLORS.get(ils, '#888888'), alpha=0.7)

                ax2.set_xticks(x + width)
                ax2.set_xticklabels(methods, rotation=45, ha='right')

        ax2.set_ylabel('Success Rate (%)', fontsize=11, fontweight='bold')
        ax2.set_title('B) Robustness: Success Rate by Method', fontsize=12, fontweight='bold')
        ax2.set_ylim(0, 105)
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.legend(frameon=True)

        # Panel C: Success vs Accuracy Trade-off
        ax3 = fig.add_subplot(gs[0, 2])

        if not sr_df.empty and not edit_data.empty:
            # Calculate mean edit distance per method per ILS
            mean_edit = edit_data.groupby(['method', 'ils_level'])['value'].mean().reset_index()
            mean_edit.columns = ['method', 'ils_level', 'mean_edit_distance']

            # Merge
            trade_off = sr_df.merge(mean_edit, on=['method', 'ils_level'], how='inner')

            for ils in ils_levels:
                ils_data = trade_off[trade_off['ils_level'] == ils]
                ax3.scatter(ils_data['success_rate'], ils_data['mean_edit_distance'],
                           s=150, alpha=0.7, color=ILS_COLORS.get(ils, '#888888'),
                           label=f'ILS {ils}', edgecolors='black', linewidths=1)

                # Add method labels
                for _, row in ils_data.iterrows():
                    ax3.annotate(row['method'],
                               (row['success_rate'], row['mean_edit_distance']),
                               fontsize=7, ha='center', va='bottom')

        ax3.set_xlabel('Success Rate (%)', fontsize=11, fontweight='bold')
        ax3.set_ylabel('Mean Edit Distance', fontsize=11, fontweight='bold')
        ax3.set_title('C) Trade-off: Success vs Accuracy', fontsize=12, fontweight='bold')
        ax3.grid(True, alpha=0.3)
        ax3.legend(frameon=True)

        # Ideal corner annotation
        ax3.annotate('Ideal\n(High success,\nLow distance)',
                    xy=(0.95, 0.05), xycoords='axes fraction',
                    fontsize=9, ha='right', va='bottom',
                    bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

        plt.tight_layout()
        fig.savefig(self.output_dir / "fig1_overall_method_performance.pdf", bbox_inches='tight')
        fig.savefig(self.output_dir / "fig1_overall_method_performance.png", bbox_inches='tight', dpi=300)
        print(f"[OK] Created: fig1_overall_method_performance")
        plt.close()

    # ========================================================================
    # FIGURE 2: Network Complexity Effects
    # ========================================================================

    def plot_complexity_effects(self, master_df):
        """
        Show how network complexity affects method performance
        Rows: Success rate, Edit distance
        Cols: H_Strict, Num_Auto, Total_WGD
        """
        complexity_metrics = [
            ('H_Strict', 'Number of Reticulations'),
            ('Num_Autopolyploidization_Events', 'Autopolyploidization Events'),
            ('Total_WGD', 'Total WGD Events')
        ]

        fig, axes = plt.subplots(2, 3, figsize=(18, 10))

        for col_idx, (metric, label) in enumerate(complexity_metrics):
            # Row 1: Success Rate
            ax_success = axes[0, col_idx]

            for config in self.config_names:
                if 'inventory' not in self.data[config]:
                    continue

                inv = self.data[config]['inventory'].merge(
                    self.network_stats[['network', metric]],
                    on='network', how='left'
                )

                ils = self._extract_ils_level(config)

                for method in sorted(inv['method'].unique()):
                    method_inv = inv[inv['method'] == method]

                    # Group by complexity value
                    grouped = method_inv.groupby(metric).apply(
                        lambda x: (x['status'] == 'exists').sum() / len(x) * 100
                    ).reset_index()
                    grouped.columns = [metric, 'success_rate']

                    if not grouped.empty:
                        ax_success.plot(grouped[metric], grouped['success_rate'],
                                      'o-', label=f'{method} ({ils})',
                                      color=METHOD_COLORS.get(method, '#000000'),
                                      alpha=0.6, markersize=5, linewidth=1.5)

            ax_success.set_xlabel(label, fontsize=10, fontweight='bold')
            if col_idx == 0:
                ax_success.set_ylabel('Success Rate (%)', fontsize=10, fontweight='bold')
            ax_success.set_title(f'{label}', fontsize=11, fontweight='bold')
            ax_success.grid(True, alpha=0.3)
            ax_success.set_ylim(-5, 105)

            if col_idx == 2:
                ax_success.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)

            # Row 2: Edit Distance (for successful runs)
            ax_edit = axes[1, col_idx]

            success_df = master_df[master_df['success']].copy()
            edit_data = success_df[success_df['metric'] == 'edit_distance']

            if not edit_data.empty:
                for method in sorted(edit_data['method'].unique()):
                    method_data = edit_data[edit_data['method'] == method]

                    # Group by complexity
                    grouped = method_data.groupby(metric)['value'].agg(['mean', 'std', 'count']).reset_index()

                    if not grouped.empty:
                        ax_edit.plot(grouped[metric], grouped['mean'],
                                   'o-', label=method,
                                   color=METHOD_COLORS.get(method, '#000000'),
                                   alpha=0.6, markersize=5, linewidth=1.5)

                        # Add error bars (std)
                        ax_edit.fill_between(grouped[metric],
                                           grouped['mean'] - grouped['std'],
                                           grouped['mean'] + grouped['std'],
                                           color=METHOD_COLORS.get(method, '#000000'),
                                           alpha=0.15)

            ax_edit.set_xlabel(label, fontsize=10, fontweight='bold')
            if col_idx == 0:
                ax_edit.set_ylabel('Mean Edit Distance', fontsize=10, fontweight='bold')
            ax_edit.grid(True, alpha=0.3)

            if col_idx == 2:
                ax_edit.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

        fig.suptitle('Network Complexity Effects on Method Performance',
                    fontsize=14, fontweight='bold', y=0.995)

        plt.tight_layout()
        fig.savefig(self.output_dir / "fig2_complexity_effects.pdf", bbox_inches='tight')
        fig.savefig(self.output_dir / "fig2_complexity_effects.png", bbox_inches='tight', dpi=300)
        print(f"[OK] Created: fig2_complexity_effects")
        plt.close()

    # ========================================================================
    # FIGURE 3: Method × Network Performance Heatmap
    # ========================================================================

    def plot_method_network_heatmap(self, master_df):
        """
        Heatmap showing edit distance for each (network, method) combination
        Shows which networks are universally hard and which methods excel where
        """
        for config in self.config_names:
            ils = self._extract_ils_level(config)

            # Filter to this config
            config_data = master_df[master_df['config'] == config]
            edit_data = config_data[config_data['metric'] == 'edit_distance']

            if edit_data.empty:
                continue

            # Pivot table: networks × methods
            pivot = edit_data.pivot_table(
                index='network',
                columns='method',
                values='value',
                aggfunc='mean'
            )

            if pivot.empty:
                continue

            # Calculate mean difficulty per network
            pivot['mean_difficulty'] = pivot.mean(axis=1)
            pivot = pivot.sort_values('mean_difficulty')
            pivot = pivot.drop('mean_difficulty', axis=1)

            # Plot
            fig, ax = plt.subplots(figsize=(10, 12))

            sns.heatmap(pivot, cmap='RdYlGn_r', center=0.5, vmin=0, vmax=1.5,
                       annot=False, fmt='.2f', linewidths=0.5,
                       cbar_kws={'label': 'Edit Distance'},
                       ax=ax)

            ax.set_title(f'Method × Network Performance Heatmap - ILS {ils.title()}',
                        fontsize=13, fontweight='bold', pad=15)
            ax.set_xlabel('Method', fontsize=11, fontweight='bold')
            ax.set_ylabel('Network (sorted by difficulty)', fontsize=11, fontweight='bold')

            plt.tight_layout()
            fig.savefig(self.output_dir / f"fig3_heatmap_{config}.pdf", bbox_inches='tight')
            fig.savefig(self.output_dir / f"fig3_heatmap_{config}.png", bbox_inches='tight', dpi=300)
            print(f"[OK] Created: fig3_heatmap_{config}")
            plt.close()

    # ========================================================================
    # FIGURE 4: Reticulation Inference Accuracy
    # ========================================================================

    def plot_reticulation_accuracy(self, master_df):
        """
        How accurately do methods infer the number of reticulations?
        Shows |inferred_H - true_H| distribution per method
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # Filter to num_rets_diff metric
        ret_data = master_df[master_df['metric'] == 'num_rets_diff']

        for idx, config in enumerate(self.config_names):
            ax = axes[idx]
            config_data = ret_data[ret_data['config'] == config]

            if config_data.empty:
                continue

            methods = sorted(config_data['method'].unique())
            data_for_box = []
            colors_for_box = []

            for method in methods:
                method_data = config_data[config_data['method'] == method]['value'].values
                data_for_box.append(method_data)
                colors_for_box.append(METHOD_COLORS.get(method, '#cccccc'))

            bp = ax.boxplot(data_for_box, labels=methods, patch_artist=True,
                           showfliers=True, widths=0.6)

            for patch, color in zip(bp['boxes'], colors_for_box):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            ils = self._extract_ils_level(config)
            ax.set_title(f'ILS {ils.title()}', fontsize=12, fontweight='bold')
            ax.set_ylabel('|Inferred H - True H|' if idx == 0 else '', fontsize=11, fontweight='bold')
            ax.tick_params(axis='x', rotation=45)
            ax.grid(True, alpha=0.3, axis='y')
            ax.axhline(y=0, color='green', linestyle='--', linewidth=1.5, alpha=0.5,
                      label='Perfect inference')
            if idx == 0:
                ax.legend()

        fig.suptitle('Reticulation Inference Accuracy (Lower = Better)',
                    fontsize=14, fontweight='bold')

        plt.tight_layout()
        fig.savefig(self.output_dir / "fig4_reticulation_accuracy.pdf", bbox_inches='tight')
        fig.savefig(self.output_dir / "fig4_reticulation_accuracy.png", bbox_inches='tight', dpi=300)
        print(f"[OK] Created: fig4_reticulation_accuracy")
        plt.close()

    # ========================================================================
    # FIGURE 5: What Makes Networks Hard?
    # ========================================================================

    def plot_difficulty_correlations(self, master_df):
        """
        Correlation heatmap: Network characteristics vs method performance
        Identifies what network properties make inference difficult
        """
        edit_data = master_df[master_df['metric'] == 'edit_distance']

        # Merge with network stats
        edit_with_stats = edit_data.merge(
            self.network_stats[['network', 'Num_Species', 'Num_Polyploids', 'Max_Copies',
                               'H_Strict', 'Num_Autopolyploidization_Events', 'Total_WGD']],
            on='network', how='left'
        )

        network_props = ['Num_Species', 'Num_Polyploids', 'Max_Copies', 'H_Strict',
                        'Num_Autopolyploidization_Events', 'Total_WGD']

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        for idx, config in enumerate(self.config_names):
            ax = axes[idx]
            config_data = edit_with_stats[edit_with_stats['config'] == config]

            if config_data.empty:
                continue

            methods = sorted(config_data['method'].unique())

            # Build correlation matrix
            corr_matrix = pd.DataFrame(index=methods, columns=network_props)

            for method in methods:
                method_data = config_data[config_data['method'] == method]
                for prop in network_props:
                    valid = method_data[[prop, 'value']].dropna()
                    if len(valid) > 3:
                        corr, _ = pearsonr(valid[prop], valid['value'])
                        corr_matrix.loc[method, prop] = corr
                    else:
                        corr_matrix.loc[method, prop] = np.nan

            corr_matrix = corr_matrix.astype(float)

            # Plot
            sns.heatmap(corr_matrix, annot=True, cmap='RdYlGn_r', center=0,
                       vmin=-1, vmax=1, fmt='.2f', linewidths=0.5,
                       cbar_kws={'label': 'Pearson Correlation'},
                       ax=ax, annot_kws={'size': 8})

            ils = self._extract_ils_level(config)
            ax.set_title(f'ILS {ils.title()}', fontsize=12, fontweight='bold')
            ax.set_xlabel('Network Characteristics', fontsize=10, fontweight='bold')
            ax.set_ylabel('Method' if idx == 0 else '', fontsize=10, fontweight='bold')

        fig.suptitle('Network Characteristics vs Edit Distance\n(Red = harder with more, Green = easier with more)',
                    fontsize=13, fontweight='bold')

        plt.tight_layout()
        fig.savefig(self.output_dir / "fig5_difficulty_correlations.pdf", bbox_inches='tight')
        fig.savefig(self.output_dir / "fig5_difficulty_correlations.png", bbox_inches='tight', dpi=300)
        print(f"[OK] Created: fig5_difficulty_correlations")
        plt.close()

    # ========================================================================
    # SUMMARY TABLE
    # ========================================================================

    def create_summary_table(self, master_df):
        """Create comprehensive summary table"""
        rows = []

        for config in self.config_names:
            ils = self._extract_ils_level(config)

            # Get inventory for success rates
            if 'inventory' in self.data[config]:
                inv = self.data[config]['inventory']

                for method in inv['method'].unique():
                    method_inv = inv[inv['method'] == method]
                    total = len(method_inv)
                    success = (method_inv['status'] == 'exists').sum()
                    success_rate = success / total * 100 if total > 0 else 0

                    # Get edit distance stats
                    method_edit = master_df[
                        (master_df['config'] == config) &
                        (master_df['method'] == method) &
                        (master_df['metric'] == 'edit_distance') &
                        (master_df['success'])
                    ]

                    if not method_edit.empty:
                        mean_edit = method_edit['value'].mean()
                        std_edit = method_edit['value'].std()
                        min_edit = method_edit['value'].min()
                        max_edit = method_edit['value'].max()
                    else:
                        mean_edit = std_edit = min_edit = max_edit = np.nan

                    # Get reticulation accuracy
                    method_ret = master_df[
                        (master_df['config'] == config) &
                        (master_df['method'] == method) &
                        (master_df['metric'] == 'num_rets_diff') &
                        (master_df['success'])
                    ]

                    if not method_ret.empty:
                        mean_ret_error = method_ret['value'].mean()
                    else:
                        mean_ret_error = np.nan

                    rows.append({
                        'ILS_Level': ils,
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
        summary_df.to_csv(self.output_dir / "table1_comprehensive_summary.csv",
                         index=False, float_format='%.4f')

        # LaTeX version
        with open(self.output_dir / "table1_comprehensive_summary.tex", 'w') as f:
            latex = summary_df.to_latex(index=False, float_format='%.2f', escape=False)
            f.write(latex)

        print(f"[OK] Created: table1_comprehensive_summary")

        return summary_df

    # ========================================================================
    # MAIN GENERATION FUNCTION
    # ========================================================================

    def generate_all_figures(self):
        """Generate all comprehensive figures"""
        print(f"\n{'='*80}")
        print("Comprehensive Method Evaluation - Generating All Figures")
        print(f"{'='*80}\n")

        # Create master dataframe
        print("[0/6] Creating master dataframe...")
        master_df = self.create_master_dataframe()

        if master_df.empty:
            print("ERROR: No data found. Make sure summary pipeline has been run.")
            return

        print(f"      Loaded {len(master_df)} comparison records\n")

        # Generate figures
        print("[1/6] Creating overall method performance comparison...")
        self.plot_overall_method_performance(master_df)

        print("[2/6] Creating network complexity effects analysis...")
        self.plot_complexity_effects(master_df)

        print("[3/6] Creating method × network heatmaps...")
        self.plot_method_network_heatmap(master_df)

        print("[4/6] Creating reticulation inference accuracy analysis...")
        self.plot_reticulation_accuracy(master_df)

        print("[5/6] Creating network difficulty correlation analysis...")
        self.plot_difficulty_correlations(master_df)

        print("[6/6] Creating comprehensive summary table...")
        self.create_summary_table(master_df)

        print(f"\n{'='*80}")
        print(f"All figures saved to: {self.output_dir}")
        print(f"{'='*80}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive phylogenetic network inference method evaluation'
    )

    parser.add_argument('--configs', nargs='+', required=True,
                       help='Configuration names')
    parser.add_argument('--network-stats', required=True,
                       help='Path to network characteristics CSV')
    parser.add_argument('--output', default='simulations/analysis/comprehensive_figures',
                       help='Output directory')

    args = parser.parse_args()

    analyzer = ComprehensiveAnalyzer(
        config_names=args.configs,
        network_stats_file=args.network_stats,
        output_dir=args.output
    )

    analyzer.generate_all_figures()


if __name__ == '__main__':
    main()
