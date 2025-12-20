#!/usr/bin/env python3
"""
create_publication_figures.py - Generate publication-quality figures

Creates comprehensive visualizations for phylogenetic network inference paper:
1. Success rate vs network characteristics (reticulations, polyploids, WGD, autopolyploids)
2. Method performance comparisons (edit distance, rankings)
3. Correlation analysis (network properties vs method performance)
4. Multi-ILS level comparisons

Usage:
    python create_publication_figures.py --configs conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M
    python create_publication_figures.py --configs conf_ils_low_10M --network-stats networks/mul_tree_final_stats.csv
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

# Publication-quality plot settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.titlesize'] = 13
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

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
    'low': '#2ecc71',
    'medium': '#f39c12',
    'high': '#e74c3c'
}


class PublicationFigureGenerator:
    """Generate publication-quality figures for phylogenetic network inference analysis"""

    def __init__(self, config_names: List[str], network_stats_file: Optional[str] = None,
                 output_dir: str = "simulations/analysis/publication_figures"):
        """
        Initialize figure generator

        Args:
            config_names: List of configuration names (e.g., ['conf_ils_low_10M'])
            network_stats_file: Path to network characteristics CSV
            output_dir: Directory to save figures
        """
        self.config_names = config_names
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load network statistics
        if network_stats_file and Path(network_stats_file).exists():
            self.network_stats = pd.read_csv(network_stats_file)
            # Clean filename to match network names
            self.network_stats['network'] = self.network_stats['Filename'].str.replace('.tre', '')
        else:
            self.network_stats = None
            print("Warning: Network stats file not provided or not found")

        # Load summary data for all configs
        self.data = {}
        for config in config_names:
            self.data[config] = self._load_config_data(config)

    def _load_config_data(self, config: str) -> Dict[str, pd.DataFrame]:
        """Load all summary data for a configuration"""
        summary_dir = Path(f"simulations/analysis/summary/{config}")

        if not summary_dir.exists():
            print(f"Warning: Summary directory not found: {summary_dir}")
            return {}

        data = {}

        # Load inventory (has completion status)
        inventory_file = summary_dir / "inventory.csv"
        if inventory_file.exists():
            data['inventory'] = pd.read_csv(inventory_file)

        # Load aggregated metrics
        agg_file = summary_dir / "aggregated_metrics.csv"
        if agg_file.exists():
            data['aggregated'] = pd.read_csv(agg_file)

        # Load comparisons (raw)
        comp_file = summary_dir / "comparisons_raw.csv"
        if comp_file.exists():
            data['comparisons'] = pd.read_csv(comp_file)

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

    def plot_success_vs_reticulations(self):
        """
        Plot number of successful runs vs number of reticulations

        Creates separate panels for each ILS level
        """
        fig, axes = plt.subplots(1, len(self.config_names), figsize=(5*len(self.config_names), 4),
                                sharey=True)

        if len(self.config_names) == 1:
            axes = [axes]

        for ax, config in zip(axes, self.config_names):
            if 'inventory' not in self.data[config]:
                continue

            inv = self.data[config]['inventory']

            # Merge with network stats if available
            if self.network_stats is not None:
                inv = inv.merge(self.network_stats[['network', 'H_Strict',
                                                    'Num_Autopolyploidization_Events', 'Total_WGD']],
                               on='network', how='left')

                # Group by reticulations and method
                success_by_ret = inv[inv['status'] == 'exists'].groupby(
                    ['H_Strict', 'method']
                ).size().reset_index(name='num_success')

                # Total attempts per (reticulations, method)
                total_by_ret = inv.groupby(['H_Strict', 'method']).size().reset_index(name='num_total')

                # Merge and calculate success rate
                success_rate = success_by_ret.merge(total_by_ret, on=['H_Strict', 'method'])
                success_rate['success_rate'] = success_rate['num_success'] / success_rate['num_total'] * 100

                # Plot for each method
                for method in success_rate['method'].unique():
                    method_data = success_rate[success_rate['method'] == method]
                    ax.plot(method_data['H_Strict'], method_data['success_rate'],
                           'o-', label=method, color=METHOD_COLORS.get(method, '#000000'),
                           markersize=6, linewidth=2, alpha=0.7)

            ils_level = self._extract_ils_level(config)
            ax.set_xlabel('Number of Reticulations (H_Strict)', fontsize=11)
            ax.set_title(f'ILS {ils_level.title()}', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_ylim(-5, 105)

        axes[0].set_ylabel('Success Rate (%)', fontsize=11)
        axes[-1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)

        plt.tight_layout()
        fig.savefig(self.output_dir / "fig1_success_vs_reticulations.pdf", bbox_inches='tight')
        fig.savefig(self.output_dir / "fig1_success_vs_reticulations.png", bbox_inches='tight')
        print(f"[OK] Created: fig1_success_vs_reticulations")
        plt.close()

    def plot_success_vs_polyploids(self):
        """Plot number of successful runs vs number of polyploid species"""
        fig, axes = plt.subplots(1, len(self.config_names), figsize=(5*len(self.config_names), 4),
                                sharey=True)

        if len(self.config_names) == 1:
            axes = [axes]

        for ax, config in zip(axes, self.config_names):
            if 'inventory' not in self.data[config]:
                continue

            inv = self.data[config]['inventory']

            if self.network_stats is not None:
                inv = inv.merge(self.network_stats[['network', 'Num_Polyploids']],
                               on='network', how='left')

                # Group by polyploids and method
                success_by_poly = inv[inv['status'] == 'exists'].groupby(
                    ['Num_Polyploids', 'method']
                ).size().reset_index(name='num_success')

                total_by_poly = inv.groupby(['Num_Polyploids', 'method']).size().reset_index(name='num_total')

                success_rate = success_by_poly.merge(total_by_poly, on=['Num_Polyploids', 'method'])
                success_rate['success_rate'] = success_rate['num_success'] / success_rate['num_total'] * 100

                # Plot
                for method in success_rate['method'].unique():
                    method_data = success_rate[success_rate['method'] == method]
                    ax.plot(method_data['Num_Polyploids'], method_data['success_rate'],
                           'o-', label=method, color=METHOD_COLORS.get(method, '#000000'),
                           markersize=6, linewidth=2, alpha=0.7)

            ils_level = self._extract_ils_level(config)
            ax.set_xlabel('Number of Polyploid Species', fontsize=11)
            ax.set_title(f'ILS {ils_level.title()}', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_ylim(-5, 105)

        axes[0].set_ylabel('Success Rate (%)', fontsize=11)
        axes[-1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)

        plt.tight_layout()
        fig.savefig(self.output_dir / "fig2_success_vs_polyploids.pdf", bbox_inches='tight')
        fig.savefig(self.output_dir / "fig2_success_vs_polyploids.png", bbox_inches='tight')
        print(f"[OK] Created: fig2_success_vs_polyploids")
        plt.close()

    def plot_success_vs_wgd(self):
        """Plot number of successful runs vs total WGD events"""
        fig, axes = plt.subplots(1, len(self.config_names), figsize=(5*len(self.config_names), 4),
                                sharey=True)

        if len(self.config_names) == 1:
            axes = [axes]

        for ax, config in zip(axes, self.config_names):
            if 'inventory' not in self.data[config]:
                continue

            inv = self.data[config]['inventory']

            if self.network_stats is not None:
                inv = inv.merge(self.network_stats[['network', 'Total_WGD']],
                               on='network', how='left')

                # Group by WGD and method
                success_by_wgd = inv[inv['status'] == 'exists'].groupby(
                    ['Total_WGD', 'method']
                ).size().reset_index(name='num_success')

                total_by_wgd = inv.groupby(['Total_WGD', 'method']).size().reset_index(name='num_total')

                success_rate = success_by_wgd.merge(total_by_wgd, on=['Total_WGD', 'method'])
                success_rate['success_rate'] = success_rate['num_success'] / success_rate['num_total'] * 100

                # Plot
                for method in success_rate['method'].unique():
                    method_data = success_rate[success_rate['method'] == method]
                    ax.plot(method_data['Total_WGD'], method_data['success_rate'],
                           'o-', label=method, color=METHOD_COLORS.get(method, '#000000'),
                           markersize=6, linewidth=2, alpha=0.7)

            ils_level = self._extract_ils_level(config)
            ax.set_xlabel('Total WGD Events (Auto + Allo)', fontsize=11)
            ax.set_title(f'ILS {ils_level.title()}', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_ylim(-5, 105)

        axes[0].set_ylabel('Success Rate (%)', fontsize=11)
        axes[-1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)

        plt.tight_layout()
        fig.savefig(self.output_dir / "fig3_success_vs_total_wgd.pdf", bbox_inches='tight')
        fig.savefig(self.output_dir / "fig3_success_vs_total_wgd.png", bbox_inches='tight')
        print(f"[OK] Created: fig3_success_vs_total_wgd")
        plt.close()

    def plot_success_vs_autopolyploidization_events(self):
        """Plot number of successful runs vs number of autopolyploidization events"""
        fig, axes = plt.subplots(1, len(self.config_names), figsize=(5*len(self.config_names), 4),
                                sharey=True)

        if len(self.config_names) == 1:
            axes = [axes]

        for ax, config in zip(axes, self.config_names):
            if 'inventory' not in self.data[config]:
                continue

            inv = self.data[config]['inventory']

            if self.network_stats is not None:
                inv = inv.merge(self.network_stats[['network', 'Num_Autopolyploidization_Events']],
                               on='network', how='left')

                # Group by autopolyploidization events and method
                success_by_auto = inv[inv['status'] == 'exists'].groupby(
                    ['Num_Autopolyploidization_Events', 'method']
                ).size().reset_index(name='num_success')

                total_by_auto = inv.groupby(['Num_Autopolyploidization_Events', 'method']).size().reset_index(name='num_total')

                success_rate = success_by_auto.merge(total_by_auto, on=['Num_Autopolyploidization_Events', 'method'])
                success_rate['success_rate'] = success_rate['num_success'] / success_rate['num_total'] * 100

                # Plot
                for method in success_rate['method'].unique():
                    method_data = success_rate[success_rate['method'] == method]
                    ax.plot(method_data['Num_Autopolyploidization_Events'], method_data['success_rate'],
                           'o-', label=method, color=METHOD_COLORS.get(method, '#000000'),
                           markersize=6, linewidth=2, alpha=0.7)

            ils_level = self._extract_ils_level(config)
            ax.set_xlabel('Number of Autopolyploidization Events', fontsize=11)
            ax.set_title(f'ILS {ils_level.title()}', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_ylim(-5, 105)

        axes[0].set_ylabel('Success Rate (%)', fontsize=11)
        axes[-1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)

        plt.tight_layout()
        fig.savefig(self.output_dir / "fig4_success_vs_autopolyploidization_events.pdf", bbox_inches='tight')
        fig.savefig(self.output_dir / "fig4_success_vs_autopolyploidization_events.png", bbox_inches='tight')
        print(f"[OK] Created: fig4_success_vs_autopolyploidization_events")
        plt.close()

    def plot_edit_distance_comparison(self):
        """
        Plot edit distance comparison across methods

        Boxplot showing distribution of edit distances for each method across all networks
        """
        fig, axes = plt.subplots(1, len(self.config_names), figsize=(6*len(self.config_names), 5),
                                sharey=True)

        if len(self.config_names) == 1:
            axes = [axes]

        for ax, config in zip(axes, self.config_names):
            if 'aggregated' not in self.data[config]:
                continue

            agg = self.data[config]['aggregated']

            # Filter to edit_distance metric
            edit_dist = agg[agg['metric'] == 'edit_distance']

            if len(edit_dist) == 0:
                continue

            # Prepare data for boxplot
            methods = sorted(edit_dist['method'].unique())
            data_for_box = [edit_dist[edit_dist['method'] == m]['mean'].values for m in methods]
            colors = [METHOD_COLORS.get(m, '#cccccc') for m in methods]

            # Create boxplot
            bp = ax.boxplot(data_for_box, labels=methods, patch_artist=True,
                           showfliers=True, widths=0.6)

            # Color boxes
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            # Overlay individual points
            for i, (method, data) in enumerate(zip(methods, data_for_box), 1):
                x = np.random.normal(i, 0.04, size=len(data))
                ax.scatter(x, data, alpha=0.4, s=20, color=colors[i-1], edgecolors='black', linewidths=0.5)

            ils_level = self._extract_ils_level(config)
            ax.set_ylabel('Edit Distance (Normalized)', fontsize=11)
            ax.set_xlabel('Method', fontsize=11)
            ax.set_title(f'ILS {ils_level.title()}', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--', axis='y')
            ax.tick_params(axis='x', rotation=45)

        plt.tight_layout()
        fig.savefig(self.output_dir / "fig5_edit_distance_boxplot.pdf", bbox_inches='tight')
        fig.savefig(self.output_dir / "fig5_edit_distance_boxplot.png", bbox_inches='tight')
        print(f"[OK] Created: fig5_edit_distance_boxplot")
        plt.close()

    def plot_correlation_heatmap(self):
        """
        Create correlation heatmap between network characteristics and method performance

        Shows which network properties correlate with method success/failure
        """
        if self.network_stats is None:
            print("Warning: Cannot create correlation heatmap without network stats")
            return

        for config in self.config_names:
            if 'aggregated' not in self.data[config]:
                continue

            agg = self.data[config]['aggregated']

            # Merge with network stats
            agg_with_stats = agg.merge(
                self.network_stats[['network', 'Num_Species', 'Num_Polyploids', 'Max_Copies',
                                   'H_Strict', 'Num_Autopolyploidization_Events', 'Total_WGD']],
                on='network', how='left'
            )

            # Filter to edit_distance
            edit_dist = agg_with_stats[agg_with_stats['metric'] == 'edit_distance']

            # Calculate correlations for each method
            methods = edit_dist['method'].unique()
            network_props = ['Num_Species', 'Num_Polyploids', 'Max_Copies', 'H_Strict',
                            'Num_Autopolyploidization_Events', 'Total_WGD']

            # Build correlation matrix
            corr_matrix = pd.DataFrame(index=methods, columns=network_props)

            for method in methods:
                method_data = edit_dist[edit_dist['method'] == method]
                for prop in network_props:
                    valid_data = method_data[[prop, 'mean']].dropna()
                    if len(valid_data) > 3:
                        corr, _ = pearsonr(valid_data[prop], valid_data['mean'])
                        corr_matrix.loc[method, prop] = corr
                    else:
                        corr_matrix.loc[method, prop] = np.nan

            # Convert to float
            corr_matrix = corr_matrix.astype(float)

            # Plot heatmap
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.heatmap(corr_matrix, annot=True, cmap='RdYlGn_r', center=0, vmin=-1, vmax=1,
                       fmt='.2f', linewidths=0.5, cbar_kws={'label': 'Pearson Correlation'},
                       ax=ax, annot_kws={'size': 9})

            ils_level = self._extract_ils_level(config)
            ax.set_title(f'Network Characteristics vs Edit Distance - ILS {ils_level.title()}',
                        fontsize=13, fontweight='bold', pad=15)
            ax.set_xlabel('Network Characteristics', fontsize=11)
            ax.set_ylabel('Method', fontsize=11)

            plt.tight_layout()
            fig.savefig(self.output_dir / f"fig6_correlation_heatmap_{config}.pdf", bbox_inches='tight')
            fig.savefig(self.output_dir / f"fig6_correlation_heatmap_{config}.png", bbox_inches='tight')
            print(f"[OK] Created: fig6_correlation_heatmap_{config}")
            plt.close()

    def create_summary_table(self):
        """Create comprehensive summary table for publication"""

        summary_rows = []

        for config in self.config_names:
            if 'aggregated' not in self.data[config] or 'inventory' not in self.data[config]:
                continue

            agg = self.data[config]['aggregated']
            inv = self.data[config]['inventory']

            ils_level = self._extract_ils_level(config)

            # Calculate success rates
            for method in inv['method'].unique():
                method_inv = inv[inv['method'] == method]
                num_total = len(method_inv)
                num_success = len(method_inv[method_inv['status'] == 'exists'])
                success_rate = num_success / num_total * 100 if num_total > 0 else 0

                # Get edit distance stats
                method_agg = agg[(agg['method'] == method) & (agg['metric'] == 'edit_distance')]
                if len(method_agg) > 0:
                    mean_edit_dist = method_agg['mean'].mean()
                    std_edit_dist = method_agg['std'].mean()
                else:
                    mean_edit_dist = np.nan
                    std_edit_dist = np.nan

                summary_rows.append({
                    'ILS_Level': ils_level,
                    'Method': method,
                    'Num_Networks': num_total,
                    'Num_Success': num_success,
                    'Success_Rate_%': success_rate,
                    'Mean_Edit_Distance': mean_edit_dist,
                    'Std_Edit_Distance': std_edit_dist
                })

        summary_df = pd.DataFrame(summary_rows)

        # Save
        summary_file = self.output_dir / "table1_method_summary.csv"
        summary_df.to_csv(summary_file, index=False, float_format='%.4f')
        print(f"[OK] Created: table1_method_summary.csv")

        # Also create LaTeX version
        latex_file = self.output_dir / "table1_method_summary.tex"
        with open(latex_file, 'w') as f:
            latex_table = summary_df.to_latex(index=False, float_format='%.2f', escape=False)
            f.write(latex_table)
        print(f"[OK] Created: table1_method_summary.tex")

        return summary_df

    def generate_all_figures(self):
        """Generate all publication figures"""
        print(f"\n{'='*80}")
        print("Generating Publication Figures")
        print(f"{'='*80}\n")

        # Success rate plots
        print("[1/7] Creating success vs reticulations plot...")
        self.plot_success_vs_reticulations()

        print("[2/7] Creating success vs polyploids plot...")
        self.plot_success_vs_polyploids()

        print("[3/7] Creating success vs total WGD plot...")
        self.plot_success_vs_wgd()

        print("[4/7] Creating success vs autopolyploidization events plot...")
        self.plot_success_vs_autopolyploidization_events()

        # Performance comparison
        print("[5/7] Creating edit distance boxplot...")
        self.plot_edit_distance_comparison()

        # Correlation analysis
        print("[6/7] Creating correlation heatmap...")
        self.plot_correlation_heatmap()

        # Summary table
        print("[7/7] Creating summary tables...")
        self.create_summary_table()

        print(f"\n{'='*80}")
        print(f"All figures saved to: {self.output_dir}")
        print(f"{'='*80}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate publication-quality figures for phylogenetic network inference analysis'
    )

    parser.add_argument('--configs', nargs='+', required=True,
                       help='Configuration names (e.g., conf_ils_low_10M conf_ils_medium_10M)')
    parser.add_argument('--network-stats', default='simulations/networks/mul_tree_final_stats.csv',
                       help='Path to network characteristics CSV')
    parser.add_argument('--output', default='simulations/analysis/publication_figures',
                       help='Output directory for figures')

    args = parser.parse_args()

    # Create figure generator
    generator = PublicationFigureGenerator(
        config_names=args.configs,
        network_stats_file=args.network_stats,
        output_dir=args.output
    )

    # Generate all figures
    generator.generate_all_figures()


if __name__ == '__main__':
    main()
