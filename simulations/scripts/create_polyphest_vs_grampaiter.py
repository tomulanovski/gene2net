#!/usr/bin/env python3
"""
create_polyphest_vs_grampaiter.py - Head-to-head comparison figures

Generates focused comparison figures between Polyphest and GRAMPA^Iter
across all simulation configurations.

Usage:
    python create_polyphest_vs_grampaiter.py
    python create_polyphest_vs_grampaiter.py --configs conf_ils_low_10M conf_ils_medium_10M
"""

import argparse
import gc
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Import merging utilities
from create_analysis_figures import (
    merge_polyphest_inventory, merge_polyphest_comparisons, reaggregate_metrics
)

SCRIPT_DIR = Path(__file__).resolve().parent
SUMMARY_BASE = SCRIPT_DIR.parent / "analysis" / "summary"
NETWORKS_DIR = SCRIPT_DIR.parent / "networks"

ALL_CONFIGS = [
    "conf_ils_low_10M",
    "conf_ils_medium_10M",
    "conf_ils_high_10M",
    "conf_dup_loss_low_10M",
    "conf_dup_loss_medium_10M",
    "conf_dup_loss_high_10M",
    "conf_dup_loss_low_10M_ne1M",
    "conf_dup_loss_medium_10M_ne1M",
    "conf_dup_loss_high_10M_ne1M",
]

CONFIG_FAMILIES = {
    'ILS': {
        'label': 'ILS Level',
        'configs': {
            'Low': 'conf_ils_low_10M',
            'Medium': 'conf_ils_medium_10M',
            'High': 'conf_ils_high_10M',
        }
    },
    'Dup/Loss\n(low ILS)': {
        'label': 'Dup/Loss Rate',
        'configs': {
            'Low': 'conf_dup_loss_low_10M',
            'Medium': 'conf_dup_loss_medium_10M',
            'High': 'conf_dup_loss_high_10M',
        }
    },
    'Dup/Loss\n(med ILS)': {
        'label': 'Dup/Loss Rate (Ne=1M)',
        'configs': {
            'Low': 'conf_dup_loss_low_10M_ne1M',
            'Medium': 'conf_dup_loss_medium_10M_ne1M',
            'High': 'conf_dup_loss_high_10M_ne1M',
        }
    },
}

LEVEL_ORDER = ['Low', 'Medium', 'High']

METHODS = ['polyphest', 'grandma_split']
METHOD_COLORS = {
    'polyphest': '#DE8F05',
    'grandma_split': '#56B4E9',
}
METHOD_DISPLAY = {
    'polyphest': 'Polyphest',
    'grandma_split': r'GRAMPA$^{Iter}$',
}

METRICS = {
    'edit_distance_multree': 'Edit Distance',
    'ret_leaf_jaccard.dist': 'Reticulation Leaf Distance',
    'ret_sisters_jaccard.dist': 'Sister-Taxa Distance',
    'ploidy_diff.dist': 'Ploidy Distance',
    'num_rets_diff': 'Reticulation Count Error (|diff|)',
    'num_rets_bias': 'Reticulation Count Bias (%)',
}

# Plot style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['lines.markersize'] = 8


def dn(method):
    return METHOD_DISPLAY.get(method, method)


# ============================================================================
# DATA LOADING
# ============================================================================

def load_data(configs, summary_base):
    """Load and merge inventory + comparisons from all configs."""
    inv_frames = []
    comp_frames = []

    for config in configs:
        config_dir = summary_base / config
        if not config_dir.exists():
            print(f"  WARNING: {config} not found, skipping")
            continue

        inv_file = config_dir / "inventory.csv"
        if inv_file.exists():
            df = pd.read_csv(inv_file)
            df['config'] = config
            inv_frames.append(df)

        comp_file = config_dir / "comparisons_raw.csv"
        if comp_file.exists():
            df = pd.read_csv(comp_file)
            if 'config' not in df.columns:
                df['config'] = config
            comp_frames.append(df)

    if not inv_frames:
        print("ERROR: No data loaded")
        sys.exit(1)

    inventory = pd.concat(inv_frames, ignore_index=True)
    comparisons = pd.concat(comp_frames, ignore_index=True) if comp_frames else pd.DataFrame()

    # Merge polyphest thresholds
    inventory = merge_polyphest_inventory(inventory)
    if not comparisons.empty:
        comparisons = merge_polyphest_comparisons(comparisons)

    # Filter to our two methods
    inventory = inventory[inventory['method'].isin(METHODS)].copy()
    if not comparisons.empty:
        comparisons = comparisons[comparisons['method'].isin(METHODS)].copy()

    # Tag family and level
    inventory = tag_family(inventory)
    if not comparisons.empty:
        comparisons = tag_family(comparisons)

    return inventory, comparisons


def tag_family(df):
    """Add family and level columns."""
    families = []
    levels = []
    for config in df['config']:
        found = False
        for fam_name, fam_info in CONFIG_FAMILIES.items():
            for level, cfg in fam_info['configs'].items():
                if config == cfg:
                    families.append(fam_name)
                    levels.append(level)
                    found = True
                    break
            if found:
                break
        if not found:
            families.append('Unknown')
            levels.append('Unknown')
    df['family'] = families
    df['level'] = levels
    return df


# ============================================================================
# FIGURES
# ============================================================================

class PolyphestVsGrampaIter:

    def __init__(self, inventory, comparisons, output_dir, network_stats=None):
        self.inventory = inventory
        self.comparisons = comparisons
        self.output_dir = output_dir
        self.plots_dir = output_dir / "plots"
        self.tables_dir = output_dir / "tables"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        self.network_stats = network_stats

        # Convert num_rets_bias from raw counts to percentage: (bias / H_Strict) * 100
        if network_stats is not None and not comparisons.empty:
            self._convert_bias_to_percentage()

        print(f"\nPolyphest vs GRAMPA^Iter analyzer:")
        print(f"  Inventory: {len(self.inventory)} rows")
        print(f"  Comparisons: {len(self.comparisons)} rows")
        print(f"  Output: {self.output_dir}")

    def _convert_bias_to_percentage(self):
        """Convert num_rets_bias values from raw counts to percentage of true H_Strict."""
        mask = self.comparisons['metric'] == 'num_rets_bias'
        if not mask.any():
            return

        # Build network -> Total_WGD lookup (auto + allo events)
        wgd_lookup = self.network_stats.set_index('network')['Total_WGD'].to_dict()

        idx = self.comparisons.index[mask]
        for i in idx:
            wgd = wgd_lookup.get(self.comparisons.at[i, 'network'])
            if wgd is not None and wgd > 0:
                self.comparisons.at[i, 'value'] = self.comparisons.at[i, 'value'] / wgd * 100
            # Total_WGD == 0 or missing: keep raw value

        print(f"  Converted num_rets_bias to percentage (÷ Total_WGD × 100)")

    def generate_all(self):
        print(f"\n{'='*70}")
        print("Generating Polyphest vs GRAMPA^Iter Figures")
        print(f"{'='*70}\n")

        self.plot_completion_comparison()
        self.plot_metric_across_conditions('edit_distance_multree')
        self.plot_metric_across_conditions('ret_leaf_jaccard.dist')
        self.plot_metric_across_conditions('ret_sisters_jaccard.dist')
        self.plot_metric_across_conditions('ploidy_diff.dist')
        self.plot_metric_across_conditions('num_rets_bias')
        self.plot_metric_across_conditions('num_rets_diff')
        self.plot_combined_accuracy_panel()
        self.plot_metric_boxplots()
        self.plot_per_network_comparison()
        self.plot_degradation_lines()
        self.generate_summary_table()

        print(f"\n{'='*70}")
        print(f"Done! Output: {self.output_dir}")
        print(f"{'='*70}\n")

    # ------------------------------------------------------------------
    # 1. Completion rate comparison
    # ------------------------------------------------------------------
    def plot_completion_comparison(self):
        """Side-by-side completion rates across all configs."""
        fig, axes = plt.subplots(1, len(CONFIG_FAMILIES), figsize=(5 * len(CONFIG_FAMILIES), 5),
                                 squeeze=False, sharey=True)
        axes = axes.flatten()

        for fam_idx, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
            ax = axes[fam_idx]
            x = np.arange(len(LEVEL_ORDER))
            bar_width = 0.35

            for m_idx, method in enumerate(METHODS):
                rates = []
                for level in LEVEL_ORDER:
                    cfg = fam_info['configs'].get(level)
                    data = self.inventory[
                        (self.inventory['method'] == method) &
                        (self.inventory['config'] == cfg)
                    ]
                    rates.append(data['inferred_exists'].mean() * 100 if len(data) > 0 else 0)

                offset = (m_idx - 0.5) * bar_width
                ax.bar(x + offset, rates, bar_width * 0.9,
                       label=dn(method), color=METHOD_COLORS[method],
                       edgecolor='white', linewidth=0.8, alpha=0.85)

                # Value labels
                for i, val in enumerate(rates):
                    ax.text(x[i] + offset, val + 1, f'{val:.0f}%',
                            ha='center', va='bottom', fontsize=9, fontweight='bold')

            ax.set_xlabel(fam_info['label'], fontsize=12, fontweight='bold')
            ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
            ax.set_xticks(x)
            ax.set_xticklabels(LEVEL_ORDER, fontsize=11)
            ax.set_ylim(0, 115)
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            if fam_idx == 0:
                ax.set_ylabel('Completion Rate (%)', fontsize=12, fontweight='bold')
                ax.legend(fontsize=10, framealpha=0.9)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "01_completion_comparison.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "01_completion_comparison.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()
        print("  [1] Completion comparison")

    # ------------------------------------------------------------------
    # 2. Single metric across conditions (grouped bars)
    # ------------------------------------------------------------------
    def plot_metric_across_conditions(self, metric_key):
        """Grouped bar chart for one metric across config families."""
        if self.comparisons.empty:
            return

        metric_label = METRICS.get(metric_key, metric_key)
        metric_data = self.comparisons[
            (self.comparisons['metric'] == metric_key) &
            (self.comparisons['status'] == 'SUCCESS')
        ]
        if metric_data.empty:
            print(f"  WARNING: No data for {metric_key}")
            return

        fig, axes = plt.subplots(1, len(CONFIG_FAMILIES),
                                 figsize=(5 * len(CONFIG_FAMILIES), 5),
                                 squeeze=False, sharey=True)
        axes = axes.flatten()

        for fam_idx, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
            ax = axes[fam_idx]
            x = np.arange(len(LEVEL_ORDER))
            bar_width = 0.35

            for m_idx, method in enumerate(METHODS):
                means = []
                sems = []
                for level in LEVEL_ORDER:
                    cfg = fam_info['configs'].get(level)
                    vals = metric_data[
                        (metric_data['method'] == method) &
                        (metric_data['config'] == cfg)
                    ]['value']
                    if len(vals) > 0:
                        means.append(vals.mean())
                        sems.append(vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0)
                    else:
                        means.append(np.nan)
                        sems.append(0)

                offset = (m_idx - 0.5) * bar_width
                ax.bar(x + offset, means, bar_width * 0.9,
                       yerr=sems, capsize=3,
                       label=dn(method), color=METHOD_COLORS[method],
                       edgecolor='white', linewidth=0.8, alpha=0.85)

            if metric_key == 'num_rets_bias':
                ax.axhline(y=0, color='black', linewidth=1, linestyle='-', alpha=0.5)

            ax.set_xlabel(fam_info['label'], fontsize=12, fontweight='bold')
            ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
            ax.set_xticks(x)
            ax.set_xticklabels(LEVEL_ORDER, fontsize=11)
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            if fam_idx == 0:
                ax.set_ylabel(metric_label, fontsize=12, fontweight='bold')
                ax.legend(fontsize=10, framealpha=0.9)

        plt.tight_layout()
        safe = metric_key.replace('.', '_')
        fig.savefig(self.plots_dir / f"02_{safe}_across_conditions.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"02_{safe}_across_conditions.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()
        print(f"  [2] {metric_label} across conditions")

    # ------------------------------------------------------------------
    # 3. Combined accuracy panel (edit dist + ret leaf + bias in one figure)
    # ------------------------------------------------------------------
    def plot_combined_accuracy_panel(self):
        """Multi-panel figure: edit distance, ret leaf distance, reticulation bias."""
        if self.comparisons.empty:
            return

        panel_metrics = [
            ('edit_distance_multree', 'Edit Distance', None),
            ('ret_leaf_jaccard.dist', 'Ret. Leaf Distance', None),
            ('num_rets_bias', 'Ret. Count Bias (%)', 0),  # 0 = reference line
        ]

        fig, axes = plt.subplots(len(panel_metrics), len(CONFIG_FAMILIES),
                                 figsize=(5 * len(CONFIG_FAMILIES), 4.5 * len(panel_metrics)),
                                 squeeze=False, sharey='row')

        for row_idx, (metric_key, metric_label, ref_line) in enumerate(panel_metrics):
            metric_data = self.comparisons[
                (self.comparisons['metric'] == metric_key) &
                (self.comparisons['status'] == 'SUCCESS')
            ]

            for fam_idx, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
                ax = axes[row_idx, fam_idx]
                x = np.arange(len(LEVEL_ORDER))
                bar_width = 0.35

                for m_idx, method in enumerate(METHODS):
                    means = []
                    sems = []
                    for level in LEVEL_ORDER:
                        cfg = fam_info['configs'].get(level)
                        vals = metric_data[
                            (metric_data['method'] == method) &
                            (metric_data['config'] == cfg)
                        ]['value']
                        if len(vals) > 0:
                            means.append(vals.mean())
                            sems.append(vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0)
                        else:
                            means.append(np.nan)
                            sems.append(0)

                    offset = (m_idx - 0.5) * bar_width
                    ax.bar(x + offset, means, bar_width * 0.9,
                           yerr=sems, capsize=3,
                           label=dn(method) if row_idx == 0 else None,
                           color=METHOD_COLORS[method],
                           edgecolor='white', linewidth=0.8, alpha=0.85)

                if ref_line is not None:
                    ax.axhline(y=ref_line, color='black', linewidth=1, linestyle='-', alpha=0.5)

                ax.set_xticks(x)
                ax.set_xticklabels(LEVEL_ORDER, fontsize=10)
                ax.grid(True, alpha=0.25, linestyle='--', axis='y')

                if row_idx == 0:
                    ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
                if row_idx == len(panel_metrics) - 1:
                    ax.set_xlabel(fam_info['label'], fontsize=11, fontweight='bold')
                if fam_idx == 0:
                    ax.set_ylabel(metric_label, fontsize=11, fontweight='bold')

        # Single legend at top
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', ncol=2, fontsize=12,
                   framealpha=0.9, bbox_to_anchor=(0.5, 1.02))

        plt.tight_layout(rect=[0, 0, 1, 0.97])
        fig.savefig(self.plots_dir / "03_combined_accuracy_panel.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "03_combined_accuracy_panel.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()
        print("  [3] Combined accuracy panel")

    # ------------------------------------------------------------------
    # 4. Boxplots (aggregated across all configs)
    # ------------------------------------------------------------------
    def plot_metric_boxplots(self):
        """Boxplots comparing the two methods on key metrics, all configs combined."""
        if self.comparisons.empty:
            return

        box_metrics = [
            ('edit_distance_multree', 'Edit Distance'),
            ('ret_leaf_jaccard.dist', 'Ret. Leaf Distance'),
            ('ret_sisters_jaccard.dist', 'Sister-Taxa Distance'),
            ('ploidy_diff.dist', 'Ploidy Distance'),
        ]

        fig, axes = plt.subplots(1, len(box_metrics),
                                 figsize=(4.5 * len(box_metrics), 5), squeeze=False)
        axes = axes.flatten()

        for idx, (metric_key, metric_label) in enumerate(box_metrics):
            ax = axes[idx]
            metric_data = self.comparisons[
                (self.comparisons['metric'] == metric_key) &
                (self.comparisons['status'] == 'SUCCESS')
            ]

            data_list = []
            labels = []
            colors = []
            for method in METHODS:
                vals = metric_data[metric_data['method'] == method]['value'].dropna()
                if len(vals) > 0:
                    data_list.append(vals.values)
                    labels.append(dn(method))
                    colors.append(METHOD_COLORS[method])

            if data_list:
                bp = ax.boxplot(data_list, labels=labels, patch_artist=True,
                               widths=0.6, showmeans=True,
                               meanprops=dict(marker='D', markerfacecolor='white',
                                             markeredgecolor='black', markersize=6))
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)

            ax.set_title(metric_label, fontsize=13, fontweight='bold')
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            ax.set_ylim(bottom=0)

        fig.suptitle('Overall Distribution (All Configurations)',
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        fig.savefig(self.plots_dir / "04_metric_boxplots.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "04_metric_boxplots.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()
        print("  [4] Metric boxplots")

    # ------------------------------------------------------------------
    # 5. Per-network comparison (paired bars or scatter)
    # ------------------------------------------------------------------
    def plot_per_network_comparison(self):
        """Per-network edit distance under low ILS: paired bars sorted by network complexity."""
        if self.comparisons.empty or self.network_stats is None:
            return

        # Use low ILS as baseline
        baseline_config = 'conf_ils_low_10M'
        metric_key = 'edit_distance_multree'

        metric_data = self.comparisons[
            (self.comparisons['metric'] == metric_key) &
            (self.comparisons['status'] == 'SUCCESS') &
            (self.comparisons['config'] == baseline_config)
        ]

        if metric_data.empty:
            return

        # Aggregate per network × method
        agg = metric_data.groupby(['network', 'method'])['value'].agg(['mean', 'std']).reset_index()

        # Sort networks by H_Strict
        net_order = self.network_stats.sort_values('H_Strict')['network'].tolist()
        net_order = [n for n in net_order if n in agg['network'].unique()]

        fig, ax = plt.subplots(figsize=(max(14, len(net_order) * 0.8), 6))
        x = np.arange(len(net_order))
        bar_width = 0.38

        for m_idx, method in enumerate(METHODS):
            method_agg = agg[agg['method'] == method].set_index('network')
            means = [method_agg.loc[n, 'mean'] if n in method_agg.index else np.nan for n in net_order]
            stds = [method_agg.loc[n, 'std'] if n in method_agg.index else 0 for n in net_order]

            offset = (m_idx - 0.5) * bar_width
            ax.bar(x + offset, means, bar_width * 0.9,
                   yerr=stds, capsize=2,
                   label=dn(method), color=METHOD_COLORS[method],
                   edgecolor='white', linewidth=0.5, alpha=0.85)

        # Add H_Strict values as secondary labels
        h_strict_vals = []
        for n in net_order:
            row = self.network_stats[self.network_stats['network'] == n]
            h_strict_vals.append(int(row['H_Strict'].values[0]) if len(row) > 0 else '?')

        ax.set_xticks(x)
        labels = [f'{n}\n(H={h})' for n, h in zip(net_order, h_strict_vals)]
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Edit Distance', fontsize=12, fontweight='bold')
        ax.set_title(f'Per-Network Edit Distance (Low ILS)', fontsize=14, fontweight='bold', pad=15)
        ax.legend(fontsize=11, framealpha=0.9, loc='upper left')
        ax.grid(True, alpha=0.25, linestyle='--', axis='y')
        ax.set_ylim(0, 1.05)

        plt.tight_layout()
        fig.savefig(self.plots_dir / "05_per_network_edit_distance.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "05_per_network_edit_distance.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()
        print("  [5] Per-network edit distance")

    # ------------------------------------------------------------------
    # 6. Degradation lines (how accuracy changes across conditions)
    # ------------------------------------------------------------------
    def plot_degradation_lines(self):
        """Line plots showing how each metric degrades across condition levels."""
        if self.comparisons.empty:
            return

        line_metrics = [
            ('edit_distance_multree', 'Edit Distance'),
            ('ret_leaf_jaccard.dist', 'Ret. Leaf Distance'),
            ('num_rets_bias', 'Ret. Bias (%)'),
        ]

        fig, axes = plt.subplots(len(line_metrics), len(CONFIG_FAMILIES),
                                 figsize=(5 * len(CONFIG_FAMILIES), 4 * len(line_metrics)),
                                 squeeze=False, sharey='row')

        for row_idx, (metric_key, metric_label) in enumerate(line_metrics):
            metric_data = self.comparisons[
                (self.comparisons['metric'] == metric_key) &
                (self.comparisons['status'] == 'SUCCESS')
            ]

            for fam_idx, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
                ax = axes[row_idx, fam_idx]

                for method in METHODS:
                    means = []
                    sems = []
                    for level in LEVEL_ORDER:
                        cfg = fam_info['configs'].get(level)
                        vals = metric_data[
                            (metric_data['method'] == method) &
                            (metric_data['config'] == cfg)
                        ]['value']
                        if len(vals) > 0:
                            means.append(vals.mean())
                            sems.append(vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0)
                        else:
                            means.append(np.nan)
                            sems.append(0)

                    ax.errorbar(LEVEL_ORDER, means, yerr=sems,
                               marker='o', markersize=8, capsize=4,
                               label=dn(method) if row_idx == 0 and fam_idx == 0 else None,
                               color=METHOD_COLORS[method], linewidth=2.5)

                if metric_key == 'num_rets_bias':
                    ax.axhline(y=0, color='black', linewidth=1, linestyle='--', alpha=0.5)

                ax.grid(True, alpha=0.25, linestyle='--')
                if row_idx == 0:
                    ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
                if row_idx == len(line_metrics) - 1:
                    ax.set_xlabel(fam_info['label'], fontsize=11, fontweight='bold')
                if fam_idx == 0:
                    ax.set_ylabel(metric_label, fontsize=11, fontweight='bold')

        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', ncol=2, fontsize=12,
                   framealpha=0.9, bbox_to_anchor=(0.5, 1.02))

        plt.tight_layout(rect=[0, 0, 1, 0.97])
        fig.savefig(self.plots_dir / "06_degradation_lines.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "06_degradation_lines.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()
        print("  [6] Degradation lines")

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------
    def generate_summary_table(self):
        """CSV table: method × config × key metrics."""
        if self.comparisons.empty:
            return

        rows = []
        for method in METHODS:
            for config in sorted(self.inventory['config'].unique()):
                inv = self.inventory[
                    (self.inventory['method'] == method) &
                    (self.inventory['config'] == config)
                ]
                row = {
                    'method': dn(method),
                    'config': config,
                    'completion_rate': inv['inferred_exists'].mean() * 100 if len(inv) > 0 else 0,
                }

                for metric_key in METRICS:
                    vals = self.comparisons[
                        (self.comparisons['method'] == method) &
                        (self.comparisons['config'] == config) &
                        (self.comparisons['metric'] == metric_key) &
                        (self.comparisons['status'] == 'SUCCESS')
                    ]['value']
                    row[f'{metric_key}_mean'] = vals.mean() if len(vals) > 0 else np.nan
                    row[f'{metric_key}_median'] = vals.median() if len(vals) > 0 else np.nan
                    row[f'{metric_key}_std'] = vals.std() if len(vals) > 1 else np.nan

                rows.append(row)

        summary = pd.DataFrame(rows)
        summary.to_csv(self.tables_dir / "01_polyphest_vs_grampaiter.csv", index=False)
        print(f"  [T] Summary table saved")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='Polyphest vs GRAMPA^Iter comparison figures')
    parser.add_argument('--configs', nargs='+', default=None,
                        help='Configs to include (default: all available)')
    parser.add_argument('--output', default=None,
                        help='Output directory')
    parser.add_argument('--network-stats', default=None,
                        help='Network stats CSV')
    args = parser.parse_args()

    configs = args.configs or ALL_CONFIGS

    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = SUMMARY_BASE / "polyphest_vs_grampaiter"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Network stats
    ns_path = Path(args.network_stats) if args.network_stats else NETWORKS_DIR / "mul_tree_final_stats.csv"
    network_stats = None
    if ns_path.exists():
        network_stats = pd.read_csv(ns_path)
        network_stats['network'] = network_stats['Filename'].str.replace('.tre', '')

    # Load data
    print(f"Loading data from: {SUMMARY_BASE}")
    inventory, comparisons = load_data(configs, SUMMARY_BASE)

    # Generate figures
    analyzer = PolyphestVsGrampaIter(inventory, comparisons, output_dir, network_stats)
    analyzer.generate_all()


if __name__ == '__main__':
    main()
