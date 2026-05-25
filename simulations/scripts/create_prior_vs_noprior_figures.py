#!/usr/bin/env python3
"""
create_prior_vs_noprior_figures.py - Ploidy Prior Impact Comparison

Head-to-head comparison of GRAMPA^Iter with and without ploidy prior,
alongside Polyphest as a reference, across all simulation configurations.

Usage:
    python create_prior_vs_noprior_figures.py
    python create_prior_vs_noprior_figures.py --configs conf_ils_low_10M conf_ils_medium_10M
    python create_prior_vs_noprior_figures.py --output analysis/summary/prior_comparison
"""

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from create_analysis_figures import (
    merge_polyphest_inventory, merge_polyphest_comparisons
)

SCRIPT_DIR   = Path(__file__).resolve().parent
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
    "conf_dup_loss_low_10M_ne2M",
    "conf_dup_loss_medium_10M_ne2M",
    "conf_dup_loss_high_10M_ne2M",
]

CONFIG_FAMILIES = {
    'ILS': {
        'label': 'ILS Level',
        'configs': {
            'Low':    'conf_ils_low_10M',
            'Medium': 'conf_ils_medium_10M',
            'High':   'conf_ils_high_10M',
        }
    },
    'Dup/Loss\n(low ILS)': {
        'label': 'Dup/Loss Rate',
        'configs': {
            'Low':    'conf_dup_loss_low_10M',
            'Medium': 'conf_dup_loss_medium_10M',
            'High':   'conf_dup_loss_high_10M',
        }
    },
    'Dup/Loss\n(med ILS)': {
        'label': 'Dup/Loss Rate (Ne=1M)',
        'configs': {
            'Low':    'conf_dup_loss_low_10M_ne1M',
            'Medium': 'conf_dup_loss_medium_10M_ne1M',
            'High':   'conf_dup_loss_high_10M_ne1M',
        }
    },
    'Dup/Loss\n(high ILS)': {
        'label': 'Dup/Loss Rate (Ne=2M)',
        'configs': {
            'Low':    'conf_dup_loss_low_10M_ne2M',
            'Medium': 'conf_dup_loss_medium_10M_ne2M',
            'High':   'conf_dup_loss_high_10M_ne2M',
        }
    },
}

LEVEL_ORDER = ['Low', 'Medium', 'High']

METHODS = ['polyphest', 'grandma_split', 'grandma_split_prior']

METHOD_COLORS = {
    'polyphest':           '#DE8F05',  # Orange
    'grandma_split':       '#56B4E9',  # Light Blue
    'grandma_split_prior': '#9467BD',  # Purple
}

METHOD_DISPLAY = {
    'polyphest':           'Polyphest',
    'grandma_split':       r'GRAMPA$^{Iter}$',
    'grandma_split_prior': r'GRAMPA$^{Iter}$ (prior)',
}

METRICS = {
    'edit_distance_multree':  'Edit Distance',
    'ret_leaf_jaccard.dist':  'Reticulation Descendants Measure',
    'ret_sisters_jaccard.dist': 'Reticulation Sister Measure',
    'ploidy_diff.dist':       'Ploidy Distance',
    'num_rets_diff':          'Reticulation Count Error (|diff|)',
    'num_rets_bias':          'Reticulation Count Bias (%)',
}

plt.rcParams['font.family']   = 'sans-serif'
plt.rcParams['font.size']     = 11
plt.rcParams['figure.dpi']    = 100
plt.rcParams['savefig.dpi']   = 300
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['lines.markersize'] = 8


def dn(method):
    return METHOD_DISPLAY.get(method, method)


# ── data loading ──────────────────────────────────────────────────────────────

def load_data(configs, summary_base):
    inv_frames  = []
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
        print("ERROR: No data loaded.")
        sys.exit(1)

    inventory   = pd.concat(inv_frames,  ignore_index=True)
    comparisons = pd.concat(comp_frames, ignore_index=True) if comp_frames else pd.DataFrame()

    inventory = merge_polyphest_inventory(inventory)
    if not comparisons.empty:
        comparisons = merge_polyphest_comparisons(comparisons)

    inventory   = inventory[inventory['method'].isin(METHODS)].copy()
    if not comparisons.empty:
        comparisons = comparisons[comparisons['method'].isin(METHODS)].copy()

    inventory   = _tag_family(inventory)
    if not comparisons.empty:
        comparisons = _tag_family(comparisons)

    return inventory, comparisons


def _tag_family(df):
    families, levels = [], []
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
    df['level']  = levels
    return df


# ── helper: aggregate one config/method/metric ───────────────────────────────

def _agg(metric_data, method, cfg, use_median=False):
    vals = metric_data[
        (metric_data['method'] == method) &
        (metric_data['config'] == cfg)
    ]['value']
    if len(vals) == 0:
        return np.nan, 0
    if use_median:
        return vals.median(), vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0
    return vals.mean(), vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0


# ── analyzer ─────────────────────────────────────────────────────────────────

class PriorComparison:

    def __init__(self, inventory, comparisons, output_dir, network_stats=None):
        self.inventory    = inventory
        self.comparisons  = comparisons
        self.output_dir   = output_dir
        self.plots_dir    = output_dir / "plots"
        self.tables_dir   = output_dir / "tables"
        self.network_stats = network_stats

        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.tables_dir.mkdir(parents=True, exist_ok=True)

        if network_stats is not None and not comparisons.empty:
            self._convert_bias_to_percentage()

        print(f"\nPrior comparison analyzer:")
        print(f"  Inventory:   {len(self.inventory)} rows")
        print(f"  Comparisons: {len(self.comparisons)} rows")
        print(f"  Output:      {self.output_dir}")

    def _convert_bias_to_percentage(self):
        mask = self.comparisons['metric'] == 'num_rets_bias'
        if not mask.any():
            return
        wgd_lookup = self.network_stats.set_index('network')['Total_WGD'].to_dict()
        for i in self.comparisons.index[mask]:
            wgd = wgd_lookup.get(self.comparisons.at[i, 'network'])
            if wgd and wgd > 0:
                self.comparisons.at[i, 'value'] /= wgd / 100

    def _save(self, fig, stem):
        fig.savefig(self.plots_dir / f"{stem}.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"{stem}.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()

    def generate_all(self):
        print(f"\n{'='*70}")
        print("Generating Prior vs No-Prior Figures")
        print(f"{'='*70}\n")

        self.plot_completion()
        for metric in METRICS:
            self.plot_metric_across_conditions(metric)
        self.plot_combined_panel()
        self.plot_boxplots()
        self.plot_per_network()
        self.plot_degradation_lines()
        self.generate_summary_table()

        print(f"\nDone! Output: {self.output_dir}\n")

    # ── 1. completion ──────────────────────────────────────────────────────────
    def plot_completion(self):
        n_fam   = len(CONFIG_FAMILIES)
        n_meth  = len(METHODS)
        bar_w   = 0.7 / n_meth

        fig, axes = plt.subplots(1, n_fam, figsize=(5 * n_fam, 5),
                                 squeeze=False, sharey=True)
        axes = axes.flatten()

        for fi, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
            ax = axes[fi]
            x  = np.arange(len(LEVEL_ORDER))

            for mi, method in enumerate(METHODS):
                rates = []
                for level in LEVEL_ORDER:
                    cfg  = fam_info['configs'].get(level)
                    data = self.inventory[
                        (self.inventory['method'] == method) &
                        (self.inventory['config'] == cfg)
                    ]
                    rates.append(data['inferred_exists'].mean() * 100 if len(data) > 0 else 0)

                offset = (mi - n_meth / 2 + 0.5) * bar_w
                ax.bar(x + offset, rates, bar_w * 0.9,
                       label=dn(method), color=METHOD_COLORS[method],
                       edgecolor='white', linewidth=0.8, alpha=0.85)
                for i, val in enumerate(rates):
                    ax.text(x[i] + offset, val + 1, f'{val:.0f}%',
                            ha='center', va='bottom', fontsize=8, fontweight='bold')

            ax.set_xlabel(fam_info['label'], fontsize=12, fontweight='bold')
            ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
            ax.set_xticks(x)
            ax.set_xticklabels(LEVEL_ORDER, fontsize=11)
            ax.set_ylim(0, 120)
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            if fi == 0:
                ax.set_ylabel('Completion Rate (%)', fontsize=12, fontweight='bold')
                ax.legend(fontsize=9, framealpha=0.9)

        plt.tight_layout()
        self._save(fig, "01_completion_comparison")
        print("  [1] Completion comparison")

    # ── 2. single metric across conditions ────────────────────────────────────
    def plot_metric_across_conditions(self, metric_key):
        if self.comparisons.empty:
            return

        metric_label = METRICS.get(metric_key, metric_key)
        metric_data  = self.comparisons[
            (self.comparisons['metric'] == metric_key) &
            (self.comparisons['status'] == 'SUCCESS')
        ]
        if metric_data.empty:
            return

        use_median = (metric_key == 'num_rets_bias')
        n_fam  = len(CONFIG_FAMILIES)
        n_meth = len(METHODS)
        bar_w  = 0.7 / n_meth

        fig, axes = plt.subplots(1, n_fam, figsize=(5 * n_fam, 5),
                                 squeeze=False, sharey=True)
        axes = axes.flatten()

        for fi, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
            ax = axes[fi]
            x  = np.arange(len(LEVEL_ORDER))

            for mi, method in enumerate(METHODS):
                centers, errs = [], []
                for level in LEVEL_ORDER:
                    cfg = fam_info['configs'].get(level)
                    c, e = _agg(metric_data, method, cfg, use_median)
                    centers.append(c)
                    errs.append(e)

                offset = (mi - n_meth / 2 + 0.5) * bar_w
                ax.bar(x + offset, centers, bar_w * 0.9,
                       yerr=errs, capsize=3,
                       label=dn(method), color=METHOD_COLORS[method],
                       edgecolor='white', linewidth=0.8, alpha=0.85)

            if metric_key == 'num_rets_bias':
                ax.axhline(0, color='black', linewidth=1, linestyle='-', alpha=0.5)

            ax.set_xlabel(fam_info['label'], fontsize=12, fontweight='bold')
            ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
            ax.set_xticks(x)
            ax.set_xticklabels(LEVEL_ORDER, fontsize=11)
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            if fi == 0:
                ax.set_ylabel(metric_label, fontsize=12, fontweight='bold')
                ax.legend(fontsize=9, framealpha=0.9)

        plt.tight_layout()
        safe = metric_key.replace('.', '_')
        self._save(fig, f"02_{safe}_across_conditions")
        print(f"  [2] {metric_label}")

    # ── 3. combined panel ─────────────────────────────────────────────────────
    def plot_combined_panel(self):
        if self.comparisons.empty:
            return

        panel_metrics = [
            ('edit_distance_multree',  'Edit Distance',      None),
            ('ret_leaf_jaccard.dist',  'Ret. Descendants Measure', None),
            ('num_rets_bias',          'Ret. Bias (%)',       0),
        ]
        n_fam  = len(CONFIG_FAMILIES)
        n_meth = len(METHODS)
        bar_w  = 0.7 / n_meth

        fig, axes = plt.subplots(len(panel_metrics), n_fam,
                                 figsize=(5 * n_fam, 4.5 * len(panel_metrics)),
                                 squeeze=False, sharey='row')

        for ri, (metric_key, metric_label, ref_line) in enumerate(panel_metrics):
            metric_data = self.comparisons[
                (self.comparisons['metric'] == metric_key) &
                (self.comparisons['status'] == 'SUCCESS')
            ]
            use_median = (metric_key == 'num_rets_bias')

            for fi, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
                ax = axes[ri, fi]
                x  = np.arange(len(LEVEL_ORDER))

                for mi, method in enumerate(METHODS):
                    centers, errs = [], []
                    for level in LEVEL_ORDER:
                        cfg = fam_info['configs'].get(level)
                        c, e = _agg(metric_data, method, cfg, use_median)
                        centers.append(c)
                        errs.append(e)

                    offset = (mi - n_meth / 2 + 0.5) * bar_w
                    ax.bar(x + offset, centers, bar_w * 0.9,
                           yerr=errs, capsize=3,
                           label=dn(method) if ri == 0 else None,
                           color=METHOD_COLORS[method],
                           edgecolor='white', linewidth=0.8, alpha=0.85)

                if ref_line is not None:
                    ax.axhline(ref_line, color='black', linewidth=1, linestyle='-', alpha=0.5)

                ax.set_xticks(x)
                ax.set_xticklabels(LEVEL_ORDER, fontsize=10)
                ax.grid(True, alpha=0.25, linestyle='--', axis='y')
                if ri == 0:
                    ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
                if ri == len(panel_metrics) - 1:
                    ax.set_xlabel(fam_info['label'], fontsize=11, fontweight='bold')
                if fi == 0:
                    ax.set_ylabel(metric_label, fontsize=11, fontweight='bold')

        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', ncol=n_meth, fontsize=11,
                   framealpha=0.9, bbox_to_anchor=(0.5, 1.02))
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        self._save(fig, "03_combined_accuracy_panel")
        print("  [3] Combined accuracy panel")

    # ── 4. boxplots ───────────────────────────────────────────────────────────
    def plot_boxplots(self):
        if self.comparisons.empty:
            return

        box_metrics = [
            ('edit_distance_multree',   'Edit Distance'),
            ('ret_leaf_jaccard.dist',   'Ret. Descendants Measure'),
            ('ret_sisters_jaccard.dist','Ret. Sister Measure'),
            ('ploidy_diff.dist',        'Ploidy Distance'),
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
            data_list, labels, colors = [], [], []
            for method in METHODS:
                vals = metric_data[metric_data['method'] == method]['value'].dropna()
                if len(vals) > 0:
                    data_list.append(vals.values)
                    labels.append(dn(method))
                    colors.append(METHOD_COLORS[method])

            if data_list:
                bp = ax.boxplot(data_list, labels=labels, patch_artist=True,
                                widths=0.5, showmeans=True,
                                meanprops=dict(marker='D', markerfacecolor='white',
                                               markeredgecolor='black', markersize=6))
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)

            ax.set_title(metric_label, fontsize=13, fontweight='bold')
            ax.set_xticklabels(labels, rotation=15, ha='right')
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            ax.set_ylim(bottom=0)

        fig.suptitle('Overall Distribution (All Configurations)',
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        self._save(fig, "04_metric_boxplots")
        print("  [4] Metric boxplots")

    # ── 5. per-network ────────────────────────────────────────────────────────
    def plot_per_network(self):
        if self.comparisons.empty or self.network_stats is None:
            return

        baseline = 'conf_ils_low_10M'
        metric_key = 'edit_distance_multree'

        metric_data = self.comparisons[
            (self.comparisons['metric'] == metric_key) &
            (self.comparisons['status'] == 'SUCCESS') &
            (self.comparisons['config'] == baseline)
        ]
        if metric_data.empty:
            return

        agg = metric_data.groupby(['network', 'method'])['value'].agg(['mean', 'std']).reset_index()
        net_order = self.network_stats.sort_values('H_Strict')['network'].tolist()
        net_order = [n for n in net_order if n in agg['network'].unique()]

        n_meth = len(METHODS)
        bar_w  = 0.7 / n_meth

        fig, ax = plt.subplots(figsize=(max(14, len(net_order) * 0.9), 6))
        x = np.arange(len(net_order))

        for mi, method in enumerate(METHODS):
            method_agg = agg[agg['method'] == method].set_index('network')
            means = [method_agg.loc[n, 'mean'] if n in method_agg.index else np.nan for n in net_order]
            stds  = [method_agg.loc[n, 'std']  if n in method_agg.index else 0       for n in net_order]

            offset = (mi - n_meth / 2 + 0.5) * bar_w
            ax.bar(x + offset, means, bar_w * 0.9,
                   yerr=stds, capsize=2,
                   label=dn(method), color=METHOD_COLORS[method],
                   edgecolor='white', linewidth=0.5, alpha=0.85)

        h_vals = []
        for n in net_order:
            row = self.network_stats[self.network_stats['network'] == n]
            h_vals.append(int(row['H_Strict'].values[0]) if len(row) > 0 else '?')

        ax.set_xticks(x)
        ax.set_xticklabels([f'{n}\n(H={h})' for n, h in zip(net_order, h_vals)],
                           rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Edit Distance', fontsize=12, fontweight='bold')
        ax.set_title('Per-Network Edit Distance (Low ILS)', fontsize=14, fontweight='bold', pad=15)
        ax.legend(fontsize=10, framealpha=0.9, loc='upper left')
        ax.grid(True, alpha=0.25, linestyle='--', axis='y')
        ax.set_ylim(0, 1.05)

        plt.tight_layout()
        self._save(fig, "05_per_network_edit_distance")
        print("  [5] Per-network edit distance")

    # ── 6. degradation lines ──────────────────────────────────────────────────
    def plot_degradation_lines(self):
        if self.comparisons.empty:
            return

        line_metrics = [
            ('edit_distance_multree',    'Edit Distance'),
            ('ret_leaf_jaccard.dist',    'Reticulation Descendants Measure'),
            ('ret_sisters_jaccard.dist', 'Reticulation Sister Measure'),
        ]
        n_fam = len(CONFIG_FAMILIES)

        fig, axes = plt.subplots(len(line_metrics), n_fam,
                                 figsize=(5 * n_fam, 4 * len(line_metrics)),
                                 squeeze=False, sharey='row')

        for ri, (metric_key, metric_label) in enumerate(line_metrics):
            metric_data = self.comparisons[
                (self.comparisons['metric'] == metric_key) &
                (self.comparisons['status'] == 'SUCCESS')
            ]
            use_median = (metric_key == 'num_rets_bias')

            for fi, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
                ax = axes[ri, fi]

                for method in METHODS:
                    centers, errs = [], []
                    for level in LEVEL_ORDER:
                        cfg = fam_info['configs'].get(level)
                        c, e = _agg(metric_data, method, cfg, use_median)
                        centers.append(c)
                        errs.append(e)

                    ax.errorbar(LEVEL_ORDER, centers, yerr=errs,
                                marker='o', markersize=8, capsize=4,
                                label=dn(method) if ri == 0 and fi == 0 else None,
                                color=METHOD_COLORS[method], linewidth=2.5)

                if metric_key == 'num_rets_bias':
                    ax.axhline(0, color='black', linewidth=1, linestyle='--', alpha=0.5)

                ax.grid(True, alpha=0.25, linestyle='--')
                if ri == 0:
                    ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
                if ri == len(line_metrics) - 1:
                    ax.set_xlabel(fam_info['label'], fontsize=11, fontweight='bold')
                if fi == 0:
                    ax.set_ylabel(metric_label, fontsize=11, fontweight='bold')

        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', ncol=n_fam, fontsize=11,
                   framealpha=0.9, bbox_to_anchor=(0.5, 1.02))
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        self._save(fig, "06_degradation_lines")
        print("  [6] Degradation lines")

    # ── summary table ─────────────────────────────────────────────────────────
    def generate_summary_table(self):
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
                    row[f'{metric_key}_mean']   = vals.mean()   if len(vals) > 0 else np.nan
                    row[f'{metric_key}_median'] = vals.median() if len(vals) > 0 else np.nan
                    row[f'{metric_key}_std']    = vals.std()    if len(vals) > 1 else np.nan
                rows.append(row)

        pd.DataFrame(rows).to_csv(self.tables_dir / "01_prior_comparison_summary.csv", index=False)
        print("  [T] Summary table saved")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--configs', nargs='+', default=None)
    parser.add_argument('--output',  default=None)
    parser.add_argument('--network-stats', default=None)
    args = parser.parse_args()

    configs    = args.configs or ALL_CONFIGS
    output_dir = Path(args.output) if args.output else SUMMARY_BASE / "prior_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)

    ns_path = Path(args.network_stats) if args.network_stats else NETWORKS_DIR / "mul_tree_final_stats.csv"
    network_stats = None
    if ns_path.exists():
        network_stats = pd.read_csv(ns_path)
        network_stats['network'] = network_stats['Filename'].str.replace('.tre', '', regex=False)

    print(f"Loading data from: {SUMMARY_BASE}")
    inventory, comparisons = load_data(configs, SUMMARY_BASE)

    analyzer = PriorComparison(inventory, comparisons, output_dir, network_stats)
    analyzer.generate_all()


if __name__ == '__main__':
    main()
