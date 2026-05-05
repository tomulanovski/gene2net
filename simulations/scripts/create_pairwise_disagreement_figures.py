#!/usr/bin/env python3
"""
create_pairwise_disagreement_figures.py

Plots pairwise disagreement between Polyphest and grandma_split across all
simulation configurations. Uses Option-B aggregation: replicate-mean per
network, then mean across networks (so every network has equal weight).

Input:
    simulations/analysis/summary/{config}/pairwise_polyphest_vs_grandma.csv

Output:
    simulations/analysis/summary/pairwise_disagreement/
        pairwise_disagreement.pdf / .png
        pairwise_disagreement_per_metric.pdf / .png
        pairwise_aggregated.csv

Usage:
    python create_pairwise_disagreement_figures.py
    python create_pairwise_disagreement_figures.py --empirical 0.72 0.55 0.61
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

SCRIPT_DIR = Path(__file__).resolve().parent
SUMMARY_BASE = SCRIPT_DIR.parent / "analysis" / "summary"
OUTPUT_DIR = SUMMARY_BASE / "pairwise_disagreement"

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

TARGET_METRICS = {
    'edit_distance_multree': 'Edit Distance',
    'ret_leaf_jaccard.dist': 'Reticulation Leaf Distance',
    'ret_sisters_jaccard.dist': 'Sister-Taxa Distance',
}

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.2


# ─── Data loading ────────────────────────────────────────────────────────────

def load_all_configs(configs):
    frames = []
    for config in configs:
        path = SUMMARY_BASE / config / "pairwise_polyphest_vs_grandma.csv"
        if not path.exists():
            print(f"  WARNING: {config} — file not found, skipping")
            continue
        df = pd.read_csv(path)
        frames.append(df)
    if not frames:
        print("ERROR: No data loaded. Run compute_pairwise_methods.py first.")
        sys.exit(1)
    return pd.concat(frames, ignore_index=True)


def aggregate_option_b(df):
    """
    Option B: for each (config, network, metric) average replicates,
    then aggregate across networks → per-config mean ± SEM.
    Returns a DataFrame with columns: config, metric, mean, sem, n_networks.
    """
    df = df[(df['status'] == 'SUCCESS') & (df['metric'].isin(TARGET_METRICS))].copy()

    # Step 1: replicate mean per network
    per_net = (
        df.groupby(['config', 'network', 'metric'])['value']
        .mean()
        .reset_index()
        .rename(columns={'value': 'net_mean'})
    )

    # Step 2: mean ± SEM across networks
    def summarise(g):
        n = len(g)
        m = g['net_mean'].mean()
        s = g['net_mean'].std(ddof=1) / np.sqrt(n) if n > 1 else 0.0
        return pd.Series({'mean': m, 'sem': s, 'n_networks': n})

    agg = (
        per_net.groupby(['config', 'metric'])
        .apply(summarise)
        .reset_index()
    )
    return agg


# ─── Figures ─────────────────────────────────────────────────────────────────

def plot_combined(agg, empirical, output_dir):
    """
    3 rows (metrics) × 3 columns (config families).
    Each panel: 3 bars (Low/Medium/High) + optional empirical reference line.
    """
    n_metrics = len(TARGET_METRICS)
    n_families = len(CONFIG_FAMILIES)
    bar_color = '#4878CF'

    fig, axes = plt.subplots(
        n_metrics, n_families,
        figsize=(5 * n_families, 4.5 * n_metrics),
        squeeze=False,
        sharey='row',
    )

    for row_idx, (metric_key, metric_label) in enumerate(TARGET_METRICS.items()):
        metric_data = agg[agg['metric'] == metric_key]

        for col_idx, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
            ax = axes[row_idx, col_idx]
            x = np.arange(len(LEVEL_ORDER))

            means, sems = [], []
            for level in LEVEL_ORDER:
                cfg = fam_info['configs'].get(level)
                row = metric_data[metric_data['config'] == cfg]
                if len(row) > 0:
                    means.append(row.iloc[0]['mean'])
                    sems.append(row.iloc[0]['sem'])
                else:
                    means.append(np.nan)
                    sems.append(0.0)

            ax.bar(x, means, yerr=sems, capsize=4,
                   color=bar_color, alpha=0.82,
                   edgecolor='white', linewidth=0.8, width=0.6)

            # Empirical reference line
            if empirical and metric_key in empirical:
                emp_val = empirical[metric_key]
                ax.axhline(emp_val, color='#D62728', linewidth=1.8,
                           linestyle='--', label='Empirical', zorder=5)
                if row_idx == 0 and col_idx == 0:
                    ax.legend(fontsize=9, framealpha=0.9)

            ax.set_xticks(x)
            ax.set_xticklabels(LEVEL_ORDER, fontsize=10)
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            ax.set_ylim(bottom=0)

            if row_idx == 0:
                ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
            if row_idx == n_metrics - 1:
                ax.set_xlabel(fam_info['label'], fontsize=11, fontweight='bold')
            if col_idx == 0:
                ax.set_ylabel(metric_label, fontsize=11, fontweight='bold')

    fig.suptitle(
        r'Polyphest vs GRAMPA$^{Iter}$ — Pairwise Disagreement by Simulation Config',
        fontsize=13, fontweight='bold', y=1.01,
    )
    plt.tight_layout()
    fig.savefig(output_dir / "pairwise_disagreement.pdf", bbox_inches='tight')
    fig.savefig(output_dir / "pairwise_disagreement.png", bbox_inches='tight', dpi=300)
    plt.close('all')
    gc.collect()
    print("  [1] Combined figure saved")


def plot_per_metric(agg, empirical, output_dir):
    """One figure per metric — larger, easier to read."""
    for metric_key, metric_label in TARGET_METRICS.items():
        metric_data = agg[agg['metric'] == metric_key]
        n_families = len(CONFIG_FAMILIES)
        bar_color = '#4878CF'

        fig, axes = plt.subplots(1, n_families,
                                 figsize=(5 * n_families, 5),
                                 squeeze=False, sharey=True)
        axes = axes.flatten()

        for col_idx, (fam_name, fam_info) in enumerate(CONFIG_FAMILIES.items()):
            ax = axes[col_idx]
            x = np.arange(len(LEVEL_ORDER))

            means, sems, n_nets = [], [], []
            for level in LEVEL_ORDER:
                cfg = fam_info['configs'].get(level)
                row = metric_data[metric_data['config'] == cfg]
                if len(row) > 0:
                    means.append(row.iloc[0]['mean'])
                    sems.append(row.iloc[0]['sem'])
                    n_nets.append(int(row.iloc[0]['n_networks']))
                else:
                    means.append(np.nan)
                    sems.append(0.0)
                    n_nets.append(0)

            bars = ax.bar(x, means, yerr=sems, capsize=4,
                          color=bar_color, alpha=0.82,
                          edgecolor='white', linewidth=0.8, width=0.6)

            # Value labels
            for i, (m, n) in enumerate(zip(means, n_nets)):
                if not np.isnan(m):
                    ax.text(x[i], m + sems[i] + 0.01, f'{m:.3f}',
                            ha='center', va='bottom', fontsize=8)

            if empirical and metric_key in empirical:
                emp_val = empirical[metric_key]
                ax.axhline(emp_val, color='#D62728', linewidth=2,
                           linestyle='--', label=f'Empirical ({emp_val:.3f})', zorder=5)
                ax.legend(fontsize=9, framealpha=0.9)

            ax.set_xticks(x)
            ax.set_xticklabels(LEVEL_ORDER, fontsize=11)
            ax.set_title(fam_name, fontsize=13, fontweight='bold', pad=10)
            ax.set_xlabel(fam_info['label'], fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            ax.set_ylim(bottom=0)
            if col_idx == 0:
                ax.set_ylabel(metric_label, fontsize=12, fontweight='bold')

        fig.suptitle(
            fr'Polyphest vs GRAMPA${{^{{Iter}}}}$ — {metric_label}',
            fontsize=13, fontweight='bold', y=1.02,
        )
        plt.tight_layout()
        safe = metric_key.replace('.', '_')
        fig.savefig(output_dir / f"per_metric_{safe}.pdf", bbox_inches='tight')
        fig.savefig(output_dir / f"per_metric_{safe}.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()
        print(f"  [2] Per-metric figure: {metric_label}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Pairwise Polyphest vs grandma_split disagreement figures'
    )
    parser.add_argument('--configs', nargs='+', default=None,
                        help='Configs to include (default: all 9)')
    parser.add_argument(
        '--empirical', nargs=3, type=float, metavar=('EDIT', 'RET_LEAF', 'RET_SISTERS'),
        default=None,
        help='Empirical pairwise values for edit_distance, ret_leaf, ret_sisters '
             '(drawn as red dashed reference lines)'
    )
    args = parser.parse_args()

    configs = args.configs or ALL_CONFIGS
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    empirical = None
    if args.empirical:
        keys = list(TARGET_METRICS.keys())
        empirical = dict(zip(keys, args.empirical))
        print(f"Empirical reference: {empirical}")

    print(f"Loading data for {len(configs)} configs...")
    raw = load_all_configs(configs)
    print(f"  Loaded {len(raw)} rows")

    print("Aggregating (Option B)...")
    agg = aggregate_option_b(raw)

    # Save aggregated table
    agg.to_csv(OUTPUT_DIR / "pairwise_aggregated.csv", index=False)
    print(f"  Aggregated table: {OUTPUT_DIR / 'pairwise_aggregated.csv'}")

    print("\nGenerating figures...")
    plot_combined(agg, empirical, OUTPUT_DIR)
    plot_per_metric(agg, empirical, OUTPUT_DIR)

    print(f"\nDone. Output: {OUTPUT_DIR}")

    # Print quick summary table to stdout
    print("\nAggregated means per config × metric:")
    pivot = agg.pivot_table(index='config', columns='metric', values='mean')
    pivot = pivot[[m for m in TARGET_METRICS if m in pivot.columns]]
    pivot.columns = [TARGET_METRICS[m] for m in pivot.columns]
    print(pivot.round(4).to_string())


if __name__ == '__main__':
    main()
