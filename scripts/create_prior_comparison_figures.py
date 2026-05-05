#!/usr/bin/env python3
"""
create_prior_comparison_figures.py - Comparison of GRAMPA^Iter with and without ploidy prior

Focused analysis for the thesis subsection comparing:
  - grandma_split_prior  (GRAMPA^Iter with prior)
  - grandma_split        (GRAMPA^Iter, no prior)
  - polyphest            (reference method)

Scans the papers/ directory directly (no papers_config.yaml entry needed).

Usage:
    python scripts/create_prior_comparison_figures.py
    python scripts/create_prior_comparison_figures.py --force-recompute
    python scripts/create_prior_comparison_figures.py --output analysis/prior_comparison
"""

import argparse
import sys
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# ── paths ──────────────────────────────────────────────────────────────────────
SCRIPTS_DIR = Path(__file__).parent.resolve()
ROOT_DIR    = SCRIPTS_DIR.parent
PAPERS_DIR  = ROOT_DIR / "papers"

sys.path.insert(0, str(SCRIPTS_DIR))
sys.path.insert(0, str(ROOT_DIR / "simulations" / "scripts"))


def _import(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


compute_comparisons = _import('compute_comparisons', SCRIPTS_DIR / 'compute_comparisons.py')
ComparisonEngine    = compute_comparisons.ComparisonEngine

# ── constants ──────────────────────────────────────────────────────────────────
METHODS = {
    'grandma_split_prior': 'final_multree.tre',
    'grandma_split':       'final_multree.tre',
    'polyphest':           'network.tre',
}

METHOD_LABELS = {
    'grandma_split_prior': r'GRAMPA$^{Iter}$ (prior)',
    'grandma_split':       r'GRAMPA$^{Iter}$',
    'polyphest':           'Polyphest',
}

METHOD_COLORS = {
    'grandma_split_prior': '#9467BD',  # deep purple
    'grandma_split':       '#CC78BC',  # muted purple
    'polyphest':           '#DE8F05',  # orange
}

METRICS = [
    'edit_distance_multree',
    'num_rets_diff',
    'polyploid_species_jaccard',
]

METRIC_LABELS = {
    'edit_distance_multree':    'Edit Distance',
    'num_rets_diff':            'Reticulation Count Difference (abs)',
    'polyploid_species_jaccard':'Polyploid Species Distance (Jaccard)',
}

plt.rcParams.update({
    'figure.dpi':        300,
    'savefig.dpi':       300,
    'font.size':         11,
    'axes.labelsize':    13,
    'axes.titlesize':    14,
    'xtick.labelsize':   11,
    'ytick.labelsize':   11,
    'legend.fontsize':   11,
    'font.family':       'sans-serif',
    'axes.linewidth':    1.2,
    'grid.linewidth':    0.8,
    'lines.linewidth':   2.5,
    'lines.markersize':  8,
})


# ── inventory ──────────────────────────────────────────────────────────────────

def build_inventory(papers_dir: Path) -> pd.DataFrame:
    """Scan papers/ and build a 3-method inventory."""
    rows = []
    for paper_dir in sorted(papers_dir.iterdir()):
        if not paper_dir.is_dir():
            continue
        network = paper_dir.name
        for method, output_file in METHODS.items():
            result_path = paper_dir / "networks" / method / output_file
            size = result_path.stat().st_size if result_path.exists() else 0
            rows.append({
                'network':      network,
                'method':       method,
                'network_path': str(result_path),
                'exists':       result_path.exists() and size > 0,
                'file_size':    size,
            })
    return pd.DataFrame(rows)


# ── helpers ────────────────────────────────────────────────────────────────────

def _label(method: str) -> str:
    return METHOD_LABELS.get(method, method)


def _color(method: str) -> str:
    return METHOD_COLORS.get(method, '#CCCCCC')


def _save(fig, plots_dir: Path, stem: str):
    fig.savefig(plots_dir / f"{stem}.pdf", bbox_inches='tight')
    fig.savefig(plots_dir / f"{stem}.png", bbox_inches='tight', dpi=300)
    plt.close(fig)


# ── figures ────────────────────────────────────────────────────────────────────

def plot_completion(inventory: pd.DataFrame, plots_dir: Path):
    """Bar chart: how many networks each method completed."""
    total = inventory['network'].nunique()
    agg = (
        inventory.groupby('method')['exists']
        .sum()
        .reset_index()
        .rename(columns={'exists': 'available'})
    )
    agg['pct'] = agg['available'] / total * 100
    agg = agg.set_index('method').reindex(METHODS.keys()).reset_index()

    fig, ax = plt.subplots(figsize=(7, 5))
    colors = [_color(m) for m in agg['method']]
    bars   = ax.bar(
        [_label(m) for m in agg['method']],
        agg['pct'],
        color=colors, alpha=0.85, edgecolor='black', linewidth=1.2,
    )
    for bar, row in zip(bars, agg.itertuples()):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 1,
            f"{int(row.available)}/{total}\n({row.pct:.0f}%)",
            ha='center', va='bottom', fontsize=10, fontweight='bold',
        )
    ax.set_ylim(0, 115)
    ax.set_ylabel('Completion Rate (%)', fontweight='bold')
    ax.set_title('Completion Rate by Method', fontweight='bold', pad=12)
    ax.grid(True, axis='y', alpha=0.25, linestyle='--')
    plt.tight_layout()
    _save(fig, plots_dir, "01_completion_rate")


def plot_availability_heatmap(inventory: pd.DataFrame, plots_dir: Path):
    """Network × method availability matrix."""
    pivot = inventory.pivot_table(
        index='network', columns='method', values='exists', aggfunc='first'
    ).astype(float)
    # reorder columns
    pivot = pivot.reindex(columns=[m for m in METHODS if m in pivot.columns])
    pivot.columns = [_label(c) for c in pivot.columns]

    fig, ax = plt.subplots(figsize=(8, max(8, len(pivot) * 0.35)))
    sns.heatmap(pivot, annot=True, fmt='.0f', cmap='RdYlGn', vmin=0, vmax=1,
                linewidths=0.5, linecolor='gray', ax=ax,
                cbar_kws={'label': 'Available (1) / Missing (0)'})
    ax.set_xlabel('Method', fontweight='bold')
    ax.set_ylabel('Network', fontweight='bold')
    ax.set_title('Result Availability Matrix', fontweight='bold', pad=12)
    plt.tight_layout()
    _save(fig, plots_dir, "02_availability_heatmap")


def plot_pairwise_heatmap(comparisons: pd.DataFrame, metric: str,
                          plots_dir: Path, fig_prefix: str):
    """Symmetric method × method heatmap averaged over networks."""
    data = comparisons[comparisons['metric'] == metric]
    if data.empty:
        print(f"  No data for metric {metric}, skipping heatmap.")
        return

    methods = sorted(METHODS.keys())
    n = len(methods)
    mat = np.full((n, n), np.nan)
    idx = {m: i for i, m in enumerate(methods)}

    for _, row in data.iterrows():
        m1, m2 = row['method1'], row['method2']
        if m1 in idx and m2 in idx:
            i, j = idx[m1], idx[m2]
            current = mat[i, j]
            mat[i, j] = row['value'] if np.isnan(current) else (current + row['value']) / 2
            mat[j, i] = mat[i, j]
    np.fill_diagonal(mat, 0)

    labels = [_label(m) for m in methods]
    df_mat  = pd.DataFrame(mat, index=labels, columns=labels)

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(df_mat, annot=True, fmt='.3f', cmap='YlGnBu_r', vmin=0,
                linewidths=1, linecolor='black', square=True, ax=ax,
                cbar_kws={'label': METRIC_LABELS.get(metric, metric)})
    ax.set_title(
        f"{METRIC_LABELS.get(metric, metric)}\n(lower = more similar)",
        fontweight='bold', pad=12,
    )
    plt.tight_layout()
    _save(fig, plots_dir, f"{fig_prefix}_heatmap_{metric.replace('.', '_')}")


def plot_per_network_metric(comparisons: pd.DataFrame, metric: str,
                            plots_dir: Path, fig_prefix: str):
    """Grouped bar chart: metric value per network, one bar per method pair."""
    data = comparisons[comparisons['metric'] == metric].copy()
    if data.empty:
        return

    data['pair'] = data.apply(
        lambda r: f"{_label(r['method1'])} vs {_label(r['method2'])}", axis=1
    )
    networks = sorted(data['network'].unique())
    pairs    = sorted(data['pair'].unique())
    x        = np.arange(len(networks))
    width    = 0.8 / max(len(pairs), 1)

    pair_colors = plt.cm.tab10(np.linspace(0, 0.6, len(pairs)))

    fig, ax = plt.subplots(figsize=(max(12, len(networks) * 0.55), 5))
    for i, (pair, color) in enumerate(zip(pairs, pair_colors)):
        vals = []
        for net in networks:
            row = data[(data['network'] == net) & (data['pair'] == pair)]
            vals.append(row['value'].values[0] if not row.empty else np.nan)
        offset = (i - len(pairs) / 2 + 0.5) * width
        ax.bar(x + offset, vals, width * 0.9, label=pair, color=color,
               alpha=0.8, edgecolor='black', linewidth=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(networks, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel(METRIC_LABELS.get(metric, metric), fontweight='bold')
    ax.set_title(f"{METRIC_LABELS.get(metric, metric)} per Dataset", fontweight='bold', pad=12)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, axis='y', alpha=0.25, linestyle='--')
    plt.tight_layout()
    _save(fig, plots_dir, f"{fig_prefix}_per_network_{metric.replace('.', '_')}")


def plot_reticulation_counts(inventory: pd.DataFrame, plots_dir: Path):
    """
    Dot/strip plot showing how many reticulations each method inferred per network.
    Reads reticulation counts from the MUL-tree files via ReticulateTree.
    """
    try:
        from reticulate_tree import ReticulateTree
    except ImportError:
        print("  reticulate_tree not importable; skipping reticulation count plot.")
        return

    records = []
    for _, row in inventory[inventory['exists']].iterrows():
        try:
            with open(row['network_path']) as f:
                newick = f.read().strip()
            if newick and not newick.endswith(';'):
                newick += ';'
            tree = ReticulateTree(newick)
            n_rets = tree.count_reticulations() if hasattr(tree, 'count_reticulations') else None
            if n_rets is None:
                # fall back to leaf count proxy
                counts = tree.get_leaf_counts()
                n_rets = sum(1 for c in counts.values() if c > 1)
            records.append({'network': row['network'], 'method': row['method'], 'n_rets': n_rets})
        except Exception:
            continue

    if not records:
        print("  Could not extract reticulation counts; skipping plot.")
        return

    df = pd.DataFrame(records)
    networks = sorted(df['network'].unique())

    fig, ax = plt.subplots(figsize=(max(14, len(networks) * 0.6), 5))
    x     = np.arange(len(networks))
    width = 0.25

    for i, method in enumerate(METHODS.keys()):
        mdf  = df[df['method'] == method].set_index('network')
        vals = [mdf.loc[n, 'n_rets'] if n in mdf.index else np.nan for n in networks]
        offset = (i - 1) * width
        ax.bar(x + offset, vals, width * 0.9, label=_label(method),
               color=_color(method), alpha=0.85, edgecolor='black', linewidth=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(networks, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Reticulations Inferred', fontweight='bold')
    ax.set_title('Reticulation Count per Dataset', fontweight='bold', pad=12)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.legend(fontsize=10)
    ax.grid(True, axis='y', alpha=0.25, linestyle='--')
    plt.tight_layout()
    _save(fig, plots_dir, "04_reticulation_counts")


def plot_summary_boxplot(comparisons: pd.DataFrame, metric: str,
                         plots_dir: Path, fig_prefix: str):
    """Boxplot of metric values (one box per method pair, across networks)."""
    data = comparisons[comparisons['metric'] == metric].copy()
    if data.empty:
        return

    data['pair'] = data.apply(
        lambda r: f"{_label(r['method1'])}\nvs\n{_label(r['method2'])}", axis=1
    )
    pairs = sorted(data['pair'].unique())
    vals  = [data[data['pair'] == p]['value'].dropna().values for p in pairs]

    fig, ax = plt.subplots(figsize=(max(8, len(pairs) * 2.5), 5))
    bp = ax.boxplot(vals, labels=pairs, patch_artist=True,
                    widths=0.5, showfliers=True,
                    boxprops=dict(linewidth=1.5),
                    whiskerprops=dict(linewidth=1.5),
                    capprops=dict(linewidth=1.5),
                    medianprops=dict(linewidth=2.5, color='red'))
    colors = plt.cm.tab10(np.linspace(0, 0.6, len(pairs)))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel(METRIC_LABELS.get(metric, metric), fontweight='bold')
    ax.set_title(f"Distribution of {METRIC_LABELS.get(metric, metric)}", fontweight='bold', pad=12)
    ax.grid(True, axis='y', alpha=0.25, linestyle='--')
    plt.tight_layout()
    _save(fig, plots_dir, f"{fig_prefix}_boxplot_{metric.replace('.', '_')}")


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--papers-dir', default=str(PAPERS_DIR),
                        help='Path to papers/ directory')
    parser.add_argument('--output', default=str(ROOT_DIR / "analysis" / "real_data_summary" / "prior_comparison"),
                        help='Output directory for figures and tables')
    parser.add_argument('--force-recompute', action='store_true',
                        help='Ignore cache and recompute all comparisons')
    args = parser.parse_args()

    papers_dir = Path(args.papers_dir)
    output_dir = Path(args.output)
    plots_dir  = output_dir / "plots"
    tables_dir = output_dir / "tables"
    cache_dir  = output_dir / "cache"
    for d in (plots_dir, tables_dir, cache_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*70}")
    print(f"  Prior Comparison Analysis")
    print(f"  Papers dir : {papers_dir}")
    print(f"  Output dir : {output_dir}")
    print(f"{'='*70}\n")

    # ── 1. inventory ──────────────────────────────────────────────────────────
    print("Building inventory...")
    inventory = build_inventory(papers_dir)
    inventory.to_csv(tables_dir / "inventory.csv", index=False)

    print("\nCompletion summary:")
    for method in METHODS:
        sub = inventory[inventory['method'] == method]
        n   = sub['exists'].sum()
        tot = len(sub)
        print(f"  {_label(method):35s}: {n}/{tot}")

    # ── 2. pairwise comparisons ───────────────────────────────────────────────
    print("\nComputing pairwise comparisons...")
    engine      = ComparisonEngine(str(cache_dir), force_recompute=args.force_recompute)
    comparisons = engine.compute_all_comparisons(inventory)

    if not comparisons.empty:
        comparisons.to_csv(tables_dir / "comparisons.csv", index=False)
        engine.print_statistics()
    else:
        print("  WARNING: No comparisons computed.")

    valid = comparisons[comparisons['status'] == 'SUCCESS'].copy() if not comparisons.empty else pd.DataFrame()

    # ── 3. figures ────────────────────────────────────────────────────────────
    print("\nGenerating figures...")

    print("  01  completion rate bar chart")
    plot_completion(inventory, plots_dir)

    print("  02  availability heatmap")
    plot_availability_heatmap(inventory, plots_dir)

    for metric in METRICS:
        prefix = "03"
        print(f"  {prefix}  pairwise heatmap — {metric}")
        plot_pairwise_heatmap(valid, metric, plots_dir, prefix)

    print("  04  reticulation counts per dataset")
    plot_reticulation_counts(inventory, plots_dir)

    for metric in METRICS:
        prefix = "05"
        print(f"  {prefix}  per-network bars — {metric}")
        plot_per_network_metric(valid, metric, plots_dir, prefix)

    for metric in METRICS:
        prefix = "06"
        print(f"  {prefix}  summary boxplot — {metric}")
        plot_summary_boxplot(valid, metric, plots_dir, prefix)

    # ── 4. summary table ──────────────────────────────────────────────────────
    if not valid.empty:
        summary = (
            valid.groupby(['method1', 'method2', 'metric'])['value']
            .agg(['mean', 'std', 'count'])
            .reset_index()
        )
        summary.to_csv(tables_dir / "pairwise_summary.csv", index=False)
        print(f"\nSummary table → {tables_dir / 'pairwise_summary.csv'}")

    print(f"\nDone. Figures saved to {plots_dir}\n")


if __name__ == '__main__':
    main()
