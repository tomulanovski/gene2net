#!/usr/bin/env python3
"""
create_real_vs_inferred_figures.py

Compare REAL (ground-truth) ploidy prior vs INFERRED ploidy prior for the two
prior-using methods, per configuration:
  - Polyphest:  polyphest (inferred)            vs  polyphest_real
  - GRANDMA:    grandma_split_prior (inferred)  vs  grandma_split_prior_real

Polyphest thresholds (p50/p70/p90) are collapsed to a single result per
(network, replicate) using the lowest completed threshold -- the same rule used
by create_analysis_figures.py -- applied independently to the inferred and the
real runs.

Metrics (lower is better for all three):
  edit_distance_multree, ret_leaf_jaccard.dist, ret_sisters_jaccard.dist

For each (method-pair, metric) the unit is one network: the mean over replicates
(SUCCESS only). A paired scatter (x = inferred, y = real, with the y=x diagonal)
shows how much the real prior helps -- points below the diagonal = real better.

Reads:  simulations/analysis/summary/<config>/comparisons_raw.csv
Writes: simulations/analysis/summary/<config>/plots/real_vs_inferred_<method>.{pdf,png}
        simulations/analysis/summary/<config>/tables/real_vs_inferred_summary.csv

Usage:
    python create_real_vs_inferred_figures.py --config conf_dup_loss_low_10M
    python create_real_vs_inferred_figures.py --config conf_dup_loss_low_10M conf_dup_loss_high_10M
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend (no X11 needed on cluster)
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

METRICS = ['edit_distance_multree', 'ret_leaf_jaccard.dist', 'ret_sisters_jaccard.dist']
METRIC_LABELS = {
    'edit_distance_multree': 'Edit Distance (MUL-tree)',
    'ret_leaf_jaccard.dist': 'Ret. Descendants (Jaccard dist.)',
    'ret_sisters_jaccard.dist': 'Ret. Sister (Jaccard dist.)',
}

POLYPHEST_THRESHOLDS = ['polyphest_p50', 'polyphest_p70', 'polyphest_p90']
POLYPHEST_REAL_THRESHOLDS = ['polyphest_real_p50', 'polyphest_real_p70', 'polyphest_real_p90']

# (figure label, inferred method, real method, threshold lists or None)
METHOD_PAIRS = [
    {
        'label': 'Polyphest',
        'slug': 'polyphest',
        'inferred': 'polyphest', 'inferred_thresholds': POLYPHEST_THRESHOLDS,
        'real': 'polyphest_real', 'real_thresholds': POLYPHEST_REAL_THRESHOLDS,
    },
    {
        'label': r'GRAMPA$^{Iter}$ (prior)',
        'slug': 'grandma',
        'inferred': 'grandma_split_prior', 'inferred_thresholds': None,
        'real': 'grandma_split_prior_real', 'real_thresholds': None,
    },
]

INFERRED_COLOR = '#0173B2'
REAL_COLOR = '#DE8F05'


# ---------------------------------------------------------------------------
# Threshold collapse (lowest completed threshold per network, replicate)
# ---------------------------------------------------------------------------

def collapse_thresholds(df: pd.DataFrame, thresholds: list, target_name: str) -> pd.DataFrame:
    """Collapse a list of threshold-methods into one, keeping the rows of the
    lowest threshold that has a SUCCESS for each (network, replicate)."""
    poly = df[df['method'].isin(thresholds)].copy()
    if poly.empty:
        return pd.DataFrame(columns=df.columns)

    group_cols = ['network', 'replicate']
    if 'config' in poly.columns:
        group_cols = ['config', 'network', 'replicate']

    kept = []
    for key, group in poly.groupby(group_cols):
        chosen = None
        for method in thresholds:
            method_rows = group[group['method'] == method]
            if len(method_rows) > 0 and (method_rows['status'] == 'SUCCESS').any():
                chosen = method
                break
        if chosen is None:
            continue
        rows = group[group['method'] == chosen].copy()
        rows['method'] = target_name
        kept.append(rows)

    if kept:
        return pd.concat(kept, ignore_index=True)
    return pd.DataFrame(columns=df.columns)


def method_rows(comparisons: pd.DataFrame, method: str, thresholds) -> pd.DataFrame:
    """Return the SUCCESS rows for a method, collapsing thresholds if given."""
    if thresholds is not None:
        df = collapse_thresholds(comparisons, thresholds, method)
    else:
        df = comparisons[comparisons['method'] == method].copy()
    return df[df['status'] == 'SUCCESS']


def per_network_mean(rows: pd.DataFrame, metric: str) -> pd.Series:
    """Mean metric value per network (over replicates), SUCCESS only."""
    sel = rows[rows['metric'] == metric]
    if sel.empty:
        return pd.Series(dtype=float)
    return sel.groupby('network')['value'].mean()


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_pair(config: str, pair: dict, comparisons: pd.DataFrame,
              plots_dir: Path) -> list:
    """Create the paired-scatter figure for one method pair. Returns summary rows."""
    inf_rows = method_rows(comparisons, pair['inferred'], pair['inferred_thresholds'])
    real_rows = method_rows(comparisons, pair['real'], pair['real_thresholds'])

    summary = []
    fig, axes = plt.subplots(1, len(METRICS), figsize=(6 * len(METRICS), 5.5))
    if len(METRICS) == 1:
        axes = [axes]

    any_data = False
    for ax, metric in zip(axes, METRICS):
        inf_mean = per_network_mean(inf_rows, metric)
        real_mean = per_network_mean(real_rows, metric)

        # Pair on networks present in BOTH
        common = inf_mean.index.intersection(real_mean.index)
        x = inf_mean.loc[common].values
        y = real_mean.loc[common].values
        n = len(common)

        n_real_better = int((y < x).sum())
        n_real_worse = int((y > x).sum())
        n_tie = int((y == x).sum())

        summary.append({
            'config': config,
            'method_pair': pair['slug'],
            'metric': metric,
            'n_networks': n,
            'mean_inferred': float(inf_mean.loc[common].mean()) if n else float('nan'),
            'mean_real': float(real_mean.loc[common].mean()) if n else float('nan'),
            'delta_inferred_minus_real': (float(inf_mean.loc[common].mean() - real_mean.loc[common].mean())
                                          if n else float('nan')),
            'n_real_better': n_real_better,
            'n_real_worse': n_real_worse,
            'n_tie': n_tie,
        })

        if n > 0:
            any_data = True
            ax.scatter(x, y, s=70, alpha=0.75, color=REAL_COLOR,
                       edgecolors='black', linewidths=0.6, zorder=3)
            lo = min(x.min(), y.min())
            hi = max(x.max(), y.max())
            pad = (hi - lo) * 0.05 if hi > lo else 0.5
            lim_lo, lim_hi = lo - pad, hi + pad
            ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], '--', color='gray',
                    linewidth=1.5, zorder=1, label='y = x (no change)')
            ax.set_xlim(lim_lo, lim_hi)
            ax.set_ylim(lim_lo, lim_hi)
            ax.set_aspect('equal', adjustable='box')
            ax.text(0.04, 0.96,
                    f"n = {n}\nreal better: {n_real_better}\nreal worse: {n_real_worse}",
                    transform=ax.transAxes, va='top', ha='left', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        else:
            ax.text(0.5, 0.5, 'no paired data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=12, color='gray')

        ax.set_xlabel(f"Inferred prior")
        ax.set_ylabel(f"Real prior")
        ax.set_title(METRIC_LABELS.get(metric, metric), fontsize=11)
        ax.legend(loc='lower right', fontsize=8)

    fig.suptitle(f"{pair['label']}: real vs inferred ploidy prior  —  {config}\n"
                 f"(below the diagonal = real prior is better; lower is better)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    if any_data:
        for ext in ('pdf', 'png'):
            out = plots_dir / f"real_vs_inferred_{pair['slug']}.{ext}"
            fig.savefig(out, dpi=200, bbox_inches='tight')
        print(f"  wrote real_vs_inferred_{pair['slug']}.pdf/png")
    else:
        print(f"  [skip] {pair['slug']}: no paired data in {config}")
    plt.close(fig)
    return summary


def process_config(config: str, summary_base: Path):
    base = summary_base / config
    comp_path = base / 'comparisons_raw.csv'
    if not comp_path.exists():
        print(f"\n{config}: comparisons_raw.csv not found ({comp_path}) -- "
              f"run run_full_summary.py first. Skipping.")
        return []

    print(f"\n{'='*60}\n{config}\n{'='*60}")
    comparisons = pd.read_csv(comp_path)

    plots_dir = base / 'plots'
    tables_dir = base / 'tables'
    plots_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    all_summary = []
    for pair in METHOD_PAIRS:
        all_summary.extend(plot_pair(config, pair, comparisons, plots_dir))

    if all_summary:
        summary_df = pd.DataFrame(all_summary)
        out_csv = tables_dir / 'real_vs_inferred_summary.csv'
        summary_df.to_csv(out_csv, index=False)
        print(f"  wrote {out_csv.name}")
    return all_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', nargs='+', required=True,
                        help='Configuration name(s), e.g. conf_dup_loss_low_10M')
    parser.add_argument('--summary-base',
                        default='simulations/analysis/summary',
                        help='Base dir holding per-config summary folders '
                             '(default: simulations/analysis/summary)')
    args = parser.parse_args()

    summary_base = Path(args.summary_base)
    for config in args.config:
        process_config(config, summary_base)


if __name__ == '__main__':
    main()
