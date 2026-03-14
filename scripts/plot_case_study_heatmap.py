#!/usr/bin/env python3
"""
Generate per-dataset pairwise heatmap for case study figures.

Usage:
    python plot_case_study_heatmap.py comparisons_raw.csv --dataset Wu_2015 --output wu_heatmap.pdf
    python plot_case_study_heatmap.py comparisons_raw.csv --dataset Wu_2015  # saves to Wu_2015_heatmap.pdf
"""

import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# Display names for methods
METHOD_DISPLAY = {
    'grampa': 'GRAMPA',
    'polyphest': 'Polyphest',
    'padre': 'PADRE',
    'mpallop': 'AlloppNET',
    'alloppnet': 'AlloppNET',
    'grandma_split': 'Grandma',
}

# Display names for metrics
METRIC_DISPLAY = {
    'edit_distance_multree': 'MUL-tree Edit Distance',
    'polyploid_species_jaccard': 'Polyploid Species\nJaccard Distance',
    'rf_distance': 'Robinson-Foulds\nDistance',
    'ret_leaf_jaccard.dist': 'Reticulation Leaf Set\nJaccard Distance',
    'ret_sisters_jaccard.dist': 'Sister Relationship\nJaccard Distance',
    'ploidy_diff.dist': 'Ploidy Difference',
    'num_rets_diff': 'Reticulation Count\nDifference (abs)',
}


def make_symmetric_matrix(df_metric, methods):
    """Build a symmetric matrix from pairwise comparison data."""
    n = len(methods)
    mat = np.full((n, n), np.nan)
    method_idx = {m: i for i, m in enumerate(methods)}

    for _, row in df_metric.iterrows():
        m1, m2 = row['method1'], row['method2']
        if m1 in method_idx and m2 in method_idx:
            i, j = method_idx[m1], method_idx[m2]
            mat[i, j] = row['value']
            mat[j, i] = row['value']

    # Diagonal = 0 (self-comparison)
    np.fill_diagonal(mat, 0.0)
    return mat


def plot_heatmap(ax, matrix, labels, title, vmin=0, vmax=1, cmap='YlOrRd', fmt='.2f'):
    """Plot a single heatmap on an axis."""
    mask_diag = np.eye(len(labels), dtype=bool)

    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal')

    # Annotate cells
    for i in range(len(labels)):
        for j in range(len(labels)):
            if i == j:
                ax.text(j, i, '—', ha='center', va='center', fontsize=12, color='gray')
            elif not np.isnan(matrix[i, j]):
                val = matrix[i, j]
                # Choose text color based on background
                text_color = 'white' if val > (vmax - vmin) * 0.65 + vmin else 'black'
                ax.text(j, i, f'{val:{fmt}}', ha='center', va='center',
                        fontsize=11, fontweight='bold', color=text_color)
            else:
                ax.text(j, i, 'N/A', ha='center', va='center', fontsize=9, color='gray')

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=11, rotation=45, ha='right')
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)

    return im


def main():
    parser = argparse.ArgumentParser(description='Generate case study pairwise heatmap')
    parser.add_argument('comparisons_csv', help='Path to comparisons_raw.csv')
    parser.add_argument('--dataset', required=True, help='Dataset name (e.g., Wu_2015)')
    parser.add_argument('--metrics', nargs='+',
                        default=['edit_distance_multree', 'polyploid_species_jaccard'],
                        help='Metrics to plot (default: edit_distance_multree polyploid_species_jaccard)')
    parser.add_argument('--output', help='Output file path (default: {dataset}_heatmap.pdf)')
    parser.add_argument('--dpi', type=int, default=300)

    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.comparisons_csv)

    # Filter to dataset
    df_dataset = df[df['network'] == args.dataset]
    if df_dataset.empty:
        print(f"Error: No data found for dataset '{args.dataset}'", file=sys.stderr)
        print(f"Available datasets: {sorted(df['network'].unique())}", file=sys.stderr)
        sys.exit(1)

    # Get methods present in this dataset
    all_methods = sorted(set(df_dataset['method1'].unique()) | set(df_dataset['method2'].unique()))
    display_labels = [METHOD_DISPLAY.get(m, m) for m in all_methods]

    # Filter to requested metrics that exist
    available_metrics = [m for m in args.metrics if m in df_dataset['metric'].unique()]
    if not available_metrics:
        print(f"Error: None of the requested metrics found for {args.dataset}", file=sys.stderr)
        print(f"Available metrics: {sorted(df_dataset['metric'].unique())}", file=sys.stderr)
        sys.exit(1)

    n_panels = len(available_metrics)
    fig_width = 5.5 * n_panels + 1.5
    fig, axes = plt.subplots(1, n_panels, figsize=(fig_width, 5.5))
    if n_panels == 1:
        axes = [axes]

    for idx, metric in enumerate(available_metrics):
        df_metric = df_dataset[df_dataset['metric'] == metric]
        matrix = make_symmetric_matrix(df_metric, all_methods)

        # Determine scale
        valid_vals = matrix[~np.isnan(matrix) & ~np.eye(len(all_methods), dtype=bool)]
        if len(valid_vals) > 0:
            vmax = max(1.0, np.nanmax(valid_vals)) if metric == 'num_rets_diff' else 1.0
            fmt = '.0f' if metric == 'num_rets_diff' else '.2f'
        else:
            vmax = 1.0
            fmt = '.2f'

        title = METRIC_DISPLAY.get(metric, metric)
        im = plot_heatmap(axes[idx], matrix, display_labels, title, vmin=0, vmax=vmax, fmt=fmt)

        # Add colorbar
        cbar = fig.colorbar(im, ax=axes[idx], fraction=0.046, pad=0.04, shrink=0.8)
        cbar.ax.tick_params(labelsize=9)

    fig.suptitle(f'{args.dataset} — Pairwise Method Comparison',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = args.output or f'{args.dataset}_heatmap.pdf'
    fig.savefig(output_path, dpi=args.dpi, bbox_inches='tight')
    print(f"Saved: {output_path}")

    # Also save PNG
    png_path = output_path.rsplit('.', 1)[0] + '.png'
    fig.savefig(png_path, dpi=args.dpi, bbox_inches='tight')
    print(f"Saved: {png_path}")


if __name__ == '__main__':
    main()
