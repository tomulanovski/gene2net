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
    'grandma_split': r'GRAMPA$^{Iter}$',
}

# Display names for metrics
METRIC_DISPLAY = {
    'edit_distance_multree': 'Edit Distance',
    'polyploid_species_jaccard': 'Polyploid Species\nDistance',
    # 'rf_distance': 'RF Distance',  # Disabled: RF not well-defined for MUL-trees
    'ret_leaf_jaccard.dist': 'Reticulation Leaf\nDistance',
    'ret_sisters_jaccard.dist': 'Sister-Taxa\nDistance',
    'ploidy_diff.dist': 'Ploidy Distance',
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


def plot_split_heatmap(ax, matrix_upper, matrix_lower, labels,
                       title_upper, title_lower,
                       vmin=0, vmax=1, cmap_upper='YlOrRd', cmap_lower='YlGnBu',
                       fmt='.2f'):
    """
    Plot a split heatmap: upper triangle shows one metric, lower triangle shows another.
    Diagonal shows '—'.
    """
    n = len(labels)
    cmap_u = plt.get_cmap(cmap_upper)
    cmap_l = plt.get_cmap(cmap_lower)

    # Normalize values for coloring
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # Draw cell-by-cell as colored rectangles
    for i in range(n):
        for j in range(n):
            if i == j:
                color = 'white'
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                             facecolor=color, edgecolor='black', linewidth=0.5))
                ax.text(j, i, '—', ha='center', va='center', fontsize=12, color='gray')
            elif i < j:
                # Upper triangle
                val = matrix_upper[i, j]
                if np.isnan(val):
                    color = '#F0F0F0'
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                 facecolor=color, edgecolor='black', linewidth=0.5))
                    ax.text(j, i, 'N/A', ha='center', va='center', fontsize=9, color='gray')
                else:
                    color = cmap_u(norm(val))
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                 facecolor=color, edgecolor='black', linewidth=0.5))
                    text_color = 'white' if val > (vmax - vmin) * 0.65 + vmin else 'black'
                    ax.text(j, i, f'{val:{fmt}}', ha='center', va='center',
                            fontsize=11, fontweight='bold', color=text_color)
            else:
                # Lower triangle
                val = matrix_lower[i, j]
                if np.isnan(val):
                    color = '#F0F0F0'
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                 facecolor=color, edgecolor='black', linewidth=0.5))
                    ax.text(j, i, 'N/A', ha='center', va='center', fontsize=9, color='gray')
                else:
                    color = cmap_l(norm(val))
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                 facecolor=color, edgecolor='black', linewidth=0.5))
                    text_color = 'white' if val > (vmax - vmin) * 0.65 + vmin else 'black'
                    ax.text(j, i, f'{val:{fmt}}', ha='center', va='center',
                            fontsize=11, fontweight='bold', color=text_color)

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, fontsize=11, rotation=45, ha='right')
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_aspect('equal')

    # Return colormaps/norm for colorbar creation
    return cmap_u, cmap_l, norm


def main():
    parser = argparse.ArgumentParser(description='Generate case study pairwise heatmap')
    parser.add_argument('comparisons_csv', help='Path to comparisons_raw.csv')
    parser.add_argument('--dataset', required=True,
                        help='Dataset name (e.g., Wu_2015)')
    parser.add_argument('--metrics', nargs='+',
                        default=['edit_distance_multree', 'ret_leaf_jaccard.dist',
                                 'ret_sisters_jaccard.dist', 'ploidy_diff.dist',
                                 'num_rets_diff'],
                        help='Metrics to plot')
    parser.add_argument('--output', help='Output directory (default: plots/{dataset}/ next to CSV)')
    parser.add_argument('--combined', action='store_true',
                        help='Also generate a single combined figure with all metrics side by side')
    parser.add_argument('--split', nargs=2, action='append', metavar=('UPPER', 'LOWER'),
                        help='Generate split heatmap with UPPER metric in upper triangle '
                             'and LOWER metric in lower triangle. Can be repeated for '
                             'multiple split heatmaps. Example: '
                             '--split ret_leaf_jaccard.dist ret_sisters_jaccard.dist '
                             '--split edit_distance_multree ploidy_diff.dist')
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

    # Default output directory
    from pathlib import Path
    if args.output:
        plots_dir = Path(args.output)
    else:
        plots_dir = Path(args.comparisons_csv).parent / 'plots' / args.dataset
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Helper to get scale and format for a metric
    def get_scale(metric, matrix):
        valid_vals = matrix[~np.isnan(matrix) & ~np.eye(len(all_methods), dtype=bool)]
        if len(valid_vals) > 0:
            vmax = max(1.0, np.nanmax(valid_vals)) if metric == 'num_rets_diff' else 1.0
            fmt = '.0f' if metric == 'num_rets_diff' else '.2f'
        else:
            vmax = 1.0
            fmt = '.2f'
        return vmax, fmt

    # Generate one figure per metric
    for metric in available_metrics:
        df_metric = df_dataset[df_dataset['metric'] == metric]
        matrix = make_symmetric_matrix(df_metric, all_methods)
        vmax, fmt = get_scale(metric, matrix)
        title = METRIC_DISPLAY.get(metric, metric)

        fig, ax = plt.subplots(1, 1, figsize=(6, 5.5))
        im = plot_heatmap(ax, matrix, display_labels, title, vmin=0, vmax=vmax, fmt=fmt)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.8)
        cbar.ax.tick_params(labelsize=9)
        plt.tight_layout()

        safe_metric = metric.replace('.', '_')
        pdf_path = plots_dir / f'{safe_metric}.pdf'
        png_path = plots_dir / f'{safe_metric}.png'
        fig.savefig(str(pdf_path), dpi=args.dpi, bbox_inches='tight')
        fig.savefig(str(png_path), dpi=args.dpi, bbox_inches='tight')
        print(f"Saved: {pdf_path}")
        plt.close(fig)

    # Generate split heatmaps (two metrics in one: upper/lower triangle)
    if args.split:
        for upper_metric, lower_metric in args.split:
            # Validate both metrics exist
            if upper_metric not in df_dataset['metric'].unique():
                print(f"WARNING: Metric '{upper_metric}' not found for {args.dataset}, skipping split")
                continue
            if lower_metric not in df_dataset['metric'].unique():
                print(f"WARNING: Metric '{lower_metric}' not found for {args.dataset}, skipping split")
                continue

            df_upper = df_dataset[df_dataset['metric'] == upper_metric]
            df_lower = df_dataset[df_dataset['metric'] == lower_metric]
            matrix_upper = make_symmetric_matrix(df_upper, all_methods)
            matrix_lower = make_symmetric_matrix(df_lower, all_methods)

            # Determine scale — use 0-1 for normalized metrics, auto for counts
            is_count = 'num_rets_diff' in (upper_metric, lower_metric)
            if is_count:
                all_vals = np.concatenate([
                    matrix_upper[~np.isnan(matrix_upper)],
                    matrix_lower[~np.isnan(matrix_lower)]
                ])
                vmax = max(1.0, np.nanmax(all_vals)) if len(all_vals) > 0 else 1.0
                fmt = '.0f'
            else:
                vmax = 1.0
                fmt = '.2f'

            title_upper = METRIC_DISPLAY.get(upper_metric, upper_metric)
            title_lower = METRIC_DISPLAY.get(lower_metric, lower_metric)

            fig, ax = plt.subplots(1, 1, figsize=(8, 7))
            cmap_u, cmap_l, norm = plot_split_heatmap(
                ax, matrix_upper, matrix_lower, display_labels,
                title_upper, title_lower,
                vmin=0, vmax=vmax, fmt=fmt,
                cmap_upper='YlOrRd', cmap_lower='YlGnBu'
            )

            # Add two thin colorbars stacked on the right, no text labels
            import matplotlib.cm as cm
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax_u = divider.append_axes("right", size="4%", pad=0.15)
            cax_l = divider.append_axes("right", size="4%", pad=0.3)
            sm_upper = cm.ScalarMappable(cmap=cmap_u, norm=norm)
            sm_lower = cm.ScalarMappable(cmap=cmap_l, norm=norm)
            fig.colorbar(sm_upper, cax=cax_u)
            fig.colorbar(sm_lower, cax=cax_l)
            cax_u.tick_params(labelsize=8)
            cax_l.tick_params(labelsize=8)

            # Title with legend for which color is which
            title_u_clean = title_upper.replace(chr(10), " ")
            title_l_clean = title_lower.replace(chr(10), " ")
            ax.set_title(f'{args.dataset.replace("_", " ")}\n'
                         f'Upper \u25b3: {title_u_clean}  |  Lower \u25bd: {title_l_clean}',
                         fontsize=11, fontweight='bold', pad=12)

            plt.tight_layout()
            safe_u = upper_metric.replace('.', '_')
            safe_l = lower_metric.replace('.', '_')
            pdf_path = plots_dir / f'split_{safe_u}_vs_{safe_l}.pdf'
            png_path = plots_dir / f'split_{safe_u}_vs_{safe_l}.png'
            fig.savefig(str(pdf_path), dpi=args.dpi, bbox_inches='tight')
            fig.savefig(str(png_path), dpi=args.dpi, bbox_inches='tight')
            print(f"Saved split heatmap: {pdf_path}")
            plt.close(fig)

    # Generate combined figure with all metrics side by side
    if args.combined and len(available_metrics) > 1:
        n_panels = len(available_metrics)
        fig_width = 5.5 * n_panels + 1.5
        fig, axes = plt.subplots(1, n_panels, figsize=(fig_width, 5.5))
        if n_panels == 1:
            axes = [axes]

        for idx, metric in enumerate(available_metrics):
            df_metric = df_dataset[df_dataset['metric'] == metric]
            matrix = make_symmetric_matrix(df_metric, all_methods)
            vmax, fmt = get_scale(metric, matrix)
            title = METRIC_DISPLAY.get(metric, metric)

            im = plot_heatmap(axes[idx], matrix, display_labels, title, vmin=0, vmax=vmax, fmt=fmt)
            cbar = fig.colorbar(im, ax=axes[idx], fraction=0.046, pad=0.04, shrink=0.8)
            cbar.ax.tick_params(labelsize=9)

        plt.tight_layout()

        pdf_path = plots_dir / 'combined_heatmaps.pdf'
        png_path = plots_dir / 'combined_heatmaps.png'
        fig.savefig(str(pdf_path), dpi=args.dpi, bbox_inches='tight')
        fig.savefig(str(png_path), dpi=args.dpi, bbox_inches='tight')
        print(f"Saved: {pdf_path}")
        plt.close(fig)


if __name__ == '__main__':
    main()
