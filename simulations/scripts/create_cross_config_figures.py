#!/usr/bin/env python3
"""
create_cross_config_figures.py - Cross-Configuration Comparative Analysis

Loads per-config summary results and generates figures comparing method
accuracy across simulation conditions (ILS levels, dup/loss rates, Ne).

Usage:
    # Auto-detect all available configs
    python create_cross_config_figures.py

    # Specify configs explicitly
    python create_cross_config_figures.py --configs conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M

    # Custom output directory
    python create_cross_config_figures.py --output /path/to/output
"""

import argparse
import gc
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

# ============================================================================
# CONSTANTS
# ============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
SUMMARY_BASE = SCRIPT_DIR.parent / "analysis" / "summary"
NETWORKS_DIR = SCRIPT_DIR.parent / "networks"

# All 9 configurations
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

# Config families — each family varies one parameter (low/medium/high)
CONFIG_FAMILIES = {
    'ILS': {
        'label': 'ILS Level',
        'description': 'Incomplete Lineage Sorting',
        'configs': {
            'Low':    'conf_ils_low_10M',
            'Medium': 'conf_ils_medium_10M',
            'High':   'conf_ils_high_10M',
        }
    },
    'Dup/Loss': {
        'label': 'Dup/Loss Rate',
        'description': 'Gene Duplication and Loss',
        'configs': {
            'Low':    'conf_dup_loss_low_10M',
            'Medium': 'conf_dup_loss_medium_10M',
            'High':   'conf_dup_loss_high_10M',
        }
    },
    'Dup/Loss (Ne=1M)': {
        'label': 'Dup/Loss Rate (Ne=1M)',
        'description': 'Gene Duplication and Loss with Ne=1M',
        'configs': {
            'Low':    'conf_dup_loss_low_10M_ne1M',
            'Medium': 'conf_dup_loss_medium_10M_ne1M',
            'High':   'conf_dup_loss_high_10M_ne1M',
        }
    },
}

# Visual style — consistent with create_analysis_figures.py
METHOD_COLORS = {
    'grampa': '#0173B2',
    'polyphest': '#DE8F05',
    'polyphest_p50': '#DE8F05',
    'polyphest_p70': '#029E73',
    'polyphest_p90': '#CC78BC',
    'mpsugar': '#8B4513',
    'padre': '#ECE133',
    'alloppnet': '#DC143C',
    'grandma_split': '#56B4E9',
}

METHOD_MARKERS = {
    'grampa': 'o',
    'polyphest': 's',
    'polyphest_p50': 's',
    'polyphest_p70': '^',
    'polyphest_p90': 'D',
    'mpsugar': 'v',
    'padre': 'P',
    'alloppnet': 'X',
    'grandma_split': 'd',
}

LEVEL_ORDER = ['Low', 'Medium', 'High']

# Display names for figures
METHOD_DISPLAY = {
    'grampa': 'GRAMPA',
    'polyphest': 'Polyphest',
    'polyphest_p50': 'Polyphest (p50)',
    'polyphest_p70': 'Polyphest (p70)',
    'polyphest_p90': 'Polyphest (p90)',
    'mpsugar': 'MPAllopp',
    'padre': 'PADRE',
    'alloppnet': 'AlloppNET',
    'grandma_split': r'GRAMPA$^{Iter}$',
}


def display_name(method: str) -> str:
    """Return publication-ready display name for a method."""
    return METHOD_DISPLAY.get(method, method)

# Key metrics for analysis
KEY_METRICS = {
    'edit_distance_multree': 'Edit Distance',
    # 'rf_distance': 'Robinson-Foulds Distance',  # Disabled: RF not well-defined for MUL-trees
    'ret_leaf_jaccard.dist': 'Reticulation Leaf Distance',
    'ret_sisters_jaccard.dist': 'Sister-Taxa Distance',
    'ploidy_diff.dist': 'Ploidy Distance',
    'num_rets_bias': 'Reticulation Bias',
}

# Plot style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['lines.markersize'] = 8


# ============================================================================
# DATA LOADING
# ============================================================================

def load_all_data(configs: List[str], summary_base: Path) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Load inventory, comparisons_raw, and aggregated_metrics for each config.

    Returns:
        dict: config_name -> {'inventory': df, 'comparisons': df, 'aggregated': df}
    """
    data = {}
    for config in configs:
        config_dir = summary_base / config
        if not config_dir.exists():
            print(f"  WARNING: No summary directory for {config}, skipping")
            continue

        entry = {}

        inv_file = config_dir / "inventory.csv"
        if inv_file.exists():
            entry['inventory'] = pd.read_csv(inv_file)
        else:
            print(f"  WARNING: No inventory.csv for {config}")
            continue

        comp_file = config_dir / "comparisons_raw.csv"
        if comp_file.exists():
            entry['comparisons'] = pd.read_csv(comp_file)

        agg_file = config_dir / "aggregated_metrics.csv"
        if agg_file.exists():
            entry['aggregated'] = pd.read_csv(agg_file)

        data[config] = entry

    return data


def build_combined_inventory(data: Dict) -> pd.DataFrame:
    """Combine inventory DataFrames from all configs, tagging each with config name."""
    from create_analysis_figures import merge_polyphest_inventory
    frames = []
    for config, dfs in data.items():
        df = dfs['inventory'].copy()
        df['config'] = config
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    return merge_polyphest_inventory(combined)


def build_combined_comparisons(data: Dict) -> pd.DataFrame:
    """Combine comparisons DataFrames from all configs."""
    from create_analysis_figures import merge_polyphest_comparisons
    frames = []
    for config, dfs in data.items():
        if 'comparisons' not in dfs:
            continue
        df = dfs['comparisons'].copy()
        if 'config' not in df.columns:
            df['config'] = config
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    return merge_polyphest_comparisons(combined)


def build_combined_aggregated(data: Dict) -> pd.DataFrame:
    """Combine aggregated metrics from all configs."""
    from create_analysis_figures import merge_polyphest_comparisons, reaggregate_metrics
    # Build from comparisons if available, to get properly merged polyphest
    comparisons = build_combined_comparisons(data)
    if not comparisons.empty:
        return reaggregate_metrics(comparisons)

    # Fallback: use pre-aggregated files
    frames = []
    for config, dfs in data.items():
        if 'aggregated' not in dfs:
            continue
        df = dfs['aggregated'].copy()
        if 'config' not in df.columns:
            df['config'] = config
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def tag_config_family(df: pd.DataFrame) -> pd.DataFrame:
    """Add 'family' and 'level' columns based on config name."""
    families = []
    levels = []
    for _, row in df.iterrows():
        config = row['config']
        family = 'Unknown'
        level = 'Unknown'
        for fam_name, fam_info in CONFIG_FAMILIES.items():
            for lev, cfg in fam_info['configs'].items():
                if config == cfg:
                    family = fam_name
                    level = lev
                    break
            if family != 'Unknown':
                break
        families.append(family)
        levels.append(level)
    df['family'] = families
    df['level'] = levels
    return df


# ============================================================================
# FIGURE GENERATION
# ============================================================================

class CrossConfigAnalyzer:
    """Generate cross-configuration comparison figures."""

    def __init__(self, data: Dict, output_dir: Path, network_stats: Optional[pd.DataFrame] = None):
        self.data = data
        self.output_dir = output_dir
        self.plots_dir = output_dir / "plots"
        self.tables_dir = output_dir / "tables"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.tables_dir.mkdir(parents=True, exist_ok=True)

        self.network_stats = network_stats

        # Build combined datasets (polyphest thresholds merged into single 'polyphest')
        self.inventory = build_combined_inventory(data)
        self.comparisons = build_combined_comparisons(data)
        # Re-aggregate from merged comparisons
        from create_analysis_figures import reaggregate_metrics
        if not self.comparisons.empty:
            self.aggregated = reaggregate_metrics(self.comparisons)
        else:
            self.aggregated = build_combined_aggregated(data)

        # Tag families
        if not self.inventory.empty:
            self.inventory = tag_config_family(self.inventory)
        if not self.comparisons.empty:
            self.comparisons = tag_config_family(self.comparisons)
        if not self.aggregated.empty:
            self.aggregated = tag_config_family(self.aggregated)

        available = list(data.keys())
        print(f"\nCross-config analyzer initialized:")
        print(f"  Configs loaded: {len(available)}")
        print(f"  Inventory rows: {len(self.inventory)}")
        print(f"  Comparison rows: {len(self.comparisons)}")
        print(f"  Aggregated rows: {len(self.aggregated)}")

    def generate_all(self):
        """Generate all cross-config figures and tables."""
        print(f"\n{'='*80}")
        print(f"Generating Cross-Configuration Analysis")
        print(f"Output: {self.output_dir}")
        print(f"{'='*80}\n")

        plot_num = 0
        total_plots = 11

        # --- 1. Completion rate heatmap ---
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Completion rate heatmap...")
        self.plot_completion_heatmap()

        # --- 2. Accuracy across conditions (one per key metric, per family) ---
        for metric_key, metric_label in [
            ('edit_distance_multree', 'Edit Distance'),
            # ('rf_distance', 'Robinson-Foulds Distance'),  # Disabled: RF not well-defined for MUL-trees
            ('ret_leaf_jaccard.dist', 'Reticulation Leaf Distance'),
            ('ploidy_diff.dist', 'Ploidy Distance'),
        ]:
            plot_num += 1
            print(f"[{plot_num}/{total_plots}] Accuracy across conditions: {metric_label}...")
            self.plot_accuracy_across_conditions(metric_key, metric_label)

        # --- 3. Accuracy heatmaps ---
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Accuracy heatmap (edit distance)...")
        self.plot_accuracy_heatmap('edit_distance_multree', 'Mean Edit Distance')

        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Accuracy heatmap (ret leaf Jaccard)...")
        self.plot_accuracy_heatmap('ret_leaf_jaccard.dist', 'Mean Reticulation Leaf Distance')

        # --- 4. Reticulation count bias across conditions ---
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Reticulation count bias across conditions...")
        self.plot_reticulation_bias_across_conditions()

        # --- 5. Network complexity interaction ---
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Network complexity × condition interaction...")
        self.plot_complexity_interaction()

        # --- 6. Method ranking across metrics ---
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Method ranking across metrics...")
        self.plot_method_ranking()

        # --- 7. Summary table ---
        plot_num += 1
        print(f"[{plot_num}/{total_plots}] Summary table...")
        self.generate_summary_table()

        print(f"\n{'='*80}")
        print(f"Cross-config analysis complete!")
        print(f"  Plots: {self.plots_dir}")
        print(f"  Tables: {self.tables_dir}")
        print(f"{'='*80}\n")

    # ========================================================================
    # FIGURE 1: Completion Rate Heatmap
    # ========================================================================

    def plot_completion_heatmap(self):
        """Heatmap: rows=methods, columns=configs, cells=completion %."""
        if self.inventory.empty:
            print("  WARNING: No inventory data, skipping")
            return

        # Compute completion rate per method × config
        completion = (
            self.inventory
            .groupby(['method', 'config'])['inferred_exists']
            .mean()
            .unstack(fill_value=0) * 100
        )

        # Order configs by family then level
        config_order = []
        for fam_name, fam_info in CONFIG_FAMILIES.items():
            for level in LEVEL_ORDER:
                cfg = fam_info['configs'].get(level)
                if cfg and cfg in completion.columns:
                    config_order.append(cfg)
        # Add any remaining
        for c in completion.columns:
            if c not in config_order:
                config_order.append(c)
        completion = completion[config_order]

        # Short config labels
        short_labels = [c.replace('conf_', '').replace('_10M', '') for c in config_order]

        fig, ax = plt.subplots(figsize=(max(14, len(config_order) * 1.5), max(5, len(completion) * 0.7)))

        im = ax.imshow(completion.values, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)

        # Annotate cells
        for i in range(completion.shape[0]):
            for j in range(completion.shape[1]):
                val = completion.values[i, j]
                text_color = 'white' if val < 40 or val > 85 else 'black'
                ax.text(j, i, f'{val:.0f}%', ha='center', va='center',
                        fontsize=10, fontweight='bold', color=text_color)

        ax.set_xticks(range(len(short_labels)))
        ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=10)
        ax.set_yticks(range(len(completion.index)))
        ax.set_yticklabels([display_name(m) for m in completion.index], fontsize=11)

        # Add family separators
        family_boundaries = []
        col_idx = 0
        for fam_name, fam_info in CONFIG_FAMILIES.items():
            n_configs = sum(1 for l in LEVEL_ORDER if fam_info['configs'].get(l) in config_order)
            if n_configs > 0:
                family_boundaries.append((col_idx, col_idx + n_configs, fam_name))
                col_idx += n_configs

        for start, end, name in family_boundaries:
            if end < len(config_order):
                ax.axvline(x=end - 0.5, color='black', linewidth=2)
            # Add family label at top
            mid = (start + end - 1) / 2
            ax.text(mid, -1.2, name, ha='center', va='bottom', fontsize=10, fontweight='bold')

        cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label('Completion Rate (%)', fontsize=12)

        ax.set_title('Method Completion Rate Across Configurations', fontsize=15, fontweight='bold', pad=40)
        plt.tight_layout()
        fig.savefig(self.plots_dir / "01_completion_rate_heatmap.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "01_completion_rate_heatmap.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()

    # ========================================================================
    # FIGURE 2: Accuracy Across Conditions (grouped bar/line per family)
    # ========================================================================

    def plot_accuracy_across_conditions(self, metric_key: str, metric_label: str):
        """
        For each config family, plot mean metric across levels, grouped by method.
        One subplot per family, methods as grouped bars.
        """
        if self.comparisons.empty:
            print(f"  WARNING: No comparison data, skipping {metric_label}")
            return

        # Filter to this metric and successful comparisons
        metric_df = self.comparisons[
            (self.comparisons['metric'] == metric_key) &
            (self.comparisons['status'] == 'SUCCESS')
        ].copy()

        if metric_df.empty:
            print(f"  WARNING: No data for metric {metric_key}, skipping")
            return

        # One subplot per family
        families_with_data = []
        for fam_name, fam_info in CONFIG_FAMILIES.items():
            fam_configs = list(fam_info['configs'].values())
            if metric_df[metric_df['config'].isin(fam_configs)].shape[0] > 0:
                families_with_data.append(fam_name)

        if not families_with_data:
            return

        n_fam = len(families_with_data)
        fig, axes = plt.subplots(1, n_fam, figsize=(7 * n_fam, 6), squeeze=False)
        axes = axes.flatten()

        for fam_idx, fam_name in enumerate(families_with_data):
            ax = axes[fam_idx]
            fam_info = CONFIG_FAMILIES[fam_name]

            fam_data = metric_df[metric_df['family'] == fam_name].copy()
            methods = sorted(fam_data['method'].unique())

            # Compute mean ± SE per method per level
            x_positions = np.arange(len(LEVEL_ORDER))
            bar_width = 0.8 / max(len(methods), 1)

            for m_idx, method in enumerate(methods):
                method_data = fam_data[fam_data['method'] == method]
                means = []
                sems = []
                for level in LEVEL_ORDER:
                    level_data = method_data[method_data['level'] == level]['value']
                    if len(level_data) > 0:
                        means.append(level_data.mean())
                        sems.append(level_data.std() / np.sqrt(len(level_data)) if len(level_data) > 1 else 0)
                    else:
                        means.append(np.nan)
                        sems.append(0)

                offset = (m_idx - len(methods) / 2 + 0.5) * bar_width
                bars = ax.bar(x_positions + offset, means, bar_width * 0.9,
                              yerr=sems, capsize=3,
                              label=display_name(method),
                              color=METHOD_COLORS.get(method, '#888888'),
                              edgecolor='white', linewidth=0.8,
                              alpha=0.85)

            ax.set_xlabel(fam_info['label'], fontsize=13, fontweight='bold')
            ax.set_ylabel(metric_label, fontsize=13, fontweight='bold')
            ax.set_title(fam_info['description'], fontsize=13, fontweight='bold', pad=10)
            ax.set_xticks(x_positions)
            ax.set_xticklabels(LEVEL_ORDER, fontsize=12)
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            ax.legend(fontsize=9, framealpha=0.9, loc='best')

        fig.suptitle(f'{metric_label} Across Simulation Conditions',
                     fontsize=15, fontweight='bold', y=1.02)
        plt.tight_layout()

        safe_name = metric_key.replace('.', '_')
        fig.savefig(self.plots_dir / f"02_{safe_name}_across_conditions.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"02_{safe_name}_across_conditions.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()

    # ========================================================================
    # FIGURE 3: Accuracy Heatmap
    # ========================================================================

    def plot_accuracy_heatmap(self, metric_key: str, title: str):
        """Heatmap: rows=methods, columns=configs, cells=mean metric value."""
        if self.comparisons.empty:
            return

        metric_df = self.comparisons[
            (self.comparisons['metric'] == metric_key) &
            (self.comparisons['status'] == 'SUCCESS')
        ]

        if metric_df.empty:
            print(f"  WARNING: No data for {metric_key}, skipping heatmap")
            return

        pivot = metric_df.groupby(['method', 'config'])['value'].mean().unstack()

        # Order columns
        config_order = []
        for fam_name, fam_info in CONFIG_FAMILIES.items():
            for level in LEVEL_ORDER:
                cfg = fam_info['configs'].get(level)
                if cfg and cfg in pivot.columns:
                    config_order.append(cfg)
        for c in pivot.columns:
            if c not in config_order:
                config_order.append(c)
        pivot = pivot.reindex(columns=config_order)

        short_labels = [c.replace('conf_', '').replace('_10M', '') for c in config_order]

        fig, ax = plt.subplots(figsize=(max(14, len(config_order) * 1.5), max(5, len(pivot) * 0.7)))

        # Use reversed colormap: lower distance = better = green
        im = ax.imshow(pivot.values, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=1)

        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                val = pivot.values[i, j]
                if np.isnan(val):
                    ax.text(j, i, '—', ha='center', va='center', fontsize=10, color='gray')
                else:
                    text_color = 'white' if val > 0.7 or val < 0.15 else 'black'
                    ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                            fontsize=10, fontweight='bold', color=text_color)

        ax.set_xticks(range(len(short_labels)))
        ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=10)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels([display_name(m) for m in pivot.index], fontsize=11)

        # Family separators
        col_idx = 0
        for fam_name, fam_info in CONFIG_FAMILIES.items():
            n_configs = sum(1 for l in LEVEL_ORDER if fam_info['configs'].get(l) in config_order)
            if n_configs > 0 and col_idx + n_configs < len(config_order):
                ax.axvline(x=col_idx + n_configs - 0.5, color='black', linewidth=2)
                mid = (col_idx + col_idx + n_configs - 1) / 2
                ax.text(mid, -1.2, fam_name, ha='center', va='bottom', fontsize=10, fontweight='bold')
            col_idx += n_configs

        cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label(title, fontsize=12)

        ax.set_title(f'{title} Across Configurations', fontsize=15, fontweight='bold', pad=40)
        plt.tight_layout()

        safe_name = metric_key.replace('.', '_')
        fig.savefig(self.plots_dir / f"03_{safe_name}_heatmap.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / f"03_{safe_name}_heatmap.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()

    # ========================================================================
    # FIGURE 4: Reticulation Bias Across Conditions
    # ========================================================================

    def plot_reticulation_bias_across_conditions(self):
        """
        Show whether methods over/under-predict reticulation count across conditions.
        Positive bias = over-predicts, negative = under-predicts.
        """
        if self.comparisons.empty:
            return

        bias_df = self.comparisons[
            (self.comparisons['metric'] == 'num_rets_bias') &
            (self.comparisons['status'] == 'SUCCESS')
        ].copy()

        if bias_df.empty:
            print("  WARNING: No num_rets_bias data, skipping")
            return

        families_with_data = []
        for fam_name, fam_info in CONFIG_FAMILIES.items():
            fam_configs = list(fam_info['configs'].values())
            if bias_df[bias_df['config'].isin(fam_configs)].shape[0] > 0:
                families_with_data.append(fam_name)

        if not families_with_data:
            return

        n_fam = len(families_with_data)
        fig, axes = plt.subplots(1, n_fam, figsize=(7 * n_fam, 6), squeeze=False)
        axes = axes.flatten()

        for fam_idx, fam_name in enumerate(families_with_data):
            ax = axes[fam_idx]
            fam_info = CONFIG_FAMILIES[fam_name]
            fam_data = bias_df[bias_df['family'] == fam_name]

            methods = sorted(fam_data['method'].unique())
            x_positions = np.arange(len(LEVEL_ORDER))
            bar_width = 0.8 / max(len(methods), 1)

            for m_idx, method in enumerate(methods):
                method_data = fam_data[fam_data['method'] == method]
                means = []
                sems = []
                for level in LEVEL_ORDER:
                    level_data = method_data[method_data['level'] == level]['value']
                    if len(level_data) > 0:
                        means.append(level_data.mean())
                        sems.append(level_data.std() / np.sqrt(len(level_data)) if len(level_data) > 1 else 0)
                    else:
                        means.append(np.nan)
                        sems.append(0)

                offset = (m_idx - len(methods) / 2 + 0.5) * bar_width
                ax.bar(x_positions + offset, means, bar_width * 0.9,
                       yerr=sems, capsize=3,
                       label=display_name(method),
                       color=METHOD_COLORS.get(method, '#888888'),
                       edgecolor='white', linewidth=0.8,
                       alpha=0.85)

            ax.axhline(y=0, color='black', linewidth=1, linestyle='-')
            ax.set_xlabel(fam_info['label'], fontsize=13, fontweight='bold')
            ax.set_ylabel('Reticulation Bias\n(positive = over-predicts)', fontsize=12, fontweight='bold')
            ax.set_title(fam_info['description'], fontsize=13, fontweight='bold', pad=10)
            ax.set_xticks(x_positions)
            ax.set_xticklabels(LEVEL_ORDER, fontsize=12)
            ax.grid(True, alpha=0.25, linestyle='--', axis='y')
            ax.legend(fontsize=9, framealpha=0.9, loc='best')

        fig.suptitle('Reticulation Bias Across Simulation Conditions',
                     fontsize=15, fontweight='bold', y=1.02)
        plt.tight_layout()
        fig.savefig(self.plots_dir / "04_reticulation_bias_across_conditions.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "04_reticulation_bias_across_conditions.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()

    # ========================================================================
    # FIGURE 5: Network Complexity × Condition Interaction
    # ========================================================================

    def plot_complexity_interaction(self):
        """
        For each config family, split networks into simple (H_Strict <= median)
        vs complex (H_Strict > median) and show how accuracy degrades differently.
        """
        if self.comparisons.empty or self.network_stats is None:
            print("  WARNING: Missing data for complexity interaction, skipping")
            return

        metric_key = 'edit_distance_multree'
        metric_label = 'Edit Distance'

        metric_df = self.comparisons[
            (self.comparisons['metric'] == metric_key) &
            (self.comparisons['status'] == 'SUCCESS')
        ].copy()

        if metric_df.empty:
            return

        # Tag networks as simple/complex based on H_Strict median
        stats = self.network_stats[['network', 'H_Strict']].copy()
        median_h = stats['H_Strict'].median()
        stats['complexity'] = stats['H_Strict'].apply(
            lambda x: 'Simple' if x <= median_h else 'Complex'
        )
        metric_df = metric_df.merge(stats[['network', 'complexity']], on='network', how='left')
        metric_df = metric_df.dropna(subset=['complexity'])

        families_with_data = []
        for fam_name in CONFIG_FAMILIES:
            fam_configs = list(CONFIG_FAMILIES[fam_name]['configs'].values())
            if metric_df[metric_df['config'].isin(fam_configs)].shape[0] > 0:
                families_with_data.append(fam_name)

        if not families_with_data:
            return

        n_fam = len(families_with_data)
        fig, axes = plt.subplots(1, n_fam, figsize=(7 * n_fam, 6), squeeze=False)
        axes = axes.flatten()

        for fam_idx, fam_name in enumerate(families_with_data):
            ax = axes[fam_idx]
            fam_info = CONFIG_FAMILIES[fam_name]
            fam_data = metric_df[metric_df['family'] == fam_name]
            methods = sorted(fam_data['method'].unique())

            for method in methods:
                method_data = fam_data[fam_data['method'] == method]
                color = METHOD_COLORS.get(method, '#888888')

                for complexity, ls in [('Simple', '-'), ('Complex', '--')]:
                    comp_data = method_data[method_data['complexity'] == complexity]
                    means = []
                    for level in LEVEL_ORDER:
                        level_data = comp_data[comp_data['level'] == level]['value']
                        means.append(level_data.mean() if len(level_data) > 0 else np.nan)

                    label = f"{method} ({complexity})"
                    ax.plot(LEVEL_ORDER, means,
                            marker=METHOD_MARKERS.get(method, 'o'),
                            color=color, linestyle=ls,
                            label=label, markersize=8,
                            markeredgewidth=1.5, markeredgecolor='white',
                            alpha=0.85)

            ax.set_xlabel(fam_info['label'], fontsize=13, fontweight='bold')
            ax.set_ylabel(metric_label, fontsize=13, fontweight='bold')
            ax.set_title(fam_info['description'], fontsize=13, fontweight='bold', pad=10)
            ax.grid(True, alpha=0.25, linestyle='--')
            ax.legend(fontsize=7, framealpha=0.9, loc='best', ncol=2)

        fig.suptitle(f'Network Complexity × Condition Interaction\n'
                     f'(Simple: H_Strict ≤ {median_h:.0f}, Complex: H_Strict > {median_h:.0f})',
                     fontsize=14, fontweight='bold', y=1.04)
        plt.tight_layout()
        fig.savefig(self.plots_dir / "05_complexity_condition_interaction.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "05_complexity_condition_interaction.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()

    # ========================================================================
    # FIGURE 6: Method Ranking Across Metrics
    # ========================================================================

    def plot_method_ranking(self):
        """
        Rank methods across key metrics (averaged over all configs).

        Produces two panels:
        - Left: heatmap of mean ranks per method per metric (lower = better)
        - Right: bar chart of overall average rank across all metrics
        """
        if self.comparisons.empty:
            print("  WARNING: No comparison data for method ranking, skipping")
            return

        # Metrics to rank on (lower is better for distance metrics,
        # but we handle them uniformly since all are distance-like)
        rank_metrics = {
            'edit_distance_multree': 'Edit\nDistance',
            # 'rf_distance': 'Robinson-Foulds\nDistance',  # Disabled: RF not well-defined for MUL-trees
            'ret_leaf_jaccard.dist': 'Reticulation\nLeaf Distance',
            'ret_sisters_jaccard.dist': 'Sister-Taxa\nDistance',
            'ploidy_diff.dist': 'Ploidy\nDistance',
        }

        # Also include completion rate (higher is better → we'll invert rank)
        # Compute per-method mean value across all configs for each metric
        method_metric_means = {}
        for metric_key in rank_metrics:
            metric_df = self.comparisons[
                (self.comparisons['metric'] == metric_key) &
                (self.comparisons['status'] == 'SUCCESS')
            ]
            if metric_df.empty:
                continue
            means = metric_df.groupby('method')['value'].mean()
            method_metric_means[metric_key] = means

        # Add completion rate (invert: 100 - rate so lower = better for ranking)
        if not self.inventory.empty:
            completion = self.inventory.groupby('method')['inferred_exists'].mean() * 100
            # Invert: methods with higher completion should rank better (lower rank)
            method_metric_means['completion_rate'] = -completion  # negative so lower = better
            rank_metrics_full = dict(rank_metrics)
            rank_metrics_full['completion_rate'] = 'Completion\nRate'
        else:
            rank_metrics_full = dict(rank_metrics)

        if not method_metric_means:
            print("  WARNING: No metric data for ranking, skipping")
            return

        # Build DataFrame: methods × metrics (mean values)
        all_methods = set()
        for means in method_metric_means.values():
            all_methods.update(means.index)
        all_methods = sorted(all_methods)

        values_df = pd.DataFrame(index=all_methods)
        for metric_key, means in method_metric_means.items():
            values_df[metric_key] = means

        # Rank within each metric (lower value = rank 1 = best)
        ranks_df = values_df.rank(axis=0, method='average', na_option='bottom')

        # Rename columns for display
        display_cols = {k: v for k, v in rank_metrics_full.items() if k in ranks_df.columns}
        ranks_display = ranks_df.rename(columns=display_cols)

        # Compute overall average rank
        avg_rank = ranks_df.mean(axis=1).sort_values()
        ranks_display = ranks_display.loc[avg_rank.index]

        # --- Create figure with two panels ---
        fig, (ax_heat, ax_bar) = plt.subplots(
            1, 2, figsize=(14, max(5, len(all_methods) * 0.7)),
            gridspec_kw={'width_ratios': [3, 1.2], 'wspace': 0.05}
        )

        # Left: rank heatmap
        n_methods = len(all_methods)
        import matplotlib.cm as cm
        cmap = cm.get_cmap('RdYlGn_r')  # Red=bad (high rank), Green=good (low rank)

        heat_data = ranks_display.values
        im = ax_heat.imshow(heat_data, cmap=cmap, aspect='auto',
                            vmin=1, vmax=n_methods)

        # Annotate cells with rank numbers
        for i in range(heat_data.shape[0]):
            for j in range(heat_data.shape[1]):
                val = heat_data[i, j]
                if not np.isnan(val):
                    text_color = 'white' if val > n_methods * 0.7 or val < n_methods * 0.3 else 'black'
                    ax_heat.text(j, i, f'{val:.1f}', ha='center', va='center',
                                fontsize=12, fontweight='bold', color=text_color)

        ax_heat.set_xticks(range(len(ranks_display.columns)))
        ax_heat.set_xticklabels(ranks_display.columns, fontsize=10, ha='center')
        ax_heat.set_yticks(range(len(ranks_display.index)))
        ax_heat.set_yticklabels([display_name(m) for m in ranks_display.index], fontsize=12, fontweight='bold')
        ax_heat.set_title('Rank per Metric (1 = best)', fontsize=13, fontweight='bold', pad=10)

        # Right: overall average rank bar chart
        colors = [METHOD_COLORS.get(m, '#888888') for m in avg_rank.index]
        bars = ax_bar.barh(range(len(avg_rank)), avg_rank.values, color=colors,
                           edgecolor='white', linewidth=1.5)

        # Add value labels on bars
        for i, (val, method) in enumerate(zip(avg_rank.values, avg_rank.index)):
            ax_bar.text(val + 0.05, i, f'{val:.2f}', va='center', fontsize=11, fontweight='bold')

        ax_bar.set_yticks(range(len(avg_rank)))
        ax_bar.set_yticklabels([])  # shared y-axis with heatmap
        ax_bar.set_xlabel('Average Rank', fontsize=12, fontweight='bold')
        ax_bar.set_title('Overall Ranking', fontsize=13, fontweight='bold', pad=10)
        ax_bar.set_xlim(0, n_methods + 0.5)
        ax_bar.invert_xaxis()  # best (lowest) on the right
        ax_bar.grid(True, axis='x', alpha=0.25, linestyle='--')

        fig.suptitle('Method Ranking Across All Metrics and Configurations',
                     fontsize=15, fontweight='bold', y=1.02)
        plt.tight_layout()
        fig.savefig(self.plots_dir / "06_method_ranking.pdf", bbox_inches='tight')
        fig.savefig(self.plots_dir / "06_method_ranking.png", bbox_inches='tight', dpi=300)
        plt.close('all')
        gc.collect()

        # Also save ranking table as CSV
        rank_table = ranks_df.copy()
        rank_table['avg_rank'] = avg_rank
        rank_table = rank_table.sort_values('avg_rank')
        rank_table.to_csv(self.tables_dir / "03_method_ranking.csv")
        print(f"  Saved: {self.tables_dir / '03_method_ranking.csv'}")

    # ========================================================================
    # TABLE: Summary across all configs and methods
    # ========================================================================

    def generate_summary_table(self):
        """Generate CSV: methods × configs × key metrics."""
        if self.comparisons.empty:
            print("  WARNING: No comparison data for summary table")
            return

        # Completion rate from inventory
        completion = (
            self.inventory
            .groupby(['method', 'config'])['inferred_exists']
            .mean() * 100
        ).reset_index()
        completion.columns = ['method', 'config', 'completion_rate']

        # Metric summaries
        rows = []
        for metric_key, metric_label in KEY_METRICS.items():
            metric_df = self.comparisons[
                (self.comparisons['metric'] == metric_key) &
                (self.comparisons['status'] == 'SUCCESS')
            ]
            if metric_df.empty:
                continue

            summary = metric_df.groupby(['method', 'config'])['value'].agg(
                ['mean', 'std', 'count']
            ).reset_index()
            summary['metric'] = metric_key
            summary['metric_label'] = metric_label
            rows.append(summary)

        if not rows:
            print("  WARNING: No metric data for summary table")
            return

        metrics_summary = pd.concat(rows, ignore_index=True)

        # Merge with completion
        metrics_summary = metrics_summary.merge(completion, on=['method', 'config'], how='left')

        # Tag families
        metrics_summary = tag_config_family(metrics_summary)

        # Save
        metrics_summary.to_csv(self.tables_dir / "01_cross_config_summary.csv", index=False)
        print(f"  Saved: {self.tables_dir / '01_cross_config_summary.csv'}")

        # Also create a wide-format pivot for easy reading
        for metric_key, metric_label in KEY_METRICS.items():
            mdf = metrics_summary[metrics_summary['metric'] == metric_key]
            if mdf.empty:
                continue
            pivot = mdf.pivot_table(index='method', columns='config', values='mean')

            # Order columns
            config_order = []
            for fam_name, fam_info in CONFIG_FAMILIES.items():
                for level in LEVEL_ORDER:
                    cfg = fam_info['configs'].get(level)
                    if cfg and cfg in pivot.columns:
                        config_order.append(cfg)
            for c in pivot.columns:
                if c not in config_order:
                    config_order.append(c)
            pivot = pivot.reindex(columns=config_order)

            safe_name = metric_key.replace('.', '_')
            pivot.to_csv(self.tables_dir / f"02_pivot_{safe_name}.csv")

        print(f"  Saved pivot tables to {self.tables_dir}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Cross-configuration comparative analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect all available configs
  %(prog)s

  # Specific configs
  %(prog)s --configs conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M

  # Custom output
  %(prog)s --output /path/to/cross_config_analysis
        """
    )
    parser.add_argument('--configs', nargs='+', default=None,
                        help='Configs to include (default: all available)')
    parser.add_argument('--output', default=None,
                        help='Output directory (default: simulations/analysis/cross_config)')
    parser.add_argument('--network-stats', default=None,
                        help='Network stats CSV (default: auto-detect)')

    args = parser.parse_args()

    # Determine configs
    configs = args.configs if args.configs else ALL_CONFIGS

    # Output directory
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = SCRIPT_DIR.parent / "analysis" / "cross_config"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Network stats
    if args.network_stats:
        network_stats_path = Path(args.network_stats)
    else:
        network_stats_path = NETWORKS_DIR / "mul_tree_final_stats.csv"

    network_stats = None
    if network_stats_path.exists():
        network_stats = pd.read_csv(network_stats_path)
        network_stats['network'] = network_stats['Filename'].str.replace('.tre', '')
        print(f"Loaded network stats: {len(network_stats)} networks")
    else:
        print(f"WARNING: Network stats not found at {network_stats_path}")

    # Load data
    print(f"\nLoading summary data from: {SUMMARY_BASE}")
    data = load_all_data(configs, SUMMARY_BASE)

    if not data:
        print("ERROR: No data loaded. Run run_full_summary.py first.", file=sys.stderr)
        sys.exit(1)

    # Generate figures
    analyzer = CrossConfigAnalyzer(data, output_dir, network_stats)
    analyzer.generate_all()


if __name__ == '__main__':
    main()
