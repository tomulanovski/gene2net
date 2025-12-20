#!/usr/bin/env python3
"""
Create a simplified summary for meetings - focuses on key findings only.
"""

import pandas as pd
import argparse
from pathlib import Path


def create_meeting_summary(config_name, output_dir):
    """Create simplified meeting-friendly summary from existing analysis."""

    summary_dir = Path(f"simulations/analysis/summary/{config_name}")

    if not summary_dir.exists():
        print(f"Error: Summary directory not found: {summary_dir}")
        print("Run 'python simulations/scripts/run_full_summary.py {config_name}' first")
        return

    # Create output directory
    meeting_dir = Path(output_dir) / config_name
    meeting_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*80}")
    print(f"Creating Meeting Summary for {config_name}")
    print(f"{'='*80}\n")

    # 1. Read aggregated metrics
    agg_file = summary_dir / "aggregated_metrics.csv"
    if not agg_file.exists():
        print(f"Error: {agg_file} not found")
        return

    agg = pd.read_csv(agg_file)

    # 2. Create simplified overall rankings - PRIMARY METRICS ONLY
    primary_metrics = [
        'edit_distance',
        'num_rets_diff',
        'ploidy_diff.dist',
        'ret_leaf_jaccard.dist',
        'ret_sisters_jaccard.dist'
    ]

    # Filter to primary metrics
    primary_data = agg[agg['metric'].isin(primary_metrics)].copy()

    # Calculate rankings per metric
    rankings = []
    for metric in primary_metrics:
        metric_data = primary_data[primary_data['metric'] == metric]

        # Group by method and calculate average
        method_avg = metric_data.groupby('method').agg({
            'mean': 'mean',
            'std': 'mean',
            'network': 'count'
        }).reset_index()
        method_avg.columns = ['method', 'avg_value', 'avg_std', 'num_networks']

        # Rank methods (lower is better for all these metrics)
        method_avg = method_avg.sort_values('avg_value')
        method_avg['rank'] = range(1, len(method_avg) + 1)
        method_avg['metric'] = metric

        rankings.append(method_avg)

    rankings_df = pd.concat(rankings, ignore_index=True)
    rankings_df = rankings_df[['metric', 'rank', 'method', 'avg_value', 'avg_std', 'num_networks']]
    rankings_df = rankings_df.sort_values(['metric', 'rank'])

    # Save
    rankings_file = meeting_dir / "01_method_rankings.csv"
    rankings_df.to_csv(rankings_file, index=False)
    print(f"[OK] Created: {rankings_file}")

    # 3. Create "best method per metric" summary
    best_per_metric = rankings_df[rankings_df['rank'] == 1][['metric', 'method', 'avg_value', 'num_networks']]
    best_file = meeting_dir / "02_best_methods.csv"
    best_per_metric.to_csv(best_file, index=False)
    print(f"[OK] Created: {best_file}")

    # 4. Create simplified per-network summary (edit_distance only)
    edit_dist_data = primary_data[primary_data['metric'] == 'edit_distance'].copy()

    # Pivot: networks as rows, methods as columns
    pivot = edit_dist_data.pivot_table(
        index='network',
        columns='method',
        values='mean',
        aggfunc='first'
    )

    # Add best method column (only use numeric columns)
    numeric_cols = pivot.select_dtypes(include='number').columns
    pivot['best_method'] = pivot[numeric_cols].idxmin(axis=1)
    pivot['best_value'] = pivot[numeric_cols].min(axis=1)

    # Sort by difficulty (higher edit_distance = harder network)
    pivot = pivot.sort_values('best_value')

    network_file = meeting_dir / "03_per_network_edit_distance.csv"
    pivot.to_csv(network_file)
    print(f"[OK] Created: {network_file}")

    # 5. Create key findings text summary
    findings_file = meeting_dir / "00_KEY_FINDINGS.txt"

    with open(findings_file, 'w') as f:
        f.write(f"{'='*80}\n")
        f.write(f"KEY FINDINGS - {config_name}\n")
        f.write(f"{'='*80}\n\n")

        f.write("BEST METHOD PER METRIC:\n")
        f.write("-" * 80 + "\n")
        for _, row in best_per_metric.iterrows():
            f.write(f"  {row['metric']:30s} -> {row['method']:20s} (value: {row['avg_value']:.4f})\n")

        f.write("\n" + "="*80 + "\n")
        f.write("OVERALL METHOD PERFORMANCE (avg rank across all metrics):\n")
        f.write("-" * 80 + "\n")

        # Calculate average rank per method
        avg_ranks = rankings_df.groupby('method')['rank'].mean().sort_values()
        for method, avg_rank in avg_ranks.items():
            f.write(f"  {method:20s} -> avg rank: {avg_rank:.2f}\n")

        f.write("\n" + "="*80 + "\n")
        f.write("EASIEST NETWORKS (lowest edit_distance):\n")
        f.write("-" * 80 + "\n")
        for network in pivot.head(5).index:
            best_method = pivot.loc[network, 'best_method']
            best_value = pivot.loc[network, 'best_value']
            f.write(f"  {network:30s} -> {best_value:.4f} ({best_method})\n")

        f.write("\n" + "="*80 + "\n")
        f.write("HARDEST NETWORKS (highest edit_distance):\n")
        f.write("-" * 80 + "\n")
        for network in pivot.tail(5).index:
            best_method = pivot.loc[network, 'best_method']
            best_value = pivot.loc[network, 'best_value']
            f.write(f"  {network:30s} -> {best_value:.4f} ({best_method})\n")

        # Add correlation insights if available
        corr_file = summary_dir / "level3_network_correlations.csv"
        if corr_file.exists():
            corr = pd.read_csv(corr_file)
            sig_corr = corr[corr['significant'] == True].sort_values('correlation', ascending=False)

            f.write("\n" + "="*80 + "\n")
            f.write("SIGNIFICANT CORRELATIONS (what makes networks harder?):\n")
            f.write("-" * 80 + "\n")

            if len(sig_corr) > 0:
                f.write("Strong positive correlations (property increases -> worse performance):\n")
                for _, row in sig_corr.head(10).iterrows():
                    f.write(f"  {row['method']:15s} | {row['metric']:25s} | "
                           f"{row['network_property']:20s} -> r={row['correlation']:.3f} (p={row['p_value']:.4f})\n")
            else:
                f.write("  No significant correlations found\n")

    print(f"[OK] Created: {findings_file}")

    # 6. Copy comparison report
    report_src = summary_dir / "comparison_report.txt"
    if report_src.exists():
        report_dst = meeting_dir / "04_comparison_report.txt"
        import shutil
        shutil.copy(report_src, report_dst)
        print(f"[OK] Created: {report_dst}")

    print(f"\n{'='*80}")
    print(f"Meeting summary created in: {meeting_dir}")
    print(f"{'='*80}")
    print("\nFILES TO REVIEW (in order):")
    print("  1. 00_KEY_FINDINGS.txt       - Start here! High-level summary")
    print("  2. 02_best_methods.csv       - Which method wins for each metric")
    print("  3. 01_method_rankings.csv    - Detailed rankings per metric")
    print("  4. 03_per_network_edit_distance.csv - Which networks are easy/hard")
    print("  5. 04_comparison_report.txt  - Technical details (failures, etc.)")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Create simplified meeting summary from full analysis results"
    )
    parser.add_argument(
        'config',
        help='Configuration name (e.g., conf_ils_low_10M)'
    )
    parser.add_argument(
        '--output',
        default='simulations/analysis/meeting_summaries',
        help='Output directory (default: simulations/analysis/meeting_summaries)'
    )

    args = parser.parse_args()

    create_meeting_summary(args.config, args.output)


if __name__ == '__main__':
    main()
