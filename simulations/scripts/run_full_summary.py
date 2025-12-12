#!/usr/bin/env python3
"""
run_full_summary.py - Summary Pipeline Orchestrator

One-command execution of the complete phylogenetic network inference evaluation pipeline.

Usage:
    python run_full_summary.py CONFIG [options]
    python run_full_summary.py conf_ils_low_10M
    python run_full_summary.py conf_ils_low_10M --methods grampa polyphest_p50 --force-recompute
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime
import yaml

# Import pipeline modules
from collect_results import ResultInventory, load_config
from compute_comparisons import ComparisonEngine
from aggregate_and_summarize import MetricAggregator, MultiLevelSummary


def run_pipeline(config_name: str, config_dict: dict, args):
    """
    Run the complete summary pipeline

    Args:
        config_name: Configuration name (e.g., 'conf_ils_low_10M')
        config_dict: Configuration dictionary from YAML
        args: Command-line arguments

    Returns:
        Tuple of (inventory_df, comparisons_df, aggregated_df)
    """
    print(f"\n{'='*80}")
    print(f"Summary Pipeline: {config_name}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*80}\n")

    # Setup output directory
    output_base = Path(args.output) if args.output else Path(config_dict['output_dir'])
    output_dir = output_base / config_name
    output_dir.mkdir(parents=True, exist_ok=True)

    cache_dir = output_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    print(f"Output directory: {output_dir}")
    print(f"Cache directory:  {cache_dir}\n")

    # ========================================================================
    # PHASE 1: Data Collection
    # ========================================================================
    print(f"{'='*80}")
    print(f"PHASE 1: Data Collection")
    print(f"{'='*80}\n")

    # Filter methods if specified
    methods_to_use = config_dict['methods']
    if args.methods:
        methods_to_use = {k: v for k, v in methods_to_use.items() if k in args.methods}
        print(f"Filtering to methods: {', '.join(methods_to_use.keys())}\n")

    # Create modified config with filtered methods
    filtered_config = config_dict.copy()
    filtered_config['methods'] = methods_to_use

    # Scan and inventory results
    inventory_scanner = ResultInventory(config_name, filtered_config)
    inventory_df = inventory_scanner.create_full_inventory()

    # Generate report
    completion_report = inventory_scanner.get_completion_report(inventory_df)
    inventory_scanner.print_completion_report(completion_report)

    # Export inventory
    inventory_file = output_dir / "inventory.csv"
    inventory_df.to_csv(inventory_file, index=False)
    print(f"Inventory saved to: {inventory_file}")

    # Dry run check
    if args.dry_run:
        print("\n[DRY RUN] Stopping before comparisons")
        print(f"Would process {len(inventory_df)} combinations")
        return inventory_df, None, None

    # ========================================================================
    # PHASE 2: Compute Comparisons
    # ========================================================================
    print(f"\n{'='*80}")
    print(f"PHASE 2: Compute Comparisons")
    print(f"{'='*80}\n")

    engine = ComparisonEngine(str(cache_dir), force_recompute=args.force_recompute)
    comparisons_df = engine.compute_all_comparisons(inventory_df)

    # Print statistics
    engine.print_statistics()

    # Export comparisons
    comparisons_file = output_dir / "comparisons_raw.csv"
    comparisons_df.to_csv(comparisons_file, index=False)
    print(f"Comparisons saved to: {comparisons_file}")

    # Write comparison report
    report_file = output_dir / "comparison_report.txt"
    engine.write_comparison_report(str(report_file))
    print(f"Comparison report saved to: {report_file}")

    # ========================================================================
    # PHASE 3: Aggregate and Summarize
    # ========================================================================
    print(f"\n{'='*80}")
    print(f"PHASE 3: Aggregate and Summarize")
    print(f"{'='*80}\n")

    # Aggregate replicates
    print("Aggregating replicates...")
    aggregator = MetricAggregator(comparisons_df)
    aggregated_df = aggregator.aggregate_replicates()

    print(f"Aggregated to {len(aggregated_df)} unique (network, method, metric) combinations\n")

    # Export aggregated metrics
    aggregated_file = output_dir / "aggregated_metrics.csv"
    aggregated_df.to_csv(aggregated_file, index=False)
    print(f"Aggregated metrics saved to: {aggregated_file}\n")

    # Generate multi-level summaries
    summary_gen = MultiLevelSummary(aggregated_df, config_name)

    # Determine network stats path
    network_stats_csv = args.network_stats
    if not network_stats_csv:
        # Try default location
        default_stats = Path(config_dict['networks_dir']) / "mul_tree_final_stats.csv"
        if default_stats.exists():
            network_stats_csv = str(default_stats)

    summary_gen.generate_all_summaries(output_dir, network_stats_csv)

    # ========================================================================
    # Pipeline Complete
    # ========================================================================
    print(f"\n{'='*80}")
    print(f"Pipeline Complete!")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*80}\n")

    print(f"Results summary:")
    print(f"  - Inventory:          {inventory_file}")
    print(f"  - Raw comparisons:    {comparisons_file}")
    print(f"  - Aggregated metrics: {aggregated_file}")
    print(f"  - Summary levels:     {output_dir}/level1_*, level2_*, level3_*")
    print(f"  - Reports:            {report_file}")
    print()

    return inventory_df, comparisons_df, aggregated_df


def main():
    parser = argparse.ArgumentParser(
        description='Run complete summary pipeline for phylogenetic network inference evaluation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process one configuration
  %(prog)s conf_ils_low_10M

  # Process with custom config file
  %(prog)s conf_ils_low_10M --config custom_config.yaml

  # Process specific methods only
  %(prog)s conf_ils_low_10M --methods grampa polyphest_p50 mpsugar

  # Force recompute (ignore cache)
  %(prog)s conf_ils_low_10M --force-recompute

  # Dry run (show what would be processed)
  %(prog)s conf_ils_low_10M --dry-run

  # Custom output location
  %(prog)s conf_ils_low_10M --output /path/to/summary/

  # Include network characteristics correlation
  %(prog)s conf_ils_low_10M --network-stats networks/mul_tree_final_stats.csv
        """
    )

    parser.add_argument('configuration',
                       help='Configuration name (e.g., conf_ils_low_10M)')
    parser.add_argument('--config', default='simulations/summary_config.yaml',
                       help='Path to configuration YAML file (default: simulations/summary_config.yaml)')
    parser.add_argument('--methods', nargs='+',
                       help='Methods to include (default: all methods in config)')
    parser.add_argument('--force-recompute', action='store_true',
                       help='Force recompute all comparisons (ignore cache)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Dry run: show what would be processed without computing')
    parser.add_argument('--output', metavar='DIR',
                       help='Output directory (default: from config file)')
    parser.add_argument('--network-stats', metavar='CSV',
                       help='Network characteristics CSV for correlation analysis')

    args = parser.parse_args()

    # Load configuration
    try:
        config_dict = load_config(args.config)
    except FileNotFoundError:
        print(f"Error: Configuration file not found: {args.config}", file=sys.stderr)
        print(f"Current directory: {os.getcwd()}", file=sys.stderr)
        print(f"\nTip: Make sure you're running from the gene2net root directory", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading configuration: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate methods if specified
    if args.methods:
        available_methods = set(config_dict['methods'].keys())
        requested_methods = set(args.methods)
        invalid_methods = requested_methods - available_methods

        if invalid_methods:
            print(f"Error: Invalid methods: {', '.join(invalid_methods)}", file=sys.stderr)
            print(f"Available methods: {', '.join(available_methods)}", file=sys.stderr)
            sys.exit(1)

    # Run pipeline
    try:
        inventory_df, comparisons_df, aggregated_df = run_pipeline(
            args.configuration,
            config_dict,
            args
        )
    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\nError running pipeline: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
