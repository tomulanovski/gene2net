#!/usr/bin/env python3
"""
run_analysis.py - Real Data Analysis Pipeline Orchestrator

One-command execution of the complete real data analysis pipeline.

Usage:
    python run_analysis.py [options]
    python run_analysis.py --force-recompute
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime
import yaml

# Add scripts directory to path for imports (must be first to avoid conflicts)
scripts_dir = Path(__file__).parent.resolve()
# Ensure scripts directory is at the very front of sys.path
if str(scripts_dir) in sys.path:
    sys.path.remove(str(scripts_dir))
sys.path.insert(0, str(scripts_dir))

# Import pipeline modules (from current scripts directory)
# Use explicit imports to avoid conflicts with simulations/scripts modules
import importlib.util

def import_from_file(module_name, file_path):
    """Import a module from a specific file path"""
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

# Import from scripts directory explicitly
collect_results = import_from_file('collect_results', scripts_dir / 'collect_results.py')
compute_comparisons = import_from_file('compute_comparisons', scripts_dir / 'compute_comparisons.py')
summarize_results = import_from_file('summarize_results', scripts_dir / 'summarize_results.py')
create_analysis_figures = import_from_file('create_analysis_figures', scripts_dir / 'create_analysis_figures.py')

ResultInventory = collect_results.ResultInventory
load_config = collect_results.load_config
ComparisonEngine = compute_comparisons.ComparisonEngine
ResultSummarizer = summarize_results.ResultSummarizer
RealDataAnalyzer = create_analysis_figures.RealDataAnalyzer


def run_pipeline(config_dict: dict, args):
    """
    Run the complete real data analysis pipeline

    Args:
        config_dict: Configuration dictionary from YAML
        args: Command-line arguments

    Returns:
        Tuple of (inventory_df, comparisons_df)
    """
    print(f"\n{'='*80}")
    print(f"Real Data Analysis Pipeline")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*80}\n")

    # Setup output directory
    output_base = Path(args.output) if args.output else Path(config_dict['output_dir'])
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = output_base / timestamp
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

    # Scan and inventory results
    inventory_scanner = ResultInventory(config_dict)
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
        print(f"Would process {len(inventory_df[inventory_df['exists']])} available networks")
        return inventory_df, None

    # ========================================================================
    # PHASE 2: Compute Pairwise Comparisons
    # ========================================================================
    print(f"\n{'='*80}")
    print(f"PHASE 2: Compute Pairwise Comparisons")
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
    # PHASE 3: Summarize Results
    # ========================================================================
    print(f"\n{'='*80}")
    print(f"PHASE 3: Summarize Results")
    print(f"{'='*80}\n")

    summarizer = ResultSummarizer(comparisons_df, inventory_df)
    summarizer.generate_all_summaries(output_dir)

    # ========================================================================
    # PHASE 4: Generate Figures
    # ========================================================================
    if not args.no_figures:
        print(f"\n{'='*80}")
        print(f"PHASE 4: Generate Figures")
        print(f"{'='*80}\n")

        analyzer = RealDataAnalyzer(
            comparisons_file=str(comparisons_file),
            inventory_file=str(inventory_file),
            output_dir=str(output_dir)
        )
        analyzer.generate_all_figures()

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
    print(f"  - Summary tables:     {output_dir}/method_availability.csv, pairwise_summary.csv, etc.")
    if not args.no_figures:
        print(f"  - Plots:              {output_dir}/plots/")
    print(f"  - Reports:            {report_file}")
    print()

    return inventory_df, comparisons_df


def main():
    parser = argparse.ArgumentParser(
        description='Run complete real data analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  %(prog)s

  # Force recompute (ignore cache)
  %(prog)s --force-recompute

  # Custom output location
  %(prog)s --output /path/to/analysis/

  # Dry run (show what would be processed)
  %(prog)s --dry-run

  # Skip figure generation
  %(prog)s --no-figures
        """
    )

    parser.add_argument('--config', default='scripts/papers_config.yaml',
                       help='Path to configuration YAML file (default: scripts/papers_config.yaml)')
    parser.add_argument('--force-recompute', action='store_true',
                       help='Force recompute all comparisons (ignore cache)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Dry run: show what would be processed without computing')
    parser.add_argument('--output', metavar='DIR',
                       help='Output directory (default: from config file)')
    parser.add_argument('--no-figures', action='store_true',
                       help='Skip figure generation')

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

    # Run pipeline
    try:
        inventory_df, comparisons_df = run_pipeline(config_dict, args)
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

