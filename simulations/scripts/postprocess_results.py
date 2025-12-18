#!/usr/bin/env python3
"""
postprocess_results.py - Post-process Program Outputs to Clean MUL-trees

Extracts clean MUL-tree strings from program output files and writes to
standardized {program}_result.tre files for the summary pipeline.

Usage:
    python postprocess_results.py CONFIG [--methods METHOD1 METHOD2 ...] [--dry-run]
    python postprocess_results.py conf_ils_low_10M
    python postprocess_results.py conf_ils_low_10M --methods grampa polyphest_p50
"""

import os
import sys
import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional
import yaml


class ResultPostprocessor:
    """Post-process program outputs to extract clean MUL-trees"""

    def __init__(self, config: str, config_dict: Dict, dry_run: bool = False):
        """
        Initialize post-processor

        Args:
            config: Configuration name (e.g., 'conf_ils_low_10M')
            config_dict: Configuration dictionary from YAML
            dry_run: If True, show what would be done without writing files
        """
        self.config = config
        self.base_dir = Path(config_dict['base_dir'])
        self.networks = config_dict['networks']
        self.methods = config_dict['methods']
        self.num_replicates = config_dict['num_replicates']
        self.dry_run = dry_run

        # Statistics
        self.stats = {
            'total': 0,
            'success': 0,
            'skipped': 0,
            'failed': 0,
            'errors': []
        }

    def extract_polyphest_multree(self, input_file: Path) -> Optional[str]:
        """
        Extract MUL-tree from Polyphest output

        Format:
            multree: (tree string here);
            network: (network string here);

        Returns:
            Clean MUL-tree string or None if failed
        """
        try:
            with open(input_file, 'r') as f:
                content = f.read()

            # Find line starting with "multree:"
            for line in content.split('\n'):
                if line.startswith('multree:'):
                    # Extract everything after "multree: "
                    multree = line.split('multree:', 1)[1].strip()
                    return multree

            return None

        except Exception as e:
            return None

    def extract_mpsugar_newick(self, input_file: Path) -> Optional[str]:
        """
        Extract Newick tree from MPSUGAR output

        Format:
            Network 1:
            Newick: (tree string here);
            Score: -12345

        Returns:
            Clean Newick string or None if failed
        """
        try:
            with open(input_file, 'r') as f:
                content = f.read()

            # Find line starting with "Newick:"
            for line in content.split('\n'):
                if line.strip().startswith('Newick:'):
                    # Extract everything after "Newick: "
                    newick = line.split('Newick:', 1)[1].strip()
                    return newick

            return None

        except Exception as e:
            return None

    def extract_grampa_best_tree(self, input_file: Path) -> Optional[str]:
        """
        Extract best MUL-tree from GRAMPA output and clean formatting

        Format (TSV):
            mul.tree    h1.node    h2.node    score    labeled.tree
            5137    <29>    <39>    3309    (best tree here)
            5116    <28>    <39>    3326    (second best tree)
            ...

        Cleans:
            - Removes * and + from taxon names (e.g., Lamium* → Lamium)
            - Removes internal node labels (e.g., <5>)

        Returns:
            Best tree (first data row, last column) or None if failed
        """
        try:
            with open(input_file, 'r') as f:
                lines = f.readlines()

            if len(lines) < 2:
                return None

            # Skip header line, get first data line
            first_data_line = lines[1].strip()

            if not first_data_line:
                return None

            # Split by tab and get last column (labeled.tree)
            columns = first_data_line.split('\t')

            if len(columns) < 5:
                return None

            best_tree = columns[-1].strip()

            # Clean the tree format
            best_tree = self._clean_grampa_tree(best_tree)

            # Ensure it ends with semicolon
            if not best_tree.endswith(';'):
                best_tree += ';'

            return best_tree

        except Exception as e:
            return None

    def _clean_grampa_tree(self, tree: str) -> str:
        """
        Clean GRAMPA tree format.

        Removes:
            - * and + from taxon names (e.g., Lamiumalbum+ → Lamiumalbum)
            - Internal node labels like <5>, <29>, etc.

        Args:
            tree: Raw GRAMPA tree string

        Returns:
            Cleaned tree string
        """
        import re

        # Remove * and + from taxon names
        # Pattern: taxon name followed by * or +
        tree = re.sub(r'([A-Za-z0-9_]+)[\*\+]', r'\1', tree)

        # Remove internal node labels like <5>, <29>, etc.
        # Pattern: <digits>
        tree = re.sub(r'<\d+>', '', tree)

        return tree

    def copy_padre_result(self, input_file: Path) -> Optional[str]:
        """
        Copy PADRE result (already clean MUL-tree)

        Args:
            input_file: Path to padre_trees-result.tre

        Returns:
            Tree string or None if failed
        """
        try:
            with open(input_file, 'r') as f:
                tree = f.read().strip()

            # Ensure it ends with semicolon
            if tree and not tree.endswith(';'):
                tree += ';'

            return tree

        except Exception as e:
            return None

    def process_file(self, method: str, network: str, replicate: int) -> bool:
        """
        Process one result file

        Args:
            method: Method name (e.g., 'grampa', 'polyphest_p50')
            network: Network name
            replicate: Replicate number (1-5)

        Returns:
            True if successful, False otherwise
        """
        # Determine input and output paths
        method_config = self.methods[method]
        method_dir = method_config['directory']

        # Determine input file location based on method
        if method == 'padre':
            # PADRE writes output to input directory (processed/)
            input_dir = (self.base_dir / network / "processed" / self.config /
                        "padre_input" / f"replicate_{replicate}")
            input_file = input_dir / "padre_trees-result.tre"

            # Output goes to results directory (for summary pipeline)
            output_dir = (self.base_dir / network / "results" / self.config /
                         method_dir / f"replicate_{replicate}")
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = output_dir / "padre_result.tre"
        else:
            # All other methods: input and output in results directory
            result_dir = (self.base_dir / network / "results" / self.config /
                         method_dir / f"replicate_{replicate}")

            # Determine input filename based on method
            if method.startswith('polyphest'):
                input_file = result_dir / "polyphest_trees-polyphest.txt"
            elif method == 'grampa':
                input_file = result_dir / "grampa-scores.txt"
            elif method == 'mpsugar':
                input_file = result_dir / "mpsugar_results.txt"
            else:
                print(f"  ERROR: Unknown method '{method}'")
                return False

            # Output file (standardized naming)
            output_file = result_dir / f"{method.split('_')[0]}_result.tre"

        # Check if input exists
        if not input_file.exists():
            self.stats['skipped'] += 1
            return False

        # Check if output already exists and is non-empty
        if output_file.exists() and output_file.stat().st_size > 0:
            self.stats['skipped'] += 1
            return True  # Already processed

        # Extract tree based on method
        tree_string = None

        if method.startswith('polyphest'):
            tree_string = self.extract_polyphest_multree(input_file)
        elif method == 'grampa':
            tree_string = self.extract_grampa_best_tree(input_file)
        elif method == 'mpsugar':
            tree_string = self.extract_mpsugar_newick(input_file)
        elif method == 'padre':
            tree_string = self.copy_padre_result(input_file)

        if tree_string is None:
            self.stats['failed'] += 1
            self.stats['errors'].append({
                'network': network,
                'method': method,
                'replicate': replicate,
                'error': f'Failed to extract tree from {input_file}'
            })
            return False

        # Write output
        if not self.dry_run:
            try:
                with open(output_file, 'w') as f:
                    f.write(tree_string)
                    if not tree_string.endswith('\n'):
                        f.write('\n')

                self.stats['success'] += 1
                return True

            except Exception as e:
                self.stats['failed'] += 1
                self.stats['errors'].append({
                    'network': network,
                    'method': method,
                    'replicate': replicate,
                    'error': f'Failed to write {output_file}: {e}'
                })
                return False
        else:
            # Dry run - just report
            print(f"  [DRY RUN] Would write: {output_file}")
            self.stats['success'] += 1
            return True

    def process_all(self, methods_filter: Optional[List[str]] = None):
        """
        Process all result files

        Args:
            methods_filter: List of methods to process (None = all)
        """
        methods_to_process = self.methods.keys()
        if methods_filter:
            methods_to_process = [m for m in methods_to_process if m in methods_filter]

        print(f"\n{'='*80}")
        print(f"Post-processing Results: {self.config}")
        print(f"{'='*80}")
        print(f"Methods: {', '.join(methods_to_process)}")
        print(f"Networks: {len(self.networks)}")
        print(f"Replicates: {self.num_replicates}")
        print(f"Dry run: {self.dry_run}")
        print(f"{'='*80}\n")

        # Process each combination
        for method in methods_to_process:
            print(f"\nProcessing {method}...")

            for network in self.networks:
                for replicate in range(1, self.num_replicates + 1):
                    self.stats['total'] += 1
                    success = self.process_file(method, network, replicate)

                    # Progress update
                    if self.stats['total'] % 50 == 0:
                        print(f"  Processed {self.stats['total']} files "
                              f"(Success: {self.stats['success']}, "
                              f"Skipped: {self.stats['skipped']}, "
                              f"Failed: {self.stats['failed']})")

        # Final statistics
        self.print_statistics()

    def print_statistics(self):
        """Print processing statistics"""
        print(f"\n{'='*80}")
        print(f"Post-processing Statistics")
        print(f"{'='*80}")
        print(f"Total files:        {self.stats['total']}")
        print(f"Successfully processed: {self.stats['success']} "
              f"({self.stats['success']/self.stats['total']*100:.1f}%)")
        print(f"Skipped (missing or already done): {self.stats['skipped']} "
              f"({self.stats['skipped']/self.stats['total']*100:.1f}%)")
        print(f"Failed:             {self.stats['failed']} "
              f"({self.stats['failed']/self.stats['total']*100:.1f}%)")

        if self.stats['errors']:
            print(f"\nErrors ({len(self.stats['errors'])}):")
            for i, error in enumerate(self.stats['errors'][:10], 1):
                print(f"  {i}. {error['network']} / {error['method']} / rep_{error['replicate']}")
                print(f"     {error['error']}")

            if len(self.stats['errors']) > 10:
                print(f"  ... and {len(self.stats['errors']) - 10} more")

        print(f"{'='*80}\n")


def load_config(config_path: str) -> Dict:
    """Load configuration from YAML file"""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(
        description='Post-process program outputs to extract clean MUL-trees',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all methods for a configuration
  %(prog)s conf_ils_low_10M

  # Process specific methods only
  %(prog)s conf_ils_low_10M --methods grampa polyphest_p50

  # Dry run (show what would be done)
  %(prog)s conf_ils_low_10M --dry-run

  # Process all configurations
  for config in conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M; do
      %(prog)s $config
  done

Output files:
  - Polyphest: polyphest_result.tre (extracted from polyphest_trees-polyphest.txt)
  - GRAMPA: grampa_result.tre (best tree from grampa-scores.txt)
  - MPSUGAR: mpsugar_result.tre (extracted from mpsugar_results.txt)
  - PADRE: padre_result.tre (copied from padre_tree-result.tre)
        """
    )

    parser.add_argument('configuration',
                       help='Configuration name (e.g., conf_ils_low_10M)')
    parser.add_argument('--config', default='simulations/summary_config.yaml',
                       help='Path to configuration YAML file (default: simulations/summary_config.yaml)')
    parser.add_argument('--methods', nargs='+',
                       help='Methods to process (default: all methods)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Dry run: show what would be done without writing files')

    args = parser.parse_args()

    # Load configuration
    try:
        config_dict = load_config(args.config)
    except FileNotFoundError:
        print(f"Error: Configuration file not found: {args.config}", file=sys.stderr)
        print(f"Current directory: {os.getcwd()}", file=sys.stderr)
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

    # Create post-processor and run
    processor = ResultPostprocessor(args.configuration, config_dict, dry_run=args.dry_run)
    processor.process_all(methods_filter=args.methods)


if __name__ == '__main__':
    main()
