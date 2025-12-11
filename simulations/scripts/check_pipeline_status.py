#!/usr/bin/env python3
"""
check_pipeline_status.py - Comprehensive pipeline status checker

Validates SimPhy simulations and phylogenetic network inference pipeline steps:
- SimPhy simulations (1000 trees per replicate)
- Input preparation for GRAMPA, Polyphest, MPSUGAR, PADRE
- Method outputs and results

Usage:
    # Check everything for a configuration
    python check_pipeline_status.py conf_ils_low_10M

    # Check only SimPhy simulations
    python check_pipeline_status.py conf_ils_low_10M --step simphy

    # Check only GRAMPA inputs
    python check_pipeline_status.py conf_ils_low_10M --method grampa --step prep

    # Check only GRAMPA outputs
    python check_pipeline_status.py conf_ils_low_10M --method grampa --step run

    # Check all methods, only outputs
    python check_pipeline_status.py conf_ils_low_10M --step run

    # Export results to CSV
    python check_pipeline_status.py conf_ils_low_10M --export results.csv

    # Verbose mode with detailed errors
    python check_pipeline_status.py conf_ils_low_10M --verbose
"""

import argparse
import os
import sys
import glob
import re
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import csv

# ============================================================================
# CONSTANTS
# ============================================================================

BASE_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
LOG_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs"

NETWORKS = [
    "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
    "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019",
    "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011",
    "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014",
    "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014"
]

NUM_REPLICATES = 5
EXPECTED_TREES_PER_REPLICATE = 1000

# Expected files for each method
METHOD_FILES = {
    'grampa': {
        'prep': ['grampa_trees.tre', 'clean_trees.tre', 'species.tre'],  # taxa_map.txt is optional
        'prep_optional': ['taxa_map.txt'],  # Only created if substring fixes needed
        'run': ['analysis/grampa-scores.txt']  # GRAMPA success marker
    },
    'polyphest': {
        'prep': ['polyphest_trees.tre', 'multi_set.txt'],
        'run': ['polyphest_trees-polyphest.txt']  # Polyphest output file
    },
    'mpsugar': {
        'prep': ['mpsugar_trees.nex', 'taxon_map.json'],  # NEXUS format and JSON
        'run': ['mpsugar_results.txt']  # MPSUGAR output file
    },
    'padre': {
        'prep': ['padre_trees.tre'],
        'run': ['padre_tree-result.tre']  # PADRE output file
    }
}

# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class ValidationResult:
    """Result of validating a single network-replicate combination"""
    network: str
    replicate: int
    status: str  # 'SUCCESS', 'PARTIAL', 'FAILED', 'MISSING'
    missing_files: List[str]
    error_messages: List[str]
    details: str

@dataclass
class MethodSummary:
    """Summary statistics for a method"""
    method: str
    step: str
    total: int
    success: int
    failed: int
    missing: int
    success_rate: float

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def count_trees_in_file(filepath: str) -> int:
    """Count number of trees in a Newick file"""
    if not os.path.exists(filepath):
        return 0

    tree_count = 0
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and line.endswith(';'):
                    tree_count += 1
    except Exception as e:
        return -1  # Error reading file

    return tree_count

def count_gene_trees_in_simphy_output(data_dir: str, replicate: int) -> Tuple[int, str]:
    """
    Count gene trees in SimPhy output directory
    Returns: (count, status_message)
    """
    replicate_dir = os.path.join(data_dir, f"replicate_{replicate}")

    if not os.path.exists(replicate_dir):
        return 0, f"Replicate directory not found: {replicate_dir}"

    # Check for batched structure (batch_1/1/g_*, batch_2/1/g_*, etc.)
    batch_dirs = sorted(glob.glob(os.path.join(replicate_dir, "batch_*")))

    if batch_dirs:
        # Batched structure
        total_trees = 0
        for batch_dir in batch_dirs:
            inner_dir = os.path.join(batch_dir, "1")
            if os.path.exists(inner_dir):
                gene_files = glob.glob(os.path.join(inner_dir, "g_*"))
                total_trees += len(gene_files)

        return total_trees, f"Found {len(batch_dirs)} batches"
    else:
        # Single batch structure (1/g_*, 1/g_tree*, etc.)
        inner_dir = os.path.join(replicate_dir, "1")
        if os.path.exists(inner_dir):
            gene_files = glob.glob(os.path.join(inner_dir, "g_*"))
            return len(gene_files), "Single batch"
        else:
            return 0, "No gene trees directory found"

def parse_slurm_log_for_errors(log_pattern: str) -> Dict[int, List[str]]:
    """
    Parse SLURM error logs and extract error messages
    Returns: dict mapping array_task_id -> list of error messages
    """
    errors = defaultdict(list)

    log_files = glob.glob(log_pattern)

    for log_file in log_files:
        # Extract array task ID from filename
        match = re.search(r'_(\d+)\.err$', log_file)
        if not match:
            continue

        task_id = int(match.group(1))

        try:
            with open(log_file, 'r') as f:
                content = f.read()

            if content.strip():  # Non-empty error log
                # Extract key error lines
                error_lines = []
                for line in content.split('\n'):
                    line = line.strip()
                    if any(keyword in line.upper() for keyword in ['ERROR', 'FAILED', 'EXCEPTION', 'TRACEBACK']):
                        error_lines.append(line[:200])  # Limit line length

                if error_lines:
                    errors[task_id] = error_lines[:5]  # Limit to first 5 errors
                else:
                    errors[task_id] = ["Non-empty error log (check file for details)"]
        except Exception as e:
            errors[task_id] = [f"Could not parse log: {e}"]

    return errors

# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

def check_simphy_simulations(config: str, verbose: bool = False) -> List[ValidationResult]:
    """Check if SimPhy simulations completed successfully"""
    results = []

    for network in NETWORKS:
        data_dir = os.path.join(BASE_DIR, network, "data", config)

        for replicate in range(1, NUM_REPLICATES + 1):
            tree_count, status_msg = count_gene_trees_in_simphy_output(data_dir, replicate)

            if tree_count == EXPECTED_TREES_PER_REPLICATE:
                status = 'SUCCESS'
                details = f"{tree_count} trees ({status_msg})"
                missing = []
                errors = []
            elif tree_count > 0:
                status = 'PARTIAL'
                details = f"{tree_count}/{EXPECTED_TREES_PER_REPLICATE} trees ({status_msg})"
                missing = [f"Missing {EXPECTED_TREES_PER_REPLICATE - tree_count} trees"]
                errors = []
            else:
                status = 'MISSING'
                details = status_msg
                missing = [f"Expected {EXPECTED_TREES_PER_REPLICATE} trees"]
                errors = []

            results.append(ValidationResult(
                network=network,
                replicate=replicate,
                status=status,
                missing_files=missing,
                error_messages=errors,
                details=details
            ))

    return results

def check_method_inputs(config: str, method: str, verbose: bool = False) -> List[ValidationResult]:
    """Check if input files for a method are ready"""
    results = []
    expected_files = METHOD_FILES[method]['prep']
    optional_files = METHOD_FILES[method].get('prep_optional', [])

    for network in NETWORKS:
        for replicate in range(1, NUM_REPLICATES + 1):
            if method == 'grampa':
                input_dir = os.path.join(BASE_DIR, network, "processed", config, "grampa_input", f"replicate_{replicate}")
            elif method == 'polyphest':
                input_dir = os.path.join(BASE_DIR, network, "processed", config, "polyphest_input", f"replicate_{replicate}")
            elif method == 'mpsugar':
                input_dir = os.path.join(BASE_DIR, network, "processed", config, "mpsugar_input", f"replicate_{replicate}")
            elif method == 'padre':
                input_dir = os.path.join(BASE_DIR, network, "processed", config, "padre_input", f"replicate_{replicate}")
            else:
                continue

            missing = []
            errors = []
            optional_present = []

            # Check required files
            for expected_file in expected_files:
                filepath = os.path.join(input_dir, expected_file)
                if not os.path.exists(filepath):
                    missing.append(expected_file)
                elif os.path.getsize(filepath) == 0:
                    errors.append(f"{expected_file} is empty")

            # Check optional files (just note if present, don't mark as missing)
            for optional_file in optional_files:
                filepath = os.path.join(input_dir, optional_file)
                if os.path.exists(filepath):
                    optional_present.append(optional_file)

            if not missing and not errors:
                status = 'SUCCESS'
                details = f"All {len(expected_files)} required files present"
                if optional_present:
                    details += f" (+{len(optional_present)} optional)"
            elif missing:
                status = 'MISSING'
                details = f"Missing {len(missing)}/{len(expected_files)} required files"
            else:
                status = 'FAILED'
                details = f"Files exist but have errors"

            results.append(ValidationResult(
                network=network,
                replicate=replicate,
                status=status,
                missing_files=missing,
                error_messages=errors,
                details=details
            ))

    return results

def check_method_outputs(config: str, method: str, percentile: int = 60,
                        iterations: int = 500, chains: int = 1,
                        verbose: bool = False) -> List[ValidationResult]:
    """Check if output files for a method exist"""
    results = []
    expected_files = METHOD_FILES[method]['run']

    for network in NETWORKS:
        for replicate in range(1, NUM_REPLICATES + 1):
            if method == 'grampa':
                output_dir = os.path.join(BASE_DIR, network, "results", config, "grampa", f"replicate_{replicate}")
            elif method == 'polyphest':
                output_dir = os.path.join(BASE_DIR, network, "results", config, f"polyphest_p{percentile}", f"replicate_{replicate}")
            elif method == 'mpsugar':
                # Note: MPSUGAR doesn't include parameters in directory name
                output_dir = os.path.join(BASE_DIR, network, "results", config, "mpsugar", f"replicate_{replicate}")
            elif method == 'padre':
                output_dir = os.path.join(BASE_DIR, network, "results", config, "padre", f"replicate_{replicate}")
            else:
                continue

            missing = []
            errors = []

            if not os.path.exists(output_dir):
                status = 'MISSING'
                details = "Output directory does not exist"
                missing = ["Output directory"]
            else:
                # Check for specific expected files
                if expected_files:
                    for expected_file in expected_files:
                        filepath = os.path.join(output_dir, expected_file)
                        if not os.path.exists(filepath):
                            missing.append(expected_file)
                        elif os.path.getsize(filepath) == 0:
                            errors.append(f"{expected_file} is empty")

                    if not missing and not errors:
                        status = 'SUCCESS'
                        details = f"All {len(expected_files)} output file(s) present"
                    elif missing:
                        status = 'FAILED'
                        details = f"Missing {len(missing)}/{len(expected_files)} output file(s)"
                    else:
                        status = 'FAILED'
                        details = "Output files exist but have errors"
                else:
                    # Fallback: check for any result files (shouldn't reach here now)
                    result_files = glob.glob(os.path.join(output_dir, "*.nwk")) + \
                                  glob.glob(os.path.join(output_dir, "*.tre")) + \
                                  glob.glob(os.path.join(output_dir, "*.txt"))

                    if result_files:
                        status = 'SUCCESS'
                        details = f"Found {len(result_files)} output file(s)"
                    else:
                        status = 'FAILED'
                        details = "Output directory exists but no result files found"
                        missing = ["Result files"]

            results.append(ValidationResult(
                network=network,
                replicate=replicate,
                status=status,
                missing_files=missing,
                error_messages=errors,
                details=details
            ))

    return results

# ============================================================================
# REPORTING FUNCTIONS
# ============================================================================

def print_summary(results: List[ValidationResult], title: str):
    """Print summary statistics"""
    total = len(results)
    success = sum(1 for r in results if r.status == 'SUCCESS')
    failed = sum(1 for r in results if r.status == 'FAILED')
    missing = sum(1 for r in results if r.status == 'MISSING')
    partial = sum(1 for r in results if r.status == 'PARTIAL')

    success_rate = (success / total * 100) if total > 0 else 0

    # Calculate what we're checking
    num_networks = len(set(r.network for r in results))
    num_replicates = len(set(r.replicate for r in results))

    print(f"\n{'='*80}")
    print(f"{title}")
    print(f"{'='*80}")
    print(f"Checking: {num_networks} networks × {num_replicates} replicates = {total} combinations")
    print(f"")
    print(f"Success:  {success:3d} / {total:3d} ({success_rate:5.1f}%)")
    if partial > 0:
        print(f"Partial:  {partial:3d} / {total:3d} ({partial/total*100:5.1f}%)")
    if failed > 0:
        print(f"Failed:   {failed:3d} / {total:3d} ({failed/total*100:5.1f}%)")
    if missing > 0:
        print(f"Missing:  {missing:3d} / {total:3d} ({missing/total*100:5.1f}%)")
    print(f"{'='*80}\n")

def print_detailed_results(results: List[ValidationResult], show_success: bool = False):
    """Print detailed results"""

    # Group by network
    by_network = defaultdict(list)
    for result in results:
        by_network[result.network].append(result)

    for network in NETWORKS:
        network_results = by_network.get(network, [])
        if not network_results:
            continue

        # Check if all replicates succeeded
        all_success = all(r.status == 'SUCCESS' for r in network_results)

        if all_success and not show_success:
            continue  # Skip successful networks unless requested

        print(f"\n{network}:")
        print(f"{'-'*80}")

        for result in sorted(network_results, key=lambda x: x.replicate):
            status_symbol = {
                'SUCCESS': '✓',
                'PARTIAL': '~',
                'FAILED': '✗',
                'MISSING': '?'
            }.get(result.status, '?')

            print(f"  Replicate {result.replicate}: [{status_symbol}] {result.status:8s} - {result.details}")

            if result.missing_files:
                print(f"    Missing: {', '.join(result.missing_files)}")

            if result.error_messages:
                for error in result.error_messages:
                    print(f"    Error: {error}")

def export_to_csv(results: List[ValidationResult], output_file: str, title: str):
    """Export results to CSV file"""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Check', 'Network', 'Replicate', 'Status', 'Details', 'Missing Files', 'Errors'])

        for result in results:
            writer.writerow([
                title,
                result.network,
                result.replicate,
                result.status,
                result.details,
                '; '.join(result.missing_files) if result.missing_files else '',
                '; '.join(result.error_messages) if result.error_messages else ''
            ])

    print(f"Results exported to: {output_file}")

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Check pipeline status for phylogenetic network inference',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check everything
  %(prog)s conf_ils_low_10M

  # Check only SimPhy simulations
  %(prog)s conf_ils_low_10M --step simphy

  # Check GRAMPA inputs only
  %(prog)s conf_ils_low_10M --method grampa --step prep

  # Check all method outputs
  %(prog)s conf_ils_low_10M --step run

  # Verbose output with all details
  %(prog)s conf_ils_low_10M --verbose

  # Export to CSV
  %(prog)s conf_ils_low_10M --export results.csv
        """
    )

    parser.add_argument('config', help='Configuration name (e.g., conf_ils_low_10M)')
    parser.add_argument('--step', choices=['simphy', 'prep', 'run', 'all'], default='all',
                       help='Which step to check (default: all)')
    parser.add_argument('--method', choices=['grampa', 'polyphest', 'mpsugar', 'padre', 'all'], default='all',
                       help='Which method to check (default: all)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Show detailed results including successful runs')
    parser.add_argument('--export', metavar='FILE',
                       help='Export results to CSV file')
    parser.add_argument('--percentile', type=int, default=60,
                       help='Polyphest percentile parameter (default: 60)')
    parser.add_argument('--iterations', type=int, default=500,
                       help='MPSUGAR iterations parameter (default: 500)')
    parser.add_argument('--chains', type=int, default=1,
                       help='MPSUGAR chains parameter (default: 1)')

    args = parser.parse_args()

    print(f"\n{'='*80}")
    print(f"Pipeline Status Check: {args.config}")
    print(f"{'='*80}")

    all_results = []

    # Check SimPhy simulations
    if args.step in ['simphy', 'all']:
        print("\n[1/3] Checking SimPhy simulations...")
        simphy_results = check_simphy_simulations(args.config, args.verbose)
        print_summary(simphy_results, f"SimPhy Simulations - {args.config}")
        if args.verbose or any(r.status != 'SUCCESS' for r in simphy_results):
            print_detailed_results(simphy_results, show_success=args.verbose)
        all_results.extend([(f"SimPhy", r) for r in simphy_results])

    # Check method inputs (prep step)
    if args.step in ['prep', 'all']:
        methods = ['grampa', 'polyphest', 'mpsugar', 'padre'] if args.method == 'all' else [args.method]

        for i, method in enumerate(methods, start=1):
            print(f"\n[{i+1}/{len(methods)+1}] Checking {method.upper()} inputs...")
            input_results = check_method_inputs(args.config, method, args.verbose)
            print_summary(input_results, f"{method.upper()} Input Files - {args.config}")
            if args.verbose or any(r.status != 'SUCCESS' for r in input_results):
                print_detailed_results(input_results, show_success=args.verbose)
            all_results.extend([(f"{method.upper()}_prep", r) for r in input_results])

    # Check method outputs (run step)
    if args.step in ['run', 'all']:
        methods = ['grampa', 'polyphest', 'mpsugar', 'padre'] if args.method == 'all' else [args.method]

        for i, method in enumerate(methods, start=1):
            print(f"\n[{i+1}/{len(methods)+1}] Checking {method.upper()} outputs...")
            output_results = check_method_outputs(
                args.config, method,
                percentile=args.percentile,
                iterations=args.iterations,
                chains=args.chains,
                verbose=args.verbose
            )
            print_summary(output_results, f"{method.upper()} Output Files - {args.config}")
            if args.verbose or any(r.status != 'SUCCESS' for r in output_results):
                print_detailed_results(output_results, show_success=args.verbose)
            all_results.extend([(f"{method.upper()}_run", r) for r in output_results])

    # Export to CSV if requested
    if args.export and all_results:
        with open(args.export, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Step', 'Network', 'Replicate', 'Status', 'Details', 'Missing Files', 'Errors'])

            for step, result in all_results:
                writer.writerow([
                    step,
                    result.network,
                    result.replicate,
                    result.status,
                    result.details,
                    '; '.join(result.missing_files) if result.missing_files else '',
                    '; '.join(result.error_messages) if result.error_messages else ''
                ])

        print(f"\n{'='*80}")
        print(f"Results exported to: {args.export}")
        print(f"{'='*80}")

    print("\n" + "="*80)
    print("Check complete!")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
