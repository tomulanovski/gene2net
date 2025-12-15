#!/usr/bin/env python3
"""
Sample GTR+Gamma parameters and alignment length from empirical distributions.
This script samples from the 2,709 empirical gene parameter sets and outputs
parameters in a format suitable for iqtree --alisim.

Usage:
    python sample_gtr_parameters.py <pickle_file> [--output-format json|iqtree|bash]

Output formats:
    - json: JSON object with all parameters
    - iqtree: IQ-TREE model string for --alisim
    - bash: Bash-compatible variable assignments (default)
"""

import pickle
import random
import sys
import json
import argparse


def load_parameters(pickle_file):
    """Load GTR parameter distributions from pickle file."""
    try:
        with open(pickle_file, 'rb') as f:
            params_dict = pickle.load(f)

        # Combine all parameters from all datasets into one list
        all_params = []
        for dataset_name, params_list in params_dict.items():
            all_params.extend(params_list)

        return all_params

    except Exception as e:
        print(f"ERROR loading pickle file: {e}", file=sys.stderr)
        sys.exit(1)


def sample_parameters(all_params):
    """Randomly sample one parameter set from the empirical distribution."""
    if not all_params:
        print("ERROR: No parameters available to sample from", file=sys.stderr)
        sys.exit(1)

    return random.choice(all_params)


def format_iqtree_model(params):
    """
    Format parameters as IQ-TREE model string for --alisim.

    Format: GTR{rAC/rAG/rAT/rCG/rCT/rGT}+FU{piA/piC/piG/piT}+G{alpha}

    For alpha > 3.0, use uniform rate distribution (no +G) as per user's decision.
    """
    # GTR rates (all relative to GT which is 1.0)
    gtr_rates = f"{params['AC']:.6f}/{params['AG']:.6f}/{params['AT']:.6f}/" \
                f"{params['CG']:.6f}/{params['CT']:.6f}/{params['GT']:.6f}"

    # Base frequencies
    freqs = f"{params['pi_A']:.6f}/{params['pi_C']:.6f}/" \
            f"{params['pi_G']:.6f}/{params['pi_T']:.6f}"

    # Alpha parameter - if > 3.0, treat as uniform (no rate heterogeneity)
    alpha = params['alpha']
    if alpha > 3.0:
        # Use GTR+FU without +G (uniform rate distribution)
        model = f"GTR{{{gtr_rates}}}+FU{{{freqs}}}"
    else:
        # Use GTR+FU+G with the alpha parameter
        model = f"GTR{{{gtr_rates}}}+FU{{{freqs}}}+G{{{alpha:.6f}}}"

    return model


def format_bash_vars(params):
    """Format parameters as bash variable assignments."""
    model = format_iqtree_model(params)
    length = params['alignment_length']

    output = f"IQTREE_MODEL='{model}'\n"
    output += f"ALIGNMENT_LENGTH={length}\n"
    output += f"ALPHA={params['alpha']:.6f}\n"

    return output


def format_json(params):
    """Format parameters as JSON."""
    output = {
        'model_string': format_iqtree_model(params),
        'alignment_length': params['alignment_length'],
        'gtr_rates': {
            'AC': params['AC'],
            'AG': params['AG'],
            'AT': params['AT'],
            'CG': params['CG'],
            'CT': params['CT'],
            'GT': params['GT']
        },
        'base_frequencies': {
            'pi_A': params['pi_A'],
            'pi_C': params['pi_C'],
            'pi_G': params['pi_G'],
            'pi_T': params['pi_T']
        },
        'alpha': params['alpha'],
        'source_dataset': params.get('gene_id', 'unknown')
    }
    return json.dumps(output, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description='Sample GTR+Gamma parameters from empirical distributions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Get bash variables (default)
  python sample_gtr_parameters.py gtr_parameters_all.pkl

  # Get IQ-TREE model string only
  python sample_gtr_parameters.py gtr_parameters_all.pkl --output-format iqtree

  # Get JSON output
  python sample_gtr_parameters.py gtr_parameters_all.pkl --output-format json

  # Use in bash script
  eval $(python sample_gtr_parameters.py gtr_parameters_all.pkl)
  iqtree --alisim output -m "$IQTREE_MODEL" -t tree.nwk --length $ALIGNMENT_LENGTH
        """
    )

    parser.add_argument('pickle_file', type=str,
                       help='Path to GTR parameters pickle file')
    parser.add_argument('--output-format', choices=['bash', 'iqtree', 'json'],
                       default='bash',
                       help='Output format (default: bash)')
    parser.add_argument('--seed', type=int, default=None,
                       help='Random seed for reproducibility')

    args = parser.parse_args()

    # Set random seed if provided
    if args.seed is not None:
        random.seed(args.seed)

    # Load and sample parameters
    all_params = load_parameters(args.pickle_file)
    sampled = sample_parameters(all_params)

    # Output in requested format
    if args.output_format == 'iqtree':
        print(format_iqtree_model(sampled))
    elif args.output_format == 'json':
        print(format_json(sampled))
    else:  # bash
        print(format_bash_vars(sampled), end='')


if __name__ == '__main__':
    main()
