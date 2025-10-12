#!/usr/bin/env python3
"""
Run MP-SUGAR on gene trees with taxon mapping

Usage:
    python mpsugar.py -t <trees.nex> -m <taxon_map.json> -o <output_file> [options]
"""

import argparse
import json
import random
import sys
from datetime import datetime
from PhyNetPy.MPSugar import MP_SUGAR


def load_taxon_map(map_file):
    """
    Load taxon map from JSON file.
    
    Args:
        map_file (str): Path to JSON file containing taxon map
        
    Returns:
        dict: Taxon mapping dictionary
    """
    try:
        with open(map_file, 'r') as f:
            taxon_map = json.load(f)
        return taxon_map
    except Exception as e:
        print(f"ERROR: Failed to load taxon map from {map_file}: {e}")
        sys.exit(1)


def run_mpsugar(trees_file, taxon_map, output_file, iter_ct=500, num_chains=1, seed=None):
    """
    Run MP-SUGAR algorithm on gene trees.
    
    Args:
        trees_file (str): Path to NEXUS file with gene trees
        taxon_map (dict): Taxon mapping dictionary
        output_file (str): Path to output file for results
        iter_ct (int): Number of iterations for hill climbing
        num_chains (int): Number of independent hill climbing chains
        seed (int): Random seed (if None, will be random)
    """
    results = []
    
    print(f"Running MP-SUGAR with {num_chains} chain(s), {iter_ct} iterations each")
    print(f"Input trees: {trees_file}")
    print(f"Output file: {output_file}")
    print("=" * 60)
    
    for chain_num in range(1, num_chains + 1):
        # Generate or use provided seed
        chain_seed = seed if seed is not None else random.randint(0, 1000000)
        
        print(f"\nChain {chain_num}/{num_chains} (seed: {chain_seed})")
        print("-" * 60)
        
        try:
            # Run MP-SUGAR
            start_time = datetime.now()
            output_networks = MP_SUGAR(
                trees_file, 
                taxon_map, 
                iter_ct=iter_ct, 
                seed=chain_seed
            )
            end_time = datetime.now()
            elapsed = (end_time - start_time).total_seconds()
            
            # Store results
            chain_results = {
                'chain': chain_num,
                'seed': chain_seed,
                'elapsed_seconds': elapsed,
                'networks': []
            }
            
            # Process each network found
            for net, score in output_networks.items():
                network_newick = net.to_newick()
                print(f"Network: {network_newick}")
                print(f"Score: {score}")
                
                chain_results['networks'].append({
                    'newick': network_newick,
                    'score': score
                })
            
            results.append(chain_results)
            print(f"Chain completed in {elapsed:.2f} seconds")
            
        except Exception as e:
            print(f"ERROR in chain {chain_num}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Write results to output file
    try:
        with open(output_file, 'w') as f:
            f.write(f"MP-SUGAR Results\n")
            f.write(f"{'=' * 80}\n")
            f.write(f"Input file: {trees_file}\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Iterations: {iter_ct}\n")
            f.write(f"Number of chains: {num_chains}\n")
            f.write(f"{'=' * 80}\n\n")
            
            for chain_result in results:
                f.write(f"Chain {chain_result['chain']} (seed: {chain_result['seed']})\n")
                f.write(f"Elapsed time: {chain_result['elapsed_seconds']:.2f} seconds\n")
                f.write(f"Networks found: {len(chain_result['networks'])}\n")
                f.write(f"{'-' * 80}\n")
                
                for i, net_info in enumerate(chain_result['networks'], 1):
                    f.write(f"\nNetwork {i}:\n")
                    f.write(f"Newick: {net_info['newick']}\n")
                    f.write(f"Score: {net_info['score']}\n")
                
                f.write(f"\n{'=' * 80}\n\n")
        
        print(f"\nResults saved to: {output_file}")
        
    except Exception as e:
        print(f"ERROR: Failed to write output file: {e}")
        sys.exit(1)


def main():
    """
    Main function to parse arguments and run MP-SUGAR.
    """
    parser = argparse.ArgumentParser(
        description="Run MP-SUGAR on gene trees with taxon mapping"
    )
    
    parser.add_argument(
        "-t", "--trees",
        required=True,
        help="Path to NEXUS file containing gene trees"
    )
    
    parser.add_argument(
        "-m", "--map",
        required=True,
        help="Path to JSON file containing taxon map"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output file for results"
    )
    
    parser.add_argument(
        "-i", "--iterations",
        type=int,
        default=500,
        help="Number of iterations for hill climbing (default: 500)"
    )
    
    parser.add_argument(
        "-c", "--chains",
        type=int,
        default=1,
        help="Number of independent hill climbing chains (default: 1)"
    )
    
    parser.add_argument(
        "-s", "--seed",
        type=int,
        default=None,
        help="Random seed (default: random)"
    )
    
    args = parser.parse_args()
    
    # Load taxon map
    print(f"Loading taxon map from: {args.map}")
    taxon_map = load_taxon_map(args.map)
    print(f"Loaded {len(taxon_map)} taxa")
    
    # Count total copies
    total_copies = sum(len(copies) for copies in taxon_map.values())
    print(f"Total copies across all taxa: {total_copies}")
    
    # Run MP-SUGAR
    run_mpsugar(
        args.trees,
        taxon_map,
        args.output,
        iter_ct=args.iterations,
        num_chains=args.chains,
        seed=args.seed
    )
    
    print("\nMP-SUGAR run completed successfully!")
    return 0


if __name__ == "__main__":
    sys.exit(main())