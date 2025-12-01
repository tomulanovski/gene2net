#!/usr/bin/env python3
"""
Summarize Notung reconciliation statistics for simulated datasets.
Calculates average duplications and losses per network.

Usage:
    python summarize_notung_stats.py <simulation_condition>
    
Example:
    python summarize_notung_stats.py low_dup_ne_100000_height_1_million_trees
    python summarize_notung_stats.py med_dup_ne_100000_height_50_million_trees
"""

import os
import re
import csv
import sys
import statistics
from pathlib import Path

# Define all networks
NETWORKS = [
    "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
    "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019", 
    "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011", 
    "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014", 
    "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014"
]

# Fixed paths
BASE_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
OUTPUT_BASE_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis"


def parse_stats_file(stats_file):
    """
    Parse a single .stats file and extract duplications and losses.
    
    Returns:
        tuple: (duplications, losses) or (None, None) if parsing fails
    """
    duplications = None
    losses = None
    
    try:
        with open(stats_file, 'r') as f:
            content = f.read()
            
            # Look for duplications - try multiple patterns
            dup_match = re.search(r'[-\s]*Duplications:\s+(\d+)', content, re.IGNORECASE)
            if dup_match:
                duplications = int(dup_match.group(1))
            
            # Look for losses
            loss_match = re.search(r'[-\s]*Losses:\s+(\d+)', content, re.IGNORECASE)
            if loss_match:
                losses = int(loss_match.group(1))
                
    except Exception as e:
        print(f"Error parsing {stats_file}: {e}")
        return None, None
    
    return duplications, losses


def summarize_network(network_name, simulation_condition):
    """
    Summarize statistics for a single network.
    
    Args:
        network_name: Name of the network
        simulation_condition: Name of simulation condition folder
    
    Returns:
        dict: Summary statistics including averages and counts
    """
    reconciled_dir = os.path.join(
        BASE_DIR, 
        network_name, 
        f"data/{simulation_condition}/reconciled_for_dup"
    )
    
    if not os.path.exists(reconciled_dir):
        print(f"WARNING: Directory not found for {network_name}: {reconciled_dir}")
        return None
    
    # Find all .stats files (Notung creates files without .txt extension)
    stats_files = []
    for pattern in ["*.stats", "*.stats.txt"]:
        stats_files.extend(list(Path(reconciled_dir).glob(pattern)))
    
    if not stats_files:
        print(f"WARNING: No stats files found for {network_name}")
        return None
    
    # Collect statistics
    duplications_list = []
    losses_list = []
    
    for stats_file in stats_files:
        dups, losses = parse_stats_file(stats_file)
        
        if dups is not None and losses is not None:
            duplications_list.append(dups)
            losses_list.append(losses)
    
    if not duplications_list:
        print(f"WARNING: No valid statistics found for {network_name}")
        return None
    
    # Calculate statistics
    avg_duplications = sum(duplications_list) / len(duplications_list)
    avg_losses = sum(losses_list) / len(losses_list)
    num_trees = len(duplications_list)
    
    summary = {
        'network': network_name,
        'num_trees': num_trees,
        'avg_duplications': round(avg_duplications, 2),
        'avg_losses': round(avg_losses, 2),
        'min_duplications': min(duplications_list),
        'max_duplications': max(duplications_list),
        'min_losses': min(losses_list),
        'max_losses': max(losses_list),
        'median_duplications': round(statistics.median(duplications_list), 2),
        'median_losses': round(statistics.median(losses_list), 2)
    }
    
    return summary


def main():
    """Main function to process all networks and write summary CSV."""
    
    # Parse command-line arguments
    if len(sys.argv) != 2:
        print("Usage: python summarize_notung_stats.py <simulation_condition>")
        print()
        print("Example:")
        print("  python summarize_notung_stats.py low_dup_ne_100000_height_1_million_trees")
        print("  python summarize_notung_stats.py med_dup_ne_100000_height_50_million_trees")
        return 1
    
    simulation_condition = sys.argv[1]
    output_file = os.path.join(OUTPUT_BASE_DIR, f"notung_{simulation_condition}_summary.csv")
    
    print("=" * 80)
    print("Notung Statistics Summary - Simulated Data")
    print("=" * 80)
    print(f"Simulation condition: {simulation_condition}")
    print(f"Output file: {output_file}")
    print("=" * 80)
    print()
    
    results = []
    
    # Process each network
    for network in NETWORKS:
        print(f"Processing {network}...")
        summary = summarize_network(network, simulation_condition)
        
        if summary:
            results.append(summary)
            print(f"  Trees: {summary['num_trees']}")
            print(f"  Avg Duplications: {summary['avg_duplications']}")
            print(f"  Avg Losses: {summary['avg_losses']}")
            print()
        else:
            print(f"  No data found or error occurred")
            print()
    
    # Write results to CSV
    if results:
        print(f"Writing results to {output_file}")
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Calculate mean and std across all networks
        avg_dups_list = [r['avg_duplications'] for r in results]
        avg_losses_list = [r['avg_losses'] for r in results]
        
        mean_duplications = statistics.mean(avg_dups_list)
        std_duplications = statistics.stdev(avg_dups_list) if len(avg_dups_list) > 1 else 0
        mean_losses = statistics.mean(avg_losses_list)
        std_losses = statistics.stdev(avg_losses_list) if len(avg_losses_list) > 1 else 0
        
        with open(output_file, 'w', newline='') as csvfile:
            fieldnames = [
                'network', 'num_trees',
                'avg_duplications', 'avg_losses',
                'median_duplications', 'median_losses',
                'min_duplications', 'max_duplications',
                'min_losses', 'max_losses'
            ]
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            # Write individual network results
            for result in results:
                writer.writerow(result)
            
            # Add empty row for separation
            writer.writerow({})
            
            # Add mean row
            writer.writerow({
                'network': 'MEAN',
                'num_trees': '',
                'avg_duplications': round(mean_duplications, 2),
                'avg_losses': round(mean_losses, 2),
                'median_duplications': '',
                'median_losses': '',
                'min_duplications': '',
                'max_duplications': '',
                'min_losses': '',
                'max_losses': ''
            })
            
            # Add std row
            writer.writerow({
                'network': 'STD',
                'num_trees': '',
                'avg_duplications': round(std_duplications, 2),
                'avg_losses': round(std_losses, 2),
                'median_duplications': '',
                'median_losses': '',
                'min_duplications': '',
                'max_duplications': '',
                'min_losses': '',
                'max_losses': ''
            })
        
        print()
        print("Summary Statistics Across All Networks:")
        print(f"Mean Duplications: {mean_duplications:.2f} ± {std_duplications:.2f}")
        print(f"Mean Losses: {mean_losses:.2f} ± {std_losses:.2f}")
        print()
        print("=" * 80)
        print(f"SUCCESS: Summary written to {output_file}")
        print(f"Total networks processed: {len(results)}/{len(NETWORKS)}")
        print("=" * 80)
    else:
        print("ERROR: No results to write")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
