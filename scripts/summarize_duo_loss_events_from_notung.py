#!/usr/bin/env python3
"""
Summarize Notung reconciliation statistics across multiple datasets.
Calculates average duplications and losses per dataset.
"""

import os
import re
import csv
import statistics
from pathlib import Path

# Define all datasets
DATASETS = [
    "Koenen_2020", "Lawrence_2016", "Diaz-Perez_2018", "Wisecaver_2023",
    "Ding_2023", "Popp_2005", "Wu_2015", "Ren_2024", "Zhao_2021",
    "Marcussen_2015", "Morales-Briones_2021"
]

# Base directory
BASE_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/papers"

# Output file
OUTPUT_FILE = "/groups/itay_mayrose/tomulanovski/gene2net/dup_loss_summary.csv"


def parse_stats_file(stats_file):
    """
    Parse a single .stats.txt file and extract duplications and losses.
    
    Returns:
        tuple: (duplications, losses) or (None, None) if parsing fails
    """
    duplications = None
    losses = None
    
    try:
        with open(stats_file, 'r') as f:
            content = f.read()
            
            # Look for duplications
            dup_match = re.search(r'- Duplications:\s+(\d+)', content)
            if dup_match:
                duplications = int(dup_match.group(1))
            
            # Look for losses
            loss_match = re.search(r'- Losses:\s+(\d+)', content)
            if loss_match:
                losses = int(loss_match.group(1))
                
    except Exception as e:
        print(f"Error parsing {stats_file}: {e}")
        return None, None
    
    return duplications, losses


def summarize_dataset(dataset_name):
    """
    Summarize statistics for a single dataset.
    
    Returns:
        dict: Summary statistics including averages and counts
    """
    reconciled_dir = os.path.join(BASE_DIR, dataset_name, "reconciled_for_dup")
    
    if not os.path.exists(reconciled_dir):
        print(f"WARNING: Directory not found for {dataset_name}: {reconciled_dir}")
        return None
    
    # Find all .stats.txt files
    stats_files = list(Path(reconciled_dir).glob("*.stats.txt"))
    
    if not stats_files:
        print(f"WARNING: No stats files found for {dataset_name}")
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
        print(f"WARNING: No valid statistics found for {dataset_name}")
        return None
    
    # Calculate averages
    avg_duplications = sum(duplications_list) / len(duplications_list)
    avg_losses = sum(losses_list) / len(losses_list)
    num_trees = len(duplications_list)
    
    summary = {
        'dataset': dataset_name,
        'num_trees': num_trees,
        'avg_duplications': round(avg_duplications, 2),
        'avg_losses': round(avg_losses, 2),
        'min_duplications': min(duplications_list),
        'max_duplications': max(duplications_list),
        'min_losses': min(losses_list),
        'max_losses': max(losses_list)
    }
    
    return summary


def main():
    """Main function to process all datasets and write summary CSV."""
    
    print("=" * 60)
    print("Notung Statistics Summary")
    print("=" * 60)
    print()
    
    results = []
    
    # Process each dataset
    for dataset in DATASETS:
        print(f"Processing {dataset}...")
        summary = summarize_dataset(dataset)
        
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
        print(f"Writing results to {OUTPUT_FILE}")
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        
        # Calculate mean and std across all datasets
        avg_dups_list = [r['avg_duplications'] for r in results]
        avg_losses_list = [r['avg_losses'] for r in results]
        
        mean_duplications = statistics.mean(avg_dups_list)
        std_duplications = statistics.stdev(avg_dups_list) if len(avg_dups_list) > 1 else 0
        mean_losses = statistics.mean(avg_losses_list)
        std_losses = statistics.stdev(avg_losses_list) if len(avg_losses_list) > 1 else 0
        
        with open(OUTPUT_FILE, 'w', newline='') as csvfile:
            fieldnames = [
                'dataset', 'num_trees',
                'avg_duplications', 'avg_losses',
                'min_duplications', 'max_duplications',
                'min_losses', 'max_losses'
            ]
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            # Write individual dataset results
            for result in results:
                writer.writerow(result)
            
            # Add empty row for separation
            writer.writerow({})
            
            # Add mean row
            writer.writerow({
                'dataset': 'MEAN',
                'num_trees': '',
                'avg_duplications': round(mean_duplications, 2),
                'avg_losses': round(mean_losses, 2),
                'min_duplications': '',
                'max_duplications': '',
                'min_losses': '',
                'max_losses': ''
            })
            
            # Add std row
            writer.writerow({
                'dataset': 'STD',
                'num_trees': '',
                'avg_duplications': round(std_duplications, 2),
                'avg_losses': round(std_losses, 2),
                'min_duplications': '',
                'max_duplications': '',
                'min_losses': '',
                'max_losses': ''
            })
        
        print()
        print("Summary Statistics Across All Datasets:")
        print()
        print("=" * 60)
        print(f"SUCCESS: Summary written to {OUTPUT_FILE}")
        print(f"Total datasets processed: {len(results)}/{len(DATASETS)}")
        print("=" * 60)
    else:
        print("ERROR: No results to write")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())