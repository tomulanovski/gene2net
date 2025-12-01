#!/usr/bin/env python3
"""
Extract alignment lengths from FASTA files across multiple datasets.
Outputs: CSV statistics, PNG visualization, and pickle file for sampling.
"""

import os
import pickle
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# Dataset paths
DATASETS = {
    'Zhao_2021': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Zhao_et_al_2021',
    'Ren_2024': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Ren_et_al_2024',
    'Morales_Briones_2021': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Morales-Briones_et_al_2021'
}

OUTPUT_DIR = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions')


def get_alignment_length(fasta_file):
    """Get alignment length from a FASTA file."""
    try:
        records = list(SeqIO.parse(fasta_file, 'fasta'))
        if records:
            return len(records[0].seq)
    except Exception as e:
        print(f"Error reading {fasta_file}: {e}")
    return None


def collect_alignment_lengths(base_dir):
    """Collect alignment lengths from all FASTA files in directory."""
    lengths = []
    fasta_extensions = ['.fasta', '.fa', '.fna', '.fas', '.aln-cln']
    
    base_path = Path(base_dir)
    for ext in fasta_extensions:
        for fasta_file in base_path.glob(f'*{ext}'):
            length = get_alignment_length(fasta_file)
            if length:
                lengths.append(length)
    
    return lengths


def calculate_statistics(lengths, dataset_name):
    """Calculate summary statistics for alignment lengths."""
    if not lengths:
        return None
    
    return {
        'Dataset': dataset_name,
        'N_genes': len(lengths),
        'Min': np.min(lengths),
        'Max': np.max(lengths),
        'Mean': np.mean(lengths),
        'Median': np.median(lengths),
        'Std': np.std(lengths),
        'Q1': np.percentile(lengths, 25),
        'Q3': np.percentile(lengths, 75)
    }


def create_visualizations(all_lengths, output_path):
    """Create visualization of alignment length distributions."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Prepare data for plotting
    datasets = []
    lengths = []
    for dataset_name, dataset_lengths in all_lengths.items():
        datasets.extend([dataset_name] * len(dataset_lengths))
        lengths.extend(dataset_lengths)
    
    df = pd.DataFrame({'Dataset': datasets, 'Length': lengths})
    
    # 1. Histograms for each dataset
    ax1 = axes[0, 0]
    for dataset_name, dataset_lengths in all_lengths.items():
        ax1.hist(dataset_lengths, alpha=0.6, label=dataset_name, bins=30)
    ax1.set_xlabel('Alignment Length (bp)', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Alignment Length Distribution by Dataset', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Box plots
    ax2 = axes[0, 1]
    df.boxplot(column='Length', by='Dataset', ax=ax2)
    ax2.set_xlabel('Dataset', fontsize=12)
    ax2.set_ylabel('Alignment Length (bp)', fontsize=12)
    ax2.set_title('Alignment Length Distributions', fontsize=14, fontweight='bold')
    plt.sca(ax2)
    plt.xticks(rotation=45, ha='right')
    
    # 3. Combined histogram
    ax3 = axes[1, 0]
    all_vals = [l for lengths in all_lengths.values() for l in lengths]
    ax3.hist(all_vals, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax3.axvline(np.mean(all_vals), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(all_vals):.1f}')
    ax3.axvline(np.median(all_vals), color='green', linestyle='--', linewidth=2, label=f'Median: {np.median(all_vals):.1f}')
    ax3.set_xlabel('Alignment Length (bp)', fontsize=12)
    ax3.set_ylabel('Frequency', fontsize=12)
    ax3.set_title('Combined Alignment Length Distribution', fontsize=14, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Violin plot
    ax4 = axes[1, 1]
    parts = ax4.violinplot([all_lengths[ds] for ds in all_lengths.keys()],
                           positions=range(len(all_lengths)),
                           showmeans=True, showmedians=True)
    ax4.set_xticks(range(len(all_lengths)))
    ax4.set_xticklabels(all_lengths.keys(), rotation=45, ha='right')
    ax4.set_ylabel('Alignment Length (bp)', fontsize=12)
    ax4.set_title('Alignment Length Distributions (Violin Plot)', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_path}")
    plt.close()


def main():
    """Main execution function."""
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Collect alignment lengths from all datasets
    print("Collecting alignment lengths from datasets...")
    all_lengths = {}
    stats_list = []
    
    for dataset_name, base_dir in DATASETS.items():
        print(f"\nProcessing {dataset_name}...")
        lengths = collect_alignment_lengths(base_dir)
        
        if lengths:
            all_lengths[dataset_name] = lengths
            stats = calculate_statistics(lengths, dataset_name)
            stats_list.append(stats)
            print(f"  Found {len(lengths)} alignments")
            print(f"  Length range: {min(lengths)} - {max(lengths)} bp")
            print(f"  Mean: {np.mean(lengths):.1f} bp, Median: {np.median(lengths):.1f} bp")
        else:
            print(f"  No FASTA files found in {base_dir}")
    
    # Calculate combined statistics
    combined_lengths = [l for lengths in all_lengths.values() for l in lengths]
    combined_stats = calculate_statistics(combined_lengths, 'Combined')
    stats_list.append(combined_stats)
    
    # Save statistics to CSV
    csv_path = OUTPUT_DIR / 'alignment_lengths.csv'
    df_stats = pd.DataFrame(stats_list)
    df_stats.to_csv(csv_path, index=False)
    print(f"\nStatistics saved to {csv_path}")
    
    # Save raw data for sampling
    pkl_path = OUTPUT_DIR / 'alignment_lengths.pkl'
    with open(pkl_path, 'wb') as f:
        pickle.dump(all_lengths, f)
    print(f"Raw data saved to {pkl_path}")
    
    # Create visualizations
    png_path = OUTPUT_DIR / 'alignment_lengths.png'
    create_visualizations(all_lengths, png_path)
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(df_stats.to_string(index=False))
    print("\nTo sample from the distribution in your code:")
    print("  import pickle, random")
    print(f"  with open('{pkl_path}', 'rb') as f:")
    print("      lengths = pickle.load(f)")
    print("  # Sample from combined distribution:")
    print("  combined = [l for ds in lengths.values() for l in ds]")
    print("  sampled_length = random.choice(combined)")


if __name__ == '__main__':
    main()
