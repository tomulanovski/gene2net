#!/usr/bin/env python3
"""
Summarize duplications and losses from Locus_Trees.csv
Usage: python summarize_dup_loss.py /path/to/directory/with/Locus_Trees.csv
"""

import csv
import sys
from pathlib import Path

def summarize_dup_loss(csv_path):
    """Calculate average n_dup and n_loss from Locus_Trees.csv"""
    
    if not Path(csv_path).exists():
        print(f"ERROR: File not found: {csv_path}")
        return None, None, None
    
    total_dup = 0
    total_loss = 0
    count = 0
    
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            try:
                total_dup += float(row['n_dup'])
                total_loss += float(row['n_loss'])
                count += 1
            except (KeyError, ValueError) as e:
                print(f"WARNING: Could not parse row: {e}")
                continue
    
    if count == 0:
        print(f"WARNING: No valid data found in {csv_path}")
        return None, None, None
    
    avg_dup = total_dup / count
    avg_loss = total_loss / count
    
    return avg_dup, avg_loss, count

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python summarize_dup_loss.py /path/to/directory")
        sys.exit(1)
    
    data_dir = Path(sys.argv[1])
    csv_path = data_dir / "Locus_Trees.csv"
    
    print(f"Processing: {csv_path}")
    
    avg_dup, avg_loss, count = summarize_dup_loss(csv_path)
    
    if avg_dup is None:
        sys.exit(1)
    
    # Write summary to file
    summary_path = data_dir / "dup_loss_summary.txt"
    with open(summary_path, 'w') as f:
        f.write(f"Duplication and Loss Summary\n")
        f.write(f"={'='*40}\n")
        f.write(f"Number of trees: {count}\n")
        f.write(f"Average duplications: {avg_dup:.4f}\n")
        f.write(f"Average losses: {avg_loss:.4f}\n")
    
    print(f"Summary written to: {summary_path}")
    print(f"  Trees: {count}")
    print(f"  Avg duplications: {avg_dup:.4f}")
    print(f"  Avg losses: {avg_loss:.4f}")