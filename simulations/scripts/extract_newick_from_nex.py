#!/usr/bin/env python3
"""
Extract and properly format Newick trees from NEXUS files for Notung.
"""

import sys
import re
from pathlib import Path

# Array of networks
NETWORKS = [
    "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
    "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019", 
    "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011", 
    "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014", 
    "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014"
]

BASE_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"


def extract_newick_from_nexus(nexus_file):
    """Extract Newick tree string from NEXUS file."""
    with open(nexus_file, 'r') as f:
        content = f.read()
    
    # Find the tree line (starts with 'tree' keyword)
    tree_match = re.search(r'^\s*tree\s+[^=]+=\s*(.+?);', content, re.MULTILINE | re.DOTALL)
    
    if not tree_match:
        return None
    
    newick = tree_match.group(1).strip()
    
    # Remove any internal newlines or extra whitespace
    newick = re.sub(r'\s+', '', newick)
    
    # Ensure it ends with semicolon
    if not newick.endswith(';'):
        newick += ';'
    
    return newick


def main():
    print("=" * 60)
    print("Extracting Newick trees from NEXUS files")
    print("=" * 60)
    print()
    
    success_count = 0
    fail_count = 0
    
    for network in NETWORKS:
        nex_file = Path(BASE_DIR) / network / "species_tree_1_million.nex"
        tre_file = Path(BASE_DIR) / network / "species_tree_1_million.tre"
        
        print(f"Processing {network}...")
        
        if not nex_file.exists():
            print(f"  WARNING: NEXUS file not found: {nex_file}")
            fail_count += 1
            continue
        
        try:
            newick = extract_newick_from_nexus(nex_file)
            
            if newick is None:
                print(f"  ERROR: Could not extract tree from {nex_file}")
                fail_count += 1
                continue
            
            # Write to file
            with open(tre_file, 'w') as f:
                f.write(newick + '\n')
            
            print(f"  ✓ Created: {tre_file}")
            print(f"    Tree length: {len(newick)} characters")
            success_count += 1
            
        except Exception as e:
            print(f"  ERROR: {e}")
            fail_count += 1
    
    print()
    print("=" * 60)
    print(f"Summary: {success_count} succeeded, {fail_count} failed")
    print("=" * 60)
    
    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    exit(main())

