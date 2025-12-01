#!/usr/bin/env python3
"""
Normalize all ultrametric species trees to have height = 1.0
(distance from root to tips = 1.0)

Usage:
    python normalize_tree_heights.py <base_dir>
    
Example:
    python normalize_tree_heights.py /groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations
"""
import sys
import os
from pathlib import Path
from io import StringIO
from Bio import Phylo

def get_root_to_tip_distance(tree):
    """
    Calculate distance from root to any tip (assumes ultrametric tree).
    Returns the maximum distance from root to any leaf.
    """
    def get_distance_to_root(clade, distance=0):
        if clade.branch_length is not None:
            distance += clade.branch_length
        return distance
    
    # Get all terminal nodes (tips)
    terminals = tree.get_terminals()
    if not terminals:
        return None
    
    # Calculate distance from root to first tip
    # (in ultrametric trees, all tips should be equidistant from root)
    distances = []
    for terminal in terminals:
        path = tree.get_path(terminal)
        distance = sum(clade.branch_length for clade in path if clade.branch_length is not None)
        distances.append(distance)
    
    # Return the maximum distance (should all be equal in ultrametric tree)
    return max(distances) if distances else None

def scale_tree(tree, factor):
    """Multiply all branch lengths by factor."""
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length *= factor
    return tree

def write_simphy_nexus(tree, output_file):
    """Write tree in SimPhy Nexus format with no root branch length."""
    handle = StringIO()
    Phylo.write([tree], handle, "newick")
    newick_str = handle.getvalue().strip().replace("\n", "")
    
    # Ensure it ends with a semicolon
    if not newick_str.endswith(";"):
        newick_str += ";"
    
    # Remove root branch length if present (SimPhy requirement)
    newick_str = newick_str.replace("):0.0;", ");")
    newick_str = newick_str.replace("):0.00000;", ");")
    newick_str = newick_str.replace("):0.000000;", ");")
    
    with open(output_file, "w") as out:
        out.write("#NEXUS\n")
        out.write("begin trees;\n")
        out.write(f"tree 1={newick_str}\n")
        out.write("end;\n")

def process_species_tree(input_file, output_file):
    """Process a single species tree file and normalize to height = 1.0"""
    try:
        # Read the tree
        trees = list(Phylo.parse(input_file, "nexus"))
        if not trees:
            print(f"  ??  No trees found in {input_file}")
            return False, None, None
        
        tree = trees[0]
        
        # Calculate current root-to-tip distance
        current_height = get_root_to_tip_distance(tree)
        if current_height is None or current_height == 0:
            print(f"  ??  Could not calculate tree height")
            return False, None, None
        
        # Calculate scaling factor to make height = 1.0
        scale_factor = 1.0 / current_height
        
        # Scale the tree
        scaled_tree = scale_tree(tree, scale_factor)
        
        # Verify the new height
        new_height = get_root_to_tip_distance(scaled_tree)
        
        # Write to output file
        write_simphy_nexus(scaled_tree, output_file)
        
        return True, current_height, new_height
        
    except Exception as e:
        print(f"  ?? Error processing {input_file}: {e}")
        return False, None, None

def process_all_networks(base_dir):
    """Process all ultrametric species trees in the base directory."""
    networks = [
        "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
        "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019", 
        "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011", 
        "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014", 
        "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014"
    ]
    
    base_path = Path(base_dir)
    networks_processed = 0
    results = []
    
    for network in networks:
        input_file = base_path / network / "species_tree_ultrametric.nex"
        output_file = base_path / network / "species_tree_ultrametric_normalized.nex"
        
        if not input_file.exists():
            print(f"??  {network}: Species tree not found: {input_file}")
            continue
        
        print(f"\nProcessing: {network}")
        print(f"  Input:  {input_file}")
        print(f"  Output: {output_file}")
        
        # Process the species tree
        success, old_height, new_height = process_species_tree(input_file, output_file)
        
        if success:
            print(f"  ? Original height: {old_height:.6f}")
            print(f"  ? New height: {new_height:.6f}")
            print(f"  ? Scale factor: {1.0/old_height:.6f}")
            networks_processed += 1
            results.append((network, old_height, new_height))
        else:
            print(f"  ? Failed to process species tree")
    
    print(f"\n{'='*60}")
    print(f"SUMMARY:")
    print(f"  Networks processed: {networks_processed}/{len(networks)}")
    print(f"  Output files: species_tree_ultrametric_normalized.nex")
    print(f"\nTree heights:")
    for network, old_h, new_h in results:
        print(f"  {network:30s}: {old_h:10.6f} ? {new_h:.6f}")
    print(f"{'='*60}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python normalize_tree_heights.py <base_dir>")
        print("\nExample:")
        print("  python normalize_tree_heights.py /groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations")
        print("\nThis will create species_tree_ultrametric_normalized.nex for each network")
        print("with tree height normalized to 1.0 (root-to-tip distance = 1.0)")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    
    if not os.path.isdir(base_dir):
        print(f"ERROR: Base directory not found: {base_dir}")
        sys.exit(1)
    
    print(f"Base directory: {base_dir}")
    print(f"Normalizing ultrametric species trees to height = 1.0...")
    
    process_all_networks(base_dir)

if __name__ == "__main__":
    main()