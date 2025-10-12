#!/usr/bin/env python3
"""
Scale branch lengths in all ultrametric species trees across all networks.
Usage:
    python scale_ultrametric_species_trees.py <base_dir> <scale_factor>
    
Example:
    python scale_ultrametric_species_trees.py /groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations 0.001
"""
import sys
import os
from pathlib import Path
from io import StringIO
from Bio import Phylo

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

def process_species_tree(input_file, output_file, scale_factor):
    """Process a single species tree file."""
    try:
        # Read the tree
        trees = list(Phylo.parse(input_file, "nexus"))
        if not trees:
            print(f"  ??  No trees found in {input_file}")
            return False
        
        # Scale the first tree
        tree = trees[0]
        scaled_tree = scale_tree(tree, scale_factor)
        
        # Write to output file
        write_simphy_nexus(scaled_tree, output_file)
        
        return True
        
    except Exception as e:
        print(f"  ? Error processing {input_file}: {e}")
        return False

def process_all_networks(base_dir, scale_factor):
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
    
    for network in networks:
        input_file = base_path / network / "species_tree_ultrametric.nex"
        output_file = base_path / network / "species_tree_ultrametric_scaled.nex"
        
        if not input_file.exists():
            print(f"??  {network}: Species tree not found: {input_file}")
            continue
        
        print(f"\nProcessing: {network}")
        print(f"  Input:  {input_file}")
        print(f"  Output: {output_file}")
        
        # Process the species tree
        success = process_species_tree(input_file, output_file, scale_factor)
        
        if success:
            print(f"  ? Successfully scaled species tree")
            networks_processed += 1
        else:
            print(f"  ? Failed to process species tree")
    
    print(f"\n{'='*60}")
    print(f"SUMMARY:")
    print(f"  Networks processed: {networks_processed}/{len(networks)}")
    print(f"  Scale factor: {scale_factor}")
    print(f"  Output files: species_tree_ultrametric_scaled.nex")
    print(f"{'='*60}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python scale_ultrametric_species_trees.py <base_dir> <scale_factor>")
        print("\nExample:")
        print("  python scale_ultrametric_species_trees.py /groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations 0.001")
        print("\nThis will create species_tree_ultrametric_scaled.nex for each network")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    scale_factor = float(sys.argv[2])
    
    if not os.path.isdir(base_dir):
        print(f"ERROR: Base directory not found: {base_dir}")
        sys.exit(1)
    
    print(f"Base directory: {base_dir}")
    print(f"Scale factor: {scale_factor}")
    print(f"Processing ultrametric species trees...")
    
    process_all_networks(base_dir, scale_factor)

if __name__ == "__main__":
    main()