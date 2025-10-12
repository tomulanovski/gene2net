#!/usr/bin/env python3
"""
Process species trees: rescale all branches to 1, convert to ultrametric, save as NEXUS.
Usage: python rescale_and_ultrametric.py <base_dir>
"""
from ete3 import Tree
import sys
import os
from pathlib import Path

def read_nexus_tree(filepath):
    """Read a tree from NEXUS file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Extract Newick string from NEXUS format
        # Look for the tree definition line
        for line in content.split('\n'):
            line = line.strip()
            if line.startswith('tree ') and '=' in line:
                # Extract the Newick part after '='
                newick_str = line.split('=', 1)[1].strip()
                # Remove any trailing semicolon and whitespace
                newick_str = newick_str.rstrip(';').strip()
                # Parse with ETE3
                tree = Tree(newick_str + ';', format=1)
                return tree
        
        raise ValueError("No tree definition found in NEXUS file")
    except Exception as e:
        print(f"Error reading NEXUS file {filepath}: {e}")
        sys.exit(1)

def rescale_branches_to_one(tree):
    """Set all branch lengths to 1."""
    for node in tree.traverse():
        if not node.is_root():
            node.dist = 1.0
    return tree

def convert_to_ultrametric(tree):
    """Convert tree to ultrametric (all leaves equidistant from root)."""
    tree.convert_to_ultrametric()
    return tree

def tree_to_nexus_newick(tree):
    """Convert ETE3 tree to Newick string with only taxa names and branch lengths."""
    def format_node(node):
        # Get node name (only for leaves)
        name = node.name if node.is_leaf() and node.name else ""
        
        # Build the node string
        if node.is_leaf():
            result = name
        else:
            children = ",".join(format_node(child) for child in node.children)
            result = f"({children})"
        
        # Add branch length
        if not node.is_root():
            result += f":{node.dist:.6g}"
        
        return result
    
    return format_node(tree) + ";"

def write_nexus(tree, output_filepath, tree_name="1"):
    """Write tree to NEXUS format file."""
    try:
        newick_str = tree_to_nexus_newick(tree)
        
        with open(output_filepath, 'w') as f:
            f.write("#NEXUS\n")
            f.write("begin trees;\n")
            f.write(f"tree {tree_name}={newick_str}\n")
            f.write("end;\n")
        
        print(f"Successfully wrote ultrametric tree to {output_filepath}")
    except Exception as e:
        print(f"Error writing NEXUS file: {e}")
        sys.exit(1)

def process_all_trees(base_dir):
    """Process all species trees in the base directory."""
    networks = [
        "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
        "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019", 
        "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011", 
        "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014", 
        "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014"
    ]
    
    base_path = Path(base_dir)
    processed_count = 0
    
    for network in networks:
        species_tree_path = base_path / network / "species_tree.nex"
        
        if not species_tree_path.exists():
            print(f"WARNING: Species tree not found: {species_tree_path}")
            continue
        
        print(f"\nProcessing: {network}")
        print(f"  Reading from: {species_tree_path}")
        
        # Read the original tree
        tree = read_nexus_tree(species_tree_path)
        print(f"  Original tree loaded")
        
        # Rescale all branches to 1
        tree = rescale_branches_to_one(tree)
        print(f"  All branches rescaled to 1")
        
        # Convert to ultrametric
        tree = convert_to_ultrametric(tree)
        print(f"  Converted to ultrametric")
        
        # Write the processed tree
        output_path = base_path / network / "species_tree_ultrametric.nex"
        write_nexus(tree, output_path, tree_name="1")
        
        processed_count += 1
    
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"Successfully processed {processed_count}/{len(networks)} trees")
    print(f"{'='*60}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python rescale_and_ultrametric.py <base_dir>")
        print("Example: python rescale_and_ultrametric.py /groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    
    if not os.path.isdir(base_dir):
        print(f"ERROR: Base directory not found: {base_dir}")
        sys.exit(1)
    
    print(f"Base directory: {base_dir}")
    print(f"Processing all species trees...")
    
    process_all_trees(base_dir)

if __name__ == "__main__":
    main()