#!/usr/bin/env python3
"""
Scale branch lengths in all ultrametric species trees across all networks.
Uses high-precision decimal arithmetic to maintain ultrametricity.
Applies ETE3 convert_to_ultrametric() after scaling to ensure ultrametricity is preserved.

Usage:
    python scale_ultrametric_species_trees.py <base_dir> <scale_factor>
    
Example:
    python scale_ultrametric_species_trees.py /groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations 1000000
"""
import sys
import os
from pathlib import Path
from io import StringIO
from Bio import Phylo
from decimal import Decimal, getcontext
from ete3 import Tree

# Set high precision for decimal calculations
getcontext().prec = 50

def scale_tree(tree, factor):
    """Multiply all branch lengths by factor using high-precision arithmetic."""
    factor_decimal = Decimal(str(factor))
    
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            # Convert to Decimal for precision, scale, then back to float
            original = Decimal(str(clade.branch_length))
            scaled = original * factor_decimal
            clade.branch_length = float(scaled)
    
    return tree

def ensure_ultrametric_with_ete3(newick_str):
    """
    Convert tree to ultrametric using ETE3.
    Takes a Newick string, applies convert_to_ultrametric(), and returns the result.
    """
    try:
        # Load tree with ETE3
        ete_tree = Tree(newick_str, format=1)
        
        # Convert to ultrametric
        ete_tree.convert_to_ultrametric()
        
        # Return the ultrametric newick string
        return ete_tree.write(format=1)
    except Exception as e:
        print(f"    Warning: ETE3 ultrametric conversion failed: {e}")
        return newick_str

def write_simphy_nexus(tree, output_file):
    """Write tree in SimPhy Nexus format with high-precision branch lengths."""
    handle = StringIO()
    Phylo.write([tree], handle, "newick")
    newick_str = handle.getvalue().strip().replace("\n", "")
    
    # Ensure it ends with a semicolon
    if not newick_str.endswith(";"):
        newick_str += ";"
    
    # Apply ETE3 ultrametric conversion
    print("    Applying ETE3 convert_to_ultrametric()...")
    newick_str = ensure_ultrametric_with_ete3(newick_str)
    
    # Ensure it ends with a semicolon after ETE3 processing
    if not newick_str.endswith(";"):
        newick_str += ";"
    
    # Remove root branch length if present (SimPhy requirement)
    import re
    newick_str = re.sub(r'\):0\.0+;', ');', newick_str)
    newick_str = re.sub(r'\):[0-9]+\.0+;', ');', newick_str)
    newick_str = re.sub(r'\):[0-9.]+;$', ');', newick_str)
    
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
        
        # Write to output file (includes ETE3 ultrametric conversion)
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
        input_file = base_path / network / "species_tree_ultrametric_normalized.nex"
        output_file = base_path / network / "species_tree_ultrametric_height_1_million.nex"
        
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
    print(f"  Output files: species_tree_ultrametric_height_1_million.nex")
    print(f"  Precision: {getcontext().prec} decimal places")
    print(f"  Ultrametric conversion: ETE3 convert_to_ultrametric()")
    print(f"{'='*60}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python scale_ultrametric_species_trees.py <base_dir> <scale_factor>")
        print("\nExample:")
        print("  python scale_ultrametric_species_trees.py /groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations 1000000")
        print("\nThis will create species_tree_ultrametric_height_1_million.nex for each network")
        print("Uses high-precision decimal arithmetic and ETE3 convert_to_ultrametric().")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    scale_factor = float(sys.argv[2])
    
    if not os.path.isdir(base_dir):
        print(f"ERROR: Base directory not found: {base_dir}")
        sys.exit(1)
    
    print(f"Base directory: {base_dir}")
    print(f"Scale factor: {scale_factor}")
    print(f"Decimal precision: {getcontext().prec} places")
    print(f"Processing ultrametric species trees with ETE3 conversion...")
    
    process_all_networks(base_dir, scale_factor)

if __name__ == "__main__":
    main()