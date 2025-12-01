#!/usr/bin/env python3
import sys
import re
from pathlib import Path

def extract_species_name(gene_name, species_list):
    """
    Match gene name to species by finding longest matching prefix from species list.
    For example: Albiziajulibrissin_1_0_0 matches Albiziajulibrissin_1
    """
    # Sort species by length (longest first) to match longest prefix
    sorted_species = sorted(species_list, key=len, reverse=True)
    
    for species in sorted_species:
        if gene_name.startswith(species + "_"):
            return species
    
    # If no match found, print warning
    print(f"WARNING: Could not find species for gene: {gene_name}", file=sys.stderr)
    return None

def get_species_from_tree(tree_string):
    """Extract species names from species tree."""
    # Find all leaf names (text before colons)
    species = re.findall(r'([A-Za-z0-9_]+):', tree_string)
    return species

def add_nhx_tags(gene_tree_string, species_list):
    """Add NHX species tags to gene tree leaf names."""
    
    def replace_leaf(match):
        gene_name = match.group(1)
        species = extract_species_name(gene_name, species_list)
        
        if species:
            # Add NHX tag
            return f"{gene_name}[&&NHX:S={species}]:"
        else:
            # Keep original if no match
            return match.group(0)
    
    # Pattern to match leaf names (everything before a colon that's not after a closing paren)
    pattern = r'([A-Za-z0-9_]+):'
    
    modified_tree = re.sub(pattern, replace_leaf, gene_tree_string)
    return modified_tree

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python add_nhx_tags.py <species_tree.tre> <input_gene_tree.trees> <output_gene_tree.trees>")
        sys.exit(1)
    
    species_tree_file = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Read species tree and extract species names
    with open(species_tree_file, 'r') as f:
        species_tree = f.read().strip()
    
    species_list = get_species_from_tree(species_tree)
    print(f"Found {len(species_list)} species: {species_list[:5]}...", file=sys.stderr)
    
    # Read and convert gene tree
    with open(input_file, 'r') as f:
        gene_tree = f.read().strip()
    
    # Add NHX tags
    modified_tree = add_nhx_tags(gene_tree, species_list)
    
    # Write output
    with open(output_file, 'w') as f:
        f.write(modified_tree)
        if not modified_tree.endswith('\n'):
            f.write('\n')
    
    print(f"Converted {input_file} -> {output_file}", file=sys.stderr)
