#!/usr/bin/env python3
"""
Script to generate taxa mapping as Python dictionary from NEXUS files.
Fixed to correctly handle species names containing periods.

Usage: python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/python_taxa_map.py" -i input_trees.nex -o taxamap.py
"""

import re
import sys
import os
import argparse

def extract_taxa_from_nexus(nexus_file):
    """Extract all taxa names from a NEXUS file."""
    try:
        with open(nexus_file, 'r') as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading file {nexus_file}: {e}")
        return set()
    
    # Find the TREES block
    trees_match = re.search(r'BEGIN TREES;(.+?)END;', content, re.DOTALL | re.IGNORECASE)
    if not trees_match:
        print(f"Warning: Could not find TREES block in {nexus_file}, treating as raw Newick")
        return extract_taxa_from_newick_content(content)
    
    trees_block = trees_match.group(1)
    
    # Extract all tree definitions
    tree_defs = re.findall(r'Tree\s+\w+\s*=\s*([^;]+);', trees_block, re.IGNORECASE)
    
    # Extract taxa names from Newick strings
    all_taxa = set()
    for tree in tree_defs:
        tree_taxa = extract_taxa_from_newick_content(tree)
        all_taxa.update(tree_taxa)
    
    # Remove potential keywords that are not taxa
    keywords = {'tree', 'end', 'begin', 'translate', 'nexus', 'trees'}
    all_taxa = {t for t in all_taxa if t.lower() not in keywords}
    
    return all_taxa

def extract_taxa_from_newick_content(content):
    """Extract taxa from Newick string content."""
    # Remove comments and quoted labels
    clean_content = re.sub(r'\[[^\]]*\]', '', content)
    
    # Multiple regex patterns to catch all possible taxa formats
    patterns = [
        r'(?<=[(:,])([A-Za-z][A-Za-z0-9_.]*)(?=[:),])',  # Standard pattern, including periods
        r'^([A-Za-z][A-Za-z0-9_.]*)(?=[:),])',           # At beginning, including periods
        r'[\'"]([A-Za-z][A-Za-z0-9_.]*)[\'"]',           # With quotes, including periods
        r'(?<=[(:,])([A-Za-z][A-Za-z0-9_.+-]*)(?=[:),])' # With special chars, including periods
    ]
    
    all_taxa = set()
    for pattern in patterns:
        taxa = re.findall(pattern, clean_content)
        all_taxa.update(taxa)
    
    # Filter out any completely numeric strings
    all_taxa = {t for t in all_taxa if not re.match(r'^\d+$', t)}
    
    return all_taxa

def group_taxa_by_root_name(taxa_list):
    """
    Group taxa by their root name, properly handling species names with periods.
    The pattern specifically looks for underscore followed by number at the end.
    """
    taxa_dict = {}
    base_names = set()
    
    # Pattern specifically looks for underscore followed by number at the end
    # This preserves periods in species names
    pattern = re.compile(r'^(.+)_(\d+)$')
    
    # First pass: identify all base names and collect all taxa
    for taxon in taxa_list:
        match = pattern.match(taxon)
        if match:
            base_name, _ = match.groups()
            base_names.add(base_name)
            
            if base_name not in taxa_dict:
                taxa_dict[base_name] = []
            taxa_dict[base_name].append(taxon)
        else:
            # Regular taxon without copy number
            if taxon not in taxa_dict:
                taxa_dict[taxon] = [taxon]
            else:
                if taxon not in taxa_dict[taxon]:
                    taxa_dict[taxon].append(taxon)
    
    # Second pass: ensure base names are included in their own mapping
    for base_name in base_names:
        # If the base name itself is also a taxon, make sure it's in its mapping
        if base_name in taxa_list and base_name not in taxa_dict[base_name]:
            taxa_dict[base_name].append(base_name)
        # If the base name doesn't exist as a taxon but has copies, add it anyway
        elif base_name not in taxa_list:
            taxa_dict[base_name].append(base_name)
    
    return taxa_dict

def generate_python_dict_string(taxa_dict):
    """Generate the taxa mapping as a Python dictionary string."""
    if not taxa_dict:
        print("Error: No taxa found to create mapping")
        return "taxa_map = {}"
    
    lines = ["taxa_map = {"]
    
    for species, copies in sorted(taxa_dict.items()):
        # Ensure no empty copies and sort them
        copies = sorted([c for c in copies if c])
        if not copies:
            continue
        
        # Format the list of taxa for this species
        taxa_list_str = "[" + ", ".join([f"'{taxon}'" for taxon in copies]) + "]"
        lines.append(f"    '{species}': {taxa_list_str},")
    
    lines.append("}")
    
    return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser(description='Generate taxa mapping as Python dictionary from NEXUS or Newick files')
    parser.add_argument('-i', '--input', required=True, nargs='+', help='Input files containing trees')
    parser.add_argument('-o', '--output', required=True, help='Output file for the taxa mapping (Python format)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show detailed output')
    
    args = parser.parse_args()
    
    all_taxa = set()
    for tree_file in args.input:
        if not os.path.isfile(tree_file):
            print(f"Error: File {tree_file} does not exist")
            continue
        
        file_taxa = extract_taxa_from_nexus(tree_file)
        if args.verbose:
            print(f"Found {len(file_taxa)} taxa in {tree_file}")
        
        all_taxa.update(file_taxa)
    
    if not all_taxa:
        print("Error: No taxa found in the provided files")
        sys.exit(1)
    
    # Group taxa ensuring base names are included
    taxa_dict = group_taxa_by_root_name(all_taxa)
    
    # Generate the Python dictionary string
    dict_string = generate_python_dict_string(taxa_dict)
    
    # Write to output file
    try:
        with open(args.output, 'w') as f:
            f.write(dict_string + '\n')
        print(f"Taxa mapping (Python dict) successfully written to {args.output}")
    except Exception as e:
        print(f"Error writing to {args.output}: {e}")
        sys.exit(1)
    
    print(f"Found {len(taxa_dict)} species with a total of {len(all_taxa)} unique taxa")
    
    # Print the dictionary to console as well
    print("\nGenerated taxa mapping:")
    print(dict_string)
    
    # Print details of the mapping
    if args.verbose:
        multi_copy_species = {sp: copies for sp, copies in taxa_dict.items() if len(copies) > 1}
        if multi_copy_species:
            print(f"\nFound {len(multi_copy_species)} potential polyploid species:")
            for species, copies in sorted(multi_copy_species.items()):
                print(f"  {species}: {', '.join(sorted(copies))}")
    
    # Write full taxa list to a separate file for reference
    taxa_list_file = os.path.splitext(args.output)[0] + "_full_taxa_list.txt"
    try:
        with open(taxa_list_file, 'w') as f:
            for taxon in sorted(all_taxa):
                f.write(f"{taxon}\n")
        print(f"Full taxa list written to {taxa_list_file}")
    except Exception as e:
        print(f"Error writing full taxa list: {e}")

if __name__ == "__main__":
    main()