#!/usr/bin/env python3
"""
Simple AlloppNET Taxa Table Generator

Creates taxa tables with a straightforward strategy:
- Diploids: All copies as different individuals with genome A
- Polyploids/Unknown with multiple copies: Random homeolog pairs (A,B)
- if odd number of sequences for polyploids then it will create another sequence with _missing suffix to pair a homeolog

Usage:
    python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/create_taxa_table.py" --nexus-dir path_to_nexus_files --ploidy-file ploidy.json(optional) --copy-number copy_numbers.tsv(optional) path_to/taxa_table.txt
"""

import sys
import argparse
import os
import json
import random
from collections import defaultdict

def parse_nexus_taxa(nexus_file):
    """Parse NEXUS file to extract taxa names from either TAXA block or DATA matrix"""
    taxa = []
    in_taxa_block = False
    in_data_block = False
    in_matrix = False
    
    try:
        with open(nexus_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Look for TAXA block (traditional format)
                if line.upper().startswith('BEGIN TAXA'):
                    in_taxa_block = True
                    continue
                    
                # Look for DATA block (your format)
                if line.upper().startswith('BEGIN DATA'):
                    in_data_block = True
                    continue
                    
                # End of blocks
                if line.upper().startswith('END') and (in_taxa_block or in_data_block):
                    break
                
                # Parse TAXLABELS in TAXA block
                if in_taxa_block and line.upper().startswith('TAXLABELS'):
                    taxa_line = line[9:].strip()
                    if taxa_line.endswith(';'):
                        taxa_line = taxa_line[:-1]
                    if taxa_line:
                        taxa.extend(taxa_line.split())
                    continue
                        
                # Handle multi-line taxa list in TAXA block
                if in_taxa_block and line and not line.startswith('DIMENSIONS'):
                    line = line.replace(';', '')
                    if line:
                        taxa.extend(line.split())
                
                # Look for MATRIX in DATA block
                if in_data_block and line.upper().startswith('MATRIX'):
                    in_matrix = True
                    continue
                
                # Parse taxa names from matrix (your format)
                if in_matrix and line and not line.startswith(';'):
                    # Skip comment lines and empty lines
                    if line.startswith('[') or not line.strip():
                        continue
                    
                    # Extract taxon name (everything before first space/tab)
                    parts = line.split()
                    if parts:
                        taxon_name = parts[0]
                        # Remove any trailing punctuation
                        taxon_name = taxon_name.rstrip('.,;:')
                        if taxon_name:
                            taxa.append(taxon_name)
        
        return taxa
        
    except FileNotFoundError:
        print(f"Warning: Could not find NEXUS file: {nexus_file}")
        return []

def parse_copy_numbers(copy_file):
    """Parse copy number TSV file (optional)"""
    if not copy_file:
        return {}
    
    copy_counts = {}
    
    try:
        with open(copy_file, 'r') as f:
            header = f.readline().strip()
            
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        species_full = parts[0]
                        try:
                            copy_count = int(parts[1])
                            copy_counts[species_full] = copy_count
                        except ValueError:
                            print(f"Warning: Could not parse copy count for {species_full}")
                            continue
        
        print(f"Loaded copy numbers for {len(copy_counts)} taxa")
        return copy_counts
        
    except FileNotFoundError:
        print(f"Warning: Could not find copy numbers file: {copy_file}")
        return {}

def parse_ploidy_file(ploidy_file):
    """Parse ploidy JSON file"""
    if not ploidy_file:
        return {}
    
    try:
        with open(ploidy_file, 'r') as f:
            ploidy_dict = json.loads(f.read())
        
        print(f"Loaded known ploidy levels for {len(ploidy_dict)} species")
        return ploidy_dict
        
    except FileNotFoundError:
        print(f"Warning: Could not find ploidy file: {ploidy_file}")
        return {}
    except json.JSONDecodeError as e:
        print(f"Error: Could not parse ploidy JSON file: {e}")
        return {}

def extract_species_info(taxa_name):
    """
    Extract species name from taxa name - everything before the first underscore
    """
    # Split by underscore and take the first part
    parts = taxa_name.split('_')
    species = parts[1]
    return species

def detect_naming_pattern(taxon_names, separator='_'):
    """Auto-detect naming pattern from taxon names"""
    separators = [separator, '_', '-', '.', ' ']
    
    for sep in separators:
        for name in taxon_names[:10]:
            parts = name.split(sep)
            if len(parts) >= 2:
                return sep, 0, 1, 2  # genus_pos, species_pos, accession_pos
    
    return separator, None, None, None

def find_nexus_files(directory):
    """Find all NEXUS files in a directory"""
    nexus_files = []
    nexus_extensions = ['.nex', '.nexus', '.nxs']
    
    if not os.path.isdir(directory):
        print(f"Error: Directory not found: {directory}")
        return []
    
    for filename in os.listdir(directory):
        file_lower = filename.lower()
        if any(file_lower.endswith(ext) for ext in nexus_extensions):
            nexus_files.append(os.path.join(directory, filename))
    
    nexus_files.sort()  # Sort for consistent order
    return nexus_files

def get_all_taxa_from_nexus_files(nexus_files):
    """Get all unique taxa from multiple NEXUS files"""
    all_taxa = set()
    
    print("Reading taxa from NEXUS files:")
    for nexus_file in nexus_files:
        taxa = parse_nexus_taxa(nexus_file)
        all_taxa.update(taxa)
        print(f"  {os.path.basename(nexus_file)}: {len(taxa)} taxa")
    
    print(f"Total unique taxa across all files: {len(all_taxa)}")
    return all_taxa

def assign_homeolog_pairs_random_with_missing(sequences, species):
    """
    Randomly assign sequences to homeolog pairs (A,B)
    If odd number of sequences, create a missing sequence for the last individual
    Returns: list of (sequence_id, individual_id, genome) tuples
    """
    assignments = []
    
    # Shuffle sequences for random assignment
    shuffled_seqs = sequences.copy()
    random.shuffle(shuffled_seqs)
    
    # Create pairs
    individual_num = 1
    for i in range(0, len(shuffled_seqs), 2):
        individual_id = f"ind{individual_num}"
        
        # First sequence gets genome A
        assignments.append((shuffled_seqs[i], individual_id, 'A'))
        
        # Second sequence gets genome B (if it exists)
        if i + 1 < len(shuffled_seqs):
            assignments.append((shuffled_seqs[i + 1], individual_id, 'B'))
        else:
            # Odd number - create missing sequence with _miss suffix
            missing_seq_id = f"{species}_miss"
            assignments.append((missing_seq_id, individual_id, 'B'))
        
        individual_num += 1
    
    return assignments

def generate_taxa_table(args):
    """Main function to generate taxa table"""
    
    # Set random seed for reproducibility
    if args.seed:
        random.seed(args.seed)
        print(f"Random seed set to: {args.seed}")
    
    # Parse inputs
    copy_counts = parse_copy_numbers(args.copy_numbers)
    # Note: copy_counts is optional and not used in main logic
    
    known_ploidy = parse_ploidy_file(args.ploidy_file)
    
    # Get NEXUS files from directory or file list
    if args.nexus_directory:
        nexus_files = find_nexus_files(args.nexus_directory)
        if not nexus_files:
            print(f"Error: No NEXUS files found in directory: {args.nexus_directory}")
            return False
        print(f"Found {len(nexus_files)} NEXUS files in directory:")
        for nf in nexus_files:
            print(f"  {os.path.basename(nf)}")
    else:
        nexus_files = args.nexus_files
    
    # Get all taxa
    all_taxa = get_all_taxa_from_nexus_files(nexus_files)
    if not all_taxa:
        print("Error: No taxa found")
        return False
    
    print(f"\nUsing species extraction method: everything before first underscore")
    
    # Group taxa by species
    species_groups = defaultdict(list)
    
    for taxon in all_taxa:
        species = extract_species_info(taxon)
        species_groups[species].append(taxon)
    
    # Create taxa table
    taxa_entries = []
    
    print(f"\n{'='*80}")
    print("GENERATING TAXA TABLE - MODIFIED STRATEGY")
    print(f"{'='*80}")
    print(f"{'Species':<20} {'Seqs':<5} {'Known':<8} {'Strategy':<30} {'Status'}")
    print("-" * 80)
    
    for species, sequences in sorted(species_groups.items()):
        num_sequences = len(sequences)
        
        # Determine ploidy
        known_ploidy_str = str(known_ploidy.get(species, '-'))
        
        if species in known_ploidy:
            ploidy = known_ploidy[species]
        else:
            # Default: treat unknown as polyploid if multiple sequences
            ploidy = 2 if num_sequences == 1 else 4
        
        # Apply strategy
        if ploidy == 2:
            # DIPLOID: All sequences as different individuals with genome A
            strategy = f"All {num_sequences} seqs as diploid indivs"
            status = "INCLUDED"
            
            for i, seq in enumerate(sequences):
                taxa_entries.append({
                    'ID': seq,
                    'species': species,
                    'individual': f"{species}_ind{i+1}",
                    'genome': 'A'
                })
                
        elif ploidy == 4:
            # POLYPLOID: Random homeolog pairs (A,B) with missing sequences if odd
            if num_sequences % 2 == 1:
                strategy = f"Random pairs from {num_sequences} seqs + missing"
                status = "INCLUDED (with missing)"
            else:
                strategy = f"Random pairs from {num_sequences} seqs"
                status = "INCLUDED"
            
            assignments = assign_homeolog_pairs_random_with_missing(sequences, species)
            
            for seq_id, individual_id, genome in assignments:
                taxa_entries.append({
                    'ID': seq_id,
                    'species': species,
                    'individual': f"{species}_{individual_id}",
                    'genome': genome
                })
        
        else:
            # Unsupported ploidy
            strategy = "-"
            status = f"SKIPPED (ploidy {ploidy}x unsupported)"
        
        print(f"{species:<20} {num_sequences:<5} {known_ploidy_str:<8} {strategy:<30} {status}")
    
    # Write output
    print(f"\nWriting taxa table to {args.output}...")
    with open(args.output, 'w') as f:
        f.write("ID species individual genome\n")
        for entry in taxa_entries:
            f.write(f"{entry['ID']} {entry['species']} {entry['individual']} {entry['genome']}\n")
    
    # Summary statistics
    species_count = len(set(e['species'] for e in taxa_entries))
    diploid_species = set()
    polyploid_species = set()
    missing_sequences = sum(1 for e in taxa_entries if '_miss' in e['ID'])
    
    for entry in taxa_entries:
        species = entry['species']
        if entry['genome'] == 'A' and not any(e['species'] == species and e['genome'] == 'B' for e in taxa_entries):
            diploid_species.add(species)
        elif any(e['species'] == species and e['genome'] == 'B' for e in taxa_entries):
            polyploid_species.add(species)
    
    print(f"\nSUCCESS: Created taxa table with {len(taxa_entries)} entries")
    print(f"  Species included: {species_count}")
    print(f"  Diploid species: {len(diploid_species)}")
    print(f"  Polyploid species: {len(polyploid_species)}")
    print(f"  Missing sequences created: {missing_sequences}")
    
    if args.seed:
        print(f"  Random seed used: {args.seed} (use same seed for reproducible results)")
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description="Generate AlloppNET taxa table with modified strategy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Strategy:
  Diploids (2x):    All copies as different individuals with genome A
  Polyploids (4x):  Random homeolog pairs with genomes A and B
                   If odd number of copies, create missing sequence with _miss suffix
  Unknown ploidy:   1 copy = diploid, 2+ copies = polyploid

Species name extraction: Everything before the first underscore

Examples:
  # Using directory with all NEXUS files (minimal)
  python taxa_table.py --nexus-dir /path/to/nexus/files output.txt
  
  # With optional copy numbers file
  python taxa_table.py --nexus-dir /path/to/nexus/files --copy-numbers copy_numbers.tsv output.txt
  
  # With known ploidy file and random seed
  python taxa_table.py --nexus-dir nexus_files --ploidy-file ploidy.json --seed 42 output.txt
  
  # Using individual files
  python taxa_table.py --nexus-files gene1.nex gene2.nex output.txt

File formats:
  Copy numbers (TSV, optional):
    Species	RepresentativeCopyNumber	Distribution
    Species1_acc1	1	1:1
    
  Ploidy file (JSON, optional):
    {"species1": 2, "species2": 4}
        """
    )
    
    parser.add_argument('output', help='Output taxa table file')
    
    # Either nexus directory OR nexus files (mutually exclusive)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--nexus-dir', dest='nexus_directory', help='Directory containing NEXUS files (.nex, .nexus, .nxs)')
    group.add_argument('--nexus-files', nargs='+', help='Individual NEXUS alignment files')
    
    parser.add_argument('--copy-numbers', dest='copy_numbers', help='Copy numbers TSV file (optional)')
    parser.add_argument('--ploidy-file', help='JSON file with known ploidy levels (optional)')
    parser.add_argument('--seed', type=int, help='Random seed for reproducible homeolog assignment')
    
    args = parser.parse_args()
    
    # Validate inputs - copy_numbers is now optional
    if args.copy_numbers and not os.path.exists(args.copy_numbers):
        print(f"Error: Copy numbers file not found: {args.copy_numbers}")
        return 1
    
    # Check NEXUS files
    if args.nexus_directory:
        if not os.path.isdir(args.nexus_directory):
            print(f"Error: Directory not found: {args.nexus_directory}")
            return 1
        nexus_files = find_nexus_files(args.nexus_directory)
        if not nexus_files:
            print(f"Error: No NEXUS files found in directory: {args.nexus_directory}")
            return 1
    else:
        nexus_files = args.nexus_files or []
        for nf in nexus_files:
            if not os.path.exists(nf):
                print(f"Error: NEXUS file not found: {nf}")
                return 1
    
    print("Modified AlloppNET Taxa Table Generator")
    print("=" * 40)
    if args.copy_numbers:
        print(f"Copy numbers: {args.copy_numbers}")
    else:
        print("Copy numbers: Not provided (optional)")
    print(f"Output: {args.output}")
    
    if args.nexus_directory:
        print(f"NEXUS directory: {args.nexus_directory}")
        print(f"Found NEXUS files: {len(find_nexus_files(args.nexus_directory))}")
    else:
        print(f"NEXUS files: {len(nexus_files)}")
    
    if args.ploidy_file:
        print(f"Known ploidy: {args.ploidy_file}")
    
    success = generate_taxa_table(args)
    
    if success:
        print(f"\nNext steps:")
        print(f"1. Review the taxa table: {args.output}")
        print(f"2. Note any missing sequences (marked with _miss suffix)")
        print(f"3. Use with AlloppNET R script and your NEXUS files")
        print(f"4. Check if AlloppNET converges with this assignment")
        return 0
    else:
        return 1

if __name__ == "__main__":
    sys.exit(main())