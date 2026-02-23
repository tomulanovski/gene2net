#!/usr/bin/env python3
"""
AlloppNET Taxa Table Generator

Creates taxa tables following the Jones manual (Section 2.4) strategy:
- Diploids: All copies as different individuals with genome A
- Polyploids with exactly 2 sequences: Pair as A/B homeologs in one individual
- Polyploids with 1 sequence: Individual with genome A + _miss placeholder for genome B
- Polyploids with 3+ sequences: Each sequence gets its own individual with genome A
  and a _miss placeholder for genome B (AlloppNET assigns Parent1/Parent2 independently)

Usage:
    python scripts/create_taxa_table.py --nexus-dir path_to_nexus_files --ploidy-file ploidy.json(optional) --copy-number copy_numbers.tsv(optional) path_to/taxa_table.txt
"""

import sys
import argparse
import os
import json
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

def extract_species_info(taxa_name, species_field=0, known_species=None):
    """
    Extract species name from taxa name.

    If known_species is provided (list of known species names), uses longest-prefix
    matching to handle tricky cases like 'saxatilis' vs 'saxatilis_var_mairei'.

    Otherwise falls back to underscore-separated field extraction.

    Args:
        taxa_name: Full taxon name (e.g., "Ephedra_sinica_KT033298")
        species_field: Which underscore-separated field to use as species name (default: 0)
        known_species: Optional list of known species names for prefix matching
    """
    # If we have known species, try longest-prefix match
    if known_species:
        for species in sorted(known_species, key=len, reverse=True):
            if taxa_name.startswith(species + '_') or taxa_name == species:
                return species

    parts = taxa_name.split('_')
    if species_field < len(parts):
        return parts[species_field]
    return parts[0]

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

def assign_polyploid_sequences(sequences, species):
    """
    Assign polyploid sequences following Jones manual Section 2.4 strategy.

    - 1 sequence: single individual with genome A + _miss placeholder for genome B
    - 2 sequences: pair as A/B homeologs in one individual
    - 3+ sequences: each sequence gets its own individual with genome A + _miss for B
      (AlloppNET assigns each to Parent1 or Parent2 independently via MCMC)

    Returns: list of (sequence_id, individual_id, genome) tuples
    """
    assignments = []
    num_seqs = len(sequences)

    if num_seqs == 2:
        # Exactly 2 sequences: pair as homeologs A/B in one individual
        assignments.append((sequences[0], 'ind1', 'A'))
        assignments.append((sequences[1], 'ind1', 'B'))

    elif num_seqs == 1:
        # Single sequence: genome A + missing B
        assignments.append((sequences[0], 'ind1', 'A'))
        assignments.append((f"{species}_miss", 'ind1', 'B'))

    else:
        # 3+ sequences: each gets its own individual with A + missing B
        for i, seq in enumerate(sorted(sequences)):
            ind_id = f"ind{i + 1}"
            assignments.append((seq, ind_id, 'A'))
            assignments.append((f"{species}_miss{i + 1}", ind_id, 'B'))

    return assignments

def generate_taxa_table(args):
    """Main function to generate taxa table"""
    
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
    
    # Use known species names from ploidy file for prefix matching (handles tricky names)
    known_species = list(known_ploidy.keys()) if known_ploidy else None
    if known_species:
        print(f"\nUsing known species from ploidy file for name matching ({len(known_species)} species)")
    else:
        print(f"\nUsing species extraction: underscore-separated field {args.species_field}")

    # Group taxa by species
    species_groups = defaultdict(list)

    for taxon in all_taxa:
        species = extract_species_info(taxon, args.species_field, known_species)
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
            # POLYPLOID: Strategy depends on number of sequences
            if num_sequences == 1:
                strategy = "1 seq: A + missing B"
                status = "INCLUDED (with missing)"
            elif num_sequences == 2:
                strategy = "2 seqs: paired as A/B"
                status = "INCLUDED"
            else:
                strategy = f"{num_sequences} seqs: each ind A + miss B"
                status = "INCLUDED (with missing)"

            assignments = assign_polyploid_sequences(sequences, species)

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
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description="Generate AlloppNET taxa table (Jones manual Section 2.4 strategy)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Strategy (per Jones manual Section 2.4):
  Diploids (2x):    All copies as different individuals with genome A
  Polyploids (4x):
    1 sequence:     Individual with genome A + _miss placeholder for genome B
    2 sequences:    Paired as A/B homeologs in one individual
    3+ sequences:   Each sequence gets its own individual (A + _miss B),
                    AlloppNET assigns each to Parent1/Parent2 via MCMC
  Unknown ploidy:   1 copy = diploid, 2+ copies = polyploid

Species name extraction:
  Uses --species-field to select which underscore-separated field is the species name.
  Default: 0 (first field, i.e., everything before the first underscore)
  Example: --species-field 1 for names like "Genus_species_accession"

Examples:
  # Using directory with all NEXUS files (minimal)
  python create_taxa_table.py --nexus-dir /path/to/nexus/files output.txt

  # With known ploidy file
  python create_taxa_table.py --nexus-dir nexus_files --ploidy-file ploidy.json output.txt

  # Species name is the second field (e.g., Genus_species_acc)
  python create_taxa_table.py --nexus-dir nexus_files --species-field 1 output.txt

  # Using individual files
  python create_taxa_table.py --nexus-files gene1.nex gene2.nex output.txt

File formats:
  Copy numbers (TSV, optional):
    Species\tRepresentativeCopyNumber\tDistribution
    Species1_acc1\t1\t1:1

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
    parser.add_argument('--species-field', dest='species_field', type=int, default=0,
                        help='Which underscore-separated field is the species name (default: 0, first field)')
    
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