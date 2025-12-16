#!/usr/bin/env python3
"""
prepare_alloppnet_input.py - Prepare AlloppNET input from simulation alignments

This script:
1. Reads PHY format sequence alignments from simulations
2. Parses taxon names (SimPhy format: species_locusid_individualid)
3. Analyzes copy number distributions per taxon across all alignments
4. Generates ploidy_level.json using kernel smoothing (robust to outliers)
5. Converts PHY alignments to NEXUS format
6. Creates taxa_table.txt for AlloppNET

Handles:
- Gene duplication (>2 copies per taxon)
- Gene loss (<2 copies per taxon)
- Complex taxon names with underscores

Usage:
    python prepare_alloppnet_input.py \\
        --network Ding_2023 \\
        --config ils_low_10M \\
        --replicate 1 \\
        --alignment-dir /path/to/alignments/ \\
        --output-dir /path/to/output/
"""

import argparse
import os
import sys
import json
import glob
from collections import defaultdict, Counter
from Bio import AlignIO, Phylo


def extract_taxon_name(seq_id):
    """
    Extract species name from sequence ID (WITHOUT network copy suffix).

    SimPhy format: species_locus_individual (3 parts)
    Modified format: species_copy_locus_individual (4+ parts)

    The network copy suffix (_0, _1, etc.) was added to simulate polyploidy
    since SimPhy can't have the same species appear multiple times.

    Examples:
        Galeopsisladanum_2_0 → Galeopsisladanum (diploid, 3 parts)
        Lamiumalbum_0_3_1 → Lamiumalbum (polyploid copy 0, 4 parts)
        Lamiumalbum_1_2_0 → Lamiumalbum (polyploid copy 1, 4 parts)

    Args:
        seq_id (str): Sequence ID

    Returns:
        str: Species name (without copy suffix)
    """
    parts = seq_id.split('_')

    # First part is always the species name
    return parts[0]


def extract_gene_number(filename):
    """
    Extract gene number from alignment filename.

    Example: alignment_0001.phy → 0001

    Args:
        filename (str): Alignment filename

    Returns:
        str: Gene number (4-digit string)
    """
    import re
    match = re.search(r'alignment_(\d+)\.phy', filename)
    if match:
        return match.group(1)
    else:
        # Fallback: return filename
        return os.path.basename(filename).replace('.phy', '')


def extract_copy_number(seq_id):
    """
    Extract network copy number from sequence ID.

    This extracts the user-added copy suffix (_0, _1, etc.) that was added
    to species names to simulate polyploidy in SimPhy.

    Format detection:
    - 3 parts (species_locus_individual): Diploid, return '0'
    - 4+ parts (species_copy_locus_individual): Polyploid, return copy number

    Examples:
        Galeopsisladanum_2_0 → '0' (diploid, 3 parts, no copy suffix)
        Lamiumalbum_0_3_1 → '0' (polyploid, 4 parts, network copy 0)
        Lamiumalbum_1_2_0 → '1' (polyploid, 4 parts, network copy 1)
        Lamiumalbum_2_5_1 → '2' (polyploid, 4 parts, network copy 2)

    Args:
        seq_id (str): Sequence ID

    Returns:
        str: Network copy number (the _i suffix added for polyploidy)
    """
    parts = seq_id.split('_')

    if len(parts) == 3:
        # Format: species_locus_individual (diploid, no copy suffix)
        return '0'
    elif len(parts) >= 4:
        # Format: species_copy_locus_individual (polyploid)
        return parts[1]  # Return the network copy number
    else:
        # Unexpected format
        return '0'


def get_representative_copy_number(copy_distribution, kernel_width=2):
    """
    Determine the most representative copy number using kernel smoothing.

    This approach is robust to outliers caused by occasional gene duplication or loss.
    It finds the peak of the distribution using a triangular kernel.

    Args:
        copy_distribution (Counter): Distribution of copy numbers across genes
        kernel_width (int): How far the kernel extends in each direction (default: 2)

    Returns:
        int: Representative copy number
    """
    # Get sorted list of copy numbers and their frequencies
    sorted_counts = sorted(copy_distribution.items())

    # If there's only one copy number, use it
    if len(sorted_counts) == 1:
        copy_number, _ = sorted_counts[0]
        return copy_number

    # If there are just two copy numbers, use the more frequent one
    if len(sorted_counts) == 2:
        if sorted_counts[0][1] >= sorted_counts[1][1]:
            return sorted_counts[0][0]
        else:
            return sorted_counts[1][0]

    # For 3 or more distinct copy numbers, use kernel smoothing
    copy_numbers = [num for num, _ in sorted_counts]
    frequencies = [freq for _, freq in sorted_counts]

    # Calculate smoothed distribution - use a triangular kernel
    smoothed_values = {}

    for i, center in enumerate(copy_numbers):
        # Apply kernel centered at each observed copy number
        for j, copy_num in enumerate(copy_numbers):
            # Calculate kernel weight based on distance
            distance = abs(center - copy_num)
            if distance <= kernel_width:
                weight = 1 - (distance / (kernel_width + 1))

                # Add weighted frequency to the smoothed value
                if center not in smoothed_values:
                    smoothed_values[center] = 0
                smoothed_values[center] += frequencies[j] * weight

    # Find the copy number with the highest smoothed value
    max_smoothed = 0
    representative = copy_numbers[0]

    for copy_num, smoothed_val in smoothed_values.items():
        if smoothed_val > max_smoothed:
            max_smoothed = smoothed_val
            representative = copy_num

    return representative


def read_phylip_manual(phy_file):
    """
    Manually parse PHYLIP file (sequential format).

    Returns:
        dict: {seq_id: sequence_string}
    """
    sequences = {}

    with open(phy_file, 'r') as f:
        lines = f.readlines()

    if not lines:
        return sequences

    # Skip header line (ntax nchar)
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue

        # Split on whitespace - first part is taxon name, rest is sequence
        parts = line.split(None, 1)
        if len(parts) == 2:
            taxon_name, sequence = parts
            sequences[taxon_name] = sequence.replace(' ', '')  # Remove any spaces in sequence

    return sequences


def count_copies_from_alignments(alignment_dir, kernel_width=2, verbose=False):
    """
    Count copies per taxon across all alignment files and determine representative copy numbers.

    Uses kernel smoothing to find the most representative copy number for each taxon,
    which is robust to outliers from occasional gene duplication or loss.

    Args:
        alignment_dir (str): Directory containing alignment_*.phy files
        kernel_width (int): Width of smoothing kernel (default: 2)
        verbose (bool): Print progress messages

    Returns:
        tuple: (taxon_representative_copies, taxon_copy_distributions)
            - taxon_representative_copies: dict mapping taxon -> representative copy count
            - taxon_copy_distributions: dict mapping taxon -> Counter of copy numbers
    """
    if verbose:
        print(f"\nCounting copies from alignments in: {alignment_dir}")

    # Get all PHY files
    phy_files = sorted(glob.glob(os.path.join(alignment_dir, "alignment_*.phy")))

    if not phy_files:
        print(f"ERROR: No alignment files found in {alignment_dir}")
        sys.exit(1)

    if verbose:
        print(f"Found {len(phy_files)} alignment files")

    # Track copies: taxon -> gene -> set of copy numbers
    taxon_gene_copies = defaultdict(lambda: defaultdict(set))

    # Process each alignment file
    for i, phy_file in enumerate(phy_files, 1):
        if verbose and i % 100 == 0:
            print(f"  Processing alignment {i}/{len(phy_files)}...")

        try:
            # Manually parse PHYLIP file
            sequences = read_phylip_manual(phy_file)
        except Exception as e:
            print(f"WARNING: Could not read {phy_file}: {e}")
            continue

        if not sequences:
            continue

        gene_num = extract_gene_number(os.path.basename(phy_file))

        # Count copies in this alignment
        for seq_id in sequences.keys():
            taxon = extract_taxon_name(seq_id)
            copy = extract_copy_number(seq_id)
            taxon_gene_copies[taxon][gene_num].add(copy)

    # Calculate copy number distributions for each taxon
    taxon_copy_distributions = {}
    taxon_representative_copies = {}

    for taxon, gene_dict in taxon_gene_copies.items():
        # Count how many genes have each copy number
        copy_distribution = Counter()
        for gene_num, copy_set in gene_dict.items():
            num_copies = len(copy_set)
            copy_distribution[num_copies] += 1

        taxon_copy_distributions[taxon] = copy_distribution

        # Find representative copy number using kernel smoothing
        representative = get_representative_copy_number(copy_distribution, kernel_width)
        taxon_representative_copies[taxon] = representative

    if verbose:
        print(f"\nFound {len(taxon_representative_copies)} unique taxa:")
        for taxon in sorted(taxon_representative_copies.keys()):
            rep_copies = taxon_representative_copies[taxon]
            dist = taxon_copy_distributions[taxon]
            dist_str = ", ".join([f"{count}:{freq}" for count, freq in dist.most_common()])
            print(f"  {taxon}: {rep_copies} copies (distribution: {dist_str})")

    return taxon_representative_copies, taxon_copy_distributions


def generate_ploidy_json(taxon_representative_copies, output_file, verbose=False):
    """
    Generate ploidy_level.json based on representative copy counts (from kernel smoothing).

    1 copy → diploid (ploidy=2)
    2+ copies → tetraploid (ploidy=4)

    Args:
        taxon_representative_copies (dict): taxon -> representative_copies mapping
        output_file (str): Output JSON file path
        verbose (bool): Print progress messages
    """
    if verbose:
        print(f"\nGenerating ploidy_level.json...")

    ploidy = {}

    for taxon, rep_copies in taxon_representative_copies.items():
        if rep_copies == 1:
            ploidy[taxon] = 2  # Diploid
        else:  # rep_copies >= 2
            ploidy[taxon] = 4  # Tetraploid

    with open(output_file, 'w') as f:
        json.dump(ploidy, f, indent=2)

    if verbose:
        diploid_count = sum(1 for p in ploidy.values() if p == 2)
        tetraploid_count = sum(1 for p in ploidy.values() if p == 4)
        print(f"  Diploid (ploidy=2): {diploid_count} taxa")
        print(f"  Tetraploid (ploidy=4): {tetraploid_count} taxa")
        print(f"  Saved to: {output_file}")


def convert_phy_to_nexus(alignment_dir, output_dir, taxon_representative_copies, verbose=False):
    """
    Convert PHY alignments to NEXUS format for AlloppNET.

    Strips locusid from sequence IDs (species_locusid_individualid → species_individualid).

    Args:
        alignment_dir (str): Input directory with *.phy files
        output_dir (str): Output directory for *.nex files
        taxon_representative_copies (dict): taxon -> representative_copies mapping
        verbose (bool): Print progress messages
    """
    if verbose:
        print(f"\nConverting PHY to NEXUS...")

    os.makedirs(output_dir, exist_ok=True)

    phy_files = sorted(glob.glob(os.path.join(alignment_dir, "alignment_*.phy")))

    for i, phy_file in enumerate(phy_files, 1):
        if verbose and i % 100 == 0:
            print(f"  Converting {i}/{len(phy_files)}...")

        try:
            # Manually parse PHYLIP file
            sequences = read_phylip_manual(phy_file)

            if not sequences:
                print(f"WARNING: No sequences found in {phy_file}")
                continue

            # Get sequence length
            seq_length = len(next(iter(sequences.values())))
            ntax = len(sequences)

            # Write NEXUS file manually
            nex_file = os.path.join(output_dir, os.path.basename(phy_file).replace('.phy', '.nex'))

            with open(nex_file, 'w') as f:
                f.write("#NEXUS\n")
                f.write("BEGIN DATA;\n")
                f.write(f"  DIMENSIONS NTAX={ntax} NCHAR={seq_length};\n")
                f.write("  FORMAT DATATYPE=DNA MISSING=N GAP=-;\n")
                f.write("  MATRIX\n")

                # Write sequences with modified names
                for seq_id, sequence in sequences.items():
                    taxon = extract_taxon_name(seq_id)
                    copy = extract_copy_number(seq_id)
                    new_id = f"{taxon}_{copy}"
                    f.write(f"    {new_id}  {sequence}\n")

                f.write("  ;\n")
                f.write("END;\n")

        except Exception as e:
            print(f"WARNING: Could not convert {phy_file}: {e}")
            continue

    if verbose:
        print(f"  Converted {len(phy_files)} alignments to NEXUS")
        print(f"  Saved to: {output_dir}")


def generate_taxa_table(taxon_representative_copies, ploidy, output_file, verbose=False):
    """
    Generate taxa_table.txt for AlloppNET.

    Format: ID species individual genome

    - Diploid taxa (1 copy): Single individual with genome A
    - Tetraploid taxa (2+ copies): Pair homeologs (A/B)
      - _0 with _1 (individual 1, genomes A/B)
      - _2 with _3 (individual 2, genomes A/B)
      - If odd number, pair last with missing data

    Args:
        taxon_representative_copies (dict): taxon -> representative_copies mapping
        ploidy (dict): taxon -> ploidy_level mapping
        output_file (str): Output taxa_table.txt path
        verbose (bool): Print progress messages
    """
    if verbose:
        print(f"\nGenerating taxa_table.txt...")

    lines = ["ID species individual genome"]

    for taxon in sorted(taxon_representative_copies.keys()):
        rep_copies = taxon_representative_copies[taxon]
        ploidy_level = ploidy[taxon]

        if ploidy_level == 2:  # Diploid (1 copy)
            # Single sequence, genome A
            lines.append(f"{taxon}_0 {taxon} {taxon}_ind1 A")

        else:  # Tetraploid (2+ copies)
            # Pair homeologs: _0 with _1, _2 with _3, etc.
            for i in range(0, rep_copies, 2):
                ind_num = (i // 2) + 1

                if i + 1 < rep_copies:
                    # Pair exists
                    lines.append(f"{taxon}_{i} {taxon} {taxon}_ind{ind_num} A")
                    lines.append(f"{taxon}_{i+1} {taxon} {taxon}_ind{ind_num} B")
                else:
                    # Odd one out - pair with missing data
                    lines.append(f"{taxon}_{i} {taxon} {taxon}_ind{ind_num} A")
                    lines.append(f"{taxon}_miss {taxon} {taxon}_ind{ind_num} B")

    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))

    if verbose:
        print(f"  Generated taxa table with {len(lines)-1} entries")
        print(f"  Saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Prepare AlloppNET input from simulation alignments',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--network', required=True,
                       help='Network name (e.g., Ding_2023)')
    parser.add_argument('--config', required=True,
                       help='Configuration name (e.g., ils_low_10M)')
    parser.add_argument('--replicate', required=True, type=int,
                       help='Replicate number (1-5)')
    parser.add_argument('--alignment-dir', required=True,
                       help='Directory containing alignment_*.phy files')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for AlloppNET input files')
    parser.add_argument('--kernel-width', type=int, default=2,
                       help='Width of smoothing kernel for ploidy detection (default: 2)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Print progress messages')

    args = parser.parse_args()

    # Verify alignment directory exists
    if not os.path.exists(args.alignment_dir):
        print(f"ERROR: Alignment directory not found: {args.alignment_dir}")
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("="*80)
    print("AlloppNET Input Preparation (with kernel smoothing)")
    print("="*80)
    print(f"Network:       {args.network}")
    print(f"Config:        {args.config}")
    print(f"Replicate:     {args.replicate}")
    print(f"Kernel width:  {args.kernel_width}")
    print(f"Input dir:     {args.alignment_dir}")
    print(f"Output dir:    {args.output_dir}")
    print("="*80)

    # Step 1: Count copies from alignments and determine representative copy numbers
    print("\n[Step 1/4] Analyzing copy number distributions (using kernel smoothing)...")
    taxon_representative_copies, taxon_copy_distributions = count_copies_from_alignments(
        args.alignment_dir, args.kernel_width, args.verbose
    )

    if not taxon_representative_copies:
        print("ERROR: No taxa found in alignments")
        sys.exit(1)

    # Step 2: Generate ploidy_level.json based on representative copy numbers
    print("\n[Step 2/4] Generating ploidy_level.json...")
    ploidy_file = os.path.join(args.output_dir, "ploidy_level.json")

    ploidy = {}
    for taxon, rep_copies in taxon_representative_copies.items():
        ploidy[taxon] = 2 if rep_copies == 1 else 4

    generate_ploidy_json(taxon_representative_copies, ploidy_file, args.verbose)

    # Step 3: Convert PHY to NEXUS
    print("\n[Step 3/4] Converting PHY to NEXUS...")
    convert_phy_to_nexus(args.alignment_dir, args.output_dir, taxon_representative_copies, args.verbose)

    # Step 4: Generate taxa_table.txt
    print("\n[Step 4/4] Generating taxa_table.txt...")
    taxa_table_file = os.path.join(args.output_dir, "taxa_table.txt")
    generate_taxa_table(taxon_representative_copies, ploidy, taxa_table_file, args.verbose)

    print("\n" + "="*80)
    print("AlloppNET input preparation complete!")
    print("="*80)
    print(f"\nOutput files:")
    print(f"  - {len(glob.glob(os.path.join(args.output_dir, '*.nex')))} NEXUS alignments")
    print(f"  - ploidy_level.json")
    print(f"  - taxa_table.txt")
    print(f"\nOutput directory: {args.output_dir}")
    print("="*80)


if __name__ == '__main__':
    main()
