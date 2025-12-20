import os
import glob
import pandas as pd
from collections import defaultdict, Counter
# Ensure reticulate_tree.py is in the same directory
from reticulate_tree import ReticulateTree

# Path to your files
DIR_PATH = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks/"

def count_autopolyploid_events(tree):
    """
    Count autopolyploidization events in a MUL-tree.

    An autopolyploidization EVENT is a single WGD that duplicates an entire clade.
    For example, in Ding_2023, there's ONE event that duplicated the Rch clade,
    which resulted in 7 autopolyploid species (Rch, Pla, Jmi, CiP, Jre, Jma, Aro).

    Strategy: Find top-level duplicated clades (identical siblings), excluding
    duplications that are nested within higher-level duplications.

    Args:
        tree: ete3 Tree object

    Returns:
        Number of autopolyploidization events
    """
    # Build canonical form for each node
    canonical_map = {}

    def compute_canonical(node):
        '''Recursively compute canonical string representation of subtree'''
        if node.is_leaf():
            canon = node.name
        else:
            child_encodings = sorted([compute_canonical(c) for c in node.children])
            canon = '(' + ','.join(child_encodings) + ')'
        canonical_map[node] = canon
        return canon

    compute_canonical(tree)

    # Track nodes that are part of already-counted duplication events
    counted_descendants = set()

    def is_descendant_of_counted(node):
        """Check if this node is a descendant of an already-counted duplication"""
        current = node
        while current:
            if id(current) in counted_descendants:
                return True
            current = current.up
        return False

    def mark_descendants(node):
        """Mark all descendants of this node as part of a counted duplication"""
        for desc in node.traverse():
            counted_descendants.add(id(desc))

    # Count autopolyploidization events
    # Process nodes top-down (root to leaves) to catch highest-level duplications first
    num_auto_events = 0

    for node in tree.traverse('preorder'):  # Top-down traversal
        # Skip if this node is already part of a counted duplication
        if id(node) in counted_descendants:
            continue

        # Skip root (no siblings)
        if not node.up:
            continue

        # Get siblings (other children of the same parent)
        siblings = [child for child in node.up.children]

        # Group siblings by canonical form
        canon_groups = defaultdict(list)
        for sibling in siblings:
            if id(sibling) not in counted_descendants:  # Only consider uncounted siblings
                canon = canonical_map[sibling]
                canon_groups[canon].append(sibling)

        # Find the group this node belongs to
        node_canon = canonical_map[node]
        group = canon_groups[node_canon]

        # If there are multiple identical siblings, this is an autopolyploidization event
        if len(group) > 1:
            # Count this as ONE event (regardless of how many copies)
            # Only count once per group (when we encounter the first member)
            if group[0] == node:  # First occurrence of this group
                num_auto_events += 1

                # Mark all members of this group and their descendants as counted
                for member in group:
                    mark_descendants(member)

    return num_auto_events

def analyze_mul_trees(directory, relaxed_threshold=0.2):
    files = glob.glob(os.path.join(directory, "*"))
    results = []

    print(f"Found {len(files)} files. Analyzing...")

    for filepath in files:
        # Filter for valid tree files only
        if os.path.isdir(filepath) or not filepath.endswith(('.tre', '.nwk', '.tree')):
            continue

        filename = os.path.basename(filepath)

        try:
            with open(filepath, 'r') as f:
                content = f.read().strip()

            if not content or not content.endswith(';'):
                continue

            # --- 1. GET BASIC STATS (From Strict Instance) ---
            rt_strict = ReticulateTree(content, is_multree=True)
            stats_strict = rt_strict.measure(printout=False)

            # --- 2. CALCULATE LEAF STATS ---
            leaf_counts = stats_strict['leaf_counts'] # Dict like {'SpA': 2, 'SpB': 1}

            num_species = len(leaf_counts)

            # Count how many species have > 1 copy
            polyploid_species = [name for name, count in leaf_counts.items() if count > 1]
            num_polyploids = len(polyploid_species)

            # Find the maximum number of copies any single species has
            max_copies = max(leaf_counts.values()) if leaf_counts else 0

            # --- 3. COUNT AUTOPOLYPLOIDIZATION EVENTS (in MUL-tree before folding) ---
            # This counts the number of WGD events, not the number of autopolyploid species
            # E.g., Ding_2023 has 1 event (Rch clade duplication) → 7 autopolyploid species
            num_auto_events = count_autopolyploid_events(rt_strict.tree)

            # --- 4. GET RELAXED STATS ---
            rt_relaxed = ReticulateTree(content, is_multree=True, threshold=relaxed_threshold, normalize=True)
            stats_relaxed = rt_relaxed.measure(printout=False)

            # --- 5. CALCULATE POLYPLOIDIZATION EVENTS ---
            # H_Strict = reticulations in network = allopolyploid/hybridization events
            # Num_Auto_Events = autopolyploidization events (WGD creating identical clades)
            h_strict = stats_strict['reticulation_count']
            h_relaxed = stats_relaxed['reticulation_count']

            # Total WGD = autopolyploidization + allopolyploidization (reticulations)
            total_wgd = num_auto_events + h_strict

            # --- 6. COMPILE ROW ---
            row = {
                'Filename': filename,
                'Num_Species': num_species,
                'Num_Polyploids': num_polyploids,
                'Max_Copies': max_copies,
                'H_Strict': h_strict,
                'H_Relaxed': h_relaxed,
                'H_Diff': h_relaxed - h_strict,
                'Num_Autopolyploidization_Events': num_auto_events,
                'Total_WGD': total_wgd,
                'Polyploid_Names': ", ".join(sorted(polyploid_species))
            }
            results.append(row)

        except Exception as e:
            print(f"Error on {filename}: {e}")
            import traceback
            traceback.print_exc()

    return pd.DataFrame(results)

if __name__ == "__main__":
    df = analyze_mul_trees(DIR_PATH)

    if not df.empty:
        # Save to CSV
        out_file = os.path.join(DIR_PATH, "mul_tree_final_stats.csv")
        df.to_csv(out_file, index=False)

        print("\n" + "="*100)
        print("Summary Statistics - Network Characteristics")
        print("="*100)

        # Print summary with new columns
        cols = ['Filename', 'Num_Species', 'Num_Polyploids', 'Max_Copies',
                'H_Strict', 'Num_Autopolyploidization_Events', 'Total_WGD']
        print(df[cols].to_string(index=False))

        print("\n" + "="*100)
        print("Column Definitions:")
        print("  - Num_Species: Total number of unique species")
        print("  - Num_Polyploids: Number of species with >1 copy (result of WGD)")
        print("  - Max_Copies: Maximum number of copies of any species")
        print("  - H_Strict: Number of reticulations (allopolyploid/hybridization events)")
        print("  - Num_Autopolyploidization_Events: Number of WGD events creating identical clades")
        print("  - Total_WGD: Total polyploidization events (Autopolyploidization + H_Strict)")
        print("\nExample (Ding_2023.tre):")
        print("  - 1 autopolyploidization event (Rch clade duplication)")
        print("  - Creates 7 autopolyploid species (Rch, Pla, Jmi, CiP, Jre, Jma, Aro)")
        print("  - 0 reticulations (no hybridization)")
        print("  - Total_WGD = 1")
        print("="*100)

        print(f"\nStats saved to: {out_file}")
    else:
        print("No valid data found.")