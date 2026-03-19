#!/usr/bin/env python3
"""
Fix naming consistency across all network.tre output files for each dataset.

For each dataset:
1. Reads all network.tre files across methods (grampa, polyphest, padre, grandma_split, alloppnet)
2. Identifies the "canonical" species names from non-AlloppNET methods
3. Fixes AlloppNET-specific issues:
   - Removes trailing digit suffixes (0, 1, 2...) from polyploid species
   - Removes underscores to match other methods' naming conventions
4. Writes fixed files (backs up originals as network.tre.bak)

Usage:
    python fix_network_naming.py /path/to/papers [--dry-run] [--dataset DATASET]

Examples:
    # Analyze all datasets (dry run)
    python fix_network_naming.py /groups/itay_mayrose/tomulanovski/gene2net/papers --dry-run

    # Fix all datasets
    python fix_network_naming.py /groups/itay_mayrose/tomulanovski/gene2net/papers

    # Fix only one dataset
    python fix_network_naming.py /groups/itay_mayrose/tomulanovski/gene2net/papers --dataset Wu_2015
"""

import os
import re
import sys
import shutil
import argparse
from collections import defaultdict


def is_alloppnet_method(method_name):
    """Check if a method name refers to AlloppNET (handles spelling variants)."""
    return 'allop' in method_name.lower()


def extract_leaf_names(tree_str):
    """Extract all leaf names from a Newick tree string.

    Handles both standard Newick and extended Newick with #H notation.
    Returns list of leaf names (may contain duplicates for MUL-trees).
    """
    # Remove branch lengths, support values, and #H annotations for parsing
    # We want to extract just the taxon names
    leaves = []

    # Match taxon names: sequences of alphanumeric/underscore chars that are
    # preceded by '(' or ',' (start of a leaf) and followed by ':', ',', ')' or end
    # Also handle #H notation in extended Newick

    # First, let's find all tokens that look like leaf names
    # A leaf name appears after ( or , and before : or , or ) or #
    # It should not be empty and should not be a pure number (branch length)

    # Strategy: walk through removing nested structures
    # Simpler: use regex to find leaf positions

    # Pattern: after ( or , grab word characters until : , ) # or whitespace
    pattern = r'(?<=[(,])\s*([A-Za-z][A-Za-z0-9_]*)'

    for m in re.finditer(pattern, tree_str):
        name = m.group(1)
        # Skip if it's a hybrid node annotation like H1, H2 etc in #H format
        # Those appear as #H1 not as standalone leaves
        leaves.append(name)

    return leaves


def extract_leaf_names_robust(tree_str):
    """More robust leaf extraction that handles edge cases."""
    leaves = []
    i = 0
    n = len(tree_str)

    while i < n:
        # Skip whitespace
        if tree_str[i] in ' \t\n\r':
            i += 1
            continue

        # After ( or , we might have a leaf name
        if tree_str[i] in '(,':
            i += 1
            # Skip whitespace
            while i < n and tree_str[i] in ' \t\n\r':
                i += 1

            # Check if next char starts a name (letter)
            if i < n and (tree_str[i].isalpha() or tree_str[i] == '_'):
                # Read the name
                start = i
                while i < n and tree_str[i] not in ':,)(\t\n\r #;':
                    i += 1
                name = tree_str[start:i].strip()
                if name and not name.replace('.', '').replace('-', '').replace('e', '').replace('E', '').isdigit():
                    leaves.append(name)
            continue

        i += 1

    return leaves


def get_canonical_names(all_names_by_method):
    """Determine canonical species names from non-AlloppNET methods.

    Returns a set of canonical names. If no non-AlloppNET methods exist,
    returns None (AlloppNET names will be self-normalized).
    """
    canonical = set()
    for method, names in all_names_by_method.items():
        if method != 'alloppnet':
            canonical.update(set(names))

    return canonical if canonical else None


def build_alloppnet_rename_map(alloppnet_names, canonical_names):
    """Build a mapping from AlloppNET names to canonical names.

    Handles:
    1. Trailing digit suffixes: species0, species1 -> species
    2. Underscore removal: foe_1887 -> foe1887
    3. Complex name differences: cf_pinnatum -> cfpinnatum
    """
    rename_map = {}

    if canonical_names is None:
        # No other methods to compare against - just remove trailing digits
        for name in set(alloppnet_names):
            stripped = re.sub(r'(\d+)$', '', name)
            # Check if the stripped version + digit was the original
            # Only strip if it looks like an AlloppNET suffix (single digit at end)
            match = re.match(r'^(.+?)(\d)$', name)
            if match:
                base = match.group(1)
                # Only strip if the base name also exists or another variant exists
                # For safety, strip single trailing digits
                rename_map[name] = base
        return rename_map

    # Build lookup: canonical names without underscores (lowercased for matching)
    canonical_normalized = {}
    for cn in canonical_names:
        # Store normalized version -> original
        norm = cn.replace('_', '').lower()
        canonical_normalized[norm] = cn

    for name in set(alloppnet_names):
        if name in canonical_names:
            continue  # Already matches

        # Try progressively stripping trailing digits (1 digit, 2 digits, etc.)
        # This handles RSC040 -> RSC04 (strip 1 digit) vs madagascariense0 -> madagascariense
        best_match = None
        for n_strip in range(1, 5):
            if len(name) <= n_strip:
                break
            if not name[-n_strip:].isdigit():
                break
            stripped = name[:-n_strip]

            # Direct match
            if stripped in canonical_names:
                best_match = stripped
                break

            # Remove underscores
            no_us = stripped.replace('_', '')
            if no_us in canonical_names:
                best_match = no_us
                break

            # Case-insensitive + no underscores
            norm = no_us.lower()
            if norm in canonical_normalized:
                best_match = canonical_normalized[norm]
                break

        if best_match:
            rename_map[name] = best_match
            continue

        # Try with original name (no digit stripping) but remove underscores
        no_underscore_orig = name.replace('_', '')
        if no_underscore_orig in canonical_names:
            rename_map[name] = no_underscore_orig
            continue

        norm_orig = no_underscore_orig.lower()
        if norm_orig in canonical_normalized:
            rename_map[name] = canonical_normalized[norm_orig]
            continue

        # Try prefix matching: AlloppNET name (without underscores) is a prefix of a canonical name
        # e.g., cf_pinnatum -> cfpinnatum is prefix of cfpinnatumDLA2015
        no_us_name = name.replace('_', '')
        prefix_matches = [cn for cn in canonical_names if cn.startswith(no_us_name)]
        if len(prefix_matches) == 1:
            rename_map[name] = prefix_matches[0]
            continue

        # Also try prefix match after stripping trailing digits
        found_prefix = False
        for n_strip in range(1, 3):
            if len(name) <= n_strip or not name[-n_strip:].isdigit():
                break
            stripped_no_us = name[:-n_strip].replace('_', '')
            prefix_matches = [cn for cn in canonical_names if cn.startswith(stripped_no_us)]
            if len(prefix_matches) == 1:
                rename_map[name] = prefix_matches[0]
                found_prefix = True
                break

        if not found_prefix:
            # Last resort: if name ends in a single digit, strip it
            # (polyploid appearing twice with suffix 0/1)
            match = re.match(r'^(.+?)(\d)$', name)
            if match:
                rename_map[name] = match.group(1)

    return rename_map


def apply_rename_to_tree(tree_str, rename_map):
    """Apply renaming to a tree string, replacing longest names first."""
    result = tree_str
    # Sort by length descending to avoid partial replacements
    for old, new in sorted(rename_map.items(), key=lambda x: len(x[0]), reverse=True):
        # Use word boundary matching to avoid partial replacements
        # In Newick, names are bounded by ( , ) : # ; or whitespace
        pattern = r'(?<=[,(])\s*' + re.escape(old) + r'(?=\s*[,:)#;])'
        result = re.sub(pattern, new, result)

        # Also handle case where name is at the very start (rare but possible)
        if result.startswith(old):
            result = new + result[len(old):]

    return result


def analyze_dataset(dataset_path, dataset_name):
    """Analyze naming consistency for one dataset."""
    networks_dir = os.path.join(dataset_path, 'networks')
    if not os.path.isdir(networks_dir):
        return None

    methods = {}
    for method_name in os.listdir(networks_dir):
        method_dir = os.path.join(networks_dir, method_name)
        tree_file = os.path.join(method_dir, 'network.tre')

        if not os.path.isfile(tree_file):
            continue

        with open(tree_file) as f:
            content = f.read().strip()

        if not content:
            continue

        leaves = extract_leaf_names_robust(content)
        if not leaves:
            # Try the regex method as fallback
            leaves = extract_leaf_names(content)

        if leaves:
            methods[method_name] = {
                'leaves': leaves,
                'unique_leaves': sorted(set(leaves)),
                'tree_content': content,
                'tree_file': tree_file
            }

    return methods


def detect_genus_prefix(method_names, canonical_names):
    """Detect if a method's leaf names have a genus prefix not present in canonical names.

    For example, if method has ['Brachypodiumsylvaticum', 'Brachypodiumdistachyon']
    and canonical has ['sylvaticum', 'distachyon'], returns 'Brachypodium'.
    Returns None if no common prefix is detected.
    """
    if not method_names or not canonical_names:
        return None

    unique_method = set(method_names)
    # Check if any method names already match canonical - if so, no prefix issue
    if unique_method & canonical_names:
        return None

    # Find the longest common prefix among all method names
    sorted_names = sorted(unique_method)
    if len(sorted_names) < 2:
        return None

    first, last = sorted_names[0], sorted_names[-1]
    prefix = ""
    for i, c in enumerate(first):
        if i < len(last) and first[i] == last[i]:
            prefix += c
        else:
            break

    if len(prefix) < 3:
        return None

    # Check if stripping the prefix maps most method names to canonical names
    matches = 0
    for name in unique_method:
        stripped = name[len(prefix):]
        if stripped in canonical_names:
            matches += 1

    # Require at least 50% match rate
    if matches / len(unique_method) >= 0.5:
        return prefix

    return None


def build_prefix_rename_map(method_names, canonical_names, prefix):
    """Build rename map by stripping a genus prefix from method names.

    Returns dict mapping old names to canonical names.
    """
    rename_map = {}
    for name in set(method_names):
        if name.startswith(prefix):
            stripped = name[len(prefix):]
            if stripped in canonical_names:
                rename_map[name] = stripped
            else:
                # Try case-insensitive match
                canonical_lower = {cn.lower(): cn for cn in canonical_names}
                if stripped.lower() in canonical_lower:
                    rename_map[name] = canonical_lower[stripped.lower()]
    return rename_map


def find_issues(dataset_name, methods_data):
    """Find naming issues within a dataset across methods."""
    issues = []

    if not methods_data or len(methods_data) == 0:
        return issues

    # Get canonical names from non-AlloppNET methods
    # Find the actual alloppnet method name (could be 'allopnet', 'alloppnet', etc.)
    alloppnet_key = None
    for m in methods_data:
        if is_alloppnet_method(m):
            alloppnet_key = m
            break

    non_alloppnet = {m: d['unique_leaves'] for m, d in methods_data.items() if not is_alloppnet_method(m)}
    alloppnet_data = methods_data.get(alloppnet_key) if alloppnet_key else None

    # Build canonical names from all non-AlloppNET methods that don't have prefix issues
    # First pass: find methods with potential prefix issues
    # Use the method with the most overlap with others as the canonical source
    canonical = set()
    for names in non_alloppnet.values():
        canonical.update(names)

    # Check each method for genus prefix naming (e.g., grandma_split with "Brachypodium" prefix)
    for method_name, method_data in methods_data.items():
        if is_alloppnet_method(method_name):
            continue  # AlloppNET handled separately below

        method_unique = set(method_data['unique_leaves'])
        # Build canonical from OTHER methods only
        other_canonical = set()
        for m, d in methods_data.items():
            if m != method_name and not is_alloppnet_method(m):
                other_canonical.update(d['unique_leaves'])

        if not other_canonical:
            continue

        # Check if this method has a genus prefix
        prefix = detect_genus_prefix(method_data['leaves'], other_canonical)
        if prefix:
            rename_map = build_prefix_rename_map(method_data['leaves'], other_canonical, prefix)
            if rename_map:
                issues.append({
                    'type': 'prefix_naming',
                    'method': method_name,
                    'prefix': prefix,
                    'extra_names': sorted(method_unique - other_canonical),
                    'rename_map': rename_map
                })

    # Check consistency among non-AlloppNET methods (excluding already-flagged prefix issues)
    prefix_methods = {i['method'] for i in issues if i['type'] == 'prefix_naming'}
    if len(non_alloppnet) > 1:
        method_names = list(non_alloppnet.keys())
        for i in range(len(method_names)):
            for j in range(i + 1, len(method_names)):
                m1, m2 = method_names[i], method_names[j]
                if m1 in prefix_methods or m2 in prefix_methods:
                    continue  # Already handled by prefix detection
                s1 = set(non_alloppnet[m1])
                s2 = set(non_alloppnet[m2])
                only_in_1 = s1 - s2
                only_in_2 = s2 - s1
                if only_in_1 or only_in_2:
                    issues.append({
                        'type': 'cross_method_mismatch',
                        'methods': (m1, m2),
                        'only_in_first': sorted(only_in_1),
                        'only_in_second': sorted(only_in_2)
                    })

    # Check AlloppNET vs others
    if alloppnet_data and non_alloppnet:
        canonical = set()
        for names in non_alloppnet.values():
            canonical.update(names)

        alloppnet_unique = set(alloppnet_data['unique_leaves'])

        # Names in AlloppNET not in canonical
        extra = alloppnet_unique - canonical
        if extra:
            rename_map = build_alloppnet_rename_map(alloppnet_data['leaves'], canonical)
            issues.append({
                'type': 'alloppnet_naming',
                'extra_names': sorted(extra),
                'rename_map': rename_map
            })
    elif alloppnet_data and not non_alloppnet:
        # Only AlloppNET exists - check for trailing digits
        names_with_digits = [n for n in alloppnet_data['unique_leaves'] if re.search(r'\d$', n)]
        if names_with_digits:
            rename_map = build_alloppnet_rename_map(alloppnet_data['leaves'], None)
            if rename_map:
                issues.append({
                    'type': 'alloppnet_naming',
                    'extra_names': sorted(names_with_digits),
                    'rename_map': rename_map
                })

    return issues


def fix_dataset(dataset_name, methods_data, issues, dry_run=False):
    """Apply fixes to a dataset's network.tre files."""
    fixes_applied = []

    for issue in issues:
        if issue['type'] in ('alloppnet_naming', 'prefix_naming') and 'rename_map' in issue:
            rename_map = issue['rename_map']
            if not rename_map:
                continue

            # Determine which method to fix
            if issue['type'] == 'alloppnet_naming':
                target_key = None
                for m in methods_data:
                    if is_alloppnet_method(m):
                        target_key = m
                        break
            else:
                target_key = issue['method']

            target_data = methods_data.get(target_key) if target_key else None
            if not target_data:
                continue

            tree_file = target_data['tree_file']
            original_content = target_data['tree_content']
            new_content = apply_rename_to_tree(original_content, rename_map)

            if new_content != original_content:
                fixes_applied.append({
                    'file': tree_file,
                    'rename_map': rename_map,
                    'method': target_key
                })

                if not dry_run:
                    # Backup
                    bak_file = tree_file + '.bak'
                    if not os.path.exists(bak_file):
                        shutil.copy2(tree_file, bak_file)

                    with open(tree_file, 'w') as f:
                        f.write(new_content)

    return fixes_applied


def main():
    parser = argparse.ArgumentParser(description='Fix naming consistency in network.tre files')
    parser.add_argument('papers_dir', help='Path to papers directory')
    parser.add_argument('--dry-run', action='store_true', help='Only analyze, do not modify files')
    parser.add_argument('--dataset', help='Process only this dataset')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show detailed leaf names')
    args = parser.parse_args()

    if not os.path.isdir(args.papers_dir):
        print(f"Error: {args.papers_dir} is not a directory")
        sys.exit(1)

    # Find all datasets
    datasets = sorted(d for d in os.listdir(args.papers_dir)
                      if os.path.isdir(os.path.join(args.papers_dir, d, 'networks'))
                      and (not args.dataset or d == args.dataset))

    if not datasets:
        print("No datasets found with networks/ directory")
        sys.exit(1)

    total_fixes = 0

    for dataset_name in datasets:
        dataset_path = os.path.join(args.papers_dir, dataset_name)
        methods_data = analyze_dataset(dataset_path, dataset_name)

        if not methods_data:
            continue

        print(f"\n{'='*60}")
        print(f"Dataset: {dataset_name}")
        print(f"{'='*60}")
        print(f"Methods with output: {', '.join(sorted(methods_data.keys()))}")

        if args.verbose:
            for method, data in sorted(methods_data.items()):
                print(f"\n  {method}:")
                print(f"    Leaves ({len(data['unique_leaves'])}): {data['unique_leaves']}")

        issues = find_issues(dataset_name, methods_data)

        if not issues:
            print("  -> OK: All names consistent across methods")
            continue

        for issue in issues:
            if issue['type'] == 'cross_method_mismatch':
                m1, m2 = issue['methods']
                print(f"\n  WARNING: Mismatch between {m1} and {m2}")
                if issue['only_in_first']:
                    print(f"    Only in {m1}: {issue['only_in_first']}")
                if issue['only_in_second']:
                    print(f"    Only in {m2}: {issue['only_in_second']}")

            elif issue['type'] == 'prefix_naming':
                method = issue['method']
                prefix = issue['prefix']
                print(f"\n  Genus prefix naming in {method} (prefix: '{prefix}'):")
                print(f"    Names not matching other methods: {issue['extra_names']}")
                if issue.get('rename_map'):
                    print(f"    Proposed renames:")
                    for old, new in sorted(issue['rename_map'].items()):
                        print(f"      {old} -> {new}")

            elif issue['type'] == 'alloppnet_naming':
                print(f"\n  AlloppNET naming issues:")
                print(f"    Names not matching other methods: {issue['extra_names']}")
                if issue.get('rename_map'):
                    print(f"    Proposed renames:")
                    for old, new in sorted(issue['rename_map'].items()):
                        print(f"      {old} -> {new}")

        # Apply fixes
        fixes = fix_dataset(dataset_name, methods_data, issues, dry_run=args.dry_run)

        if fixes:
            total_fixes += len(fixes)
            for fix in fixes:
                action = "Would fix" if args.dry_run else "Fixed"
                print(f"\n  {action}: {fix['file']}")
                for old, new in sorted(fix['rename_map'].items()):
                    print(f"    {old} -> {new}")

    print(f"\n{'='*60}")
    if args.dry_run:
        print(f"DRY RUN: {total_fixes} files would be modified")
        print("Run without --dry-run to apply fixes")
    else:
        print(f"Done: {total_fixes} files modified (originals backed up as .bak)")


if __name__ == '__main__':
    main()
