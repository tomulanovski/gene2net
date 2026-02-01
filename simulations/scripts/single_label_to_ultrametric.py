#!/usr/bin/env python3
"""
Convert single-label Newick trees to ultrametric trees and save as NEXUS.

Reads Newick files (handles :nan branch lengths by ignoring them; uses topology only).
Resolves polytomies, assigns branch lengths from topology, and writes ultrametric NEXUS.

Usage:
  python single_label_to_ultrametric.py <tree_height> [options]
  python single_label_to_ultrametric.py 10000000 --base_dir simulations/simulations/conf_ils_low_10M
  python single_label_to_ultrametric.py 10000000 --output_filename my_single.nex
"""
from ete3 import Tree
import sys
import os
import re
import argparse
from pathlib import Path


def _newick_to_topology_only(newick_str):
    """
    Strip branch lengths and internal labels to leave topology only (leaf names and parens).
    Use with format=9 when all other ETE3 formats fail.
    """
    s = newick_str.strip().rstrip(';')
    # Remove all branch lengths :value (value = number, nan, NA, etc.)
    s = re.sub(r':[^,\))]+', '', s)
    # Remove internal node labels: )label when followed by , or )
    s = re.sub(r'\)([^,\))]+)(?=[,\)])', ')', s)
    return s + ';' if not s.endswith(';') else s


def _parse_newick_string(newick_str):
    """
    Parse a Newick string with ETE3. Tries formats 0-9, then fallback to
    topology-only + format=9 so a wide variety of Newick files are accepted.
    """
    if not newick_str.strip().endswith(';'):
        newick_str = newick_str.strip() + ';'
    for fmt in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9):
        try:
            return Tree(newick_str, format=fmt)
        except Exception:
            continue
    # Fallback: strip to topology only (we only need topology for ultrametric conversion)
    try:
        topology_str = _newick_to_topology_only(newick_str)
        return Tree(topology_str, format=9)
    except Exception as e:
        raise ValueError(
            f"Could not parse Newick with any ETE3 format (0-9) or topology-only fallback: {e}"
        )


def _normalize_invalid_branch_lengths(newick_str):
    """Replace :nan, :NA, :na, etc. with :0 so ETE3 can parse (we ignore lengths for ultrametric)."""
    s = newick_str
    s = re.sub(r':\s*nan\b', ':0', s, flags=re.IGNORECASE)
    s = re.sub(r':\s*na\b', ':0', s, flags=re.IGNORECASE)  # R-style NA, or bare 'na'
    s = re.sub(r':\s*n\/a\b', ':0', s, flags=re.IGNORECASE)
    s = re.sub(r':\s*nan\s*([,\)])', r':0\1', s, flags=re.IGNORECASE)
    return s


def _fill_empty_leaves(newick_str):
    """
    Replace empty leaf nodes with placeholder _empty so ETE3 can parse.
    ETE3 requires every leaf to have a name; some Newick files have ,: or (,. etc.
    """
    s = newick_str
    s = s.replace('(,', '(_empty,')
    s = s.replace(',)', ',_empty)')
    s = s.replace(',,', ',_empty,')
    # Empty leaf before branch length: ,: or (,:  -> ,_empty: or (_empty:
    s = re.sub(r',\s*:', ',_empty:', s)
    s = re.sub(r'\(\s*:', '(_empty:', s)
    return s


def read_newick_tree(filepath):
    """
    Read a tree from a Newick file. Handles :nan, :NA, etc. by replacing with :0.
    Tries ETE3 formats 0-9, then topology-only fallback.
    """
    try:
        with open(filepath, 'r') as f:
            content = f.read().strip()
        content = _normalize_invalid_branch_lengths(content)
        content = _fill_empty_leaves(content)
        tree = _parse_newick_string(content)
        return tree
    except Exception as e:
        print(f"Error reading Newick file {filepath}: {e}")
        sys.exit(1)


def read_nexus_tree(filepath):
    """Read a tree from NEXUS file (optional support)."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        for line in content.split('\n'):
            line = line.strip()
            if line.startswith('tree ') and '=' in line:
                newick_str = line.split('=', 1)[1].strip().rstrip(';').strip()
                newick_str = _normalize_invalid_branch_lengths(newick_str)
                newick_str = _fill_empty_leaves(newick_str)
                tree = _parse_newick_string(newick_str)
                return tree
        raise ValueError("No tree definition found in NEXUS file")
    except Exception as e:
        print(f"Error reading NEXUS file {filepath}: {e}")
        sys.exit(1)


def read_tree(filepath):
    """Read tree from Newick or NEXUS file (auto-detect)."""
    with open(filepath, 'r') as f:
        peek = f.read(20)
    if peek.strip().upper().startswith('#NEXUS'):
        return read_nexus_tree(filepath)
    return read_newick_tree(filepath)


def resolve_polytomies(tree, seed=None):
    """Resolve all polytomies randomly with temporary branch lengths."""
    if seed is not None:
        import random
        random.seed(seed)
    tree.resolve_polytomy(recursive=True, default_dist=1)
    return tree


def make_ultrametric_from_topology(tree, target_height=1000000):
    """
    Convert a topology to an ultrametric tree with EXACT ultrametricity.
    Assigns heights from topology only; all tips end at target_height.
    """
    node_depths_from_tips = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node_depths_from_tips[node] = 0
        else:
            node_depths_from_tips[node] = 1 + max(
                node_depths_from_tips[child] for child in node.children
            )
    max_depth = node_depths_from_tips[tree]
    if max_depth == 0:
        print("  WARNING: Tree has only one node")
        return tree
    node_heights = {}
    for node, depth_from_tips in node_depths_from_tips.items():
        node_heights[node] = target_height * (max_depth - depth_from_tips) / max_depth
    for node in tree.traverse():
        if not node.is_root() and not node.is_leaf():
            parent = node.up
            node.dist = node_heights[node] - node_heights[parent]
    for leaf in tree.iter_leaves():
        dist_to_parent = 0.0
        node = leaf.up
        while not node.is_root():
            dist_to_parent += node.dist
            node = node.up
        leaf.dist = max(0.0, round(target_height - dist_to_parent, 12))
    return tree


def verify_ultrametric(tree, tolerance=1.0):
    """Verify that the tree is ultrametric by checking all root-to-tip distances."""
    distances = [tree.get_distance(leaf) for leaf in tree.iter_leaves()]
    if not distances:
        return True, 0
    max_dist = max(distances)
    min_dist = min(distances)
    if abs(max_dist - min_dist) < tolerance:
        return True, max_dist
    return False, (min_dist, max_dist)


def tree_to_nexus_newick(tree):
    """Convert ETE3 tree to Newick string with taxa names and branch lengths."""
    def format_node(node):
        name = node.name if node.is_leaf() and node.name else ""
        if node.is_leaf():
            result = name
        else:
            children = ",".join(format_node(child) for child in node.children)
            result = f"({children})"
        if not node.is_root():
            result += f":{node.dist:.12f}"
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
        print(f"  Successfully wrote ultrametric tree to {output_filepath}")
    except Exception as e:
        print(f"Error writing NEXUS file: {e}")
        sys.exit(1)


# Same network list as rescale_and_keep_ultrametric.py
NETWORKS = [
    "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
    "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019",
    "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011",
    "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014",
    "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014",
]


def process_all_trees(base_dir, input_filename, output_filename, tree_height, seed=None):
    """Process all single-label trees in the base directory (one per network)."""
    base_path = Path(base_dir)
    processed_count = 0
    for network in NETWORKS:
        input_path = base_path / network / input_filename
        if not input_path.exists():
            print(f"WARNING: Input not found: {input_path}")
            continue
        print(f"\nProcessing: {network}")
        print(f"  Reading from: {input_path}")
        tree = read_tree(str(input_path))
        print(f"  Tree loaded ({tree.get_leaf_count()} tips)")
        tree = resolve_polytomies(tree, seed=seed)
        print(f"  Polytomies resolved")
        tree = make_ultrametric_from_topology(tree, target_height=tree_height)
        print(f"  Converted to ultrametric (height={tree_height:,})")
        is_ultrametric, result = verify_ultrametric(tree, tolerance=1e-6)
        if is_ultrametric:
            print(f"  Verified ultrametric: all tips at distance {result:.2f}")
        else:
            print(f"  WARNING: Not ultrametric! Distance range: {result[0]:.2f} - {result[1]:.2f}")
        output_path = base_path / network / output_filename
        write_nexus(tree, str(output_path), tree_name="1")
        processed_count += 1
    print(f"\n{'='*60}")
    print(f"Done. Processed {processed_count}/{len(NETWORKS)} trees. Target height: {tree_height:,}")
    if seed is not None:
        print(f"Random seed: {seed}")
    print(f"{'='*60}")


def main():
    script_dir = Path(__file__).resolve().parent
    default_base = script_dir.parent / "simulations"  # simulations/simulations/ (Bendiksby_2011/, etc.)
    parser = argparse.ArgumentParser(
        description='Convert single-label Newick trees to ultrametric NEXUS (topology-based).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python single_label_to_ultrametric.py 10000000
  python single_label_to_ultrametric.py 10000000 --base_dir simulations/simulations/conf_ils_low_10M
  python single_label_to_ultrametric.py 10000000 --output_filename single_species.nex
  python single_label_to_ultrametric.py 5000000 --input_filename astral_4.tre --seed 42
        """
    )
    parser.add_argument('tree_height', type=float,
                        help='Target height for ultrametric tree (e.g. 10000000)')
    parser.add_argument('--base_dir', default=None,
                        help='Base directory with network subdirs (default: simulations/ relative to script)')
    parser.add_argument('--input_filename', default='astral_4.tre',
                        help='Input Newick filename per network (default: astral_4.tre)')
    parser.add_argument('--output_filename', default='single_label.nex',
                        help='Output NEXUS filename per network (default: single_label.nex)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for polytomy resolution (optional)')
    args = parser.parse_args()
    base_dir = args.base_dir if args.base_dir is not None else str(default_base)
    if not os.path.isdir(base_dir):
        print(f"ERROR: Base directory not found: {base_dir}")
        sys.exit(1)
    if args.tree_height <= 0:
        print(f"ERROR: Tree height must be positive, got: {args.tree_height}")
        sys.exit(1)
    print(f"Base directory: {base_dir}")
    print(f"Input filename: {args.input_filename}")
    print(f"Output filename: {args.output_filename}")
    print(f"Target height: {args.tree_height:,}")
    if args.seed is not None:
        print(f"Seed: {args.seed}")
    process_all_trees(base_dir, args.input_filename, args.output_filename,
                      args.tree_height, seed=args.seed)


if __name__ == "__main__":
    main()
