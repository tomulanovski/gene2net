#!/usr/bin/env python3
"""
Process species trees: convert topology to ultrametric tree of specified height.
Usage: python "/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/rescale_and_keep_ultrametric.py" <output_filename> <tree_height> [--base_dir PATH] [--input_filename NAME] [--seed SEED]
"""
from ete3 import Tree
import sys
import os
import argparse
from pathlib import Path

def read_nexus_tree(filepath):
    """Read a tree from NEXUS file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Extract Newick string from NEXUS format
        for line in content.split('\n'):
            line = line.strip()
            if line.startswith('tree ') and '=' in line:
                newick_str = line.split('=', 1)[1].strip()
                newick_str = newick_str.rstrip(';').strip()
                tree = Tree(newick_str + ';', format=1)
                return tree
        
        raise ValueError("No tree definition found in NEXUS file")
    except Exception as e:
        print(f"Error reading NEXUS file {filepath}: {e}")
        sys.exit(1)

def resolve_polytomies(tree, seed=None):
    """Resolve all polytomies randomly with temporary branch lengths."""
    if seed is not None:
        import random
        random.seed(seed)
    
    # Use default_dist=1 just for resolution (will be recalculated)
    tree.resolve_polytomy(recursive=True, default_dist=1)
    return tree

def make_ultrametric_from_topology(tree, target_height=1000000):
    """
    Convert a topology to an ultrametric tree with EXACT ultrametricity.
    
    Strategy: Assign heights to nodes, then adjust terminal branches
    to guarantee all tips are at exactly target_height.
    """
    # Step 1: Calculate depth from tips for each node
    node_depths_from_tips = {}
    
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node_depths_from_tips[node] = 0
        else:
            node_depths_from_tips[node] = 1 + max(node_depths_from_tips[child] 
                                                    for child in node.children)
    
    # Step 2: Maximum depth is at the root
    max_depth = node_depths_from_tips[tree]
    
    if max_depth == 0:
        print("  WARNING: Tree has only one node")
        return tree
    
    # Step 3: Assign heights to nodes
    # Use simple proportion but we'll fix terminal branches later
    node_heights = {}
    for node, depth_from_tips in node_depths_from_tips.items():
        node_heights[node] = target_height * (max_depth - depth_from_tips) / max_depth
    
    # Step 4: Set branch lengths for INTERNAL branches
    for node in tree.traverse():
        if not node.is_root() and not node.is_leaf():
            parent = node.up
            node.dist = node_heights[node] - node_heights[parent]
    
    # Step 5: Set TERMINAL branches by working backwards from tips
    # This GUARANTEES each path sums to exactly target_height
    for leaf in tree.iter_leaves():
        # Calculate distance from root to parent of this leaf
        dist_to_parent = 0.0
        node = leaf.up
        while not node.is_root():
            dist_to_parent += node.dist
            node = node.up
        
        # Set terminal branch so total = target_height EXACTLY
        leaf.dist = max(0.0, round(target_height - dist_to_parent, 12))
    
    return tree

def verify_ultrametric(tree, tolerance=1.0):
    """Verify that the tree is ultrametric by checking all root-to-tip distances."""
    distances = []
    for leaf in tree.iter_leaves():
        dist = tree.get_distance(leaf)
        distances.append(dist)
    
    if not distances:
        return True, 0
    
    max_dist = max(distances)
    min_dist = min(distances)
    
    # Check if all distances are equal (within tolerance)
    if abs(max_dist - min_dist) < tolerance:
        return True, max_dist
    else:
        return False, (min_dist, max_dist)

def tree_to_nexus_newick(tree):
    """Convert ETE3 tree to Newick string with only taxa names and branch lengths."""
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

def process_all_trees(base_dir, input_filename, output_filename, tree_height, seed=None):
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
        species_tree_path = base_path / network / input_filename
        
        if not species_tree_path.exists():
            print(f"WARNING: Species tree not found: {species_tree_path}")
            continue
        
        print(f"\nProcessing: {network}")
        print(f"  Reading from: {species_tree_path}")
        
        # Read the original tree (topology only, no branch lengths needed)
        tree = read_nexus_tree(species_tree_path)
        print(f"  Original tree loaded")
        
        # Resolve polytomies if any exist
        tree = resolve_polytomies(tree, seed=seed)
        print(f"  Polytomies resolved")
        
        # Convert topology to ultrametric tree with specified height
        # This is the KEY function that guarantees exact ultrametricity
        tree = make_ultrametric_from_topology(tree, target_height=tree_height)
        print(f"  Converted to ultrametric (height={tree_height:,})")
        
        # Verify ultrametric property
        is_ultrametric, result = verify_ultrametric(tree, tolerance=1e-6)
        if is_ultrametric:
            print(f"  ? Verified ultrametric: all tips at distance {result:.2f}")
        else:
            print(f"  ? WARNING: Not ultrametric! Distance range: {result[0]:.2f} - {result[1]:.2f}")
        
        # Write the processed tree
        output_path = base_path / network / output_filename
        write_nexus(tree, output_path, tree_name="1")
        
        processed_count += 1
    
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"Successfully processed {processed_count}/{len(networks)} trees")
    print(f"Target height: {tree_height:,}")
    if seed is not None:
        print(f"Random seed used: {seed}")
    print(f"{'='*60}")

def main():
    parser = argparse.ArgumentParser(
        description='Process species trees: resolve polytomies and make ultrametric',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python rescale_and_ultrametric.py output.nex 1000000
  python rescale_and_ultrametric.py ultrametric.nex 500000 --seed 42
  python rescale_and_ultrametric.py out.nex 1000000 --base_dir /custom/path --input_filename input.nex
        """
    )
    parser.add_argument('output_filename',
                       help='Name of output tree file (e.g., species_tree_ultrametric.nex)')
    parser.add_argument('tree_height', type=float,
                       help='Target height for ultrametric tree')
    parser.add_argument('--base_dir', 
                       default='/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations/',
                       help='Base directory containing network subdirectories (default: /groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations/)')
    parser.add_argument('--input_filename',
                       default='species_tree.nex',
                       help='Name of input tree file (default: species_tree.nex)')
    parser.add_argument('--seed', type=int, default=None,
                       help='Random seed for polytomy resolution (optional)')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.base_dir):
        print(f"ERROR: Base directory not found: {args.base_dir}")
        sys.exit(1)
    
    if args.tree_height <= 0:
        print(f"ERROR: Tree height must be positive, got: {args.tree_height}")
        sys.exit(1)
    
    print(f"Base directory: {args.base_dir}")
    print(f"Input filename: {args.input_filename}")
    print(f"Output filename: {args.output_filename}")
    print(f"Target height: {args.tree_height:,}")
    if args.seed is not None:
        print(f"Using random seed: {args.seed}")
    print(f"Processing all species trees...")
    
    process_all_trees(args.base_dir, args.input_filename, args.output_filename, 
                     args.tree_height, seed=args.seed)

if __name__ == "__main__":
    main()