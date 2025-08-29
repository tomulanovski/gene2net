"""
Check if a phylogenetic tree is rooted.
Usage: python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/is_rooted.py" "newick_string"
   or: python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/is_rooted.py" tree_file
Returns: True if rooted, False if unrooted
"""
import sys
import os
from ete3 import Tree

def is_rooted(tree_input):
    """Check if tree is rooted (root has 2 children, not >2)."""
    tree = None
    
    # Try to determine if input is file or string
    if os.path.isfile(tree_input):
        try:
            # Try common Newick formats
            for fmt in [0, 1, 2, 5]:  # Most common formats
                try:
                    tree = Tree(tree_input, format=fmt)
                    break
                except:
                    continue
            if tree is None:
                raise ValueError("Could not parse tree file with any standard format")
        except Exception as e:
            raise ValueError(f"Error reading file: {e}")
    else:
        # Treat as Newick string
        try:
            # Try common formats for string input
            for fmt in [0, 1]:
                try:
                    tree = Tree(tree_input, format=fmt)
                    break
                except:
                    continue
            if tree is None:
                raise ValueError("Could not parse Newick string")
        except Exception as e:
            raise ValueError(f"Error parsing Newick string: {e}")
    
    root = tree.get_tree_root()
    n_children = len(root.children)
    
    # Rooted tree: exactly 2 children at root
    # Unrooted tree: typically 3 or more children at root
    return n_children == 2

def main():
    if len(sys.argv) != 2:
        print('Usage: python is_rooted.py "newick_string"')
        print('   or: python is_rooted.py tree_file')
        sys.exit(1)
    
    tree_input = sys.argv[1]
    
    try:
        result = is_rooted(tree_input)
        print("True" if result else "False")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()