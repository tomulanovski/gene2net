import re
import sys
'''
path: python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/simplify_tree_for_mpallop.py" input_tree.tre output_tree.tre
'''

def clean_tree(tree_str):
    # Remove branch lengths (including all numeric formats)
    # This pattern matches: :1, :0.123, :1e-06, etc.
    cleaned_tree = re.sub(r':\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', '', tree_str)
    
    # Remove bootstrap values (numbers after parentheses)
    cleaned_tree = re.sub(r'\)\d+(?:\.\d+)?', ')', cleaned_tree)
    
    # Remove underscores from all characters in the tree string
    cleaned_tree = cleaned_tree.replace('_', '')
    
    return cleaned_tree

def process_file(input_file, output_file=None):
    with open(input_file, 'r') as f:
        content = f.read()
    
    cleaned_content = clean_tree(content)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(cleaned_content)
        print(f"Cleaned tree written to {output_file}")
    else:
        print(cleaned_content)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python clean_tree.py input_file [output_file]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    process_file(input_file, output_file)