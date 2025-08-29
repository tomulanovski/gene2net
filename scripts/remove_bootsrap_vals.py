#!/usr/bin/env python3

import re
import sys

def remove_support_values(tree_string):
    """
    Remove support values from a tree string.
    Support values are numbers that appear after a closing parenthesis and before a colon.
    Example: (((A:1.0,B:1.0)100:1.0,(C:1.0,D:1.0)100:1.0)90:1.0,E:1.0) 
             becomes (((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0):1.0,E:1.0)
    """
    # Regular expression to match support values
    # Looks for numbers between a closing parenthesis and a colon
    pattern = r'\)(\d+):'
    
    # Replace with just the closing parenthesis and colon
    cleaned_tree = re.sub(pattern, '):', tree_string)
    
    return cleaned_tree

def process_file(input_file, output_file=None):
    """
    Process a file containing tree strings, removing support values from each line.
    
    Args:
        input_file: Path to the input file
        output_file: Path to the output file (optional, defaults to stdout)
    """
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        cleaned_lines = [remove_support_values(line) for line in lines]
        
        if output_file:
            with open(output_file, 'w') as f:
                f.writelines(cleaned_lines)
            print(f"Processed trees written to {output_file}")
        else:
            for line in cleaned_lines:
                print(line, end='')
            
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        sys.exit(1)

def main():
    if len(sys.argv) < 2:
        print("Usage: python remove_support_values.py input_file [output_file]")
        print("If output_file is not specified, results will be printed to stdout.")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    process_file(input_file, output_file)

if __name__ == "__main__":
    main()