
'''
use: python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/concatenate_trees.py" /path/to/rearrange/files /path/to/output/all_trees.tre
'''

import os
import sys

def concat_rearrange_files(input_dir, output_file):
    with open(output_file, "w") as outfile:
        for file in sorted(os.listdir(input_dir)):
            if file.endswith(".resolve.0"):
                with open(os.path.join(input_dir, file), "r") as infile:
                    outfile.write(infile.read().strip() + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python concat_rearrange.py <directory> <output_file>")
    else:
        input_dir = sys.argv[1]
        output_file = sys.argv[2]
        concat_rearrange_files(input_dir, output_file)
        print(f"Concatenated .rearrange.0 files from {input_dir} into {output_file}")