
'''
use: python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/create_notung_batch.py" path_to_files
'''

import os
import sys

def create_notung_batch(input_dir, species_tree):
    output_file = os.path.join(input_dir, "notung_batch.txt")
    
    with open(output_file, "w") as f:
        # Write species tree path as first line
        f.write(species_tree + "\n")
        
        # Write all gene tree paths
        for file in sorted(os.listdir(input_dir)):
            if file.endswith(".tre"):
                full_path = os.path.join(input_dir, file)
                f.write(full_path + "\n")
    
    print(f"? notung_batch.txt created in {input_dir}")
    print(f"? Species tree: {species_tree}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python create_notung_batch.py <directory> <species_tree_path>")
        sys.exit(1)
    else:
        input_dir = sys.argv[1]
        species_tree = sys.argv[2]
        create_notung_batch(input_dir, species_tree)