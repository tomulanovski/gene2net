
'''
use: python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/create_notung_batch.py" path_to_files
'''

import os
import sys

def create_notung_batch(input_dir):
    output_file = os.path.join(input_dir, "notung_batch.txt")
    with open(output_file, "w") as f:
        for file in sorted(os.listdir(input_dir)):
            if file.endswith(".tre"):
                full_path = os.path.join(input_dir, file)
                f.write(full_path + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python make_notung_batch.py <directory>")
    else:
        input_dir = sys.argv[1]
        create_notung_batch(input_dir)
        print(f"notung_batch.txt created in {input_dir}")
