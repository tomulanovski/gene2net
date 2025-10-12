#!/usr/bin/env python3
import sys
import csv
import os
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

'''
python script that takes copies tsv file and input and output dir and reduce the fasta files in input dir to the labels in the tsv file and write them with _reduce extension in the name to the output dir

usage: python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/reduce_fasta.py" labels.tsv path/input_dir path/output_dir

'''

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <labels.tsv> <input_dir> <output_dir>")
    sys.exit(1)

tsv_file = sys.argv[1]
input_dir = sys.argv[2]
output_dir = sys.argv[3]
file_format = "fasta"  # change if needed

# Read first column of TSV (skip header)
labels_to_keep = set()
with open(tsv_file, newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)  # skip header
    for row in reader:
        if row:
            labels_to_keep.add(row[0].strip())

print(f"Loaded {len(labels_to_keep)} labels from {tsv_file}")

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Process all .fasta files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        input_path = os.path.join(input_dir, filename)
        
        # Add "_reduced" before file extension
        name, ext = os.path.splitext(filename)
        output_filename = f"{name}_reduced{ext}"
        output_path = os.path.join(output_dir, output_filename)
        
        # Read alignment
        alignment = AlignIO.read(input_path, file_format)
        
        # Filter by exact ID match
        filtered_records = [rec for rec in alignment if rec.id in labels_to_keep]
        
        # Write filtered alignment
        AlignIO.write(MultipleSeqAlignment(filtered_records), output_path, file_format)
        
        print(f"{filename}: kept {len(filtered_records)} sequences out of {len(alignment)}")
        print(f"Saved filtered alignment to: {output_path}")
