#!/usr/bin/env python3
"""
Reduce a FASTA alignment to only the sequences whose IDs match the
leaf names of a (reduced) gene tree.

Usage:
    python reduce_fasta_by_tree.py <tree.tre> <alignment.fasta> <output.fasta>
"""
import sys
from ete3 import Tree
from Bio import SeqIO

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} <tree.tre> <alignment.fasta> <output.fasta>")
    sys.exit(1)

tree_file  = sys.argv[1]
fasta_file = sys.argv[2]
out_file   = sys.argv[3]

# --- Get leaf names from tree ---
try:
    tree = Tree(tree_file, format=1)
except:
    tree = Tree(tree_file, format=0)

leaf_names = {leaf.name for leaf in tree.get_leaves()}
print(f"Leaf names in tree ({len(leaf_names)}):")
for n in sorted(leaf_names):
    print(f"  {n}")

# --- Filter FASTA ---
kept, missed = [], []
for rec in SeqIO.parse(fasta_file, "fasta"):
    if rec.id in leaf_names:
        kept.append(rec)
    else:
        missed.append(rec.id)

SeqIO.write(kept, out_file, "fasta")
print(f"\nKept {len(kept)} sequences → {out_file}")

if missed:
    print(f"\nSequences in FASTA NOT in tree ({len(missed)}) - skipped:")
    for m in missed:
        print(f"  {m}")

not_found = leaf_names - {r.id for r in kept}
if not_found:
    print(f"\nWARNING: Tree leaves NOT found in FASTA ({len(not_found)}):")
    for n in sorted(not_found):
        print(f"  {n}")
