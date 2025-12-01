#!/usr/bin/env python3
import sys
import re

def extract_species(gene_name, mode):
    parts = gene_name.split("_")
    
    if mode == "gene":
        # remove last 2 underscore parts
        if len(parts) < 3:
            return gene_name   # fallback (not enough parts)
        return "_".join(parts[:-2])
    
    elif mode == "locus":
        # remove last 1 underscore part
        if len(parts) < 2:
            return gene_name
        return "_".join(parts[:-1])
    
    else:
        print(f"ERROR: mode must be gene or locus, got {mode}", file=sys.stderr)
        sys.exit(1)

def add_nhx(tree_string, mode):
    def repl(m):
        gene = m.group(1)
        species = extract_species(gene, mode)
        return f"{gene}[&&NHX:S={species}]:"
    
    return re.sub(r'([A-Za-z0-9_]+):', repl, tree_string)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python add_nhx_single.py <gene|locus> <input_tree> <output_tree>", file=sys.stderr)
        sys.exit(1)
    
    mode = sys.argv[1]
    infile = sys.argv[2]
    outfile = sys.argv[3]
    
    with open(infile) as f:
        t = f.read().strip()
    
    out = add_nhx(t, mode)
    
    with open(outfile, "w") as f:
        f.write(out)
        if not out.endswith("\n"):
            f.write("\n")
    
    print(f"done {infile} -> {outfile}", file=sys.stderr)
