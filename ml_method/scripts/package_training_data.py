#!/usr/bin/env python
"""Package SimPhy + ASTRAL output into Gene2Net training samples (.pt files).

Reads the generated MUL-trees, SimPhy gene trees, and ASTRAL species trees,
then creates Gene2NetSample objects and saves them for training.

Usage:
    # Package one example (called by SLURM array job)
    python package_training_data.py --index 42 --mul-trees-dir /path/to/mul_trees --config ils_low

    # Package a batch locally
    python package_training_data.py --index 0 --n-batch 10 --mul-trees-dir /path/to/mul_trees --config ils_low
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree
from gene2net_gnn.data.dataset import Gene2NetSample


def load_gene_trees(gene_trees_file, max_trees=500):
    """Load gene trees from concatenated file (one Newick per line)."""
    trees = []
    with open(gene_trees_file) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    trees.append(Tree(line, format=1))
                except Exception:
                    continue
            if len(trees) >= max_trees:
                break
    return trees


def load_nexus_tree(filepath):
    """Load a tree from NEXUS format file."""
    with open(filepath) as f:
        content = f.read()
    for line in content.split("\n"):
        line = line.strip()
        if line.startswith("tree ") and "=" in line:
            newick_str = line.split("=", 1)[1].strip()
            return Tree(newick_str, format=1)
    raise ValueError(f"No tree found in NEXUS file: {filepath}")


def load_astral_tree(filepath):
    """Load ASTRAL species tree (plain Newick)."""
    with open(filepath) as f:
        newick = f.read().strip()
    return Tree(newick, format=1)


def package_one(index, mul_trees_dir, config, replicate=1, max_gene_trees=500):
    """Package one training example from existing files."""
    idx_str = f"{index:04d}"

    # Paths
    mul_tree_path = os.path.join(mul_trees_dir, f"mul_tree_{idx_str}.nex")
    metadata_path = os.path.join(mul_trees_dir, f"metadata_{idx_str}.json")
    rep_dir = os.path.join(mul_trees_dir, "simphy", config, idx_str, f"replicate_{replicate}")
    gene_trees_path = os.path.join(rep_dir, "gene_trees.tre")
    astral_path = os.path.join(rep_dir, "astral_species.tre")

    # Validate all files exist
    for path, desc in [
        (mul_tree_path, "MUL-tree"),
        (gene_trees_path, "gene trees"),
        (astral_path, "ASTRAL species tree"),
    ]:
        if not os.path.exists(path):
            print(f"  SKIP [{idx_str}]: {desc} not found: {path}")
            return None

    # Load everything
    mul_tree = load_nexus_tree(mul_tree_path)
    gene_trees = load_gene_trees(gene_trees_path, max_trees=max_gene_trees)
    astral_tree = load_astral_tree(astral_path)

    if len(gene_trees) == 0:
        print(f"  SKIP [{idx_str}]: No gene trees loaded")
        return None

    # Get species list from ASTRAL tree (single-label, no duplicates)
    species_list = sorted(set(astral_tree.get_leaf_names()))

    # Build training sample
    sample = Gene2NetSample.from_trees(
        species_tree=astral_tree,
        gene_trees=gene_trees,
        species_list=species_list,
        mul_tree=mul_tree,
    )

    return sample


def main():
    parser = argparse.ArgumentParser(description="Package training data for Gene2Net-GNN")
    parser.add_argument("--index", type=int, required=True,
                        help="Starting index (or SLURM_ARRAY_TASK_ID)")
    parser.add_argument("--n-batch", type=int, default=1,
                        help="Number of examples to package (batch mode)")
    parser.add_argument("--mul-trees-dir", required=True,
                        help="Directory with mul_tree_NNNN.nex files and simphy/ output")
    parser.add_argument("--config", required=True,
                        help="SimPhy config name (e.g., ils_low)")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory for .pt samples (default: mul_trees_dir/training/config/)")
    parser.add_argument("--replicate", type=int, default=1,
                        help="Which replicate to use (default: 1)")
    parser.add_argument("--max-gene-trees", type=int, default=500,
                        help="Max gene trees per sample (default: 500)")
    args = parser.parse_args()

    output_dir = args.output_dir or os.path.join(args.mul_trees_dir, "training", args.config)
    os.makedirs(output_dir, exist_ok=True)

    packaged = 0
    skipped = 0

    for i in range(args.n_batch):
        idx = args.index + i
        sample = package_one(
            idx, args.mul_trees_dir, args.config,
            replicate=args.replicate,
            max_gene_trees=args.max_gene_trees,
        )

        if sample is None:
            skipped += 1
            continue

        # Save
        sample_dir = os.path.join(output_dir, f"sample_{idx:04d}")
        sample.save(sample_dir)
        packaged += 1

        n_wgd = 0
        if sample.labels and sample.labels.wgd_counts is not None:
            n_wgd = int((sample.labels.wgd_counts > 0).sum())

        print(f"[{idx:04d}] species={sample.n_species}, "
              f"gene_trees={len(sample.gene_tree_edge_indices)}, "
              f"wgd_edges={n_wgd} -> {sample_dir}")

    print(f"\nDone: {packaged} packaged, {skipped} skipped")


if __name__ == "__main__":
    main()
