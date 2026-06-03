#!/usr/bin/env python
"""Append the 5 WGD-detection edge features to existing packaged samples.

Avoids a full repackage. For each already-packaged sample this:
  1. Re-reads the source gene trees + ASTRAL species tree (needed because the
     copy-pair divergence feature requires faithful branch lengths, which the
     stored .pkl tensors do not preserve — packaging sanitises 0-length branches
     to 1.0).
  2. Computes ONLY the new detection features (synchrony, mirrored-sister,
     copy-pair divergence mean/cv, clade-duplicated fraction).
  3. Appends them to the sample's stored edge_features (4 -> 9 dims) in place.

The expensive node features (co-clustering) and base edge features
(concordance) are NOT recomputed — they are read from the existing .pkl
untouched. Edge ordering matches the stored features because the same ASTRAL
tree is re-read and edges are enumerated in the identical preorder.

Usage:
    # One batch (called by SLURM array job)
    python augment_edge_features.py --index 0 --n-batch 20 \
        --mul-trees-dir /path/to/mul_trees_2k --config ils_low
"""
import argparse
import os
import pickle
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
from ete3 import Tree

from gene2net_gnn.data.features import compute_species_tree_edge_detection_features

# Order must match dataset.py's edge feature assembly.
DET_KEYS = [
    "dup_synchrony",
    "mirrored_sister_frac",
    "copy_pair_div_mean",
    "copy_pair_div_cv",
    "frac_clade_duplicated",
]
BASE_EDGE_DIM = 4
FULL_EDGE_DIM = 9


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


def load_astral_tree(filepath):
    """Load ASTRAL species tree (plain Newick)."""
    with open(filepath) as f:
        newick = f.read().strip()
    return Tree(newick, format=1)


def augment_one(index, mul_trees_dir, config, replicate=1, max_gene_trees=500):
    """Append detection features to one packaged sample. Returns a status string."""
    idx_str = f"{index:04d}"

    sample_dir = os.path.join(mul_trees_dir, "training", config, f"sample_{idx_str}")
    pkl_path = os.path.join(sample_dir, "sample.pkl")
    if not os.path.exists(pkl_path):
        return f"SKIP [{idx_str}]: no packaged sample at {pkl_path}"

    with open(pkl_path, "rb") as f:
        data = pickle.load(f)

    edge_feats = data.get("species_tree_edge_features")
    if edge_feats is None:
        return f"SKIP [{idx_str}]: sample has no edge features"
    if edge_feats.shape[1] >= FULL_EDGE_DIM:
        return f"SKIP [{idx_str}]: already augmented ({edge_feats.shape[1]} dims)"
    if edge_feats.shape[1] != BASE_EDGE_DIM:
        return f"SKIP [{idx_str}]: unexpected edge dim {edge_feats.shape[1]}"
    n_edges = edge_feats.shape[0]

    rep_dir = os.path.join(mul_trees_dir, "simphy", config, idx_str, f"replicate_{replicate}")
    gene_trees_path = os.path.join(rep_dir, "gene_trees.tre")
    astral_path = os.path.join(rep_dir, "astral_species.tre")
    if not os.path.exists(gene_trees_path) or not os.path.exists(astral_path):
        return f"SKIP [{idx_str}]: source trees missing in {rep_dir}"

    astral_tree = load_astral_tree(astral_path)
    gene_trees = load_gene_trees(gene_trees_path, max_trees=max_gene_trees)
    if len(gene_trees) == 0:
        return f"SKIP [{idx_str}]: no gene trees loaded"

    det = compute_species_tree_edge_detection_features(astral_tree, gene_trees)
    if len(det) != n_edges:
        return (f"SKIP [{idx_str}]: edge count mismatch "
                f"(stored {n_edges}, recomputed {len(det)}) — ASTRAL tree differs")

    extra = [[det.get(i, {}).get(k, 0.0) for k in DET_KEYS] for i in range(n_edges)]
    extra_t = torch.tensor(extra, dtype=torch.float)
    data["species_tree_edge_features"] = torch.cat([edge_feats, extra_t], dim=1)

    with open(pkl_path, "wb") as f:
        pickle.dump(data, f)

    return f"[{idx_str}] edge feats {BASE_EDGE_DIM} -> {FULL_EDGE_DIM}"


def main():
    parser = argparse.ArgumentParser(description="Augment packaged samples with detection edge features")
    parser.add_argument("--index", type=int, required=True, help="Starting index (or SLURM_ARRAY_TASK_ID)")
    parser.add_argument("--n-batch", type=int, default=1, help="Number of samples to process")
    parser.add_argument("--mul-trees-dir", required=True, help="Directory with training/ and simphy/ output")
    parser.add_argument("--config", required=True, help="Config name (e.g., ils_low)")
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--max-gene-trees", type=int, default=500)
    args = parser.parse_args()

    done = 0
    skipped = 0
    for i in range(args.n_batch):
        msg = augment_one(
            args.index + i, args.mul_trees_dir, args.config,
            replicate=args.replicate, max_gene_trees=args.max_gene_trees,
        )
        print(msg)
        if msg.startswith("SKIP"):
            skipped += 1
        else:
            done += 1

    print(f"\nDone: {done} augmented, {skipped} skipped")


if __name__ == "__main__":
    main()
