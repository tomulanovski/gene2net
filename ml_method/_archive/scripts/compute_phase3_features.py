"""Compute Phase 3 gene tree topology features and save augmented samples.

Reads existing Phase 1 .pkl samples, computes 3 new per-edge topology
features from gene tree structure, concatenates with existing 4-dim edge
features to produce 7-dim edge features.

Unlike Phase 2 (copy-count aggregations, which GAT already learns),
these features capture WHERE paralogs sit in the gene tree — the key
WGD signature.

Usage:
    python scripts/compute_phase3_features.py \
        --input-dir data/mul_trees_2k/training/ils_low \
        --output-dir data/mul_trees_2k/training_phase3/ils_low

    # Limit gene trees for speed
    python scripts/compute_phase3_features.py \
        --input-dir data/mul_trees_2k/training/ils_low \
        --output-dir data/mul_trees_2k/training_phase3/ils_low \
        --max-gene-trees 200
"""
import argparse
import os
import sys
import time

import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.data.gene_tree_features import get_edge_clades
from gene2net_gnn.data.gene_tree_topology import compute_topology_features


def augment_sample(sample: Gene2NetSample, max_gene_trees: int = 0) -> Gene2NetSample:
    """Add Phase 3 topology features to a sample's edge features.

    Args:
        sample: Original sample with 4-dim edge features.
        max_gene_trees: Max gene trees to use (0 = all).

    Returns:
        Sample with 7-dim edge features (original 4 + 3 topology).
    """
    edge_index = sample.species_tree_edge_index
    is_leaf = sample.species_tree_is_leaf
    n_edges = edge_index.shape[1] // 2

    # Get clade membership for each edge
    edge_clades = get_edge_clades(edge_index, is_leaf)

    # Optionally subsample gene trees
    ei_list = sample.gene_tree_edge_indices
    sp_ids_list = sample.gene_tree_species_ids
    leaf_masks_list = sample.gene_tree_leaf_masks
    if max_gene_trees > 0 and len(ei_list) > max_gene_trees:
        ei_list = ei_list[:max_gene_trees]
        sp_ids_list = sp_ids_list[:max_gene_trees]
        leaf_masks_list = leaf_masks_list[:max_gene_trees]

    # Compute 3 topology features
    topo_features = compute_topology_features(
        edge_clades=edge_clades,
        gene_tree_edge_indices=ei_list,
        gene_tree_species_ids=sp_ids_list,
        gene_tree_leaf_masks=leaf_masks_list,
        species_tree_is_leaf=is_leaf,
        species_tree_node_names=sample.species_tree_node_names,
        species_list=sample.species_list,
        n_edges=n_edges,
    )

    # Use original 4-dim edge features + 3 topology = 7-dim
    existing = sample.species_tree_edge_features
    if existing is not None:
        base = existing[:, :4]  # keep only original 4 features
    else:
        base = torch.zeros(n_edges, 4, dtype=torch.float)
    sample.species_tree_edge_features = torch.cat([base, topo_features], dim=-1)

    return sample


def main():
    parser = argparse.ArgumentParser(description="Compute Phase 3 gene tree topology features")
    parser.add_argument("--input-dir", required=True,
                        help="Input directory with Phase 1 samples")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for augmented samples")
    parser.add_argument("--max-gene-trees", type=int, default=0,
                        help="Max gene trees per sample (0 = all)")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"ERROR: Input directory not found: {args.input_dir}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    sample_dirs = sorted([
        d for d in os.listdir(args.input_dir)
        if os.path.isdir(os.path.join(args.input_dir, d))
    ])

    print(f"Input: {args.input_dir}")
    print(f"Output: {args.output_dir}")
    print(f"Samples: {len(sample_dirs)}")
    if args.max_gene_trees > 0:
        print(f"Max gene trees: {args.max_gene_trees}")
    print()

    processed = 0
    skipped = 0
    start_time = time.time()

    for i, dirname in enumerate(sample_dirs):
        input_path = os.path.join(args.input_dir, dirname)
        output_path = os.path.join(args.output_dir, dirname)

        try:
            sample = Gene2NetSample.load(input_path)
            sample = augment_sample(sample, max_gene_trees=args.max_gene_trees)
            sample.save(output_path)
            processed += 1
        except Exception as e:
            print(f"  SKIP {dirname}: {e}")
            skipped += 1
            continue

        if (i + 1) % 50 == 0:
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed
            remaining = (len(sample_dirs) - i - 1) / rate
            print(f"  Processed {i + 1}/{len(sample_dirs)} "
                  f"({rate:.1f} samples/s, ~{remaining:.0f}s remaining)")

    elapsed = time.time() - start_time
    print(f"\nDone: {processed} processed, {skipped} skipped ({elapsed:.0f}s)")
    print(f"Edge features: 4 -> 7 dims (4 original + 3 topology)")


if __name__ == "__main__":
    main()
