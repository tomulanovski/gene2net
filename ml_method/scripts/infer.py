#!/usr/bin/env python
"""Inference CLI: gene trees → MUL-tree.

Usage:
    python infer.py --gene-trees input.nwk --model model.pt --output mul_tree.nwk
    python infer.py --gene-trees input.nwk --species-tree astral.nwk --model model.pt --output mul_tree.nwk
"""
import argparse
import os
import sys
import yaml
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree
from gene2net_gnn.data.tree_io import load_gene_trees_from_file, get_species_set
from gene2net_gnn.model.gene2net_model import Gene2NetModel
from gene2net_gnn.inference.predict import predict_mul_tree


def main():
    parser = argparse.ArgumentParser(description="Gene2Net-GNN Inference")
    parser.add_argument("--gene-trees", required=True, help="Path to gene trees (Newick)")
    parser.add_argument("--species-tree", help="Path to species tree (if not provided, uses ASTRAL-Pro)")
    parser.add_argument("--model", required=True, help="Path to trained model (.pt)")
    parser.add_argument("--config", default="configs/default.yaml", help="Config file")
    parser.add_argument("--output", required=True, help="Output MUL-tree path")
    parser.add_argument("--threshold", type=float, default=0.5, help="WGD confidence threshold")
    parser.add_argument("--max-events", type=int, default=None, help="Maximum number of WGD events")
    parser.add_argument("--device", default="cpu", help="Device (cpu or cuda)")
    args = parser.parse_args()

    device = torch.device(args.device)

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Load gene trees
    print(f"Loading gene trees from {args.gene_trees}...")
    gene_trees = load_gene_trees_from_file(args.gene_trees)
    print(f"Loaded {len(gene_trees)} gene trees")

    species = get_species_set(gene_trees)
    print(f"Species: {len(species)}")

    # Load or infer species tree
    if args.species_tree:
        print(f"Loading species tree from {args.species_tree}...")
        species_tree = Tree(open(args.species_tree).read(), format=1)
    else:
        print("No species tree provided. Please provide one with --species-tree.")
        print("(ASTRAL-Pro integration for automatic inference is planned for cluster deployment)")
        sys.exit(1)

    # Load model
    model_config = config.get("model", {})
    data_config = config.get("data", {})
    max_species = data_config.get("max_species", 80)

    # handcrafted_dim = 8 copy count features + 5 clustering summary stats = 13 (fixed)
    handcrafted_dim = 13

    model = Gene2NetModel(
        n_species=max_species,
        embed_dim=model_config.get("embed_dim", 64),
        n_attention_heads=model_config.get("num_attention_heads", 8),
        n_gat_layers=model_config.get("num_gat_layers", 3),
        handcrafted_dim=handcrafted_dim,
        conv_type=model_config.get("conv_type", "gat"),
        dropout=0.0,  # no dropout at inference
    )
    model.load_state_dict(torch.load(args.model, map_location=device, weights_only=True))
    model.to(device)
    print("Model loaded")

    # Run inference
    print("Running inference...")
    mul_tree, events, raw = predict_mul_tree(
        model, gene_trees, species_tree,
        device=device,
        threshold=args.threshold,
        max_events=args.max_events,
    )

    # Output
    with open(args.output, "w") as f:
        f.write(mul_tree.write(format=5) + "\n")

    print(f"\nPredicted {len(events)} WGD events:")
    for e in sorted(events, key=lambda x: -x.confidence):
        clade = ",".join(sorted(e.wgd_edge_clade))
        partner = ",".join(sorted(e.partner_edge_clade))
        etype = "auto" if e.wgd_edge_clade == e.partner_edge_clade else "allo"
        print(f"  {etype}: clade={{{clade}}} → partner={{{partner}}} (conf={e.confidence:.3f})")

    print(f"\nMUL-tree written to {args.output}")


if __name__ == "__main__":
    main()
