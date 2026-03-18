#!/usr/bin/env python
"""Training entry point for Gene2Net-GNN."""
import argparse
import os
import sys
import yaml
import random
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.model.gene2net_model import Gene2NetModel
from gene2net_gnn.training.trainer import Trainer


def main():
    parser = argparse.ArgumentParser(description="Train Gene2Net-GNN model")
    parser.add_argument("--data-dir", required=True, help="Directory with training examples")
    parser.add_argument("--config", default="configs/default.yaml", help="Config YAML file")
    parser.add_argument("--output-dir", default="output", help="Output directory for model checkpoints")
    parser.add_argument("--device", default="auto", help="Device (cpu, cuda, auto)")
    parser.add_argument("--val-split", type=float, default=0.1, help="Validation split fraction")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Device
    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    print(f"Using device: {device}")

    # Seed
    random.seed(args.seed)
    torch.manual_seed(args.seed)

    # Load dataset
    print(f"Loading dataset from {args.data_dir}...")
    dataset = Gene2NetDataset(args.data_dir)
    print(f"Found {len(dataset)} examples")

    # Split
    n_val = max(1, int(len(dataset) * args.val_split))
    indices = list(range(len(dataset)))
    random.shuffle(indices)
    val_indices = indices[:n_val]
    train_indices = indices[n_val:]

    train_samples = [dataset[i] for i in train_indices]
    val_samples = [dataset[i] for i in val_indices]
    print(f"Train: {len(train_samples)}, Val: {len(val_samples)}")

    # Determine dimensions from data and config
    data_config = config.get("data", {})
    max_species = data_config.get("max_species", 80)
    sample = train_samples[0]
    handcrafted_dim = sample.species_tree_node_features.shape[1]

    # Create model — use max_species for embedding table (supports varying species counts)
    model_config = config.get("model", {})
    model = Gene2NetModel(
        n_species=max_species,
        embed_dim=model_config.get("embed_dim", 64),
        n_attention_heads=model_config.get("num_attention_heads", 8),
        n_gat_layers=model_config.get("num_gat_layers", 3),
        handcrafted_dim=handcrafted_dim,
        conv_type=model_config.get("conv_type", "gat"),
        dropout=model_config.get("dropout", 0.2),
    )
    print(f"Model: {sum(p.numel() for p in model.parameters())} parameters")

    # Train
    training_config = config.get("training", {})
    trainer = Trainer(model, training_config, device=device)
    trainer.train(train_samples, val_samples, args.output_dir)

    print(f"\nTraining complete! Model saved to {args.output_dir}/best_model.pt")


if __name__ == "__main__":
    main()
