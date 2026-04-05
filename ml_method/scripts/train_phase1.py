"""Train Phase 1 model: features-only GNN for binary WGD detection.

Usage:
    python scripts/train_phase1.py --data-dir /path/to/training/ils_low
    python scripts/train_phase1.py --data-dir /path/to/training/ils_low --config configs/phase1.yaml
"""
import argparse
import os
import random
import sys

import torch
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2
from gene2net_gnn.training.trainer_phase1 import Phase1Trainer


def main():
    parser = argparse.ArgumentParser(description="Train Phase 1 GNN")
    parser.add_argument("--data-dir", required=True, nargs="+",
                        help="Training data directory (one or more)")
    parser.add_argument("--config", default=None, help="Config YAML (default: configs/phase1.yaml)")
    parser.add_argument("--output-dir", default=None, help="Output directory (default: output/phase1)")
    args = parser.parse_args()

    # Resolve paths
    base_dir = os.path.join(os.path.dirname(__file__), "..")
    config_path = args.config or os.path.join(base_dir, "configs", "phase1.yaml")
    output_dir = args.output_dir or os.path.join(base_dir, "output", "phase1")

    # Load config
    with open(config_path) as f:
        config = yaml.safe_load(f)

    model_config = config.get("model", {})
    train_config = config.get("training", {})

    # Device
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Print header
    print("=" * 70)
    print("Phase 1: Features-only GNN for binary WGD detection")
    print("=" * 70)
    print(f"Data dirs: {args.data_dir}")
    print(f"Config: {config_path}")
    print(f"Output: {output_dir}")
    print(f"Device: {device}")
    if device == "cuda":
        print(f"GPU: {torch.cuda.get_device_name()}")
    print("=" * 70)

    # Load dataset(s)
    all_samples = []
    skipped = 0
    for data_dir in args.data_dir:
        print(f"Loading dataset from {data_dir}...")
        dataset = Gene2NetDataset(data_dir)
        dir_loaded = 0
        for i in range(len(dataset)):
            try:
                sample = dataset[i]
                if sample.labels is not None:
                    all_samples.append(sample)
                    dir_loaded += 1
                else:
                    skipped += 1
            except Exception:
                skipped += 1
        print(f"  {dir_loaded} samples from {os.path.basename(data_dir)}")

    print(f"Total: {len(all_samples)} samples ({skipped} skipped)")

    # Train/val split
    val_split = float(train_config.get("val_split", 0.2))
    random.seed(42)
    random.shuffle(all_samples)
    n_val = int(len(all_samples) * val_split)
    val_samples = all_samples[:n_val]
    train_samples = all_samples[n_val:]
    print(f"Train: {len(train_samples)}, Val: {len(val_samples)}")

    # Count parameters
    model = SpeciesTreeGNNv2(
        node_feat_dim=int(model_config.get("node_feat_dim", 13)),
        edge_feat_dim=int(model_config.get("edge_feat_dim", 4)),
        hidden_dim=int(model_config.get("hidden_dim", 64)),
        n_gat_layers=int(model_config.get("n_gat_layers", 3)),
        n_gat_heads=int(model_config.get("n_gat_heads", 4)),
        dropout=float(model_config.get("dropout", 0.2)),
    )
    n_params = sum(p.numel() for p in model.parameters())
    print(f"Model: {n_params:,} parameters")

    # Train
    trainer = Phase1Trainer(model, train_config, device=device)
    trainer.train(train_samples, val_samples, output_dir)

    print(f"\nTraining complete! Model saved to {output_dir}")


if __name__ == "__main__":
    main()
