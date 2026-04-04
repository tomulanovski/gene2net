"""Train WGD Detector: GIN gene tree encoder + GAT species tree.

Usage:
    python scripts/train_wgd.py --data-dir /path/to/training/ils_low
    python scripts/train_wgd.py --data-dir /path/to/training/ils_low --config configs/wgd_detector.yaml
"""
import argparse
import os
import random
import sys

import torch
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.model.wgd_detector import WGDDetector
from gene2net_gnn.training.trainer_wgd import WGDTrainer


def main():
    parser = argparse.ArgumentParser(description="Train WGD Detector")
    parser.add_argument("--data-dir", required=True, help="Training data directory")
    parser.add_argument("--config", default=None, help="Config YAML")
    parser.add_argument("--output-dir", default=None, help="Output directory")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    config_path = args.config or os.path.join(base_dir, "configs", "wgd_detector.yaml")
    output_dir = args.output_dir or os.path.join(base_dir, "output", "wgd_detector")

    with open(config_path) as f:
        config = yaml.safe_load(f)

    model_config = config.get("model", {})
    train_config = config.get("training", {})

    device = "cuda" if torch.cuda.is_available() else "cpu"

    print("=" * 70)
    print("WGD Detector: GIN Gene Tree Encoder + GAT Species Tree")
    print("=" * 70)
    print(f"Data: {args.data_dir}")
    print(f"Config: {config_path}")
    print(f"Output: {output_dir}")
    print(f"Device: {device}")
    if device == "cuda":
        print(f"GPU: {torch.cuda.get_device_name()}")
    print("=" * 70)

    # Load dataset
    print(f"Loading dataset from {args.data_dir}...")
    dataset = Gene2NetDataset(args.data_dir)
    all_samples = []
    skipped = 0
    max_species = 0

    for i in range(len(dataset)):
        try:
            sample = dataset[i]
            if sample.labels is not None:
                all_samples.append(sample)
                max_species = max(max_species, sample.n_species)
            else:
                skipped += 1
        except Exception:
            skipped += 1

    print(f"Loaded {len(all_samples)} samples ({skipped} skipped)")
    print(f"Max species in dataset: {max_species}")

    # Train/val split
    val_split = float(train_config.get("val_split", 0.2))
    random.seed(42)
    random.shuffle(all_samples)
    n_val = int(len(all_samples) * val_split)
    val_samples = all_samples[:n_val]
    train_samples = all_samples[n_val:]
    print(f"Train: {len(train_samples)}, Val: {len(val_samples)}")

    # Use max_species from config or auto-detect
    n_species = int(model_config.get("n_species", max_species))
    print(f"Species embedding size: {n_species}")

    # Build model
    model = WGDDetector(
        n_species=n_species,
        node_feat_dim=int(model_config.get("node_feat_dim", 13)),
        edge_feat_dim=int(model_config.get("edge_feat_dim", 4)),
        hidden_dim=int(model_config.get("hidden_dim", 64)),
        n_gin_layers=int(model_config.get("n_gin_layers", 2)),
        n_gat_layers=int(model_config.get("n_gat_layers", 3)),
        n_gat_heads=int(model_config.get("n_gat_heads", 4)),
        dropout=float(model_config.get("dropout", 0.2)),
    )
    n_params = sum(p.numel() for p in model.parameters())
    print(f"Model: {n_params:,} parameters")

    # Count gene trees
    total_gt = sum(len(s.gene_tree_edge_indices) for s in all_samples)
    avg_gt = total_gt / len(all_samples)
    print(f"Gene trees: {total_gt:,} total, {avg_gt:.0f} avg per sample")

    # Train
    trainer = WGDTrainer(model, train_config, device=device)
    trainer.train(train_samples, val_samples, output_dir)

    print(f"\nTraining complete! Model saved to {output_dir}")


if __name__ == "__main__":
    main()
