"""Train the reconstruction model: joint WGD detection + partner-edge prediction.

Same data + species-tree GNN as Phase 1, plus a partner head. Partner prediction
subsumes auto/allo (self-partner = auto, other = allo), so its output can be fed
to build_mul_tree and scored against the existing methods with edit distance/Jaccard.

Usage:
    python scripts/train_reconstruct.py --data-dir <dir...> [--output-dir output/reconstruct]
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
from gene2net_gnn.training.trainer_reconstruct import (
    ReconstructTrainer,
    filter_compatible_state_dict,
)


def main():
    parser = argparse.ArgumentParser(description="Train reconstruction model (detection + partner)")
    parser.add_argument("--data-dir", required=True, nargs="+")
    parser.add_argument("--config", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--init-from", default=None,
                        help="Optional Phase 1 checkpoint to warm-start the shared backbone")
    parser.add_argument("--clade-labels", action="store_true",
                        help="Train on corrected clade-level labels (labels_clade.pkl)")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    config_path = args.config or os.path.join(base_dir, "configs", "reconstruct.yaml")
    output_dir = args.output_dir or os.path.join(base_dir, "output", "reconstruct")

    with open(config_path) as f:
        config = yaml.safe_load(f)
    model_config = config.get("model", {})
    train_config = config.get("training", {})

    device = "cuda" if torch.cuda.is_available() else "cpu"
    expected_edge_dim = int(model_config.get("edge_feat_dim", 9))

    print("=" * 70)
    print("Reconstruction: WGD detection + partner-edge prediction")
    print("=" * 70)
    print(f"Data dirs: {args.data_dir}")
    print(f"Config: {config_path}\nOutput: {output_dir}\nDevice: {device}")
    print("=" * 70)

    all_samples = []
    skipped = wrong_dim = 0
    for data_dir in args.data_dir:
        print(f"Loading dataset from {data_dir}...")
        dataset = Gene2NetDataset(data_dir, clade_labels=args.clade_labels)
        dir_loaded = dir_wrong = 0
        for i in range(len(dataset)):
            try:
                sample = dataset[i]
                if sample.labels is None:
                    skipped += 1
                    continue
                ef = sample.species_tree_edge_features
                if ef is None or ef.shape[1] != expected_edge_dim:
                    wrong_dim += 1
                    dir_wrong += 1
                    continue
                all_samples.append(sample)
                dir_loaded += 1
            except Exception:
                skipped += 1
        print(f"  {dir_loaded} samples from {os.path.basename(data_dir)}"
              + (f"  ({dir_wrong} wrong edge dim, skipped)" if dir_wrong else ""))

    print(f"Total: {len(all_samples)} samples ({skipped} no-label, {wrong_dim} wrong-dim skipped)")

    val_split = float(train_config.get("val_split", 0.2))
    random.seed(42)
    random.shuffle(all_samples)
    n_val = int(len(all_samples) * val_split)
    val_samples = all_samples[:n_val]
    train_samples = all_samples[n_val:]
    print(f"Train: {len(train_samples)}, Val: {len(val_samples)}")

    model = SpeciesTreeGNNv2(
        node_feat_dim=int(model_config.get("node_feat_dim", 13)),
        edge_feat_dim=expected_edge_dim,
        hidden_dim=int(model_config.get("hidden_dim", 64)),
        n_gat_layers=int(model_config.get("n_gat_layers", 3)),
        n_gat_heads=int(model_config.get("n_gat_heads", 4)),
        dropout=float(model_config.get("dropout", 0.2)),
        partner_pair_feat_dim=int(model_config.get("partner_pair_feat_dim", 2)),
        n_parents=int(model_config.get("n_parents", 2)),   # two-parent phaser head
    )
    if args.init_from and os.path.exists(args.init_from):
        # Warm-start shared weights from a prior model. Shape-incompatible layers
        # (e.g. the partner head when its pairwise feature width changed) are
        # dropped so they reinitialize fresh instead of erroring on load.
        state = torch.load(args.init_from, map_location="cpu", weights_only=True)
        state = filter_compatible_state_dict(state, model.state_dict())
        missing, unexpected = model.load_state_dict(state, strict=False)
        print(f"Warm-started from {args.init_from} (missing={len(missing)}, unexpected={len(unexpected)})")

    n_params = sum(p.numel() for p in model.parameters())
    print(f"Model: {n_params:,} parameters")

    trainer = ReconstructTrainer(model, train_config, device=device)
    trainer.train(train_samples, val_samples, output_dir)
    print(f"\nTraining complete! Model saved to {output_dir}")


if __name__ == "__main__":
    main()
