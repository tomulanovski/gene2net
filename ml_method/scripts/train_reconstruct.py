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


def select_train_subset(train_samples, fraction=None, max_samples=None):
    """Slice the already-shuffled training pool for the learning curve. `fraction`
    (0-1) takes precedence over `max_samples`; both return a nested prefix."""
    if fraction is not None:
        return train_samples[:int(len(train_samples) * fraction)]
    if max_samples is not None:
        return train_samples[:max_samples]
    return train_samples


def run_training(data_dirs, model_config, train_config, output_dir, *,
                 clade_labels=False, away_labels=False, init_from=None,
                 train_fraction=None, max_train_samples=None, device=None):
    """Load data, build the model, and train. Returns the trainer result dict
    (best_val_loss, best_partner_acc, best_allo_acc, best_f1, history, model).

    No silent error handling: sample-load errors propagate, and a specified
    init_from that does not exist raises. Only the documented no-label and
    wrong-edge-dim filters drop samples, and both are reported.
    """
    train_config = dict(train_config)
    train_config["use_n_eff"] = model_config.get("use_n_eff", False)
    train_config["coclust_condition_on_dup"] = model_config.get("coclust_condition_on_dup", False)
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    expected_edge_dim = int(model_config.get("edge_feat_dim", 9))

    all_samples = []
    no_label = wrong_dim = 0
    for data_dir in data_dirs:
        print(f"Loading dataset from {data_dir}...")
        dataset = Gene2NetDataset(data_dir, clade_labels=clade_labels, away_labels=away_labels)
        dir_loaded = 0
        for i in range(len(dataset)):
            sample = dataset[i]
            if sample.labels is None:
                no_label += 1
                continue
            ef = sample.species_tree_edge_features
            if ef is None or ef.shape[1] != expected_edge_dim:
                wrong_dim += 1
                continue
            all_samples.append(sample)
            dir_loaded += 1
        print(f"  {dir_loaded} samples from {os.path.basename(data_dir)}")

    print(f"Total: {len(all_samples)} samples "
          f"({no_label} no-label, {wrong_dim} wrong-edge-dim filtered)")
    if wrong_dim:
        print(f"WARNING: {wrong_dim} samples had edge dim != {expected_edge_dim}; "
              f"check config edge_feat_dim against the data.")

    val_split = float(train_config.get("val_split", 0.2))
    random.seed(42)
    random.shuffle(all_samples)
    n_val = int(len(all_samples) * val_split)
    val_samples = all_samples[:n_val]
    train_samples = select_train_subset(all_samples[n_val:],
                                        fraction=train_fraction,
                                        max_samples=max_train_samples)
    print(f"Train: {len(train_samples)}, Val: {len(val_samples)}"
          + (f" (fraction={train_fraction})" if train_fraction is not None else "")
          + (f" (max={max_train_samples})" if max_train_samples is not None else ""))

    model = SpeciesTreeGNNv2(
        node_feat_dim=int(model_config.get("node_feat_dim", 13)),
        edge_feat_dim=expected_edge_dim,
        hidden_dim=int(model_config.get("hidden_dim", 64)),
        n_gat_layers=int(model_config.get("n_gat_layers", 3)),
        n_gat_heads=int(model_config.get("n_gat_heads", 4)),
        dropout=float(model_config.get("dropout", 0.2)),
        partner_pair_feat_dim=int(model_config.get("partner_pair_feat_dim", 2)),
        n_parents=int(model_config.get("n_parents", 2)),
        conv_type=str(model_config.get("conv_type", "gat")),
    )
    if init_from is not None:
        if not os.path.exists(init_from):
            raise FileNotFoundError(f"--init-from path does not exist: {init_from}")
        # Warm-start shared weights; shape-incompatible layers are dropped so they
        # reinitialize fresh instead of erroring on load.
        state = torch.load(init_from, map_location="cpu", weights_only=True)
        state = filter_compatible_state_dict(state, model.state_dict())
        missing, unexpected = model.load_state_dict(state, strict=False)
        print(f"Warm-started from {init_from} (missing={len(missing)}, unexpected={len(unexpected)})")

    n_params = sum(p.numel() for p in model.parameters())
    print(f"Model: {n_params:,} parameters | conv={model_config.get('conv_type', 'gat')} "
          f"hidden={model_config.get('hidden_dim', 64)} layers={model_config.get('n_gat_layers', 3)}")

    trainer = ReconstructTrainer(model, train_config, device=device)
    return trainer.train(train_samples, val_samples, output_dir)


def main():
    parser = argparse.ArgumentParser(description="Train reconstruction model (detection + partner)")
    parser.add_argument("--data-dir", required=True, nargs="+")
    parser.add_argument("--config", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--init-from", default=None,
                        help="Optional Phase 1 checkpoint to warm-start the shared backbone")
    parser.add_argument("--clade-labels", action="store_true",
                        help="Train on corrected clade-level labels (labels_clade.pkl)")
    parser.add_argument("--away-labels", action="store_true",
                        help="Train on away-parent labels (labels_away.pkl): each allo partner "
                             "retargeted to the non-home parent. Takes precedence over --clade-labels.")
    parser.add_argument("--max-train-samples", type=int, default=None,
                        help="Cap the TRAINING set to this many samples (validation is untouched). "
                             "The seed-42 shuffle makes the subsets nested.")
    parser.add_argument("--train-fraction", type=float, default=None,
                        help="Train on this fraction (0-1) of the shuffled training pool. "
                             "Takes precedence over --max-train-samples. For the learning curve.")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    config_path = args.config or os.path.join(base_dir, "configs", "reconstruct_final.yaml")
    output_dir = args.output_dir or os.path.join(base_dir, "output", "reconstruct")

    with open(config_path) as f:
        config = yaml.safe_load(f)
    model_config = config.get("model", {})
    train_config = config.get("training", {})

    print("=" * 70)
    print("Reconstruction: WGD detection + partner-edge prediction")
    print("=" * 70)
    print(f"Data dirs: {args.data_dir}")
    print(f"Config: {config_path}\nOutput: {output_dir}")
    print("=" * 70)

    run_training(
        args.data_dir, model_config, train_config, output_dir,
        clade_labels=args.clade_labels, away_labels=args.away_labels,
        init_from=args.init_from, train_fraction=args.train_fraction,
        max_train_samples=args.max_train_samples,
    )
    print(f"\nTraining complete! Model saved to {output_dir}")


if __name__ == "__main__":
    main()
