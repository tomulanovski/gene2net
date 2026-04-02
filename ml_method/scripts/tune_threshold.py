"""Tune decision threshold on validation set using saved Phase 1 model.

Instead of argmax (implicit threshold=0.5), find the threshold on the WGD
class probability that maximizes F1.

Usage:
    python scripts/tune_threshold.py --data-dir /path/to/training/ils_low --model-dir output/phase1
"""
import argparse
import os
import random
import sys

import numpy as np
import torch
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2
from gene2net_gnn.training.trainer_phase1 import prepare_sample


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--model-dir", default=None)
    parser.add_argument("--config", default=None)
    parser.add_argument("--max-samples", type=int, default=400, help="Max val samples to use")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    model_dir = args.model_dir or os.path.join(base_dir, "output", "phase1")
    config_path = args.config or os.path.join(base_dir, "configs", "phase1.yaml")

    with open(config_path) as f:
        config = yaml.safe_load(f)
    model_config = config.get("model", {})

    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Load model
    model = SpeciesTreeGNNv2(
        node_feat_dim=int(model_config.get("node_feat_dim", 13)),
        edge_feat_dim=int(model_config.get("edge_feat_dim", 4)),
        hidden_dim=int(model_config.get("hidden_dim", 64)),
        n_gat_layers=int(model_config.get("n_gat_layers", 3)),
        n_gat_heads=int(model_config.get("n_gat_heads", 4)),
        dropout=float(model_config.get("dropout", 0.2)),
    )

    # Try best_f1_model first, then best_model
    for model_name in ["best_f1_model.pt", "best_model.pt"]:
        model_path = os.path.join(model_dir, model_name)
        if os.path.exists(model_path):
            model.load_state_dict(torch.load(model_path, map_location=device, weights_only=True))
            print(f"Loaded {model_name}")
            break
    else:
        print("ERROR: No model found")
        return

    model = model.to(device)
    model.eval()

    # Load dataset indices, use same split as training (without loading all into memory)
    dataset = Gene2NetDataset(args.data_dir)
    n_total = len(dataset)

    # Reproduce the same train/val split
    indices = list(range(n_total))
    random.seed(42)
    random.shuffle(indices)
    n_val = int(n_total * 0.2)
    val_indices = indices[:min(n_val, args.max_samples)]
    print(f"Val samples: {len(val_indices)}")

    # Collect all predictions and targets (one at a time to save memory)
    all_probs = []
    all_targets = []

    with torch.no_grad():
        for i, idx in enumerate(val_indices):
            try:
                sample = dataset[idx]
            except Exception:
                continue
            if sample.labels is None:
                continue
            prepared = prepare_sample(sample, torch.device(device))
            if prepared is None:
                continue
            inputs, targets, mask = prepared
            wgd_logits, _ = model(**inputs)

            probs = torch.softmax(wgd_logits, dim=-1)[:, 1]  # P(WGD)
            all_probs.append(probs[mask].cpu().numpy())
            all_targets.append(targets[mask].cpu().numpy())

            # Free memory
            del sample, prepared, inputs, targets, mask, wgd_logits, probs

            if (i + 1) % 100 == 0:
                print(f"  Processed {i+1}/{len(val_indices)}")

    probs = np.concatenate(all_probs)
    targets = np.concatenate(all_targets)

    print(f"\nTotal edges: {len(targets)}")
    print(f"Positive: {targets.sum()} ({100*targets.mean():.1f}%)")

    # Sweep thresholds
    print(f"\n{'Threshold':>10} | {'Prec':>6} | {'Rec':>6} | {'F1':>6} | {'TP':>6} | {'FP':>6} | {'FN':>6}")
    print("-" * 65)

    best_f1 = 0
    best_thresh = 0.5

    for thresh in np.arange(0.05, 0.95, 0.05):
        preds = (probs >= thresh).astype(int)
        tp = ((preds == 1) & (targets == 1)).sum()
        fp = ((preds == 1) & (targets == 0)).sum()
        fn = ((preds == 0) & (targets == 1)).sum()

        prec = tp / max(tp + fp, 1)
        rec = tp / max(tp + fn, 1)
        f1 = 2 * prec * rec / max(prec + rec, 1e-8)

        marker = " <-- best" if f1 > best_f1 else ""
        print(f"{thresh:>10.2f} | {prec:>6.3f} | {rec:>6.3f} | {f1:>6.3f} | {tp:>6} | {fp:>6} | {fn:>6}{marker}")

        if f1 > best_f1:
            best_f1 = f1
            best_thresh = thresh

    print(f"\nBest threshold: {best_thresh:.2f} → F1={best_f1:.3f}")

    # Save results
    out_path = os.path.join(model_dir, "threshold_tuning.txt")
    with open(out_path, "w") as f:
        f.write(f"Best threshold: {best_thresh}\n")
        f.write(f"Best F1: {best_f1}\n")
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()
