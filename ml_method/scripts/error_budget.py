"""Per-config error budget for the reconstruction model.

Runs the trained model's own evaluation (same metrics as training) on each config
separately, so you can see WHERE each error type dominates:

  precision  low  -> detection over-predicts (WGD on non-polyploids) in this config
  recall     low  -> misses true WGD
  allo_acc        -> partner (parent-2) prediction — the persistent weak link
  auto_acc        -> autopolyploidy partner (usually easy)

This is E3 (detection/over-prediction) + E4 (partner) of the error budget. Pair it with
the structural diagnostics for the full picture:
  E1 diploid backbone : backbone_polyploid_localization.py   (per config)
  E2 polyploid placement (the 81%) : backbone_polyploid_placement_accuracy.py (per config)
  E5 build floor : oracle_test.py --events labels + score    (per config)

Run in final_project (torch + ete3).
Usage:
  python scripts/error_budget.py --model-dir output/reconstruct_cladelabels_rooted \
      --model-config configs/reconstruct_base.yaml --data-root data/mul_trees_2k \
      --configs ils_low ils_medium ils_high \
                dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M
"""
import argparse
import os
import sys

import torch
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2
from gene2net_gnn.training.trainer_reconstruct import ReconstructTrainer


def load_model(model_dir, model_config, device):
    with open(model_config) as f:
        cfg = yaml.safe_load(f)
    mc = cfg.get("model", {})
    model = SpeciesTreeGNNv2(
        node_feat_dim=int(mc.get("node_feat_dim", 13)),
        edge_feat_dim=int(mc.get("edge_feat_dim", 9)),
        hidden_dim=int(mc.get("hidden_dim", 64)),
        n_gat_layers=int(mc.get("n_gat_layers", 3)),
        n_gat_heads=int(mc.get("n_gat_heads", 4)),
        dropout=float(mc.get("dropout", 0.2)),
        partner_pair_feat_dim=int(mc.get("partner_pair_feat_dim", 2)),
    )
    ckpt = None
    for name in ("best_model.pt", "best_partner_model.pt", "model.pt"):
        p = os.path.join(model_dir, name)
        if os.path.exists(p):
            ckpt = p
            break
    if ckpt is None:
        raise FileNotFoundError(f"No checkpoint (best_model.pt/...) in {model_dir}")
    state = torch.load(ckpt, map_location="cpu", weights_only=True)
    model.load_state_dict(state)
    model.to(device).eval()
    return model, cfg.get("training", {}), ckpt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", required=True)
    ap.add_argument("--model-config", required=True)
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--subdir", default="training_rooted")
    ap.add_argument("--configs", nargs="+", required=True)
    ap.add_argument("--max-samples", type=int, default=300)
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model, train_cfg, ckpt = load_model(args.model_dir, args.model_config, device)
    trainer = ReconstructTrainer(model, train_cfg, device=device)
    print(f"model: {ckpt}   device: {device}\n")

    hdr = f"{'config':<26}{'F1':>7}{'prec':>7}{'rec':>7}{'allo':>7}{'auto':>7}{'n_allo':>8}"
    print(hdr)
    print("-" * len(hdr))
    for cfg in args.configs:
        d = os.path.join(args.data_root, args.subdir, cfg)
        if not os.path.isdir(d):
            print(f"{cfg:<26}  (missing {d})")
            continue
        ds = Gene2NetDataset(d, clade_labels=True)
        n = min(len(ds), args.max_samples)

        def stream():
            # yield one sample at a time so we never hold the whole config in memory
            for i in range(n):
                try:
                    s = ds[i]
                except Exception:
                    continue
                if s.labels is not None:
                    yield s

        with torch.no_grad():
            m = trainer.evaluate(stream())
        if m.get("allo_total", 0) + m.get("auto_total", 0) == 0:
            print(f"{cfg:<26}  (no labelled samples)")
            continue
        print(f"{cfg:<26}{m['f1']:>7.3f}{m['precision']:>7.3f}{m['recall']:>7.3f}"
              f"{m['allo_acc']:>7.3f}{m['auto_acc']:>7.3f}{m['allo_total']:>8}")

    print("\nReading:")
    print("  low precision  -> detection over-predicts (fires WGD on non-polyploids)")
    print("  low allo_acc   -> partner / parent-2 weakness (the persistent gap)")
    print("Compare across configs: where does each error type worsen? Focus effort there.")


if __name__ == "__main__":
    main()
