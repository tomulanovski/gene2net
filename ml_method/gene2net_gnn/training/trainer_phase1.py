"""Phase 1 trainer: binary WGD detection with features-only GNN.

Simplified training loop:
  - Single task: binary WGD per edge (no partner, no count)
  - Focal loss with class weights for imbalance
  - Tracks precision, recall, F1 per epoch
  - Gradient accumulation for effective batch size
"""
import os
import random
from typing import List

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from torch.optim.lr_scheduler import ReduceLROnPlateau

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2


def focal_loss(logits, targets, alpha=0.25, gamma=2.0, class_weights=None):
    """Focal loss for binary classification with class imbalance."""
    if class_weights is not None:
        ce = F.cross_entropy(logits, targets, weight=class_weights, reduction="none")
    else:
        ce = F.cross_entropy(logits, targets, reduction="none")
    pt = torch.exp(-ce)
    return (alpha * (1 - pt) ** gamma * ce).mean()


def prepare_sample(sample: Gene2NetSample, device: torch.device):
    """Convert a sample into model inputs and binary WGD targets.

    Returns None if the sample has no valid labels.
    """
    labels = sample.labels
    if labels is None:
        return None

    n_edges = sample.species_tree_edge_index.shape[1] // 2

    # Binary targets: 0 = no WGD, 1 = WGD
    wgd_targets = torch.zeros(n_edges, dtype=torch.long, device=device)
    mask = torch.ones(n_edges, dtype=torch.bool, device=device)

    for edge_idx in range(min(len(labels.wgd_counts), n_edges)):
        if labels.wgd_counts[edge_idx] > 0:
            wgd_targets[edge_idx] = 1

    # Apply unmappable mask
    if labels.mask is not None:
        for i in range(len(labels.wgd_edges)):
            edge_idx = labels.wgd_edges[i]
            if edge_idx < n_edges and i < len(labels.mask) and not labels.mask[i]:
                mask[edge_idx] = False

    model_inputs = {
        "node_features": sample.species_tree_node_features.to(device),
        "edge_index": sample.species_tree_edge_index.to(device),
        "edge_features": sample.species_tree_edge_features.to(device),
        "is_leaf": sample.species_tree_is_leaf.to(device),
    }

    return model_inputs, wgd_targets, mask


class Phase1Trainer:
    """Training loop for Phase 1 binary WGD detection."""

    def __init__(self, model: SpeciesTreeGNNv2, config: dict, device: str = "cpu"):
        self.model = model.to(device)
        self.device = torch.device(device)
        self.config = config

        # Class weights for imbalance (12:1 ratio → weight positive class higher)
        pos_weight = float(config.get("pos_class_weight", 12.0))
        self.class_weights = torch.tensor([1.0, pos_weight], dtype=torch.float, device=self.device)

        self.focal_alpha = float(config.get("focal_alpha", 0.25))
        self.focal_gamma = float(config.get("focal_gamma", 2.0))

        self.optimizer = Adam(
            model.parameters(),
            lr=float(config.get("lr", 3e-4)),
            weight_decay=float(config.get("weight_decay", 1e-4)),
        )
        self.scheduler = ReduceLROnPlateau(
            self.optimizer, mode="min", factor=0.5, patience=10,
        )

        self.batch_size = max(1, int(config.get("batch_size", 8)))
        self.patience = int(config.get("patience", 30))
        self.max_epochs = int(config.get("max_epochs", 200))

    def train_epoch(self, train_samples: List[Gene2NetSample]) -> float:
        self.model.train()
        total_loss = 0.0
        n_samples = 0
        random.shuffle(train_samples)
        self.optimizer.zero_grad()

        for i, sample in enumerate(train_samples):
            prepared = prepare_sample(sample, self.device)
            if prepared is None:
                continue

            inputs, targets, mask = prepared
            wgd_logits, _ = self.model(**inputs)

            if mask.sum() == 0:
                continue

            loss = focal_loss(
                wgd_logits[mask], targets[mask],
                alpha=self.focal_alpha, gamma=self.focal_gamma,
                class_weights=self.class_weights,
            )

            (loss / self.batch_size).backward()
            total_loss += loss.item()
            n_samples += 1

            if (i + 1) % self.batch_size == 0 or (i + 1) == len(train_samples):
                nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                self.optimizer.step()
                self.optimizer.zero_grad()

        return total_loss / max(n_samples, 1)

    @torch.no_grad()
    def evaluate(self, val_samples: List[Gene2NetSample]) -> dict:
        self.model.eval()
        total_loss = 0.0
        total_tp = 0
        total_fp = 0
        total_fn = 0
        total_tn = 0
        n_samples = 0

        for sample in val_samples:
            prepared = prepare_sample(sample, self.device)
            if prepared is None:
                continue

            inputs, targets, mask = prepared
            wgd_logits, _ = self.model(**inputs)

            if mask.sum() == 0:
                continue

            loss = focal_loss(
                wgd_logits[mask], targets[mask],
                alpha=self.focal_alpha, gamma=self.focal_gamma,
                class_weights=self.class_weights,
            )
            total_loss += loss.item()
            n_samples += 1

            # Binary metrics
            preds = wgd_logits[mask].argmax(dim=-1)
            true = targets[mask]
            total_tp += (preds & true).sum().item()
            total_fp += (preds & ~true.bool()).sum().item()
            total_fn += (~preds.bool() & true.bool()).sum().item()
            total_tn += (~preds.bool() & ~true.bool()).sum().item()

        precision = total_tp / max(total_tp + total_fp, 1)
        recall = total_tp / max(total_tp + total_fn, 1)
        f1 = 2 * precision * recall / max(precision + recall, 1e-8)

        return {
            "val_loss": total_loss / max(n_samples, 1),
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "tp": total_tp,
            "fp": total_fp,
            "fn": total_fn,
            "tn": total_tn,
        }

    def train(self, train_samples, val_samples, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        best_val_loss = float("inf")
        best_f1 = 0.0
        patience_counter = 0

        for epoch in range(self.max_epochs):
            train_loss = self.train_epoch(train_samples)
            val = self.evaluate(val_samples)
            self.scheduler.step(val["val_loss"])

            lr = self.optimizer.param_groups[0]["lr"]
            print(
                f"Epoch {epoch+1}/{self.max_epochs} | "
                f"Train Loss: {train_loss:.4f} | "
                f"Val Loss: {val['val_loss']:.4f} | "
                f"F1: {val['f1']:.3f} | "
                f"Prec: {val['precision']:.3f} | "
                f"Rec: {val['recall']:.3f} | "
                f"TP:{val['tp']} FP:{val['fp']} FN:{val['fn']} | "
                f"LR: {lr:.6f}"
            )

            # Save best model by val loss
            if val["val_loss"] < best_val_loss:
                best_val_loss = val["val_loss"]
                patience_counter = 0
                torch.save(self.model.state_dict(), os.path.join(output_dir, "best_model.pt"))
            else:
                patience_counter += 1

            # Also track best F1 (saved separately)
            if val["f1"] > best_f1:
                best_f1 = val["f1"]
                torch.save(self.model.state_dict(), os.path.join(output_dir, "best_f1_model.pt"))

            if patience_counter >= self.patience:
                print(f"Early stopping at epoch {epoch+1}")
                break

        print(f"\nBest val loss: {best_val_loss:.4f}")
        print(f"Best F1: {best_f1:.3f}")

        # Load best model (by val loss)
        self.model.load_state_dict(
            torch.load(os.path.join(output_dir, "best_model.pt"), weights_only=True)
        )
        return self.model
