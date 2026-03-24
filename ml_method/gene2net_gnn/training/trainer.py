"""Training loop for Gene2Net."""
import os
import torch
import torch.nn as nn
from torch.optim import Adam
from torch.optim.lr_scheduler import ReduceLROnPlateau
import random
import yaml
from typing import List, Optional

from gene2net_gnn.model.gene2net_model import Gene2NetModel
from gene2net_gnn.training.loss import Gene2NetLoss
from gene2net_gnn.training.metrics import wgd_f1, partner_accuracy, event_count_error
from gene2net_gnn.data.dataset import Gene2NetSample


def prepare_batch(sample: Gene2NetSample, device: torch.device, max_gene_trees: int = 300):
    """Convert a Gene2NetSample into model inputs and loss targets.

    Subsamples gene trees if there are more than max_gene_trees.
    """
    # Subsample gene trees
    n_gt = len(sample.gene_tree_edge_indices)
    if n_gt > max_gene_trees:
        indices = random.sample(range(n_gt), max_gene_trees)
    else:
        indices = list(range(n_gt))

    gene_tree_data = []
    for i in indices:
        gene_tree_data.append({
            "edge_index": sample.gene_tree_edge_indices[i].to(device),
            "species_ids": sample.gene_tree_species_ids[i].to(device),
            "branch_lengths": sample.gene_tree_branch_lengths[i].to(device),
            "leaf_mask": sample.gene_tree_leaf_masks[i].to(device),
        })

    # Build species_ids for species tree nodes
    sp_to_idx = {sp: i for i, sp in enumerate(sample.species_list)}
    species_ids = []
    for name in sample.species_tree_node_names:
        species_ids.append(sp_to_idx.get(name, -1))

    species_tree_data = {
        "edge_index": sample.species_tree_edge_index.to(device),
        "node_features": sample.species_tree_node_features.to(device),
        "is_leaf": sample.species_tree_is_leaf.to(device),
        "species_ids": torch.tensor(species_ids, dtype=torch.long, device=device),
    }
    if sample.species_tree_edge_features is not None:
        species_tree_data["edge_features"] = sample.species_tree_edge_features.to(device)

    # Build targets from labels
    labels = sample.labels
    n_edges = sample.species_tree_edge_index.shape[1] // 2

    wgd_targets = torch.zeros(n_edges, dtype=torch.long, device=device)
    mask = torch.ones(n_edges, dtype=torch.bool, device=device)
    # Per-edge partner target: for edges with WGD, which edge is the partner
    partner_per_edge = torch.zeros(n_edges, dtype=torch.long, device=device)

    if labels is not None:
        # wgd_counts is always per-edge (length n_edges) from map_events_to_astral
        for edge_idx in range(min(len(labels.wgd_counts), n_edges)):
            count = labels.wgd_counts[edge_idx]
            if count > 0:
                wgd_targets[edge_idx] = min(count, 3)  # clamp to max_events

        # Map events to per-edge partner targets and mask
        # labels.wgd_edges, partner_edges, mask are parallel per-event arrays
        for i in range(len(labels.wgd_edges)):
            edge_idx = labels.wgd_edges[i]
            if edge_idx >= n_edges:
                continue
            # Set partner for this edge (last event wins if multiple)
            if i < len(labels.partner_edges):
                partner_per_edge[edge_idx] = min(labels.partner_edges[i], n_edges - 1)
            # Set mask: if any event on this edge is unmappable, mask it out
            if hasattr(labels, 'mask') and labels.mask is not None and i < len(labels.mask):
                if not labels.mask[i]:
                    mask[edge_idx] = False

    wgd_mask = wgd_targets > 0
    # partner_targets aligned with wgd_mask: one target per WGD edge
    partner_targets = partner_per_edge[wgd_mask]
    true_count = torch.tensor([float(sum(labels.wgd_counts))] if labels else [0.0], device=device)

    return gene_tree_data, species_tree_data, wgd_targets, mask, wgd_mask, partner_targets, true_count


class Trainer:
    """Training loop with validation, early stopping, and checkpointing."""

    def __init__(self, model, config, device="cpu"):
        self.model = model.to(device)
        self.device = device
        self.config = config

        self.loss_fn = Gene2NetLoss(
            lambda_wgd=float(config.get("lambda_wgd", 1.0)),
            lambda_partner=float(config.get("lambda_partner", 1.0)),
            lambda_count=float(config.get("lambda_count", 0.1)),
            focal_alpha=float(config.get("focal_alpha", 0.25)),
            focal_gamma=float(config.get("focal_gamma", 2.0)),
        )

        self.optimizer = Adam(
            model.parameters(),
            lr=float(config.get("lr", 1e-3)),
            weight_decay=float(config.get("weight_decay", 1e-4)),
        )
        self.scheduler = ReduceLROnPlateau(
            self.optimizer, mode="min", factor=0.5, patience=10
        )

        self.max_gene_trees = int(config.get("num_gene_trees", 300))
        self.patience = int(config.get("patience", 20))
        self.max_epochs = int(config.get("max_epochs", 200))
        self.accumulation_steps = max(1, int(config.get("batch_size", 1)))

    def train_epoch(self, train_samples: List[Gene2NetSample]):
        self.model.train()
        total_loss = 0.0
        random.shuffle(train_samples)
        self.optimizer.zero_grad()

        for i, sample in enumerate(train_samples):
            gt_data, sp_data, wgd_targets, mask, wgd_mask, partner_targets, true_count = \
                prepare_batch(sample, self.device, self.max_gene_trees)

            wgd_logits, partner_logits, _ = self.model(gt_data, sp_data)

            loss = self.loss_fn(
                wgd_logits, partner_logits,
                wgd_targets, mask, wgd_mask, partner_targets, true_count
            )
            # Scale loss by accumulation steps for correct gradient magnitude
            (loss / self.accumulation_steps).backward()
            total_loss += loss.item()

            if (i + 1) % self.accumulation_steps == 0 or (i + 1) == len(train_samples):
                nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                self.optimizer.step()
                self.optimizer.zero_grad()

        return total_loss / max(len(train_samples), 1)

    @torch.no_grad()
    def evaluate(self, val_samples: List[Gene2NetSample]):
        self.model.eval()
        total_loss = 0.0
        all_f1 = []
        all_partner_acc = []
        all_count_err = []

        for sample in val_samples:
            gt_data, sp_data, wgd_targets, mask, wgd_mask, partner_targets, true_count = \
                prepare_batch(sample, self.device, self.max_gene_trees)

            wgd_logits, partner_logits, _ = self.model(gt_data, sp_data)

            loss = self.loss_fn(
                wgd_logits, partner_logits,
                wgd_targets, mask, wgd_mask, partner_targets, true_count
            )
            total_loss += loss.item()

            # Metrics
            metrics = wgd_f1(wgd_logits, wgd_targets, mask)
            all_f1.append(metrics["f1"])

            if wgd_mask.sum() > 0:
                wgd_partner_logits = partner_logits[wgd_mask]
                all_partner_acc.append(partner_accuracy(wgd_partner_logits, partner_targets))

            all_count_err.append(event_count_error(wgd_logits, int(true_count.item())))

        n = max(len(val_samples), 1)
        return {
            "val_loss": total_loss / n,
            "wgd_f1": sum(all_f1) / max(len(all_f1), 1),
            "partner_acc": sum(all_partner_acc) / max(len(all_partner_acc), 1),
            "count_error": sum(all_count_err) / max(len(all_count_err), 1),
        }

    def train(self, train_samples, val_samples, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        best_val_loss = float("inf")
        patience_counter = 0

        for epoch in range(self.max_epochs):
            train_loss = self.train_epoch(train_samples)
            val_metrics = self.evaluate(val_samples)
            val_loss = val_metrics["val_loss"]

            self.scheduler.step(val_loss)

            print(f"Epoch {epoch+1}/{self.max_epochs} | "
                  f"Train Loss: {train_loss:.4f} | "
                  f"Val Loss: {val_loss:.4f} | "
                  f"WGD F1: {val_metrics['wgd_f1']:.3f} | "
                  f"Partner Acc: {val_metrics['partner_acc']:.3f} | "
                  f"Count Err: {val_metrics['count_error']:.1f}")

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0
                torch.save(self.model.state_dict(), os.path.join(output_dir, "best_model.pt"))
            else:
                patience_counter += 1

            if patience_counter >= self.patience:
                print(f"Early stopping at epoch {epoch+1}")
                break

        # Load best model
        self.model.load_state_dict(torch.load(os.path.join(output_dir, "best_model.pt"), weights_only=True))
        return self.model
