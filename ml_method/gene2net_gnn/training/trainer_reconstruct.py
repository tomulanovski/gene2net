"""Reconstruction trainer: joint WGD detection + partner-edge prediction.

Builds on the Phase 1 detection model by also training the partner head.
For every WGD edge the partner head points to the edge its duplicate attaches
to; pointing to itself means autopolyploidy, to another edge means
allopolyploidy. So auto/allo is not a separate task — it falls out of the
partner prediction.

Loss = focal detection loss + partner_weight * partner cross-entropy
(the latter only over WGD edges that have a mappable partner label).
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
from gene2net_gnn.data.features import (
    species_coclustering_matrix,
    edge_clades_species,
    pairwise_partner_features,
)
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2
from gene2net_gnn.training.trainer_phase1 import prepare_sample, focal_loss


def build_pairwise_feat(sample):
    """Compute (and cache) the [E, E, 2] clade co-clustering pairwise feature.

    Uses the preorder edge ordering (sample._edge_index_pre, set by
    prepare_sample) so it aligns with the model's edge embeddings.
    """
    cached = getattr(sample, "_pairwise_feat", None)
    if cached is not None:
        return cached
    n_species = len(sample.species_list)
    C = species_coclustering_matrix(
        sample.gene_tree_edge_indices, sample.gene_tree_species_ids,
        sample.gene_tree_leaf_masks, n_species,
    )
    sp_to_idx = {sp: i for i, sp in enumerate(sample.species_list)}
    node_species = [sp_to_idx.get(name, -1) for name in sample.species_tree_node_names]
    clades = edge_clades_species(sample._edge_index_pre, sample.species_tree_is_leaf, node_species)
    feat = pairwise_partner_features(C, clades)
    sample._pairwise_feat = feat
    return feat


def build_partner_targets(sample: Gene2NetSample, n_edges: int, device: torch.device):
    """Per-WGD-edge partner index target (-1 = no target / ignore).

    For each mappable event, target[wgd_edge] = partner_edge (== wgd_edge for
    autopolyploidy). Edges that are not WGD, or whose event was unmappable, stay
    at -1 and are excluded from the partner loss.
    """
    labels = sample.labels
    targets = torch.full((n_edges,), -1, dtype=torch.long, device=device)
    if labels is None or not labels.wgd_edges or not labels.partner_edges:
        return targets

    mask = labels.mask or []
    for k in range(len(labels.wgd_edges)):
        if k < len(mask) and not mask[k]:
            continue  # unmappable event
        w = labels.wgd_edges[k]
        p = labels.partner_edges[k] if k < len(labels.partner_edges) else None
        if p is None:
            continue
        if 0 <= w < n_edges and 0 <= p < n_edges:
            targets[w] = p
    return targets


class ReconstructTrainer:
    """Joint detection + partner-prediction training loop."""

    def __init__(self, model: SpeciesTreeGNNv2, config: dict, device: str = "cpu"):
        self.model = model.to(device)
        self.device = torch.device(device)
        self.config = config

        pos_weight = float(config.get("pos_class_weight", 12.0))
        self.class_weights = torch.tensor([1.0, pos_weight], dtype=torch.float, device=self.device)
        self.focal_alpha = float(config.get("focal_alpha", 0.25))
        self.focal_gamma = float(config.get("focal_gamma", 2.0))
        self.partner_weight = float(config.get("partner_weight", 1.0))

        self.optimizer = Adam(
            model.parameters(),
            lr=float(config.get("lr", 1e-3)),
            weight_decay=float(config.get("weight_decay", 1e-4)),
        )
        self.scheduler = ReduceLROnPlateau(self.optimizer, mode="min", factor=0.5, patience=10)

        self.batch_size = max(1, int(config.get("batch_size", 8)))
        self.patience = int(config.get("patience", 30))
        self.max_epochs = int(config.get("max_epochs", 200))

    def _losses_for_sample(self, sample):
        """Return (det_loss, partner_loss_or_None, n_edges) for one sample."""
        prepared = prepare_sample(sample, self.device)
        if prepared is None:
            return None
        inputs, wgd_targets, mask = prepared
        n_edges = wgd_targets.shape[0]

        wgd_logits, edge_emb = self.model(**inputs)
        if mask.sum() == 0:
            return None

        det_loss = focal_loss(
            wgd_logits[mask], wgd_targets[mask],
            alpha=self.focal_alpha, gamma=self.focal_gamma,
            class_weights=self.class_weights,
        )

        partner_targets = build_partner_targets(sample, n_edges, self.device)
        pmask = partner_targets >= 0
        partner_loss = None
        if pmask.any():
            pairwise_feat = build_pairwise_feat(sample).to(self.device)
            partner_scores = self.model.compute_partner_scores(edge_emb, pairwise_feat)  # [E, E]
            partner_loss = F.cross_entropy(partner_scores[pmask], partner_targets[pmask])

        return det_loss, partner_loss, wgd_logits, mask, wgd_targets

    def train_epoch(self, train_samples):
        self.model.train()
        total = 0.0
        n = 0
        random.shuffle(train_samples)
        self.optimizer.zero_grad()

        for i, sample in enumerate(train_samples):
            out = self._losses_for_sample(sample)
            if out is None:
                continue
            det_loss, partner_loss, *_ = out
            loss = det_loss
            if partner_loss is not None:
                loss = loss + self.partner_weight * partner_loss

            (loss / self.batch_size).backward()
            total += loss.item()
            n += 1

            if (i + 1) % self.batch_size == 0 or (i + 1) == len(train_samples):
                nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                self.optimizer.step()
                self.optimizer.zero_grad()

        return total / max(n, 1)

    @torch.no_grad()
    def evaluate(self, val_samples):
        self.model.eval()
        total_loss = 0.0
        n = 0
        tp = fp = fn = 0
        partner_correct = 0
        partner_total = 0
        auto_correct = allo_correct = auto_total = allo_total = 0

        for sample in val_samples:
            prepared = prepare_sample(sample, self.device)
            if prepared is None:
                continue
            inputs, wgd_targets, mask = prepared
            n_edges = wgd_targets.shape[0]
            wgd_logits, edge_emb = self.model(**inputs)
            if mask.sum() == 0:
                continue

            det_loss = focal_loss(
                wgd_logits[mask], wgd_targets[mask],
                alpha=self.focal_alpha, gamma=self.focal_gamma,
                class_weights=self.class_weights,
            )
            total_loss += det_loss.item()
            n += 1

            preds = wgd_logits[mask].argmax(dim=-1)
            true = wgd_targets[mask]
            tp += int((preds & true).sum())
            fp += int((preds & ~true.bool()).sum())
            fn += int((~preds.bool() & true.bool()).sum())

            # Partner accuracy on edges with a partner label
            partner_targets = build_partner_targets(sample, n_edges, self.device)
            pmask = partner_targets >= 0
            if pmask.any():
                pairwise_feat = build_pairwise_feat(sample).to(self.device)
                scores = self.model.compute_partner_scores(edge_emb, pairwise_feat)
                pred_partner = scores[pmask].argmax(dim=-1)
                tgt = partner_targets[pmask]
                idx = pmask.nonzero(as_tuple=True)[0]
                partner_correct += int((pred_partner == tgt).sum())
                partner_total += int(pmask.sum())
                # auto = target points to self; allo = target points elsewhere
                is_auto = (tgt == idx)
                auto_total += int(is_auto.sum())
                allo_total += int((~is_auto).sum())
                auto_correct += int(((pred_partner == tgt) & is_auto).sum())
                allo_correct += int(((pred_partner == tgt) & ~is_auto).sum())

        precision = tp / max(tp + fp, 1)
        recall = tp / max(tp + fn, 1)
        f1 = 2 * precision * recall / max(precision + recall, 1e-8)
        return {
            "val_loss": total_loss / max(n, 1),
            "precision": precision, "recall": recall, "f1": f1,
            "partner_acc": partner_correct / max(partner_total, 1),
            "auto_acc": auto_correct / max(auto_total, 1),
            "allo_acc": allo_correct / max(allo_total, 1),
            "partner_total": partner_total, "auto_total": auto_total, "allo_total": allo_total,
        }

    def train(self, train_samples, val_samples, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        best_val_loss = float("inf")
        best_partner = 0.0
        patience_counter = 0

        for epoch in range(self.max_epochs):
            train_loss = self.train_epoch(train_samples)
            val = self.evaluate(val_samples)
            self.scheduler.step(val["val_loss"])
            lr = self.optimizer.param_groups[0]["lr"]

            print(
                f"Epoch {epoch+1}/{self.max_epochs} | Train {train_loss:.4f} | "
                f"Val {val['val_loss']:.4f} | F1 {val['f1']:.3f} "
                f"(P {val['precision']:.3f} R {val['recall']:.3f}) | "
                f"PartnerAcc {val['partner_acc']:.3f} "
                f"(auto {val['auto_acc']:.3f}/{val['auto_total']}, "
                f"allo {val['allo_acc']:.3f}/{val['allo_total']}) | LR {lr:.6f}"
            )

            if val["val_loss"] < best_val_loss:
                best_val_loss = val["val_loss"]
                patience_counter = 0
                torch.save(self.model.state_dict(), os.path.join(output_dir, "best_model.pt"))
            else:
                patience_counter += 1

            if val["partner_acc"] > best_partner:
                best_partner = val["partner_acc"]
                torch.save(self.model.state_dict(), os.path.join(output_dir, "best_partner_model.pt"))

            if patience_counter >= self.patience:
                print(f"Early stopping at epoch {epoch+1}")
                break

        print(f"\nBest val loss: {best_val_loss:.4f}")
        print(f"Best partner accuracy: {best_partner:.3f}")
        self.model.load_state_dict(
            torch.load(os.path.join(output_dir, "best_model.pt"), weights_only=True)
        )
        return self.model
