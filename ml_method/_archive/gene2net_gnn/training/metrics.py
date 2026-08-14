"""Evaluation metrics for Gene2Net predictions."""
import torch
from typing import Dict


def wgd_f1(pred_logits: torch.Tensor, true_labels: torch.Tensor, mask: torch.Tensor) -> Dict[str, float]:
    """F1 score for WGD event detection (binary: event vs no event).

    Args:
        pred_logits: [E, C] raw logits
        true_labels: [E] target classes (0 = no event, >0 = event)
        mask: [E] bool mask for valid edges
    Returns:
        dict with precision, recall, f1
    """
    if mask.sum() == 0:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0}

    pred = pred_logits[mask].argmax(dim=-1) > 0  # binary: has event
    true = true_labels[mask] > 0

    tp = (pred & true).sum().float()
    fp = (pred & ~true).sum().float()
    fn = (~pred & true).sum().float()

    precision = (tp / (tp + fp)).item() if (tp + fp) > 0 else 0.0
    recall = (tp / (tp + fn)).item() if (tp + fn) > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0

    return {"precision": precision, "recall": recall, "f1": f1}


def partner_accuracy(pred_logits: torch.Tensor, true_partners: torch.Tensor) -> float:
    """Accuracy of partner prediction for detected WGD edges.

    Args:
        pred_logits: [N_wgd, E] partner scores
        true_partners: [N_wgd] true partner edge indices
    Returns:
        accuracy (float)
    """
    if len(true_partners) == 0:
        return 0.0
    pred = pred_logits.argmax(dim=-1)
    return (pred == true_partners).float().mean().item()


def event_count_error(pred_logits: torch.Tensor, true_count: int) -> float:
    """Absolute error in predicted total event count.

    Args:
        pred_logits: [E, C] WGD logits
        true_count: true total number of events
    Returns:
        absolute error
    """
    pred_count = pred_logits.argmax(dim=-1).sum().item()
    return abs(pred_count - true_count)
