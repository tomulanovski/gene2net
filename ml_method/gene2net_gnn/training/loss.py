"""Loss functions for Gene2Net training."""
import torch
import torch.nn as nn
import torch.nn.functional as F


def focal_loss(logits, targets, alpha=0.25, gamma=2.0, class_weights=None):
    """Focal loss for class imbalance with per-class weights.

    Args:
        logits: [B, C] raw logits
        targets: [B] class indices
        alpha: base weighting factor
        gamma: focusing parameter (higher = more focus on hard examples)
        class_weights: [C] per-class weights (e.g., upweight rare WGD classes)
    Returns:
        scalar loss
    """
    if class_weights is not None:
        ce_loss = F.cross_entropy(logits, targets, weight=class_weights, reduction='none')
    else:
        ce_loss = F.cross_entropy(logits, targets, reduction='none')
    pt = torch.exp(-ce_loss)  # probability of correct class
    focal_weight = alpha * (1 - pt) ** gamma
    return (focal_weight * ce_loss).mean()


class Gene2NetLoss(nn.Module):
    """Combined loss for Gene2Net model.

    L = lambda_wgd * L_wgd + lambda_partner * L_partner + lambda_count * L_count
    """
    def __init__(self, lambda_wgd=1.0, lambda_partner=1.0, lambda_count=0.1,
                 focal_alpha=0.25, focal_gamma=2.0,
                 wgd_class_weights=None):
        super().__init__()
        self.lambda_wgd = lambda_wgd
        self.lambda_partner = lambda_partner
        self.lambda_count = lambda_count
        self.focal_alpha = focal_alpha
        self.focal_gamma = focal_gamma
        # Class weights: [weight_for_0_events, weight_for_1, weight_for_2, weight_for_3]
        # Higher weight for rare classes forces the model to pay attention to WGD edges
        if wgd_class_weights is not None:
            self.register_buffer("wgd_class_weights", torch.tensor(wgd_class_weights, dtype=torch.float))
        else:
            # Default: heavily upweight WGD classes
            self.register_buffer("wgd_class_weights", torch.tensor([1.0, 10.0, 10.0, 10.0], dtype=torch.float))

    def wgd_loss(self, wgd_logits, wgd_targets, mask):
        """Focal loss on WGD count prediction, masked for unmappable edges."""
        if mask.sum() == 0:
            return torch.tensor(0.0, device=wgd_logits.device)
        weights = self.wgd_class_weights.to(wgd_logits.device)
        return focal_loss(
            wgd_logits[mask], wgd_targets[mask],
            alpha=self.focal_alpha, gamma=self.focal_gamma,
            class_weights=weights,
        )

    def partner_loss(self, partner_logits, partner_targets, wgd_mask):
        """Cross-entropy on partner prediction, only for WGD edges."""
        if wgd_mask.sum() == 0:
            return torch.tensor(0.0, device=partner_logits.device)
        # Select only rows for edges that have WGD events
        wgd_partner_logits = partner_logits[wgd_mask]  # [N_wgd, E]
        return F.cross_entropy(wgd_partner_logits, partner_targets)

    def count_loss(self, wgd_logits, true_count):
        """MSE on predicted vs true total event count."""
        # Predicted count = sum of argmax per edge (or soft sum)
        pred_probs = F.softmax(wgd_logits, dim=-1)  # [E, C]
        event_values = torch.arange(pred_probs.shape[1], device=pred_probs.device, dtype=torch.float)
        pred_counts_per_edge = (pred_probs * event_values).sum(dim=-1)  # [E]
        pred_total = pred_counts_per_edge.sum()
        return F.mse_loss(pred_total.unsqueeze(0), true_count)

    def forward(self, wgd_logits, partner_logits, wgd_targets, mask, wgd_mask, partner_targets, true_count):
        """Compute combined loss."""
        l_wgd = self.wgd_loss(wgd_logits, wgd_targets, mask)
        l_partner = self.partner_loss(partner_logits, partner_targets, wgd_mask)
        l_count = self.count_loss(wgd_logits, true_count)

        return self.lambda_wgd * l_wgd + self.lambda_partner * l_partner + self.lambda_count * l_count
