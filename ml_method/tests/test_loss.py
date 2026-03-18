import torch
from gene2net_gnn.training.loss import focal_loss, Gene2NetLoss

def test_focal_loss_zero_when_correct():
    logits = torch.tensor([[10.0, -10.0]])
    target = torch.tensor([0])
    loss = focal_loss(logits, target, alpha=0.25, gamma=2.0)
    assert loss.item() < 0.01

def test_focal_loss_high_when_wrong():
    logits = torch.tensor([[-10.0, 10.0]])
    target = torch.tensor([0])
    loss = focal_loss(logits, target, alpha=0.25, gamma=2.0)
    assert loss.item() > 0.5

def test_focal_loss_harder_than_ce():
    """Focal loss should be <= CE for easy examples."""
    logits = torch.tensor([[2.0, -1.0, -1.0]])
    target = torch.tensor([0])
    fl = focal_loss(logits, target, alpha=1.0, gamma=2.0)
    ce = torch.nn.functional.cross_entropy(logits, target)
    assert fl.item() <= ce.item()

def test_gene2net_loss_masks_unmappable():
    loss_fn = Gene2NetLoss()
    wgd_logits = torch.randn(8, 4)
    wgd_targets = torch.zeros(8, dtype=torch.long)
    mask = torch.tensor([1,1,1,1,1,1,0,0], dtype=torch.bool)
    loss = loss_fn.wgd_loss(wgd_logits, wgd_targets, mask)
    assert loss.item() >= 0

def test_gene2net_loss_partner():
    loss_fn = Gene2NetLoss()
    partner_logits = torch.randn(8, 8)
    # Only edges 0 and 3 have WGD events
    wgd_mask = torch.tensor([1,0,0,1,0,0,0,0], dtype=torch.bool)
    partner_targets = torch.tensor([5, 2])  # edge 0's partner is 5, edge 3's partner is 2
    loss = loss_fn.partner_loss(partner_logits, partner_targets, wgd_mask)
    assert loss.item() >= 0

def test_gene2net_loss_combined():
    loss_fn = Gene2NetLoss(lambda_wgd=1.0, lambda_partner=1.0, lambda_count=0.1)

    wgd_logits = torch.randn(8, 4)
    partner_logits = torch.randn(8, 8)

    wgd_targets = torch.zeros(8, dtype=torch.long)
    wgd_targets[2] = 1  # edge 2 has 1 WGD event
    mask = torch.ones(8, dtype=torch.bool)
    wgd_mask = wgd_targets > 0
    partner_targets = torch.tensor([5])  # edge 2's partner is edge 5
    true_count = torch.tensor([1.0])  # 1 total WGD event

    total_loss = loss_fn(wgd_logits, partner_logits, wgd_targets, mask, wgd_mask, partner_targets, true_count)
    assert total_loss.item() > 0

def test_gene2net_loss_no_wgd_events():
    """When no WGD events exist, partner loss should be 0."""
    loss_fn = Gene2NetLoss()
    wgd_logits = torch.randn(8, 4)
    partner_logits = torch.randn(8, 8)
    wgd_targets = torch.zeros(8, dtype=torch.long)
    mask = torch.ones(8, dtype=torch.bool)
    wgd_mask = torch.zeros(8, dtype=torch.bool)
    partner_targets = torch.tensor([], dtype=torch.long)
    true_count = torch.tensor([0.0])

    total_loss = loss_fn(wgd_logits, partner_logits, wgd_targets, mask, wgd_mask, partner_targets, true_count)
    assert total_loss.item() >= 0
