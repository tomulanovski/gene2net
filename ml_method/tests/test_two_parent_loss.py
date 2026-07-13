import torch
from gene2net_gnn.training.trainer_reconstruct import two_parent_loss


def test_permutation_invariant():
    # 1 query, 4 edges, target parents {1, 3}. Loss identical for (a,b) and (b,a).
    Q, E = 1, 4
    torch.manual_seed(0)
    scores = torch.randn(Q, E, 2)
    a = torch.tensor([1]); b = torch.tensor([3])
    l1 = two_parent_loss(scores, a, b)
    l2 = two_parent_loss(scores, b, a)
    assert torch.allclose(l1, l2)


def test_confident_correct_pair_has_low_loss():
    Q, E = 1, 4
    scores = torch.full((Q, E, 2), -5.0)
    scores[0, 1, 0] = 10.0   # slot 0 -> edge 1
    scores[0, 3, 1] = 10.0   # slot 1 -> edge 3
    a = torch.tensor([1]); b = torch.tensor([3])
    loss = two_parent_loss(scores, a, b)
    assert loss.item() < 0.05


def test_auto_both_slots_same_edge():
    # auto: a == b == self edge. A model confident on that edge in both slots is low-loss.
    Q, E = 1, 4
    scores = torch.full((Q, E, 2), -5.0)
    scores[0, 2, 0] = 10.0
    scores[0, 2, 1] = 10.0
    a = torch.tensor([2]); b = torch.tensor([2])
    assert two_parent_loss(scores, a, b).item() < 0.05
