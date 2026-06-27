"""The partner head must consume only the pairwise-feature columns it was built
for, ignoring extras.

build_pairwise_feat now emits a 4-dim feature (2 co-clustering + 2 cluster
support). A model built with partner_pair_feat_dim=2 (e.g. an older checkpoint)
must still run against that wider feature, using only its first 2 columns, rather
than crashing on a shape mismatch.
"""
import torch

from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2


def test_partner_head_ignores_extra_pairwise_dims():
    torch.manual_seed(0)
    m = SpeciesTreeGNNv2(partner_pair_feat_dim=2).eval()
    E = 5
    emb = torch.randn(E, m.hidden_dim)
    wide = torch.randn(E, E, 4)          # feature wider than the model's 2
    q = torch.arange(E)

    out_wide = m.compute_partner_scores_rows(emb, q, wide)
    out_sliced = m.compute_partner_scores_rows(emb, q, wide[..., :2])
    assert out_wide.shape == (E, E)
    assert torch.allclose(out_wide, out_sliced)

    full_wide = m.compute_partner_scores(emb, wide)
    full_sliced = m.compute_partner_scores(emb, wide[..., :2])
    assert torch.allclose(full_wide, full_sliced)
