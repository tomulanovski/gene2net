import torch
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2


def test_partner_scores_have_two_slots():
    m = SpeciesTreeGNNv2(node_feat_dim=13, edge_feat_dim=4, hidden_dim=16,
                         n_gat_layers=2, n_gat_heads=2)
    m.eval()  # disable dropout so the two code paths are deterministic
    E = 5
    edge_emb = torch.randn(E, 16)
    full = m.compute_partner_scores(edge_emb)
    assert full.shape == (E, E, 2)
    q = torch.tensor([0, 3])
    rows = m.compute_partner_scores_rows(edge_emb, q)
    assert rows.shape == (2, E, 2)
    # rows must equal the gathered full rows
    assert torch.allclose(rows, full[q], atol=1e-5)
