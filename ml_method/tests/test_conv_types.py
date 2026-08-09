import torch
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2


def _tiny_inputs():
    # root(0) with two leaves(1,2); undirected edges, even indices are parent->child
    node_features = torch.randn(3, 13)
    edge_index = torch.tensor([[0, 1, 0, 2], [1, 0, 2, 0]], dtype=torch.long)
    edge_features = torch.randn(2, 9)          # 2 undirected edges, edge_feat_dim=9
    is_leaf = torch.tensor([False, True, True])
    return node_features, edge_index, edge_features, is_leaf


def test_all_conv_types_build_and_forward():
    nf, ei, ef, il = _tiny_inputs()
    for ct in ("gat", "gin", "gcn"):
        m = SpeciesTreeGNNv2(node_feat_dim=13, edge_feat_dim=9, hidden_dim=16,
                             n_gat_layers=2, conv_type=ct,
                             partner_pair_feat_dim=4, n_parents=1)
        logits, edge_emb = m(nf, ei, ef, il)
        assert logits.shape == (2, 2)
        assert edge_emb.shape == (2, 16)


def test_default_conv_is_gat_and_param_count_matches():
    common = dict(node_feat_dim=13, edge_feat_dim=9, hidden_dim=16, n_gat_layers=2)
    m_default = SpeciesTreeGNNv2(**common)
    m_gat = SpeciesTreeGNNv2(**common, conv_type="gat")
    assert m_default.conv_type == "gat"
    assert sum(p.numel() for p in m_default.parameters()) == \
           sum(p.numel() for p in m_gat.parameters())


def test_unknown_conv_type_raises():
    import pytest
    with pytest.raises(ValueError):
        SpeciesTreeGNNv2(node_feat_dim=13, edge_feat_dim=9, hidden_dim=16,
                         n_gat_layers=1, conv_type="sage")
