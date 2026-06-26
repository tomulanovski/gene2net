"""Warm-start compatibility filtering for retraining the reconstruct model.

When the partner head's input width changes (e.g. adding the copy-aware
cluster-support feature: 2H+2 -> 2H+4), warm-starting from an older reconstruct
checkpoint must drop the now-incompatible partner_head weight so it reinitializes
fresh, instead of crashing load_state_dict on a size mismatch.
"""
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2
from gene2net_gnn.training.trainer_reconstruct import filter_compatible_state_dict


def test_filter_drops_shape_mismatched_keys():
    old = SpeciesTreeGNNv2(partner_pair_feat_dim=2)
    new = SpeciesTreeGNNv2(partner_pair_feat_dim=4)

    filtered = filter_compatible_state_dict(old.state_dict(), new.state_dict())

    # The partner head's first linear weight changed shape -> dropped.
    assert "partner_head.0.weight" not in filtered
    # Shared, unchanged weights are kept.
    assert "node_proj.0.weight" in filtered
    # Loading the filtered dict no longer raises (even non-strict).
    missing, unexpected = new.load_state_dict(filtered, strict=False)
    assert "partner_head.0.weight" in missing
    assert unexpected == []


def test_filter_keeps_everything_when_compatible():
    a = SpeciesTreeGNNv2(partner_pair_feat_dim=4)
    b = SpeciesTreeGNNv2(partner_pair_feat_dim=4)
    filtered = filter_compatible_state_dict(a.state_dict(), b.state_dict())
    assert set(filtered.keys()) == set(b.state_dict().keys())
