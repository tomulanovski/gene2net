import pytest
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.data.label_extractor import _get_edge_bipartitions, decompose_mul_tree
from gene2net_gnn.data.metadata_labels import (
    events_from_metadata,
    sample_edge_bipartitions,
    labels_from_metadata_for_sample,
)


def test_events_from_metadata_allo_and_auto():
    md = [
        {"target_clade": ["D", "E"], "partner_clade": ["A", "B"], "event_type": "allo"},
        {"target_clade": ["A", "B", "C"], "partner_clade": None, "event_type": "auto"},
    ]
    evs = events_from_metadata(md)
    assert evs[0].wgd_edge_clade == frozenset({"D", "E"})
    assert evs[0].partner_edge_clade == frozenset({"A", "B"})
    # auto: partner == target (self)
    assert evs[1].wgd_edge_clade == evs[1].partner_edge_clade == frozenset({"A", "B", "C"})


def test_events_from_metadata_rejects_bad_input():
    with pytest.raises(ValueError):
        events_from_metadata([{"target_clade": ["D"], "partner_clade": None, "event_type": "allo"}])
    with pytest.raises(ValueError):
        events_from_metadata([{"target_clade": ["D"], "event_type": "weird"}])


def _build_sample():
    sp_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    gts = [Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1) for _ in range(6)]
    gts += [Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1) for _ in range(4)]
    s = Gene2NetSample.from_trees(sp_tree, gts, ["A", "B", "C", "D", "E"])
    return s


def _sample_dict(s):
    return {
        "species_tree_edge_index": s.species_tree_edge_index,
        "species_tree_node_names": s.species_tree_node_names,
        "species_tree_is_leaf": s.species_tree_is_leaf,
        "species_list": s.species_list,
        "species_tree_edge_features": s.species_tree_edge_features,
    }


def test_clade_level_allo_labeled_on_ancestral_edge():
    s = _build_sample()
    sd = _sample_dict(s)
    md = [{"target_clade": ["D", "E"], "partner_clade": ["A", "B"], "event_type": "allo"}]
    labels = labels_from_metadata_for_sample(md, sd)

    # exactly one mappable event
    assert labels.mask == [True]
    # interpret the chosen edge indices via the sample's bipartitions
    bip = dict(sample_edge_bipartitions(sd))
    assert set(bip[labels.wgd_edges[0]]) == {"D", "E"}        # on the ancestral clade edge
    assert set(bip[labels.partner_edges[0]]) == {"A", "B"}    # true external partner
    # alignment guard: label edge-count == feature rows
    assert labels.n_edges == sd["species_tree_edge_features"].shape[0]


def test_metadata_path_does_not_fragment_like_decompose():
    # The same clade-level event, decomposed from a MUL-tree, fragments to reciprocal tips.
    mul = Tree("(((A:1,B:1):1,(D:1,E:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    frag = [set(e.wgd_edge_clade) for e in decompose_mul_tree(mul)]
    assert {"D"} in frag and {"E"} in frag          # old behavior: fragmented tips
    assert {"D", "E"} not in frag                    # clade never appears as one event
