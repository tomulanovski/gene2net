from ete3 import Tree
from gene2net_gnn.inference.mul_tree_builder import build_mul_tree, WGDEvent

def test_build_mul_tree_single_allo():
    species_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    events = [WGDEvent(
        wgd_edge_clade=frozenset({"D"}),
        partner_edge_clade=frozenset({"A", "B"}),
        confidence=0.9
    )]
    mul_tree = build_mul_tree(species_tree, events)
    d_count = sum(1 for l in mul_tree.get_leaf_names() if l == "D")
    assert d_count == 2

def test_build_mul_tree_auto():
    species_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    events = [WGDEvent(
        wgd_edge_clade=frozenset({"D", "E"}),
        partner_edge_clade=frozenset({"D", "E"}),  # auto = partner is same clade
        confidence=0.9
    )]
    mul_tree = build_mul_tree(species_tree, events)
    d_count = sum(1 for l in mul_tree.get_leaf_names() if l == "D")
    e_count = sum(1 for l in mul_tree.get_leaf_names() if l == "E")
    assert d_count == 2
    assert e_count == 2

def test_build_mul_tree_no_events():
    species_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    mul_tree = build_mul_tree(species_tree, [])
    assert sorted(mul_tree.get_leaf_names()) == sorted(species_tree.get_leaf_names())

def test_build_mul_tree_two_events():
    """Two independent allo events."""
    species_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    events = [
        WGDEvent(wgd_edge_clade=frozenset({"D"}), partner_edge_clade=frozenset({"A", "B"}), confidence=0.9),
        WGDEvent(wgd_edge_clade=frozenset({"C"}), partner_edge_clade=frozenset({"A", "B"}), confidence=0.8),
    ]
    mul_tree = build_mul_tree(species_tree, events)
    d_count = sum(1 for l in mul_tree.get_leaf_names() if l == "D")
    c_count = sum(1 for l in mul_tree.get_leaf_names() if l == "C")
    assert d_count == 2
    assert c_count == 2
