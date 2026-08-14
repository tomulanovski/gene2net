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


def test_build_mul_tree_superset_fallback_applies_contaminated_clade():
    # event1 grafts a copy of A into the (D,E) region, so the ancestor clade
    # {C,D,E} becomes {A,C,D,E}. The auto event2 on {C,D,E} no longer matches
    # exactly, but _find_node_by_leaf_set's superset fallback (added to recover
    # clade-level events after an inner graft) finds {A,C,D,E} and APPLIES event2
    # there instead of dropping it -> dropped == 0. Deliberate tradeoff on the
    # decompose->build path (NOT the reported scoring path): recovering a
    # contaminated clade is preferred over silently losing the event.
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    events = [
        WGDEvent(frozenset({"A"}), frozenset({"D", "E"}), 0.9),
        WGDEvent(frozenset({"C", "D", "E"}), frozenset({"C", "D", "E"}), 0.9),
    ]
    tree, dropped = build_mul_tree(sp, events, return_dropped=True)
    assert dropped == 0


def test_build_mul_tree_no_drops_clean_case():
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    events = [WGDEvent(frozenset({"D"}), frozenset({"A", "B"}), 0.9)]
    tree, dropped = build_mul_tree(sp, events, return_dropped=True)
    assert dropped == 0
    assert tree.get_leaf_names().count("D") == 2


def test_build_mul_tree_whole_tree_event_no_hang():
    # A whole-tree clade would alias the root in _apply_wgd_event and hang on the
    # next traversal. It must be skipped (dropped), leaving the tree unchanged.
    sp = Tree("(A:1,(B:1,C:1):1);", format=1)
    events = [WGDEvent(frozenset({"A", "B", "C"}), frozenset({"A", "B", "C"}), 0.9)]
    tree, dropped = build_mul_tree(sp, events, return_dropped=True)
    assert dropped == 1
    assert sorted(tree.get_leaf_names()) == ["A", "B", "C"]
