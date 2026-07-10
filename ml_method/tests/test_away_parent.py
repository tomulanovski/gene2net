from ete3 import Tree

from gene2net_gnn.data.metadata_labels import relabel_events_partner_as_away
from gene2net_gnn.inference.mul_tree_builder import WGDEvent


def _bip(clades):
    return [(i, frozenset(c)) for i, c in enumerate(clades)]


def test_home_is_partner_retargets_to_other_parent():
    # True tree: X's other parent (its sibling) is A.
    true_tree = Tree("((X,A),(C,D));", format=1)
    # ASTRAL placed X next to B  -> home == labelled partner B  (the ~55% bug case)
    edge_bip = _bip([{"X"}, {"B"}, {"X", "B"}, {"C"}, {"D"}])
    ev = WGDEvent(frozenset({"X"}), frozenset({"B"}), 1.0)
    out = relabel_events_partner_as_away([ev], true_tree, edge_bip)
    assert out[0].partner_edge_clade == frozenset({"A"})   # retargeted to the away parent


def test_home_is_other_parent_keeps_label():
    true_tree = Tree("((X,A),(C,D));", format=1)
    # ASTRAL placed X next to A  -> partner B is already the away parent, keep it
    edge_bip = _bip([{"X"}, {"A"}, {"X", "A"}, {"C"}, {"D"}])
    ev = WGDEvent(frozenset({"X"}), frozenset({"B"}), 1.0)
    out = relabel_events_partner_as_away([ev], true_tree, edge_bip)
    assert out[0].partner_edge_clade == frozenset({"B"})   # unchanged


def test_auto_event_unchanged():
    true_tree = Tree("((X,A),(C,D));", format=1)
    edge_bip = _bip([{"X"}, {"X", "A"}])
    ev = WGDEvent(frozenset({"X"}), frozenset({"X"}), 1.0)   # auto: partner == target
    out = relabel_events_partner_as_away([ev], true_tree, edge_bip)
    assert out[0].partner_edge_clade == frozenset({"X"})     # auto untouched


def test_clade_level_event_unchanged():
    true_tree = Tree("((X,A),(C,D));", format=1)
    edge_bip = _bip([{"X", "Y"}, {"B"}])
    ev = WGDEvent(frozenset({"X", "Y"}), frozenset({"B"}), 1.0)   # clade-level target
    out = relabel_events_partner_as_away([ev], true_tree, edge_bip)
    assert out[0].partner_edge_clade == frozenset({"B"})     # left as-is (|target|>1)
