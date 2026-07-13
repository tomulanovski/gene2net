from ete3 import Tree
from gene2net_gnn.inference.mul_tree_builder import (
    TwoParentEvent, build_mul_tree_two_parent,
)


def _count(tree, name):
    return sum(1 for l in tree.get_leaf_names() if l == name)


def test_allo_detaches_and_grafts_at_both_parents():
    # X = D. True parents A = {A,B}, B(partner) = {C}. D must end up with
    # exactly 2 copies, and NEITHER should remain at D's old position under E.
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    ev = TwoParentEvent(frozenset({"D"}), frozenset({"A", "B"}), frozenset({"C"}), 1.0)
    mul = build_mul_tree_two_parent(sp, [ev])
    assert _count(mul, "D") == 2
    # The old (D,E) cherry must be broken: E is no longer directly paired with a
    # bare D leaf (D was detached from its backbone position and re-grafted at
    # its two true parents {A,B} and {C}).
    e_leaf = mul.search_nodes(name="E")[0]
    direct_siblings = [ch for ch in e_leaf.up.get_children() if ch is not e_leaf]
    assert not any(ch.is_leaf() and ch.name == "D" for ch in direct_siblings)


def test_auto_event_uses_sibling_duplication():
    # parent_a == parent_b == target => autopolyploidy: two identical D,E siblings.
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    ev = TwoParentEvent(frozenset({"D", "E"}), frozenset({"D", "E"}), frozenset({"D", "E"}), 1.0)
    mul = build_mul_tree_two_parent(sp, [ev])
    assert _count(mul, "D") == 2 and _count(mul, "E") == 2


def test_true_backbone_roundtrip_is_faithful():
    # Building the true events on the TRUE backbone must reproduce a MUL-tree
    # with the right multiplicities (the faithfulness guard from spec Sec 7).
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    ev = TwoParentEvent(frozenset({"D"}), frozenset({"A", "B"}), frozenset({"C"}), 1.0)
    mul, dropped = build_mul_tree_two_parent(sp, [ev], return_dropped=True)
    assert dropped == 0
    assert _count(mul, "D") == 2


def test_missing_parent_is_dropped_not_crash():
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    ev = TwoParentEvent(frozenset({"D"}), frozenset({"Z"}), frozenset({"C"}), 1.0)  # Z absent
    mul, dropped = build_mul_tree_two_parent(sp, [ev], return_dropped=True)
    assert dropped == 1
    assert _count(mul, "D") == 1  # unchanged
