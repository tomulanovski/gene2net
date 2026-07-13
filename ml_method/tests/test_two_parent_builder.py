from ete3 import Tree
from gene2net_gnn.inference.mul_tree_builder import (
    TwoParentEvent, build_mul_tree_two_parent,
)


def _count(tree, name):
    return sum(1 for l in tree.get_leaf_names() if l == name)


def _sibling_leaves(tree, name):
    leaf = tree.search_nodes(name=name)[0]
    sibs = set()
    for ch in leaf.up.get_children():
        if ch is not leaf:
            sibs |= set(ch.get_leaf_names())
    return sibs


def test_graft_keeps_home_and_adds_copy_at_away():
    # X = D, home is E (its current sibling), away parent = {C}. graft mode keeps D
    # next to E and adds a copy next to C -> D appears twice, home preserved.
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    ev = TwoParentEvent(frozenset({"D"}), frozenset({"A", "B"}), frozenset({"C"}), 1.0)
    mul = build_mul_tree_two_parent(sp, [ev])            # default mode="graft"
    assert _count(mul, "D") == 2
    assert "E" in _sibling_leaves(mul, "D") or "D" in _sibling_leaves(mul, "E")  # home kept
    assert "C" in _sibling_leaves(mul, "D") or "D" in _sibling_leaves(mul, "C")  # copy at away


def test_auto_event_uses_sibling_duplication():
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    ev = TwoParentEvent(frozenset({"D", "E"}), frozenset({"D", "E"}), frozenset({"D", "E"}), 1.0)
    mul = build_mul_tree_two_parent(sp, [ev])
    assert _count(mul, "D") == 2 and _count(mul, "E") == 2


def test_nested_events_compose_under_graft():
    # sample_0057 pattern: an outer clade event {sp0,sp3} + an inner event {sp3}.
    # graft mode must apply BOTH -> sp0 x2, sp3 x3. (detach mode drops the outer one.)
    sp = Tree("(((sp0:1,sp3:1):1,(sp2:1,sp4:1):1):1,sp1:1);", format=1)
    events = [
        TwoParentEvent(frozenset({"sp0", "sp3"}), frozenset({"sp2", "sp4"}), frozenset({"sp1"}), 1.0),
        TwoParentEvent(frozenset({"sp3"}), frozenset({"sp0"}), frozenset({"sp2"}), 1.0),
    ]
    mul, dropped = build_mul_tree_two_parent(sp, events, return_dropped=True)  # graft
    assert dropped == 0
    assert _count(mul, "sp0") == 2
    assert _count(mul, "sp3") == 3


def test_detach_mode_drops_nested_outer_event():
    # The old detach mode tears sp3 out of {sp0,sp3}, so the outer event is dropped.
    sp = Tree("(((sp0:1,sp3:1):1,(sp2:1,sp4:1):1):1,sp1:1);", format=1)
    events = [
        TwoParentEvent(frozenset({"sp0", "sp3"}), frozenset({"sp2", "sp4"}), frozenset({"sp1"}), 1.0),
        TwoParentEvent(frozenset({"sp3"}), frozenset({"sp0"}), frozenset({"sp2"}), 1.0),
    ]
    mul, dropped = build_mul_tree_two_parent(sp, events, return_dropped=True, mode="detach")
    # sp0 stays at 1 copy because the {sp0,sp3} event could not be applied.
    assert _count(mul, "sp0") == 1 or dropped >= 1


def test_true_backbone_roundtrip_is_faithful():
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    ev = TwoParentEvent(frozenset({"D"}), frozenset({"A", "B"}), frozenset({"C"}), 1.0)
    mul, dropped = build_mul_tree_two_parent(sp, [ev], return_dropped=True)
    assert dropped == 0
    assert _count(mul, "D") == 2


def test_missing_parent_is_dropped_not_crash():
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    ev = TwoParentEvent(frozenset({"D"}), frozenset({"Z"}), frozenset({"W"}), 1.0)  # both absent
    mul, dropped = build_mul_tree_two_parent(sp, [ev], return_dropped=True)
    assert dropped == 1
    assert _count(mul, "D") == 1  # unchanged
