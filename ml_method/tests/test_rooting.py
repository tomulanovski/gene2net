from ete3 import Tree

from gene2net_gnn.data.rooting import robust_set_outgroup


def basal_sides(tree):
    return set(frozenset(c.get_leaf_names()) for c in tree.children)


def test_reroot_simple_clade():
    t = Tree("((A,B),(C,D));", format=1)
    assert robust_set_outgroup(t, {"A", "B"})
    assert basal_sides(t) == {frozenset({"A", "B"}), frozenset({"C", "D"})}


def test_reroot_outgroup_spanning_arbitrary_root():
    # Arbitrary root at A: the outgroup {A,B} spans the current root, so the naive
    # set_outgroup(get_common_ancestor(...)) throws. The robust version must handle it.
    t = Tree("(A,(B,(C,D)));", format=1)
    assert robust_set_outgroup(t, {"A", "B"})
    assert basal_sides(t) == {frozenset({"A", "B"}), frozenset({"C", "D"})}


def test_reroot_single_leaf():
    t = Tree("(A,(B,(C,D)));", format=1)
    assert robust_set_outgroup(t, {"C"})
    assert frozenset({"C"}) in basal_sides(t)


def test_reroot_non_monophyletic_returns_false():
    # {A,C} is not a clade in this unrooted topology -> cannot root on it.
    t = Tree("((A,B),(C,D));", format=1)
    assert not robust_set_outgroup(t, {"A", "C"})
