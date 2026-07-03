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


def test_root_at_reference_matches_reference_basal_split():
    from gene2net_gnn.data.rooting import root_at_reference
    ref = Tree("((A,B),(C,D));", format=1)     # rooted: {A,B} | {C,D}
    tree = Tree("(A,(B,(C,D)));", format=1)     # same topology, arbitrary root
    assert root_at_reference(tree, ref)
    assert basal_sides(tree) == {frozenset({"A", "B"}), frozenset({"C", "D"})}


def test_root_at_reference_with_name_map():
    # Reference uses original names, tree uses renamed taxa (as in the benchmark).
    from gene2net_gnn.data.rooting import root_at_reference
    ref = Tree("((Ao,Bo),(Co,Do));", format=1)
    tree = Tree("(A,(B,(C,D)));", format=1)
    name_map = {"Ao": "A", "Bo": "B", "Co": "C", "Do": "D"}
    assert root_at_reference(tree, ref, name_map=name_map)
    assert basal_sides(tree) == {frozenset({"A", "B"}), frozenset({"C", "D"})}
