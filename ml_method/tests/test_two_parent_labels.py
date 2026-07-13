from ete3 import Tree
from gene2net_gnn.data.metadata_labels import home_edge_for_event
from gene2net_gnn.data.label_extractor import _get_edge_bipartitions, _best_matching_edge


def test_home_edge_is_true_sibling():
    # true tree: D's sibling is E. Backbone has an edge whose clade is {E}.
    backbone = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    true_tree = backbone.copy()
    edges = _get_edge_bipartitions(backbone)
    home_clade = home_edge_for_event(frozenset({"D"}), "allo", true_tree)
    assert home_clade == frozenset({"E"})
    idx, score = _best_matching_edge(home_clade, edges)
    assert score == 1.0  # {E} is an exact backbone edge


def test_auto_home_is_self():
    true_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    assert home_edge_for_event(frozenset({"D", "E"}), "auto", true_tree) == frozenset({"D", "E"})


def test_clade_home_is_clade_sibling():
    true_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    assert home_edge_for_event(frozenset({"D", "E"}), "allo", true_tree) == frozenset({"C"})
