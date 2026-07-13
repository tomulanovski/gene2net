from ete3 import Tree
from scripts.two_parent_oracle import (
    sibling_of_clade, snap_to_backbone, two_parent_events, predicted_parents_coclust,
)


def test_sibling_of_single_species():
    t = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    assert sibling_of_clade(t, frozenset({"D"})) == frozenset({"E"})


def test_sibling_of_clade():
    t = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    assert sibling_of_clade(t, frozenset({"D", "E"})) == frozenset({"C"})


def test_snap_exact_when_clade_present():
    t = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    snapped, score = snap_to_backbone(frozenset({"D", "E"}), t)
    assert snapped == frozenset({"D", "E"})
    assert score == 1.0


def test_two_parent_events_allo_from_metadata():
    true_t = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    backbone = true_t.copy()
    md = [{"event_type": "allo", "target_clade": ["D"], "partner_clade": ["C"]}]
    events, scores = two_parent_events(md, true_t, backbone)
    assert len(events) == 1
    ev = events[0]
    assert ev.target_clade == frozenset({"D"})
    # parent A = D's true sibling (E); parent B = metadata partner (C)
    assert ev.parent_a_clade == frozenset({"E"})
    assert ev.parent_b_clade == frozenset({"C"})


def test_predicted_parents_coclust_returns_top2():
    # X sisters P1 in 2/3 trees, P2 in 1/3 -> top-2 = [P1, P2].
    gene_trees = [
        Tree("((X:1,P1:1):1,(P2:1,P3:1):1);", format=1),
        Tree("((X:1,P1:1):1,(P2:1,P3:1):1);", format=1),
        Tree("((X:1,P2:1):1,(P1:1,P3:1):1);", format=1),
    ]
    all_species = {"X", "P1", "P2", "P3"}
    top2 = predicted_parents_coclust(gene_trees, all_species, "X")
    assert top2[:2] == ["P1", "P2"]


def test_two_parent_events_coclust_mode_uses_predicted_parents():
    true_t = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    backbone = true_t.copy()
    md = [{"event_type": "allo", "target_clade": ["D"], "partner_clade": ["C"]}]
    # D sisters A in 2/3 trees, C in 1/3 -> predicted parents {A}, {C}
    gene_trees = [
        Tree("((A:1,D:1):1,(B:1,(C:1,E:1):1):1);", format=1),
        Tree("((A:1,D:1):1,(B:1,(C:1,E:1):1):1);", format=1),
        Tree("((C:1,D:1):1,(B:1,(A:1,E:1):1):1);", format=1),
    ]
    all_species = {"A", "B", "C", "D", "E"}
    events, scores = two_parent_events(md, true_t, backbone, mode="coclust",
                                       gene_trees=gene_trees, all_species=all_species)
    assert len(events) == 1
    ev = events[0]
    assert ev.target_clade == frozenset({"D"})
    # predicted top-2 parents snapped to backbone: {A} and {C}
    assert {ev.parent_a_clade, ev.parent_b_clade} == {frozenset({"A"}), frozenset({"C"})}
