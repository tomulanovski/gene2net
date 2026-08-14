import torch
from gene2net_gnn.inference.predict import events_from_two_parent_scores


def test_decode_picks_argmax_parents():
    E = 4
    scores2 = torch.full((E, E, 2), -5.0)
    scores2[0, 1, 0] = 9.0   # edge0 parent_a -> edge1
    scores2[0, 2, 1] = 9.0   # edge0 parent_b -> edge2
    clades = [frozenset({"X"}), frozenset({"A"}), frozenset({"B"}), frozenset({"C"})]
    events = events_from_two_parent_scores(scores2, [0], clades)
    assert len(events) == 1
    e = events[0]
    assert e.target_clade == frozenset({"X"})
    assert {e.parent_a_clade, e.parent_b_clade} == {frozenset({"A"}), frozenset({"B"})}


def test_decode_auto_both_slots_same_edge():
    E = 3
    scores2 = torch.full((E, E, 2), -5.0)
    scores2[1, 1, 0] = 9.0   # both slots point to self -> autopolyploidy
    scores2[1, 1, 1] = 9.0
    clades = [frozenset({"A"}), frozenset({"B"}), frozenset({"C"})]
    events = events_from_two_parent_scores(scores2, [1], clades)
    e = events[0]
    assert e.target_clade == e.parent_a_clade == e.parent_b_clade == frozenset({"B"})
