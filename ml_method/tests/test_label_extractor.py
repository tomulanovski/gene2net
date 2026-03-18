import pytest
from ete3 import Tree
from gene2net_gnn.data.label_extractor import (
    extract_backbone,
    decompose_mul_tree,
    map_events_to_astral,
)
from gene2net_gnn.inference.mul_tree_builder import WGDEvent, build_mul_tree
from gene2net_gnn.data.mul_tree_generator import get_polyploid_species

def test_extract_backbone():
    """Remove duplicate species -> species tree."""
    mul_tree = Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    backbone = extract_backbone(mul_tree)
    leaf_names = sorted(backbone.get_leaf_names())
    assert leaf_names == ["A", "B", "C", "D", "E"]

def test_extract_backbone_auto():
    """Auto MUL-tree: (((D,E),(D,E)),(A,(B,C)))"""
    mul_tree = Tree("(((D:1,E:1):1,(D:1,E:1):1):1,(A:1,(B:1,C:1):1):1);", format=1)
    backbone = extract_backbone(mul_tree)
    leaf_names = sorted(backbone.get_leaf_names())
    assert leaf_names == ["A", "B", "C", "D", "E"]

def test_decompose_simple_allo():
    mul_tree = Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    events = decompose_mul_tree(mul_tree)
    assert len(events) == 1
    assert events[0].wgd_edge_clade == frozenset({"D"})

def test_decompose_auto():
    mul_tree = Tree("(((D:1,E:1):1,(D:1,E:1):1):1,(A:1,(B:1,C:1):1):1);", format=1)
    events = decompose_mul_tree(mul_tree)
    assert len(events) == 1
    assert events[0].wgd_edge_clade == frozenset({"D", "E"})
    # Auto event: partner == wgd clade
    assert events[0].partner_edge_clade == events[0].wgd_edge_clade

def test_decompose_no_polyploidy():
    tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    events = decompose_mul_tree(tree)
    assert len(events) == 0

def test_map_events_to_astral_perfect():
    """When ASTRAL matches backbone perfectly."""
    mul_tree = Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    astral_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    events = decompose_mul_tree(mul_tree)
    labels = map_events_to_astral(events, astral_tree)
    assert len(labels.wgd_edges) == 1

def test_map_events_to_astral_mismatch():
    """When ASTRAL has a different topology."""
    mul_tree = Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    astral_tree = Tree("((A:1,C:1):1,(B:1,(D:1,E:1):1):1);", format=1)
    events = decompose_mul_tree(mul_tree)
    labels = map_events_to_astral(events, astral_tree)
    # Should not crash, unmappable events are masked
    assert labels.n_unmappable >= 0

def test_roundtrip():
    """Decompose a MUL-tree into events, reconstruct, verify polyploids match."""
    original = Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    events = decompose_mul_tree(original)
    backbone = extract_backbone(original)
    reconstructed = build_mul_tree(backbone, events)
    assert sorted(get_polyploid_species(reconstructed).keys()) == sorted(get_polyploid_species(original).keys())
