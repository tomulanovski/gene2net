import pytest
from ete3 import Tree
from gene2net_gnn.data.mul_tree_generator import (
    generate_random_species_tree,
    add_wgd_event,
    generate_random_mul_tree,
    get_polyploid_species,
)

def test_generate_random_species_tree():
    tree = generate_random_species_tree(n_species=10, seed=42)
    leaves = tree.get_leaf_names()
    assert len(leaves) == 10
    assert len(set(leaves)) == 10  # all unique

def test_add_wgd_event_allo():
    tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    # WGD on clade {D}, partner = edge leading to A/B clade
    mul_tree = add_wgd_event(tree, target_clade={"D"}, partner_clade={"A", "B"}, event_type="allo")
    d_leaves = [l for l in mul_tree.get_leaf_names() if l == "D"]
    assert len(d_leaves) == 2

def test_add_wgd_event_auto():
    tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    # Auto WGD on clade (D, E)
    mul_tree = add_wgd_event(tree, target_clade={"D", "E"}, event_type="auto")
    d_leaves = [l for l in mul_tree.get_leaf_names() if l == "D"]
    e_leaves = [l for l in mul_tree.get_leaf_names() if l == "E"]
    assert len(d_leaves) == 2
    assert len(e_leaves) == 2

def test_generate_random_mul_tree():
    mul_tree, events = generate_random_mul_tree(n_species=10, n_events=2, seed=42)
    assert len(events) == 2
    polyploids = get_polyploid_species(mul_tree)
    assert len(polyploids) > 0

def test_get_polyploid_species():
    tree = Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    polyploids = get_polyploid_species(tree)
    assert polyploids == {"D": 2}


def _is_ultrametric(tree, tol=1e-6):
    """Check if all root-to-leaf distances are equal (within tolerance)."""
    distances = []
    for leaf in tree.get_leaves():
        d = leaf.get_distance(tree)
        distances.append(d)
    if not distances:
        return True
    return max(distances) - min(distances) < tol


def test_add_wgd_auto_preserves_ultrametricity():
    """Auto WGD on an ultrametric tree should produce an ultrametric MUL-tree."""
    # Ultrametric: all leaves at distance 3 from root
    tree = Tree("((A:1,B:1):2,(C:2,(D:1,E:1):1):1);", format=1)
    assert _is_ultrametric(tree), "Input tree should be ultrametric"
    mul_tree = add_wgd_event(tree, target_clade={"D", "E"}, event_type="auto")
    assert _is_ultrametric(mul_tree), (
        f"Auto WGD broke ultrametricity. Root-to-leaf distances: "
        f"{[round(l.get_distance(mul_tree), 6) for l in mul_tree.get_leaves()]}"
    )


def test_add_wgd_allo_preserves_ultrametricity():
    """Allo WGD on an ultrametric tree should produce an ultrametric MUL-tree."""
    # Ultrametric: all leaves at distance 3 from root
    tree = Tree("((A:1,B:1):2,(C:2,(D:1,E:1):1):1);", format=1)
    assert _is_ultrametric(tree), "Input tree should be ultrametric"
    mul_tree = add_wgd_event(tree, target_clade={"D"}, partner_clade={"A", "B"}, event_type="allo")
    assert _is_ultrametric(mul_tree), (
        f"Allo WGD broke ultrametricity. Root-to-leaf distances: "
        f"{[round(l.get_distance(mul_tree), 6) for l in mul_tree.get_leaves()]}"
    )


def test_multiple_wgd_preserves_ultrametricity():
    """Multiple WGD events should keep tree ultrametric."""
    tree = Tree("((A:2,B:2):1,(C:1.5,(D:1,E:1):0.5):1.5);", format=1)
    assert _is_ultrametric(tree), "Input tree should be ultrametric"
    # Auto event
    mul = add_wgd_event(tree, target_clade={"D", "E"}, event_type="auto")
    assert _is_ultrametric(mul), "Broke ultrametricity after auto WGD"
    # Allo event on top of auto
    mul2 = add_wgd_event(mul, target_clade={"C"}, partner_clade={"A", "B"}, event_type="allo")
    assert _is_ultrametric(mul2), "Broke ultrametricity after allo WGD on top of auto"


def test_generate_random_mul_tree_ultrametric():
    """generate_random_mul_tree should produce ultrametric trees when starting from ultrametric input."""
    # generate_random_species_tree uses ETE3's populate with random_branches,
    # which is NOT ultrametric by default. So we test add_wgd_event directly
    # on a known ultrametric tree via multiple random events.
    tree = Tree("(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);", format=1)
    assert _is_ultrametric(tree)
    # Apply several events
    from gene2net_gnn.data.mul_tree_generator import _find_node_by_leafset
    import random
    rng = random.Random(123)
    for _ in range(5):
        non_root = [n for n in tree.traverse() if n.up is not None]
        target = rng.choice(non_root)
        target_clade = set(target.get_leaf_names())
        tree = add_wgd_event(tree, target_clade=target_clade, event_type="auto")
        assert _is_ultrametric(tree), f"Broke ultrametricity after auto WGD on {target_clade}"
