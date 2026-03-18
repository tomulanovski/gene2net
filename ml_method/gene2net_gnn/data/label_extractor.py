"""Decompose a MUL-tree into WGD events and map them to ASTRAL tree edges."""
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

from ete3 import Tree, TreeNode

from gene2net_gnn.inference.mul_tree_builder import WGDEvent


@dataclass
class TrainingLabels:
    wgd_edges: List[int] = field(default_factory=list)       # ASTRAL edge indices with WGD
    partner_edges: List[int] = field(default_factory=list)    # corresponding partner edge indices
    wgd_counts: List[int] = field(default_factory=list)       # events per edge (usually 0 or 1)
    mask: List[bool] = field(default_factory=list)            # True = mappable
    n_unmappable: int = 0
    n_edges: int = 0                                          # total edges in ASTRAL tree


def extract_backbone(mul_tree: Tree) -> Tree:
    """Remove duplicate species from a MUL-tree, returning a backbone species tree.

    For each species appearing more than once, keep only the first copy
    encountered in a preorder traversal and prune the rest. Then resolve
    any degree-2 internal nodes.
    """
    tree = mul_tree.copy("deepcopy")

    # Find polyploid species
    counts: Dict[str, int] = Counter(leaf.name for leaf in tree.get_leaves())
    polyploids = {sp for sp, cnt in counts.items() if cnt > 1}

    if not polyploids:
        return tree

    # Keep track of which species we've already seen
    seen: Set[str] = set()
    to_remove: List[TreeNode] = []

    for leaf in tree.iter_leaves():
        if leaf.name in polyploids:
            if leaf.name in seen:
                to_remove.append(leaf)
            else:
                seen.add(leaf.name)

    # Remove extra copies
    for leaf in to_remove:
        parent = leaf.up
        leaf.detach()
        # Resolve degree-2 nodes up the tree
        _resolve_single_child_chain(parent)

    # Also resolve the root if needed
    _resolve_single_child_chain(tree)

    # Handle case where root itself became degree-1
    while len(tree.get_children()) == 1:
        child = tree.get_children()[0]
        child.detach()
        # Copy child's children to tree (replace root content)
        for grandchild in list(child.get_children()):
            grandchild.detach()
            tree.add_child(grandchild)
        if child.is_leaf():
            tree.name = child.name
            break

    return tree


def _resolve_single_child_chain(node: Optional[TreeNode]) -> None:
    """Walk up from node, resolving any internal nodes that have exactly one child."""
    while node is not None and not node.is_root():
        children = node.get_children()
        if len(children) == 1:
            # Bypass this node: connect child directly to parent
            child = children[0]
            parent = node.up
            dist = node.dist + child.dist
            node.detach()
            child.detach()
            parent.add_child(child, dist=dist)
            node = parent
        else:
            break


def decompose_mul_tree(mul_tree: Tree) -> List[WGDEvent]:
    """Decompose a MUL-tree into a list of WGDEvent objects.

    Detects both autopolyploidy (identical sibling subtrees) and
    allopolyploidy (duplicate copies in distant locations).
    """
    counts: Dict[str, int] = Counter(leaf.name for leaf in mul_tree.get_leaves())
    polyploids = {sp for sp, cnt in counts.items() if cnt > 1}

    if not polyploids:
        return []

    events: List[WGDEvent] = []
    explained_species: Set[str] = set()

    # Step 1: detect auto events (identical sibling subtrees)
    for node in mul_tree.traverse("postorder"):
        if node.is_leaf():
            continue
        children = node.get_children()
        if len(children) != 2:
            continue
        left_set = set(children[0].get_leaf_names())
        right_set = set(children[1].get_leaf_names())
        if left_set == right_set and left_set & polyploids:
            # Check these species aren't already explained
            if not (left_set & explained_species):
                clade = frozenset(left_set)
                events.append(WGDEvent(
                    wgd_edge_clade=clade,
                    partner_edge_clade=clade,
                    confidence=1.0,
                ))
                explained_species |= left_set

    # Step 2: detect allo events for remaining unexplained polyploids
    remaining = polyploids - explained_species
    if remaining:
        backbone = extract_backbone(mul_tree)
        backbone_locations = _get_backbone_locations(backbone)

        for species in sorted(remaining):  # sorted for determinism
            _detect_allo_event(mul_tree, backbone_locations, species, events)

    return events


def _get_backbone_locations(backbone: Tree) -> Dict[str, TreeNode]:
    """Map each leaf name to its node in the backbone tree."""
    return {leaf.name: leaf for leaf in backbone.get_leaves()}


def _detect_allo_event(
    mul_tree: Tree,
    backbone_locations: Dict[str, TreeNode],
    species: str,
    events: List[WGDEvent],
) -> None:
    """Detect an allopolyploidy event for a given polyploid species.

    The "primary" copy is the one whose neighborhood best matches the backbone.
    The "extra" copy's sibling leaf set defines the partner_edge_clade.
    """
    # Find all copies of this species in the MUL-tree
    copies = [leaf for leaf in mul_tree.get_leaves() if leaf.name == species]
    if len(copies) < 2:
        return

    # Determine which copy is primary vs extra using backbone context
    # The primary copy is the one whose sibling subtree leaves overlap more
    # with what we'd expect from the backbone
    backbone_node = backbone_locations.get(species)
    if backbone_node is None:
        return

    # Get expected sibling species from backbone
    if backbone_node.up is not None:
        expected_siblings = set()
        for child in backbone_node.up.get_children():
            if child is not backbone_node:
                expected_siblings |= set(child.get_leaf_names())
    else:
        expected_siblings = set()

    # Score each copy: how well does its local neighborhood match the backbone?
    best_score = -1
    primary_copy = copies[0]
    for copy_node in copies:
        if copy_node.up is None:
            continue
        sibling_leaves = set()
        for child in copy_node.up.get_children():
            if child is not copy_node:
                sibling_leaves |= set(child.get_leaf_names())
        # Remove the species itself from sibling leaves (may appear multiple times)
        sibling_leaves.discard(species)
        score = len(sibling_leaves & expected_siblings)
        if score > best_score:
            best_score = score
            primary_copy = copy_node

    # The extra copies define events
    for copy_node in copies:
        if copy_node is primary_copy:
            continue
        # The sibling subtree of the extra copy defines the partner clade
        if copy_node.up is not None:
            sibling_leaves: Set[str] = set()
            for child in copy_node.up.get_children():
                if child is not copy_node:
                    sibling_leaves |= set(child.get_leaf_names())
            # Remove duplicate species names for clean clade definition
            partner_clade = frozenset(sibling_leaves) if sibling_leaves else frozenset({species})
        else:
            partner_clade = frozenset({species})

        wgd_clade = frozenset({species})
        events.append(WGDEvent(
            wgd_edge_clade=wgd_clade,
            partner_edge_clade=partner_clade,
            confidence=1.0,
        ))


def _jaccard(a: FrozenSet[str], b: FrozenSet[str]) -> float:
    """Compute Jaccard similarity between two sets."""
    if not a and not b:
        return 1.0
    intersection = len(a & b)
    union = len(a | b)
    return intersection / union if union > 0 else 0.0


def _get_edge_bipartitions(tree: Tree) -> List[Tuple[int, FrozenSet[str]]]:
    """Get bipartitions for each edge in the tree (preorder traversal).

    Each non-root node defines an edge to its parent. The bipartition
    is the leaf set below that node.

    Returns list of (edge_index, leaf_set_below).
    """
    all_leaves = frozenset(tree.get_leaf_names())
    edges: List[Tuple[int, FrozenSet[str]]] = []
    idx = 0
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        below = frozenset(node.get_leaf_names())
        edges.append((idx, below))
        idx += 1
    return edges


def _best_matching_edge(
    target_clade: FrozenSet[str],
    edges: List[Tuple[int, FrozenSet[str]]],
) -> Tuple[int, float]:
    """Find edge with highest Jaccard similarity to target_clade.

    Returns (edge_index, best_jaccard_score).
    """
    best_idx = -1
    best_score = -1.0
    for idx, below in edges:
        score = _jaccard(target_clade, below)
        if score > best_score:
            best_score = score
            best_idx = idx
    return best_idx, best_score


def map_events_to_astral(
    events: List[WGDEvent],
    astral_tree: Tree,
    jaccard_threshold: float = 0.5,
) -> TrainingLabels:
    """Map WGD events to edges of an ASTRAL species tree using bipartition matching.

    Parameters
    ----------
    events : list of WGDEvent
        WGD events extracted from a MUL-tree.
    astral_tree : Tree
        The ASTRAL-inferred species tree.
    jaccard_threshold : float
        Minimum Jaccard similarity for a match to be considered valid.

    Returns
    -------
    TrainingLabels
        Edge-level labels for GNN training.
    """
    edges = _get_edge_bipartitions(astral_tree)
    n_edges = len(edges)

    # Initialize per-edge counts
    edge_wgd_counts = [0] * n_edges

    wgd_edges: List[int] = []
    partner_edges: List[int] = []
    mask: List[bool] = []
    n_unmappable = 0

    for event in events:
        wgd_idx, wgd_score = _best_matching_edge(event.wgd_edge_clade, edges)
        partner_idx, partner_score = _best_matching_edge(event.partner_edge_clade, edges)

        mappable = wgd_score >= jaccard_threshold and partner_score >= jaccard_threshold

        if mappable:
            wgd_edges.append(wgd_idx)
            partner_edges.append(partner_idx)
            mask.append(True)
            edge_wgd_counts[wgd_idx] += 1
        else:
            wgd_edges.append(wgd_idx)
            partner_edges.append(partner_idx)
            mask.append(False)
            n_unmappable += 1

    return TrainingLabels(
        wgd_edges=wgd_edges,
        partner_edges=partner_edges,
        wgd_counts=edge_wgd_counts,
        mask=mask,
        n_unmappable=n_unmappable,
        n_edges=n_edges,
    )
