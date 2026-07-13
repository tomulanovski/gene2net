"""Construct a MUL-tree from a species tree + predicted WGD events."""
from dataclasses import dataclass
from typing import FrozenSet, List, Optional
from ete3 import Tree, TreeNode


@dataclass
class WGDEvent:
    wgd_edge_clade: FrozenSet[str]     # leaf set of the clade being duplicated
    partner_edge_clade: FrozenSet[str]  # leaf set below the partner edge
    confidence: float


def _find_node_by_leaf_set(tree: Tree, leaf_set: FrozenSet[str]) -> Optional[TreeNode]:
    """Find the node whose descendant leaf set exactly matches the given set."""
    for node in tree.traverse("postorder"):
        if frozenset(node.get_leaf_names()) == leaf_set:
            return node
    return None


def _apply_wgd_event(tree: Tree, event: WGDEvent) -> bool:
    """Apply a single WGD event to the tree in-place.

    Returns True if the event was applied, False if the clade was not found
    (which can happen when prior events have already modified the tree).
    """
    wgd_node = _find_node_by_leaf_set(tree, event.wgd_edge_clade)
    if wgd_node is None:
        return False

    # A whole-tree clade is the root: duplicating it aliases `tree` into its own
    # child (self-cycle) and hangs any later traversal. Skip it. Real predicted
    # WGD edges are always non-root, so this only guards degenerate/hand-built input.
    if wgd_node.up is None:
        return False

    is_auto = event.partner_edge_clade == event.wgd_edge_clade

    # Deep-copy the subtree to be duplicated
    wgd_copy = wgd_node.copy("deepcopy")

    if is_auto:
        # Autopolyploidy: duplicate the clade under a new internal node
        parent = wgd_node.up
        if parent is None:
            # wgd_node is the root — wrap it
            new_internal = Tree()
            new_internal.add_child(wgd_node.detach())
            new_internal.add_child(wgd_copy)
            # Replace tree content (can't reassign root directly with ete3)
            # Copy new_internal structure into tree
            tree.__init__()
            for child in new_internal.get_children():
                tree.add_child(child.detach())
        else:
            branch_length = wgd_node.dist
            wgd_node.detach()
            new_internal = parent.add_child(dist=branch_length)
            new_internal.add_child(wgd_node)
            new_internal.add_child(wgd_copy)
    else:
        # Allopolyploidy: graft duplicate onto partner edge
        partner_node = _find_node_by_leaf_set(tree, event.partner_edge_clade)
        if partner_node is None:
            return False

        partner_parent = partner_node.up
        if partner_parent is None:
            return False

        partner_dist = partner_node.dist
        partner_node.detach()

        # New internal node replaces partner in its parent
        new_internal = partner_parent.add_child(dist=partner_dist)
        new_internal.add_child(partner_node, dist=0)
        new_internal.add_child(wgd_copy, dist=0)

    return True


@dataclass
class TwoParentEvent:
    target_clade: FrozenSet[str]      # the polyploid clade to place
    parent_a_clade: FrozenSet[str]    # first true parent (leaf set below its edge)
    parent_b_clade: FrozenSet[str]    # second true parent
    confidence: float


def _graft_copy_at(node: TreeNode, subtree_copy: TreeNode) -> None:
    """Insert a new internal node in place of `node` holding both `node` and the copy."""
    parent = node.up
    dist = node.dist
    node.detach()
    new_internal = parent.add_child(dist=dist)
    new_internal.add_child(node, dist=0)
    new_internal.add_child(subtree_copy, dist=0)


def _apply_two_parent_event(tree: Tree, event: TwoParentEvent) -> bool:
    """Detach the target clade and graft a copy at BOTH true parents.

    Auto (parent_a == parent_b == target) delegates to the existing sibling-
    duplication path. Returns False (drop) if any clade is not found, is the
    root, or a parent overlaps the target. When both parents resolve to the same
    backbone node, both copies land near that lineage (the honest "signal could
    not separate the two parents" collapse; ploidy is still 2).
    """
    target, a_clade, b_clade = event.target_clade, event.parent_a_clade, event.parent_b_clade

    if a_clade == b_clade == target:
        return _apply_wgd_event(tree, WGDEvent(target, target, event.confidence))

    x_node = _find_node_by_leaf_set(tree, target)
    if x_node is None or x_node.up is None:
        return False
    a_node = _find_node_by_leaf_set(tree, a_clade)
    b_node = _find_node_by_leaf_set(tree, b_clade)
    if a_node is None or b_node is None or a_node.up is None or b_node.up is None:
        return False
    if a_node is x_node or b_node is x_node:
        return False

    copy_a = x_node.copy("deepcopy")
    copy_b = x_node.copy("deepcopy")

    # Detach the target from its backbone position and collapse the unary parent.
    parent = x_node.up
    x_node.detach()
    if parent is not tree and len(parent.children) == 1:
        parent.delete(preserve_branch_length=True)

    # Re-find parents (leaf sets unchanged, but the tree object moved).
    a_node = _find_node_by_leaf_set(tree, a_clade)
    b_node = _find_node_by_leaf_set(tree, b_clade)
    if a_node is None or b_node is None or a_node.up is None or b_node.up is None:
        return False

    _graft_copy_at(a_node, copy_a)
    _graft_copy_at(b_node, copy_b)
    return True


def _apply_two_parent_grafted(tree: Tree, event: TwoParentEvent) -> bool:
    """Place a two-parent allopolyploid WITHOUT detaching the target.

    The target stays at its current position (its 'home' parent) and a copy is
    grafted at the OTHER parent (the one it is not already next to). This matches
    how the ground-truth MUL-tree is written (home copy kept, partner copy added),
    so nested events compose: an inner event never tears leaves out of an outer
    event's target clade. It still fixes the sp39 collapse, because the copy goes
    to the away parent, not on top of the home.
    """
    target, a_clade, b_clade = event.target_clade, event.parent_a_clade, event.parent_b_clade

    if a_clade == b_clade == target:
        return _apply_wgd_event(tree, WGDEvent(target, target, event.confidence))

    x_node = _find_node_by_leaf_set(tree, target)
    if x_node is None or x_node.up is None:
        return False

    # Current home = the leaves X currently sits next to.
    home = set()
    for ch in x_node.up.get_children():
        if ch is not x_node:
            home |= set(ch.get_leaf_names())

    # Graft at the parent LESS overlapping with the home (the 'away' parent);
    # keep X on the home side.
    away = b_clade if len(home & a_clade) >= len(home & b_clade) else a_clade
    away_node = _find_node_by_leaf_set(tree, away)
    if away_node is None or away_node.up is None or away_node is x_node:
        return False

    _graft_copy_at(away_node, x_node.copy("deepcopy"))
    return True


def build_mul_tree_two_parent(species_tree: Tree, events: List["TwoParentEvent"],
                              return_dropped: bool = False, mode: str = "graft"):
    """Build a MUL-tree by applying two-parent events bottom-up (smallest target first).

    mode='graft' (default): keep the target at its home, add a copy at the away
        parent. Composes for nested events; reproduces the ground-truth convention.
    mode='detach': detach the target and graft a copy at BOTH parents. Cleaner for
        a single isolated event but DROPS outer events whose target clade an inner
        event has already torn apart (see _apply_two_parent_grafted docstring).
    """
    apply = _apply_two_parent_grafted if mode == "graft" else _apply_two_parent_event
    mul_tree = species_tree.copy("deepcopy")
    sorted_events = sorted(events, key=lambda e: len(e.target_clade))
    dropped = 0
    for event in sorted_events:
        if not apply(mul_tree, event):
            dropped += 1
    if return_dropped:
        return mul_tree, dropped
    return mul_tree


def build_mul_tree(species_tree: Tree, events: List[WGDEvent], return_dropped: bool = False):
    """Build MUL-tree by applying WGD events to species tree.

    Events are sorted bottom-up (smallest clades first) so nested events
    are applied before outer events that may duplicate them.

    An event whose target/partner clade cannot be found in the current tree is
    dropped (this happens when an earlier allo graft added foreign leaves into
    that clade's ancestors). With return_dropped=True, returns (tree, n_dropped)
    so callers can see how many predicted events were lost; otherwise returns the
    tree (backward-compatible).
    """
    mul_tree = species_tree.copy("deepcopy")

    # Sort by clade size ascending (deepest events first)
    sorted_events = sorted(events, key=lambda e: len(e.wgd_edge_clade))

    dropped = 0
    for event in sorted_events:
        if not _apply_wgd_event(mul_tree, event):
            dropped += 1

    if return_dropped:
        return mul_tree, dropped
    return mul_tree
