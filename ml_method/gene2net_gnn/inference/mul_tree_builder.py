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


def build_mul_tree(species_tree: Tree, events: List[WGDEvent]) -> Tree:
    """Build MUL-tree by applying WGD events to species tree.

    Events are sorted bottom-up (smallest clades first) so nested events
    are applied before outer events that may duplicate them.
    """
    mul_tree = species_tree.copy("deepcopy")

    # Sort by clade size ascending (deepest events first)
    sorted_events = sorted(events, key=lambda e: len(e.wgd_edge_clade))

    for event in sorted_events:
        _apply_wgd_event(mul_tree, event)

    return mul_tree
