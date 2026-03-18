"""MUL-tree generator: create random species trees and add WGD events."""
import random
import copy
from typing import Dict, List, Optional, Set, Tuple

from ete3 import Tree


def generate_random_species_tree(n_species: int, seed: Optional[int] = None) -> Tree:
    """Create a random binary tree with n_species unique leaf labels (sp_0, sp_1, ...).

    Uses ETE3's populate() to generate the random topology, then renames leaves.
    NOTE: This tree is NOT ultrametric. For ultrametric trees, use
    generate_birth_death_tree() instead.
    """
    rng = random.Random(seed)
    tree = Tree()
    tree.populate(n_species, random_branches=True, branch_range=(0.5, 2.0))
    leaves = tree.get_leaves()
    names = [f"sp_{i}" for i in range(n_species)]
    rng.shuffle(names)
    for leaf, name in zip(leaves, names):
        leaf.name = name
    return tree


def generate_birth_death_tree(
    n_species: int,
    birth_rate: float = 1.0,
    death_rate: float = 0.5,
    tree_height: float = 10_000_000.0,
    seed: Optional[int] = None,
) -> Tree:
    """Generate an ultrametric species tree under a birth-death process.

    Simulates a forward-time birth-death process until exactly n_species extant
    lineages exist, then assigns ultrametric branch lengths and rescales so the
    root-to-tip distance equals tree_height (in generations).

    Parameters
    ----------
    n_species : int
        Target number of extant species (leaves).
    birth_rate : float
        Speciation rate (lambda). Must be > death_rate for net growth.
    death_rate : float
        Extinction rate (mu). Set to 0 for pure-birth (Yule) process.
    tree_height : float
        Total root-to-tip distance in generations (default 10M, matching
        existing SimPhy network files).
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    Tree
        An ultrametric ETE3 Tree with n_species leaves named sp_0..sp_{n-1}.
    """
    rng = random.Random(seed)

    # Forward-time birth-death: simulate until exactly n_species extant
    max_attempts = 500
    for _ in range(max_attempts):
        lineages = [0]  # list of node IDs (just track topology)
        tree_splits = {}  # parent_id -> (child1_id, child2_id, split_time)
        next_id = 1
        t = 0.0

        while 0 < len(lineages) < n_species:
            n_alive = len(lineages)
            total_rate = n_alive * (birth_rate + death_rate)
            dt = rng.expovariate(total_rate)
            t += dt

            # Pick a random lineage
            idx = rng.randint(0, n_alive - 1)
            node_id = lineages[idx]

            if rng.random() < birth_rate / (birth_rate + death_rate):
                # Speciation
                c1, c2 = next_id, next_id + 1
                next_id += 2
                tree_splits[node_id] = (c1, c2, t)
                lineages[idx] = c1
                lineages.append(c2)
            else:
                # Extinction
                lineages.pop(idx)

        if len(lineages) == n_species:
            break
    else:
        # Fallback to Yule if birth-death keeps failing
        return _generate_yule_tree(n_species, tree_height, seed)

    # Build full tree, then prune extinct lineages
    extant = set(lineages)
    total_time = t
    root = _build_tree_from_splits(0, tree_splits, total_time)

    # Remove extinct leaves and clean up unary internal nodes
    root = _prune_extinct(root, extant)

    # Rescale so root-to-tip = tree_height
    current_height = _subtree_height(root)
    if current_height > 0:
        scale = tree_height / current_height
        for node in root.traverse():
            if node.up is not None:
                node.dist *= scale

    _force_ultrametric(root, tree_height)

    # Name leaves
    leaves = root.get_leaves()
    names = [f"sp_{i}" for i in range(len(leaves))]
    rng.shuffle(names)
    for leaf, name in zip(leaves, names):
        leaf.name = name

    return root


def _prune_extinct(root, extant_ids):
    """Remove extinct leaves from tree and clean up unary internal nodes.

    Each leaf in the tree built from splits has an implicit ID. We label leaves
    with sequential IDs during tree building, then prune those not in extant_ids.
    """
    # Label all leaves with their node IDs for identification
    # Leaves in the built tree correspond to node_ids NOT in tree_splits
    # We need to identify which leaves to keep. Since we don't store IDs on nodes,
    # re-traverse and tag them.

    # Actually, we can identify extant leaves by their position in the tree.
    # Simpler approach: tag leaves during build. Let's use node.name as the ID.
    # The build function already creates unnamed leaves. Let's fix this differently:
    # just remove leaves whose name (set to node_id during build) is not in extant.

    # Remove extinct leaves iteratively until only extant remain
    changed = True
    while changed:
        changed = False
        for leaf in root.get_leaves():
            leaf_id = int(leaf.name) if leaf.name.isdigit() else -1
            if leaf_id not in extant_ids:
                parent = leaf.up
                leaf.detach()
                changed = True
                # Clean up unary parent
                if parent and len(parent.children) == 1:
                    child = parent.children[0]
                    child.dist += parent.dist  # merge branch lengths
                    grandparent = parent.up
                    if grandparent:
                        parent.detach()
                        grandparent.add_child(child)
                    else:
                        # parent was root — child becomes new root
                        child.detach()
                        root = child
                elif parent and len(parent.children) == 0:
                    # parent became a leaf too — will be handled next iteration
                    pass
                break  # restart iteration after tree modification

    return root


def _build_tree_from_splits(node_id, splits, total_time, birth_time=0.0):
    """Recursively build ETE3 tree from birth-death split records."""
    node = Tree()
    node.name = str(node_id)  # tag with ID for pruning

    if node_id in splits:
        c1_id, c2_id, split_time = splits[node_id]
        c1 = _build_tree_from_splits(c1_id, splits, total_time, split_time)
        c2 = _build_tree_from_splits(c2_id, splits, total_time, split_time)
        node.add_child(c1)
        node.add_child(c2)
        node.dist = split_time - birth_time
    else:
        # Leaf: lives from birth_time to total_time
        node.dist = total_time - birth_time

    return node


def _generate_yule_tree(n_species: int, tree_height: float, seed: Optional[int] = None) -> Tree:
    """Fallback: pure-birth (Yule) tree. Always succeeds."""
    rng = random.Random(seed)
    root = Tree()
    leaves = [root]

    while len(leaves) < n_species:
        idx = rng.randint(0, len(leaves) - 1)
        node = leaves.pop(idx)
        c1 = node.add_child(dist=0.0)
        c2 = node.add_child(dist=0.0)
        leaves.append(c1)
        leaves.append(c2)

    _assign_ultrametric_times(root, rng)

    current_height = _subtree_height(root)
    if current_height > 0:
        scale = tree_height / current_height
        for node in root.traverse():
            if node.up is not None:
                node.dist *= scale

    _force_ultrametric(root, tree_height)

    leaf_nodes = root.get_leaves()
    names = [f"sp_{i}" for i in range(len(leaf_nodes))]
    rng.shuffle(names)
    for leaf, name in zip(leaf_nodes, names):
        leaf.name = name

    return root


def _force_ultrametric(root, target_height: float) -> None:
    """Adjust leaf branch lengths so all root-to-leaf distances equal target_height exactly.

    Walks backwards from each leaf to root, sums internal branch lengths,
    then sets leaf.dist = target - sum. Same approach as
    simulations/scripts/rescale_and_keep_ultrametric.py.
    """
    for leaf in root.iter_leaves():
        dist_to_parent = 0.0
        node = leaf.up
        while node is not None and not node.is_root():
            dist_to_parent += node.dist
            node = node.up
        leaf.dist = max(0.0, round(target_height - dist_to_parent, 12))


def _assign_ultrametric_times(root, rng) -> None:
    """Assign random ultrametric branch lengths to a topology (in-place).

    Uses a coalescent-like approach: assign random internal node times from
    root (time=1.0) to tips (time=0.0), then set dist = parent_time - child_time.
    """
    # Collect internal nodes (non-leaves, non-root)
    internals = [n for n in root.traverse("preorder") if not n.is_leaf() and n != root]

    # Assign times: root=1.0, leaves=0.0, internals=random in between
    # Sort internal nodes by depth to ensure parent time > child time
    node_times = {root: 1.0}
    for leaf in root.get_leaves():
        node_times[leaf] = 0.0

    # Assign times top-down: each internal node gets a time between
    # its parent's time and 0 (leaf time)
    for node in root.traverse("preorder"):
        if node == root or node.is_leaf():
            continue
        parent_time = node_times[node.up]
        # Random fraction of parent's remaining time
        node_times[node] = parent_time * rng.uniform(0.1, 0.9)

    # Set branch lengths
    for node in root.traverse("preorder"):
        if node == root:
            node.dist = 0.0
        else:
            node.dist = node_times[node.up] - node_times[node]


def _subtree_height(node) -> float:
    """Return the height of a subtree (max distance from node to any leaf)."""
    if node.is_leaf():
        return 0.0
    return max(child.dist + _subtree_height(child) for child in node.children)


def _find_node_by_leafset(tree: Tree, target_clade: Set[str]) -> Optional[object]:
    """Return the node whose descendant leaf set exactly matches target_clade."""
    for node in tree.traverse("postorder"):
        if set(node.get_leaf_names()) == target_clade:
            return node
    return None


def add_wgd_event(
    tree: Tree,
    target_clade: Set[str],
    partner_clade: Optional[Set[str]] = None,
    event_type: str = "allo",
) -> Tree:
    """Add a WGD event to a deep copy of tree. Never mutates the input.

    Parameters
    ----------
    tree : Tree
        Input species tree (not modified).
    target_clade : set of str
        Leaf names defining the clade to duplicate.
    partner_clade : set of str or None
        Leaf names defining the edge where the duplicate attaches (allo only).
    event_type : {"auto", "allo"}
        "auto" → duplicate placed as sister of original (same parent).
        "allo" → duplicate attached to the partner edge.

    Returns
    -------
    Tree
        A new tree with the WGD event applied.
    """
    tree = tree.copy("deepcopy")

    target_node = _find_node_by_leafset(tree, target_clade)
    if target_node is None:
        raise ValueError(f"No node found with leaf set {target_clade}")

    # Deep-copy the subtree to duplicate
    dup = target_node.copy("deepcopy")

    if event_type == "auto":
        # Place duplicate as sister of target_node under the same parent.
        # Split the edge above target_node at a fraction f to maintain ultrametricity.
        # parent --[d*f]--> new_internal --[d*(1-f)]--> target_node
        #                                --[d*(1-f)]--> dup
        parent = target_node.up
        if parent is None:
            # target is the root – make a new root with a small stem
            new_root = Tree()
            new_root.dist = 0.0
            stem = 0.1  # small stem above old root
            target_node.detach()
            new_root.add_child(target_node, dist=stem)
            new_root.add_child(dup, dist=stem)
            return new_root
        else:
            original_dist = target_node.dist
            # Split at midpoint (f=0.5) for determinism; could randomize
            split_above = original_dist * 0.5   # parent → new_internal
            split_below = original_dist * 0.5   # new_internal → each child
            target_node.detach()
            new_internal = parent.add_child(dist=split_above)
            new_internal.add_child(target_node, dist=split_below)
            new_internal.add_child(dup, dist=split_below)
            return tree

    elif event_type == "allo":
        if partner_clade is None:
            raise ValueError("partner_clade must be provided for allo events")

        partner_node = _find_node_by_leafset(tree, partner_clade)
        if partner_node is None:
            raise ValueError(f"No node found with leaf set {partner_clade}")

        # Subdivide the edge above partner_node and attach dup with correct stem
        # to maintain ultrametricity.
        #
        # Split partner edge at fraction f:
        #   partner_parent --[d*f]--> new_internal --[d*(1-f)]--> partner_node
        #                                          --[stem]-----> dup
        #
        # For ultrametricity: stem = d*(1-f) + h_partner - h_target
        # We need stem > 0, so: f < 1 - (h_target - h_partner) / d
        # If h_target > d + h_partner, placement is impossible (raise ValueError).
        partner_parent = partner_node.up
        h_target = _subtree_height(dup)
        h_partner = _subtree_height(partner_node)

        if partner_parent is None:
            # partner is root – create new root
            new_root = Tree()
            new_root.dist = 0.0
            # stem for partner side; dup_stem computed for ultrametricity
            stem = max(h_target - h_partner, 0.1)
            dup_stem = stem + h_partner - h_target
            if dup_stem <= 0:
                dup_stem = stem  # symmetric if heights match
            partner_node.detach()
            new_root.add_child(partner_node, dist=stem)
            new_root.add_child(dup, dist=dup_stem)
            return new_root
        else:
            partner_dist = partner_node.dist

            # Check feasibility: h_target must be < partner_dist + h_partner
            if h_target >= partner_dist + h_partner:
                raise ValueError(
                    f"Cannot maintain ultrametricity: target subtree too tall "
                    f"(h_target={h_target:.1f}) for partner edge "
                    f"(d={partner_dist:.1f}, h_partner={h_partner:.1f})"
                )

            # Choose split fraction so dup_stem is comfortably positive
            # dup_stem = d*(1-f) + h_partner - h_target
            # For dup_stem > 0: f < 1 - (h_target - h_partner) / d
            max_f = 1.0 - (h_target - h_partner) / partner_dist if h_target > h_partner else 1.0
            f = min(0.5, max_f * 0.9)  # use midpoint or 90% of max, whichever is smaller
            f = max(f, 0.01)  # ensure we split at least a little

            split_above = partner_dist * f
            split_below = partner_dist * (1.0 - f)
            dup_stem = split_below + h_partner - h_target

            partner_node.detach()
            new_internal = partner_parent.add_child(dist=split_above)
            new_internal.add_child(partner_node, dist=split_below)
            new_internal.add_child(dup, dist=dup_stem)
            return tree

    else:
        raise ValueError(f"Unknown event_type: {event_type!r}. Use 'auto' or 'allo'.")


def generate_random_mul_tree(
    n_species: int,
    n_events: int,
    seed: Optional[int] = None,
) -> Tuple[Tree, List[Dict]]:
    """Generate a random species tree and apply n_events random WGD events.

    Parameters
    ----------
    n_species : int
        Number of species in the base species tree.
    n_events : int
        Number of WGD events to apply sequentially.
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    mul_tree : Tree
        The resulting MUL-tree.
    events : list of dict
        Each dict has keys: target_clade, partner_clade, event_type.
    """
    rng = random.Random(seed)
    tree = generate_random_species_tree(n_species=n_species, seed=seed)
    events = []

    for _ in range(n_events):
        # Collect all internal nodes (non-root) as candidate targets
        # Exclude trivial single-leaf clades only occasionally; allow them for simplicity
        nodes = list(tree.traverse("postorder"))
        # exclude root to ensure every node has a parent for partner selection
        non_root_nodes = [n for n in nodes if n.up is not None]
        if not non_root_nodes:
            break

        target_node = rng.choice(non_root_nodes)
        target_clade = set(target_node.get_leaf_names())

        event_type = rng.choice(["auto", "allo"])
        partner_clade = None

        if event_type == "allo":
            # Pick a partner node that doesn't overlap with target
            candidate_partners = [
                n for n in non_root_nodes
                if n is not target_node and not (set(n.get_leaf_names()) & target_clade)
            ]
            if not candidate_partners:
                # Fall back to auto if no valid partner exists
                event_type = "auto"
            else:
                partner_node = rng.choice(candidate_partners)
                partner_clade = set(partner_node.get_leaf_names())

        tree = add_wgd_event(
            tree,
            target_clade=target_clade,
            partner_clade=partner_clade,
            event_type=event_type,
        )

        events.append(
            {
                "target_clade": target_clade,
                "partner_clade": partner_clade,
                "event_type": event_type,
            }
        )

    return tree, events


def get_polyploid_species(tree: Tree) -> Dict[str, int]:
    """Count occurrences of each leaf name. Return {species: count} for count > 1.

    Parameters
    ----------
    tree : Tree
        A MUL-tree (or any ETE3 tree).

    Returns
    -------
    dict
        Mapping from species name to the number of times it appears as a leaf,
        restricted to species appearing more than once.
    """
    counts: Dict[str, int] = {}
    for leaf in tree.get_leaves():
        counts[leaf.name] = counts.get(leaf.name, 0) + 1
    return {sp: cnt for sp, cnt in counts.items() if cnt > 1}
