"""MUL-tree build strategies from per-edge WGD + partner predictions.

The model gives, per species-tree edge: a WGD probability and a partner
distribution. Turning those into a MUL-tree involves choosing which edges become
events. Two strategies are kept:

  cap          - flagged edges (prob >= threshold), but cap duplications per
                 species at the inferred copy bound (fixes genuine over-ploidy).
                 Confidence-ordered.
  bound_driven - no threshold: take edges in confidence order until each species
                 reaches its copy bound. Removes the threshold knob. (default)
  detect       - detection-driven, multiset as a SOFT prior (a floor, not a cap):
                 fill to the copy bound, then ADD any edge the detection head is
                 confident about (prob >= threshold) beyond that bound. Recovers
                 events fractionation deleted from the copy count, which
                 bound_driven and cap cannot (they never exceed the bound).

Event selection returns edge indices (preorder, non-root). Partners are computed
afterward by the caller (needs the model + embeddings).
"""
import collections
from typing import Dict, FrozenSet, List, Optional


def build_parent_edge_map(tree) -> List[Optional[int]]:
    """parent_edge[i] = the edge index above edge i's parent node, or None if
    that parent is the root. Edges are preorder, non-root (matching the model)."""
    node_to_edge = {}
    order = []
    idx = 0
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        node_to_edge[id(node)] = idx
        order.append(node)
        idx += 1
    parent_edge: List[Optional[int]] = []
    for node in order:
        p = node.up
        if p is None or p.is_root():
            parent_edge.append(None)
        else:
            parent_edge.append(node_to_edge[id(p)])
    return parent_edge


def infer_copy_bound(gene_trees) -> Dict[str, int]:
    """Per-species consensus copy number: the largest k such that >= half the
    gene trees contain >= k copies of the species. (Same idea as Polyphest's
    consensus multiset.) Minimum 1."""
    n = len(gene_trees)
    per_tree = []
    species = set()
    for gt in gene_trees:
        c = collections.Counter(leaf.name for leaf in gt.get_leaves())
        per_tree.append(c)
        species.update(c.keys())
    bound = {}
    for s in species:
        max_k = max((c.get(s, 0) for c in per_tree), default=0)
        b = 1
        for k in range(max_k, 0, -1):
            if sum(1 for c in per_tree if c.get(s, 0) >= k) >= n / 2:
                b = k
                break
        bound[s] = max(b, 1)
    return bound


def _cap(candidates: List[int], wgd_probs, clades: List[FrozenSet[str]],
         copy_bound: Dict[str, int]) -> List[int]:
    """Greedily select candidate edges in descending confidence, skipping any
    that would push a species past its copy bound."""
    copies = collections.defaultdict(lambda: 1)
    selected = []
    for i in sorted(candidates, key=lambda k: -float(wgd_probs[k])):
        clade = clades[i]
        if all(copies[s] + 1 <= copy_bound.get(s, 1) for s in clade):
            selected.append(i)
            for s in clade:
                copies[s] += 1
    return selected


def select_event_edges(
    strategy: str,
    wgd_probs,
    threshold: float,
    parent_edge: List[Optional[int]],
    clades: List[FrozenSet[str]],
    copy_bound: Optional[Dict[str, int]] = None,
) -> List[int]:
    """Return the edge indices that become WGD events under the given strategy."""
    n_edges = len(clades)

    if strategy == "cap":
        if copy_bound is None:
            raise ValueError("cap strategy requires copy_bound")
        flagged = [i for i in range(n_edges) if float(wgd_probs[i]) >= threshold]
        return _cap(flagged, wgd_probs, clades, copy_bound)

    if strategy == "bound_driven":
        if copy_bound is None:
            raise ValueError("bound_driven strategy requires copy_bound")
        return _cap(list(range(n_edges)), wgd_probs, clades, copy_bound)

    if strategy == "detect":
        # Multiset as a SOFT prior: fill to the copy bound (the baseline ploidy),
        # then ADD any edge the detection head is confident about beyond that bound.
        # threshold controls how sure a detection must be to exceed the ploidy
        # estimate; tune it on the validation split. Under fractionation the bound
        # collapses so `base` is small and the confident detections do the work,
        # recovering events a ploidy-only method (Polyphest) structurally cannot.
        if copy_bound is None:
            raise ValueError("detect strategy requires copy_bound")
        base = _cap(list(range(n_edges)), wgd_probs, clades, copy_bound)
        base_set = set(base)
        extra = [i for i in range(n_edges)
                 if float(wgd_probs[i]) >= threshold and i not in base_set]
        return base + extra

    raise ValueError(
        f"Unknown strategy: {strategy!r}. Supported: 'bound_driven', 'cap', 'detect'."
    )
