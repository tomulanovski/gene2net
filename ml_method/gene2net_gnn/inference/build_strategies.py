"""MUL-tree build strategies from per-edge WGD + partner predictions.

The model gives, per species-tree edge: a WGD probability and a partner
distribution. Turning those into a MUL-tree involves choices about which edges
become events. We expose several strategies so they can be validated against the
reconstruction measures:

  raw          - every flagged edge (prob >= threshold) is its own event (baseline).
  collapse     - one event per connected block of flagged edges, at the block top
                 (flagged edge whose parent edge is not flagged). Fixes
                 over-fragmentation: a clade-spanning event becomes ONE clade
                 duplication -> folds to one reticulation.
  dedup        - drop any flagged edge whose clade is a strict subset of another
                 flagged edge's clade, keeping only the most ancestral flagged
                 edge in each containment hierarchy. Motivated by FP analysis:
                 ~87% of false positives are tips sitting INSIDE a truly-
                 duplicated (flagged) clade — redundant with the ancestral event.
                 Unlike collapse this uses clade containment, not parent-
                 adjacency, so it catches descendants even across an unflagged
                 intermediate edge. More conservative than clade_collapse: it
                 never invents an ancestral event or merges separate flagged edges.
  dedup_cap    - dedup, then cap.
  cap          - flagged edges, but cap duplications per species at the inferred
                 copy bound (fixes genuine over-ploidy). Confidence-ordered.
  collapse_cap - collapse, then cap.
  bound_driven - no threshold: take edges in confidence order until each species
                 reaches its copy bound. Removes the threshold knob.

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
    flagged = [i for i in range(n_edges) if float(wgd_probs[i]) >= threshold]
    flagged_set = set(flagged)

    if strategy == "raw":
        return flagged

    if strategy == "collapse":
        # block tops: a flagged edge whose parent edge is not flagged
        return [i for i in flagged if parent_edge[i] not in flagged_set]

    if strategy in ("dedup", "dedup_cap"):
        # Keep a flagged edge only if its clade is NOT a strict subset of another
        # flagged edge's clade (i.e. drop redundant descendants like the tip flags
        # inside a truly-duplicated clade). Containment-based, not parent-adjacency.
        kept = [i for i in flagged
                if not any(i != j and clades[i] < clades[j] for j in flagged)]
        if strategy == "dedup_cap":
            if copy_bound is None:
                raise ValueError("dedup_cap strategy requires copy_bound")
            kept = _cap(kept, wgd_probs, clades, copy_bound)
        return kept

    if strategy == "cap":
        if copy_bound is None:
            raise ValueError("cap strategy requires copy_bound")
        return _cap(flagged, wgd_probs, clades, copy_bound)

    if strategy == "collapse_cap":
        if copy_bound is None:
            raise ValueError("collapse_cap strategy requires copy_bound")
        tops = [i for i in flagged if parent_edge[i] not in flagged_set]
        return _cap(tops, wgd_probs, clades, copy_bound)

    if strategy == "bound_driven":
        if copy_bound is None:
            raise ValueError("bound_driven strategy requires copy_bound")
        return _cap(list(range(n_edges)), wgd_probs, clades, copy_bound)

    if strategy in ("clade_collapse", "clade_collapse_cap"):
        # The model flags individual polyploid species (tips). When a set of
        # flagged tips forms a clade, the true event is one duplication at their
        # common ancestor — so place a single event at the most ancestral edge
        # whose clade is entirely flagged. Greedy maximal-clade cover: take the
        # largest all-flagged clade edges first; remaining flagged tips stay as
        # their own events.
        flagged_species = set()
        for i in flagged:
            flagged_species |= clades[i]
        candidates = [i for i in range(n_edges) if clades[i] <= flagged_species]
        candidates.sort(key=lambda i: -len(clades[i]))  # largest clades first
        selected, covered = [], set()
        for i in candidates:
            if not (clades[i] & covered):
                selected.append(i)
                covered |= clades[i]
        if strategy == "clade_collapse_cap":
            if copy_bound is None:
                raise ValueError("clade_collapse_cap strategy requires copy_bound")
            selected = _cap(selected, wgd_probs, clades, copy_bound)
        return selected

    raise ValueError(f"Unknown strategy: {strategy}")
