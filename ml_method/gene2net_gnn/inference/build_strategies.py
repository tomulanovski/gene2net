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
  detect_only  - ploidy-independent: events are exactly the edges with prob >=
                 threshold, with no copy bound at all. The analogue of iterative
                 GRAMPA using reconciliation instead of an inferred ploidy. Robust to
                 fractionation, but over-predicts on clean data without a cap.

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


def _representative_copy_number(copy_distribution, kernel_width: int = 2) -> int:
    """Kernel-smoothed representative copy number, robust to occasional duplication
    and loss. Faithful reimplementation of the triangular-kernel estimator used to
    build the inferred-ploidy multiset. One distinct value returns it; two return the
    more frequent; three or more smooth with a triangular kernel and take the peak."""
    sorted_counts = sorted(copy_distribution.items())
    if len(sorted_counts) == 1:
        return sorted_counts[0][0]
    if len(sorted_counts) == 2:
        return sorted_counts[0][0] if sorted_counts[0][1] >= sorted_counts[1][1] else sorted_counts[1][0]
    copy_numbers = [num for num, _ in sorted_counts]
    frequencies = [freq for _, freq in sorted_counts]
    smoothed = {}
    for center in copy_numbers:
        for j, copy_num in enumerate(copy_numbers):
            distance = abs(center - copy_num)
            if distance <= kernel_width:
                weight = 1 - (distance / (kernel_width + 1))
                smoothed[center] = smoothed.get(center, 0) + frequencies[j] * weight
    max_smoothed = 0
    representative = copy_numbers[0]
    for copy_num, val in smoothed.items():
        if val > max_smoothed:
            max_smoothed = val
            representative = copy_num
    return representative


def infer_copy_bound_kernel(gene_trees, kernel_width: int = 2) -> Dict[str, int]:
    """Per-species copy bound via kernel-smoothed estimation of the gene-tree copy
    distribution. Self-contained ploidy inference for the method: for each species,
    count its copies in every gene tree (0 when absent), smooth that distribution with
    a triangular kernel, and take the peak. More robust to occasional duplication and
    loss than the majority-consensus infer_copy_bound. Minimum 1."""
    per_tree = [collections.Counter(leaf.name for leaf in gt.get_leaves()) for gt in gene_trees]
    species = set()
    for c in per_tree:
        species.update(c.keys())
    bound = {}
    for s in species:
        dist = collections.Counter()
        for c in per_tree:
            dist[c.get(s, 0)] += 1
        bound[s] = max(_representative_copy_number(dist, kernel_width), 1)
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

    if strategy == "detect_only":
        # Ploidy-INDEPENDENT: events are exactly the edges the detection head is
        # confident about, with no copy-number bound anywhere. This is the analogue of
        # iterative GRAMPA relying on reconciliation rather than an inferred ploidy, so
        # it is robust when fractionation destroys the copy count. It has no ploidy cap,
        # so it can over-predict on clean data. Tune threshold on the validation split.
        return [i for i in range(n_edges) if float(wgd_probs[i]) >= threshold]

    raise ValueError(
        f"Unknown strategy: {strategy!r}. Supported: 'bound_driven', 'cap', 'detect', "
        "'detect_only'."
    )
