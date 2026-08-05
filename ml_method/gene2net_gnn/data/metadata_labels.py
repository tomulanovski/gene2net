"""Build clade-level training labels from ground-truth metadata events.

The default packaging path (decompose_mul_tree) fragments clade-level
allopolyploidy into per-species reciprocal tip events. The metadata
(data/mul_trees_2k/metadata_*.json) stores the true events in clade-level form,
so we build labels from those instead. Events are mapped onto the sample's OWN
stored edge bipartitions, which guarantees the resulting label indices align
with the sample's species_tree_edge_features rows.
"""
from typing import FrozenSet, List, Tuple

from gene2net_gnn.inference.mul_tree_builder import WGDEvent
from gene2net_gnn.data.label_extractor import TrainingLabels, _best_matching_edge
from gene2net_gnn.data.tree_io import reorder_edge_index_preorder
from gene2net_gnn.data.features import edge_clades_species


def events_from_metadata(metadata_events: List[dict]) -> List[WGDEvent]:
    """Convert metadata event dicts to clade-level WGDEvents. Fails loud on bad input."""
    out: List[WGDEvent] = []
    for ev in metadata_events:
        target = frozenset(ev.get("target_clade") or [])
        if not target:
            raise ValueError(f"empty/missing target_clade in event {ev!r}")
        etype = ev.get("event_type")
        if etype == "auto":
            partner = target
        elif etype == "allo":
            pc = ev.get("partner_clade")
            if not pc:
                raise ValueError(f"allo event missing partner_clade: {ev!r}")
            partner = frozenset(pc)
        else:
            raise ValueError(f"unknown event_type {etype!r} in event {ev!r}")
        out.append(WGDEvent(wgd_edge_clade=target, partner_edge_clade=partner, confidence=1.0))
    return out


def sample_edge_bipartitions(sample_dict: dict) -> List[Tuple[int, FrozenSet[str]]]:
    """(edge_index, leaf-name set) per preorder edge, in the SAME order as the
    sample's species_tree_edge_features rows."""
    ei = reorder_edge_index_preorder(sample_dict["species_tree_edge_index"])
    species_list = sample_dict["species_list"]
    sp_to_idx = {sp: j for j, sp in enumerate(species_list)}
    idx_to_sp = {j: sp for sp, j in sp_to_idx.items()}
    node_species = [sp_to_idx.get(nm, -1) for nm in sample_dict["species_tree_node_names"]]
    clades_idx = edge_clades_species(ei, sample_dict["species_tree_is_leaf"], node_species)
    return [(i, frozenset(idx_to_sp[s] for s in clade)) for i, clade in enumerate(clades_idx)]


def map_events_to_edges(events, edge_bipartitions, jaccard_threshold: float = 0.5) -> TrainingLabels:
    """Map clade-level events to precomputed edge bipartitions by Jaccard.

    Mirrors label_extractor.map_events_to_astral but takes bipartitions directly
    so the indices align with the sample's stored features.
    """
    n_edges = len(edge_bipartitions)
    edge_wgd_counts = [0] * n_edges
    wgd_edges: List[int] = []
    partner_edges: List[int] = []
    mask: List[bool] = []
    n_unmappable = 0
    for ev in events:
        w_idx, w_score = _best_matching_edge(ev.wgd_edge_clade, edge_bipartitions)
        p_idx, p_score = _best_matching_edge(ev.partner_edge_clade, edge_bipartitions)
        mappable = w_score >= jaccard_threshold and p_score >= jaccard_threshold
        wgd_edges.append(w_idx)
        partner_edges.append(p_idx)
        mask.append(mappable)
        if mappable:
            edge_wgd_counts[w_idx] += 1
        else:
            n_unmappable += 1
    return TrainingLabels(
        wgd_edges=wgd_edges, partner_edges=partner_edges, wgd_counts=edge_wgd_counts,
        mask=mask, n_unmappable=n_unmappable, n_edges=n_edges,
    )


def _sibling_clade(true_tree, sp):
    """Leaf-name set of sp's sibling(s) in the true species tree (sp's OTHER parent)."""
    nodes = true_tree.search_nodes(name=sp)
    if not nodes or nodes[0].up is None:
        return None
    x = nodes[0]
    sib = set()
    for ch in x.up.get_children():
        if ch is not x:
            sib |= set(ch.get_leaf_names())
    sib.discard(sp)
    return frozenset(sib)


def _astral_home_set(x, clades):
    """Smallest ASTRAL edge-clade strictly containing species set x, minus x (x's ASTRAL
    home). Generalizes _astral_home from a single species to a clade."""
    x = frozenset(x)
    supersets = [c for c in clades if x < c]
    if not supersets:
        return None
    return frozenset(min(supersets, key=len) - x)


def _astral_home(sp, clades):
    """sp's sibling clade in the ASTRAL tree: smallest edge-clade strictly containing sp, minus sp."""
    return _astral_home_set({sp}, clades)


def relabel_events_partner_as_away(events, true_tree, edge_bipartitions):
    """For each ALLO event, set partner = the true parent that is NOT the target's ASTRAL
    home.

    Fixes the ~55% of allo labels whose partner is the home (verified by label_audit):
    when ASTRAL placed the target next to its labelled partner B, the old target == home,
    which gives the model a contradictory objective and makes the build graft a copy where
    the target already sits (the sp39 collapse). Here we retarget those to the OTHER true
    parent A. Auto events (partner==target) are left unchanged. Single-species targets are
    always handled. Clade-level targets are handled only when the clade is monophyletic in
    ASTRAL, so its home is a single well-defined sibling; otherwise the event is left
    unchanged rather than risk a fuzzy retarget.
    """
    clades = [c for _, c in edge_bipartitions]
    clade_set = set(clades)                                 # for the ASTRAL-monophyly guard
    leaf_map = {leaf.name: leaf for leaf in true_tree.get_leaves()}  # one traversal, not per-event

    def true_sibling(C):
        """Sibling leaf set of C's MRCA in the true tree (the other parent A). Monophyly-
        guarded: returns None unless C is exactly the leaf set of a true-tree node. For a
        single species this reduces to that leaf's sibling (identical to the old behavior)."""
        leaves = [leaf_map[s] for s in C if s in leaf_map]
        if len(leaves) < len(C):
            return None
        mrca = leaves[0] if len(leaves) == 1 else true_tree.get_common_ancestor(leaves)
        if set(mrca.get_leaf_names()) != set(C) or mrca.up is None:
            return None
        sib = set()
        for ch in mrca.up.get_children():
            if ch is not mrca:
                sib |= set(ch.get_leaf_names())
        return frozenset(sib - set(C))

    out = []
    for ev in events:
        tgt, B = ev.wgd_edge_clade, ev.partner_edge_clade
        if B == tgt:                                        # auto: leave unchanged
            out.append(ev)
            continue
        if len(tgt) > 1 and frozenset(tgt) not in clade_set:  # clade not monophyletic in ASTRAL
            out.append(ev)
            continue
        A = true_sibling(tgt)                   # the other (true-home) parent
        H = _astral_home_set(tgt, clades)       # where ASTRAL placed the target
        if A is None or H is None:
            out.append(ev)
            continue
        new_partner = A if (H & B) else B       # home == labelled partner -> retarget to A
        out.append(WGDEvent(wgd_edge_clade=tgt, partner_edge_clade=new_partner,
                            confidence=ev.confidence))
    return out


def labels_from_metadata_for_sample(metadata_events, sample_dict, jaccard_threshold: float = 0.5,
                                    true_tree=None) -> TrainingLabels:
    events = events_from_metadata(metadata_events)
    bip = sample_edge_bipartitions(sample_dict)
    if true_tree is not None:
        events = relabel_events_partner_as_away(events, true_tree, bip)
    return map_events_to_edges(events, bip, jaccard_threshold=jaccard_threshold)


def home_edge_for_event(target_clade, event_type, true_tree):
    """The 'home' parent clade of an event.

    allo: the target's sibling clade in the TRUE species tree (the other parent).
    auto: the target itself (both parents identical).
    Returns None if the target clade is not found / has no sibling.
    """
    target = frozenset(target_clade)
    if event_type == "auto":
        return target
    node = None
    for n in true_tree.traverse("postorder"):
        if frozenset(n.get_leaf_names()) == target:
            node = n
            break
    if node is None or node.up is None:
        return None
    sib = set()
    for ch in node.up.get_children():
        if ch is not node:
            sib |= set(ch.get_leaf_names())
    sib -= set(target)
    return frozenset(sib) if sib else None


def two_parent_labels_from_metadata(metadata_events, sample_dict, true_tree,
                                    jaccard_threshold: float = 0.5) -> TrainingLabels:
    """TrainingLabels carrying BOTH parents: partner_edges (B = metadata partner)
    and home_edges (A = target's true-tree sibling; == wgd edge for auto)."""
    bip = sample_edge_bipartitions(sample_dict)
    events = events_from_metadata(metadata_events)          # target + partner (B)
    labels = map_events_to_edges(events, bip, jaccard_threshold=jaccard_threshold)
    home_edges: List[int] = []
    for ev, md in zip(events, metadata_events):
        hc = home_edge_for_event(ev.wgd_edge_clade, md.get("event_type"), true_tree)
        if hc is None:
            home_edges.append(-1)
            continue
        h_idx, h_score = _best_matching_edge(hc, bip)
        home_edges.append(h_idx if h_score >= jaccard_threshold else -1)
    labels.home_edges = home_edges
    return labels
