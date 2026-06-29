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


def labels_from_metadata_for_sample(metadata_events, sample_dict, jaccard_threshold: float = 0.5) -> TrainingLabels:
    events = events_from_metadata(metadata_events)
    bip = sample_edge_bipartitions(sample_dict)
    return map_events_to_edges(events, bip, jaccard_threshold=jaccard_threshold)
