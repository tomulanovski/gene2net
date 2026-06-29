import os
import pytest
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from scripts.relabel_from_metadata import relabel_one_sample


def _write_sample(tmp_path):
    # MUL-tree has D,E duplicated as a clade -> from_trees produces (fragmented) labels
    sp_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    mul = Tree("(((A:1,B:1):1,(D:1,E:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    gts = [Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1) for _ in range(8)]
    gts += [Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1) for _ in range(2)]
    s = Gene2NetSample.from_trees(sp_tree, gts, ["A", "B", "C", "D", "E"], mul_tree=mul)
    d = tmp_path / "sample_0000"
    s.save(str(d))
    return str(d)


def test_relabel_writes_clade_labels(tmp_path):
    d = _write_sample(tmp_path)
    metadata = {
        "polyploid_species": {"D": 2, "E": 2},
        "events": [{"target_clade": ["D", "E"], "partner_clade": ["A", "B"], "event_type": "allo"}],
    }
    labels = relabel_one_sample(d, metadata)
    assert os.path.exists(f"{d}/labels_clade.pkl")
    assert labels.mask == [True]
    assert len(labels.wgd_edges) == 1  # one clade event, NOT two fragmented tips


def test_relabel_allows_subset_polyploids(tmp_path):
    # The sample's labeled polyploids ({D,E}) are a SUBSET of metadata's
    # ({C,D,E}) — the old pipeline dropped an unmappable event. This must NOT
    # raise; metadata is the fuller ground truth.
    d = _write_sample(tmp_path)
    metadata = {
        "n_species": 5,
        "polyploid_species": {"C": 2, "D": 2, "E": 2},
        "events": [{"target_clade": ["D", "E"], "partner_clade": ["A", "B"], "event_type": "allo"}],
    }
    labels = relabel_one_sample(d, metadata)
    assert labels.mask == [True]


def test_relabel_rejects_wrong_species_count(tmp_path):
    d = _write_sample(tmp_path)
    metadata = {
        "n_species": 99,  # sample has 5 species -> wrong index
        "polyploid_species": {"D": 2, "E": 2},
        "events": [{"target_clade": ["D", "E"], "partner_clade": ["A", "B"], "event_type": "allo"}],
    }
    with pytest.raises(ValueError, match="species-count"):
        relabel_one_sample(d, metadata)


def test_relabel_rejects_index_mismatch(tmp_path):
    d = _write_sample(tmp_path)
    # metadata for a DIFFERENT sample (wrong polyploids) -> must fail loud
    metadata = {
        "polyploid_species": {"A": 2, "C": 2},
        "events": [{"target_clade": ["A"], "partner_clade": ["C"], "event_type": "allo"}],
    }
    with pytest.raises(ValueError, match="polyploid"):
        relabel_one_sample(d, metadata)
