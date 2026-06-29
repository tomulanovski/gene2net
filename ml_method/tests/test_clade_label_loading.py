import pickle
import pytest
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.data.label_extractor import TrainingLabels


def _write_sample(tmp_path):
    sp_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    gts = [Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1) for _ in range(3)]
    s = Gene2NetSample.from_trees(sp_tree, gts, ["A", "B", "C", "D", "E"])
    d = tmp_path / "sample_0000"
    s.save(str(d))
    return str(d)


def test_clade_labels_loaded_from_sidecar(tmp_path):
    d = _write_sample(tmp_path)
    sentinel = TrainingLabels(wgd_edges=[3], partner_edges=[0], wgd_counts=[],
                              mask=[True], n_unmappable=0, n_edges=8)
    with open(f"{d}/labels_clade.pkl", "wb") as f:
        pickle.dump(sentinel, f)

    s_default = Gene2NetSample.load(d)
    s_clade = Gene2NetSample.load(d, clade_labels=True)
    assert s_clade.labels.wgd_edges == [3]               # came from sidecar
    assert s_default.labels != sentinel or s_default.labels is None  # default path unaffected


def test_clade_labels_missing_raises(tmp_path):
    d = _write_sample(tmp_path)
    with pytest.raises(FileNotFoundError):
        Gene2NetSample.load(d, clade_labels=True)
