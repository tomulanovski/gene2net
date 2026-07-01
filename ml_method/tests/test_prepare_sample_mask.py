import torch
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.data.label_extractor import TrainingLabels
from gene2net_gnn.training.trainer_phase1 import prepare_sample


def test_unmappable_event_does_not_mask_a_positive_edge():
    """An unmappable event whose best-match edge collides with a mappable POSITIVE
    edge must not remove that positive from the detection loss/metrics."""
    sp = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    gts = [Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1) for _ in range(3)]
    s = Gene2NetSample.from_trees(sp, gts, ["A", "B", "C", "D", "E"])
    n_edges = s.species_tree_edge_index.shape[1] // 2

    wgd_counts = [0] * n_edges
    wgd_counts[5] = 1  # edge 5 is a real positive (from a mappable event)
    # event 0: mappable, WGD on edge 5. event 1: unmappable, best-match also edge 5.
    s.labels = TrainingLabels(
        wgd_edges=[5, 5], partner_edges=[0, 3], wgd_counts=wgd_counts,
        mask=[True, False], n_unmappable=1, n_edges=n_edges,
    )

    _, wgd_targets, mask = prepare_sample(s, torch.device("cpu"))
    assert wgd_targets[5].item() == 1          # it's a positive
    assert bool(mask[5].item()) is True        # must NOT be masked out
