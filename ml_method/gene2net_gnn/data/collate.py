"""Custom collation for Gene2Net batches."""
from typing import Dict, List
from gene2net_gnn.data.dataset import Gene2NetSample


def gene2net_collate(samples: List[Gene2NetSample]) -> Dict:
    """Collate samples into a batch dict.

    Since each sample has a different number of gene trees and potentially
    different tree sizes, we don't stack tensors - we keep lists.
    """
    return {
        "species_trees": [
            {
                "edge_index": s.species_tree_edge_index,
                "node_features": s.species_tree_node_features,
                "edge_features": s.species_tree_edge_features,
                "is_leaf": s.species_tree_is_leaf,
                "node_names": s.species_tree_node_names,
                "n_species": s.n_species,
                "species_list": s.species_list,
            }
            for s in samples
        ],
        "gene_trees": [
            {
                "edge_indices": s.gene_tree_edge_indices,
                "species_ids": s.gene_tree_species_ids,
                "branch_lengths": s.gene_tree_branch_lengths,
                "leaf_masks": s.gene_tree_leaf_masks,
            }
            for s in samples
        ],
        "labels": [s.labels for s in samples],
    }
