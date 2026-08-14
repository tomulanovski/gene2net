"""WGD Detector: GIN gene tree encoder + species tree GAT.

End-to-end model that:
1. Encodes gene trees with GIN to extract per-edge gene tree features
2. Combines with hand-crafted node/edge features
3. Runs GAT on species tree
4. Predicts binary WGD per edge

This replaces the Phase 1 features-only model by adding learned gene tree
representations, while keeping the architecture simple and trainable.
"""
import torch
import torch.nn as nn
from torch_geometric.nn import GATConv
from typing import Dict, List, Set

from gene2net_gnn.model.gene_tree_gin import GeneTreeGIN


class WGDDetector(nn.Module):
    """Full model for WGD edge detection.

    Args:
        n_species: Number of distinct species in the dataset.
        node_feat_dim: Dimension of hand-crafted node features (default: 13).
        edge_feat_dim: Dimension of hand-crafted edge features (default: 4).
        hidden_dim: Hidden dimension throughout the model.
        n_gin_layers: Number of GIN layers for gene tree encoding.
        n_gat_layers: Number of GAT layers for species tree.
        n_gat_heads: Number of attention heads per GAT layer.
        dropout: Dropout rate.
    """

    def __init__(
        self,
        n_species: int,
        node_feat_dim: int = 13,
        edge_feat_dim: int = 4,
        hidden_dim: int = 64,
        n_gin_layers: int = 2,
        n_gat_layers: int = 3,
        n_gat_heads: int = 4,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.node_feat_dim = node_feat_dim
        self.hidden_dim = hidden_dim

        # Gene tree encoder
        self.gene_tree_encoder = GeneTreeGIN(
            n_species=n_species,
            embed_dim=hidden_dim,
            n_gin_layers=n_gin_layers,
            dropout=dropout,
        )

        # Node feature projection (hand-crafted → hidden_dim)
        self.node_proj = nn.Sequential(
            nn.Linear(node_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

        # GAT layers on species tree
        self.gat_layers = nn.ModuleList()
        self.gat_norms = nn.ModuleList()
        for _ in range(n_gat_layers):
            self.gat_layers.append(
                GATConv(hidden_dim, hidden_dim, heads=n_gat_heads,
                        concat=False, dropout=dropout)
            )
            self.gat_norms.append(nn.LayerNorm(hidden_dim))
        self.dropout = nn.Dropout(dropout)

        # Edge embedding MLP:
        # parent node [hidden_dim] + child node [hidden_dim]
        # + hand-crafted edge features [edge_feat_dim]
        # + gene tree features [output_dim from GIN: 2*hidden_dim + copy count stats]
        gene_tree_out_dim = self.gene_tree_encoder.output_dim
        edge_input_dim = 2 * hidden_dim + edge_feat_dim + gene_tree_out_dim
        self.edge_mlp = nn.Sequential(
            nn.Linear(edge_input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )

        # WGD prediction head
        self.wgd_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 2),
        )

    def propagate_features_to_internal(self, node_features, edge_index, is_leaf):
        """Fill internal node features by averaging descendant leaf features."""
        x = node_features.clone()
        n_nodes = x.shape[0]

        parent_to_children = {}
        child_to_parent = {}
        for i in range(0, edge_index.shape[1], 2):
            parent = edge_index[0, i].item()
            child = edge_index[1, i].item()
            parent_to_children.setdefault(parent, []).append(child)
            child_to_parent[child] = parent

        leaf_indices = is_leaf.nonzero(as_tuple=True)[0].tolist()
        internal_indices = (~is_leaf).nonzero(as_tuple=True)[0].tolist()

        descendant_sum = torch.zeros(n_nodes, self.node_feat_dim, device=x.device)
        descendant_count = torch.zeros(n_nodes, 1, device=x.device)

        for leaf_idx in leaf_indices:
            descendant_sum[leaf_idx] = x[leaf_idx]
            descendant_count[leaf_idx] = 1.0

        processed = set(leaf_indices)
        queue = list(set(child_to_parent.get(l, -1) for l in leaf_indices) - {-1})

        while queue:
            next_queue = []
            for node in queue:
                children = parent_to_children.get(node, [])
                if all(c in processed for c in children):
                    for c in children:
                        descendant_sum[node] += descendant_sum[c]
                        descendant_count[node] += descendant_count[c]
                    processed.add(node)
                    parent = child_to_parent.get(node, -1)
                    if parent != -1 and parent not in processed:
                        next_queue.append(parent)
                else:
                    next_queue.append(node)
            if set(next_queue) == set(queue):
                break
            queue = list(set(next_queue))

        for idx in internal_indices:
            if descendant_count[idx] > 0:
                x[idx] = descendant_sum[idx] / descendant_count[idx]

        return x

    def _get_edge_clades_as_species(
        self,
        edge_index: torch.Tensor,
        is_leaf: torch.Tensor,
        species_tree_node_names: List[str],
        species_list: List[str],
    ) -> Dict[int, Set[int]]:
        """Get species indices in each edge's clade."""
        sp_to_idx = {sp: i for i, sp in enumerate(species_list)}
        n_edges = edge_index.shape[1] // 2

        # Build parent→children map
        parent_to_children: Dict[int, List[int]] = {}
        for i in range(0, edge_index.shape[1], 2):
            p = edge_index[0, i].item()
            c = edge_index[1, i].item()
            parent_to_children.setdefault(p, []).append(c)

        def get_descendant_leaves(node: int) -> Set[int]:
            if is_leaf[node]:
                name = species_tree_node_names[node]
                if name in sp_to_idx:
                    return {sp_to_idx[name]}
                return set()
            result = set()
            for child in parent_to_children.get(node, []):
                result |= get_descendant_leaves(child)
            return result

        edge_clades = {}
        for edge_idx in range(n_edges):
            child_node = edge_index[1, edge_idx * 2].item()
            edge_clades[edge_idx] = get_descendant_leaves(child_node)

        return edge_clades

    def forward(
        self,
        node_features: torch.Tensor,
        edge_index: torch.Tensor,
        edge_features: torch.Tensor,
        is_leaf: torch.Tensor,
        gene_tree_edge_indices: List[torch.Tensor],
        gene_tree_species_ids: List[torch.Tensor],
        gene_tree_leaf_masks: List[torch.Tensor],
        species_tree_node_names: List[str],
        species_list: List[str],
    ):
        """Forward pass.

        Args:
            node_features: [N, node_feat_dim] species tree node features.
            edge_index: [2, 2E] species tree edges.
            edge_features: [E, edge_feat_dim] species tree edge features.
            is_leaf: [N] boolean mask.
            gene_tree_edge_indices: List of [2, 2E_gt] per gene tree.
            gene_tree_species_ids: List of [N_gt] per gene tree.
            gene_tree_leaf_masks: List of [N_gt] boolean per gene tree.
            species_tree_node_names: Node names for species tree.
            species_list: Sorted species names.

        Returns:
            wgd_logits: [E, 2] binary WGD logits.
            edge_embeddings: [E, hidden_dim].
        """
        n_edges = edge_index.shape[1] // 2

        # 1. Get edge clades as species index sets
        edge_clades = self._get_edge_clades_as_species(
            edge_index, is_leaf, species_tree_node_names, species_list
        )

        # 2. Encode gene trees → per-edge features [n_edges, 2*hidden_dim]
        gene_tree_features = self.gene_tree_encoder(
            gene_tree_edge_indices,
            gene_tree_species_ids,
            gene_tree_leaf_masks,
            edge_clades,
            n_edges,
        )

        # 3. Process species tree node features
        x = self.propagate_features_to_internal(node_features, edge_index, is_leaf)
        x = self.node_proj(x)

        # 4. GAT message passing on species tree
        for gat, norm in zip(self.gat_layers, self.gat_norms):
            x = x + gat(x, edge_index)
            x = norm(x)
            x = self.dropout(x)

        # 5. Edge embeddings: node pair + hand-crafted + gene tree features
        src = edge_index[0, 0::2]
        dst = edge_index[1, 0::2]
        edge_input = torch.cat([
            x[src], x[dst],
            edge_features,
            gene_tree_features,
        ], dim=-1)
        edge_emb = self.edge_mlp(edge_input)

        # 6. WGD prediction
        wgd_logits = self.wgd_head(edge_emb)

        return wgd_logits, edge_emb
