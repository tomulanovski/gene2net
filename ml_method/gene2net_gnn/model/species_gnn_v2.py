"""Species Tree GNN for Phase 1: hand-crafted features only, no gene tree encoder.

This is a clean, focused model for WGD edge detection using only:
  - Node features: 13-dim (copy counts + clustering summary) for leaves,
    propagated to internal nodes via descendant averaging.
  - Edge features: 4-dim (concordance, branch length, clade size, depth).
  - Graph structure: GAT message passing on the species tree.

The key insight is that WGD signal exists in the hand-crafted features
(demonstrated by baselines achieving F1 ~0.3), so this model focuses on
leveraging graph structure to improve over non-graph baselines.
"""
import torch
import torch.nn as nn
from torch_geometric.nn import GATConv


class SpeciesTreeGNNv2(nn.Module):
    """GNN on species tree with hand-crafted features only.

    Architecture:
        1. Propagate leaf features to internal nodes (average of descendant leaves)
        2. Project node features to hidden_dim
        3. GAT message passing (with residual connections + layer norm)
        4. Edge embeddings from endpoint node pairs + edge features
        5. Binary WGD classification per edge

    Args:
        node_feat_dim: Dimension of input node features (default: 13).
        edge_feat_dim: Dimension of input edge features (default: 4).
        hidden_dim: Hidden dimension throughout the model.
        n_gat_layers: Number of GAT message passing layers.
        n_gat_heads: Number of attention heads per GAT layer.
        dropout: Dropout rate.
    """

    def __init__(
        self,
        node_feat_dim: int = 13,
        edge_feat_dim: int = 4,
        hidden_dim: int = 64,
        n_gat_layers: int = 3,
        n_gat_heads: int = 4,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.node_feat_dim = node_feat_dim
        self.edge_feat_dim = edge_feat_dim
        self.hidden_dim = hidden_dim

        # Node feature projection
        self.node_proj = nn.Sequential(
            nn.Linear(node_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

        # GAT layers with residual connections and layer norm
        self.gat_layers = nn.ModuleList()
        self.gat_norms = nn.ModuleList()
        for _ in range(n_gat_layers):
            self.gat_layers.append(
                GATConv(hidden_dim, hidden_dim, heads=n_gat_heads, concat=False, dropout=dropout)
            )
            self.gat_norms.append(nn.LayerNorm(hidden_dim))
        self.dropout = nn.Dropout(dropout)

        # Edge embedding: parent node + child node + edge features → hidden_dim
        edge_input_dim = 2 * hidden_dim + edge_feat_dim
        self.edge_mlp = nn.Sequential(
            nn.Linear(edge_input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )

        # WGD prediction head: edge embedding → binary (0=no WGD, 1=WGD)
        self.wgd_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 2),
        )

    def propagate_features_to_internal(self, node_features, edge_index, is_leaf):
        """Fill internal node features by averaging their descendant leaf features.

        This is critical: without this, internal nodes have zero features and
        the GNN has no signal for internal edges (where many WGD events occur).

        Uses iterative bottom-up propagation on the tree structure.
        """
        x = node_features.clone()
        n_nodes = x.shape[0]

        # Build adjacency: for each node, find its children
        # edge_index has parent→child at even indices, child→parent at odd
        parent_to_children = {}
        child_to_parent = {}
        for i in range(0, edge_index.shape[1], 2):
            parent = edge_index[0, i].item()
            child = edge_index[1, i].item()
            parent_to_children.setdefault(parent, []).append(child)
            child_to_parent[child] = parent

        # Find leaf nodes and internal nodes
        leaf_indices = is_leaf.nonzero(as_tuple=True)[0].tolist()
        internal_indices = (~is_leaf).nonzero(as_tuple=True)[0].tolist()

        # Bottom-up: for each internal node, average all descendant leaf features
        # Use BFS from leaves upward
        descendant_sum = torch.zeros(n_nodes, self.node_feat_dim, device=x.device)
        descendant_count = torch.zeros(n_nodes, 1, device=x.device)

        # Initialize leaves
        for leaf_idx in leaf_indices:
            descendant_sum[leaf_idx] = x[leaf_idx]
            descendant_count[leaf_idx] = 1.0

        # Process bottom-up: nodes whose all children are already processed
        processed = set(leaf_indices)
        queue = list(set(child_to_parent.get(l, -1) for l in leaf_indices) - {-1})

        while queue:
            next_queue = []
            for node in queue:
                children = parent_to_children.get(node, [])
                if all(c in processed for c in children):
                    # Aggregate children
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
                # Prevent infinite loop on edge cases
                break
            queue = list(set(next_queue))

        # Fill internal nodes with average of descendant leaves
        for idx in internal_indices:
            if descendant_count[idx] > 0:
                x[idx] = descendant_sum[idx] / descendant_count[idx]

        return x

    def forward(self, node_features, edge_index, edge_features, is_leaf):
        """Forward pass.

        Args:
            node_features: [N, node_feat_dim] hand-crafted features (zeros for internal).
            edge_index: [2, 2E] undirected edges.
            edge_features: [E, edge_feat_dim] per undirected edge.
            is_leaf: [N] boolean mask.

        Returns:
            wgd_logits: [E, 2] binary WGD logits per edge.
            edge_embeddings: [E, hidden_dim] for potential downstream use.
        """
        # Step 1: Propagate leaf features to internal nodes
        x = self.propagate_features_to_internal(node_features, edge_index, is_leaf)

        # Step 2: Project to hidden dim
        x = self.node_proj(x)

        # Step 3: GAT message passing with residual connections
        for gat, norm in zip(self.gat_layers, self.gat_norms):
            x = x + gat(x, edge_index)
            x = norm(x)
            x = self.dropout(x)

        # Step 4: Edge embeddings (parent-child pairs from even-indexed edges)
        src = edge_index[0, 0::2]  # parent nodes
        dst = edge_index[1, 0::2]  # child nodes
        edge_input = torch.cat([x[src], x[dst], edge_features], dim=-1)
        edge_emb = self.edge_mlp(edge_input)

        # Step 5: WGD prediction
        wgd_logits = self.wgd_head(edge_emb)

        return wgd_logits, edge_emb
