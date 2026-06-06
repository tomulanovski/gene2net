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


def propagate_to_internal(node_features, edge_index, is_leaf, node_feat_dim):
    """Fill internal node features by averaging their descendant leaf features.

    Weight-independent: depends only on the input features and tree structure,
    so the result can be cached once per sample instead of recomputed every
    forward pass. Shared by the model's method and the training data prep.

    Uses iterative bottom-up propagation on the tree structure.
    """
    x = node_features.clone()
    n_nodes = x.shape[0]

    # edge_index has parent→child at even indices, child→parent at odd
    parent_to_children = {}
    child_to_parent = {}
    for i in range(0, edge_index.shape[1], 2):
        parent = edge_index[0, i].item()
        child = edge_index[1, i].item()
        parent_to_children.setdefault(parent, []).append(child)
        child_to_parent[child] = parent

    leaf_indices = is_leaf.nonzero(as_tuple=True)[0].tolist()
    internal_indices = (~is_leaf).nonzero(as_tuple=True)[0].tolist()

    descendant_sum = torch.zeros(n_nodes, node_feat_dim, device=x.device)
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

        # Partner head: scores edge j as the partner of WGD edge i from the
        # ordered pair of edge embeddings [emb_i || emb_j]. The diagonal (j==i)
        # represents autopolyploidy; an off-diagonal partner is allopolyploidy.
        self.partner_head = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def compute_partner_scores(self, edge_emb):
        """Score every ordered edge pair for the partner relationship.

        Args:
            edge_emb: [E, hidden_dim] per-edge embeddings (from forward).

        Returns:
            scores: [E, E] where scores[i, j] is the compatibility of edge j as
            the partner of WGD edge i. Row i softmaxed over j (including j==i,
            which means autopolyploidy) gives the partner distribution.
        """
        n = edge_emb.shape[0]
        ei = edge_emb.unsqueeze(1).expand(n, n, -1)   # [E, E, H]  (i along dim 0)
        ej = edge_emb.unsqueeze(0).expand(n, n, -1)   # [E, E, H]  (j along dim 1)
        pair = torch.cat([ei, ej], dim=-1)            # [E, E, 2H]
        return self.partner_head(pair).squeeze(-1)    # [E, E]

    def propagate_features_to_internal(self, node_features, edge_index, is_leaf):
        """Fill internal node features by averaging descendant leaves.

        Thin wrapper over the shared module-level ``propagate_to_internal`` so
        the same logic is reused by the cached data-prep path.
        """
        return propagate_to_internal(node_features, edge_index, is_leaf, self.node_feat_dim)

    def forward(self, node_features, edge_index, edge_features, is_leaf,
                node_features_propagated=None):
        """Forward pass.

        Args:
            node_features: [N, node_feat_dim] hand-crafted features (zeros for internal).
            edge_index: [2, 2E] undirected edges.
            edge_features: [E, edge_feat_dim] per undirected edge.
            is_leaf: [N] boolean mask.
            node_features_propagated: optional [N, node_feat_dim] precomputed
                internal-node propagation. When provided, skips the (weight-
                independent) BFS — this is the cached fast path. Falls back to
                computing it when None, so the model still works standalone.

        Returns:
            wgd_logits: [E, 2] binary WGD logits per edge.
            edge_embeddings: [E, hidden_dim] for potential downstream use.
        """
        # Step 1: Propagate leaf features to internal nodes (cached if provided)
        if node_features_propagated is not None:
            x = node_features_propagated
        else:
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
