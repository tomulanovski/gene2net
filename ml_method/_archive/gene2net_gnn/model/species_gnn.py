"""Module 3: GNN on the species tree producing node and edge embeddings."""
import torch
import torch.nn as nn
from torch_geometric.nn import GATConv, GCNConv, GINConv, SAGEConv


CONV_REGISTRY = {
    "gat": lambda dim, **kw: GATConv(dim, dim, heads=kw.get("n_gat_heads", 4), concat=False, dropout=kw.get("dropout", 0.2)),
    "gcn": lambda dim, **kw: GCNConv(dim, dim),
    "gin": lambda dim, **kw: GINConv(nn.Sequential(nn.Linear(dim, dim), nn.ReLU(), nn.Linear(dim, dim))),
    "sage": lambda dim, **kw: SAGEConv(dim, dim),
}


class SpeciesTreeGNN(nn.Module):
    """GNN operating on the species tree.

    Supports multiple GNN architectures (GAT, GCN, GIN, GraphSAGE) via
    the ``conv_type`` parameter. Produces both node and edge embeddings.
    Edge embeddings are computed by concatenating endpoint node embeddings
    and projecting through an MLP.
    """
    def __init__(self, input_dim, hidden_dim, num_layers=3, conv_type="gat", n_gat_heads=4, dropout=0.2, edge_feat_dim=0):
        super().__init__()
        self.input_proj = nn.Linear(input_dim, hidden_dim) if input_dim != hidden_dim else nn.Identity()

        if conv_type not in CONV_REGISTRY:
            raise ValueError(f"Unknown conv_type '{conv_type}'. Choose from: {list(CONV_REGISTRY.keys())}")

        self.layers = nn.ModuleList()
        for _ in range(num_layers):
            self.layers.append(CONV_REGISTRY[conv_type](hidden_dim, n_gat_heads=n_gat_heads, dropout=dropout))

        self.norms = nn.ModuleList([nn.LayerNorm(hidden_dim) for _ in range(num_layers)])
        self.dropout = nn.Dropout(dropout)

        # Edge embedding MLP: node pair + optional hand-crafted edge features
        self.edge_feat_dim = edge_feat_dim
        edge_input_dim = 2 * hidden_dim + edge_feat_dim
        self.edge_mlp = nn.Sequential(
            nn.Linear(edge_input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

    def forward(self, x, edge_index, edge_features=None):
        """
        Args:
            x: [N, input_dim] node features
            edge_index: [2, 2E] undirected edges
            edge_features: [E, edge_feat_dim] optional hand-crafted edge features
                           (one per undirected edge, same order as even-indexed edges)
        Returns:
            node_embeddings: [N, hidden_dim]
            edge_embeddings: [E, hidden_dim] one per undirected edge
        """
        x = self.input_proj(x)

        for layer, norm in zip(self.layers, self.norms):
            x = x + layer(x, edge_index)  # residual
            x = norm(x)
            x = self.dropout(x)

        # Compute edge embeddings from parent→child edges (even indices)
        src = edge_index[0, 0::2]  # parent nodes
        dst = edge_index[1, 0::2]  # child nodes
        edge_input = torch.cat([x[src], x[dst]], dim=-1)
        if edge_features is not None:
            edge_input = torch.cat([edge_input, edge_features], dim=-1)
        elif self.edge_feat_dim > 0:
            # Zero-pad when edge features expected but not provided
            edge_input = torch.cat([edge_input, torch.zeros(src.shape[0], self.edge_feat_dim, device=x.device)], dim=-1)
        edge_emb = self.edge_mlp(edge_input)

        return x, edge_emb
