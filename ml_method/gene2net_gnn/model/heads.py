"""Prediction heads for WGD event detection and partner edge prediction."""
import torch
import torch.nn as nn


class WGDHead(nn.Module):
    """Predict number of WGD events per edge.

    Outputs logits for classes 0, 1, ..., max_events.
    """
    def __init__(self, edge_dim, max_events=3):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(edge_dim, edge_dim),
            nn.ReLU(),
            nn.Linear(edge_dim, max_events + 1),
        )

    def forward(self, edge_embeddings):
        """
        Args:
            edge_embeddings: [E, edge_dim]
        Returns:
            logits: [E, max_events+1]
        """
        return self.mlp(edge_embeddings)


class PartnerHead(nn.Module):
    """Predict partner edge for a WGD event via bilinear attention.

    For each query (WGD) edge, scores all candidate edges as potential partners.
    """
    def __init__(self, edge_dim):
        super().__init__()
        self.query_proj = nn.Linear(edge_dim, edge_dim, bias=False)
        self.key_proj = nn.Linear(edge_dim, edge_dim, bias=False)
        self.scale = edge_dim ** 0.5

    def forward(self, query_edges, all_edges):
        """
        Args:
            query_edges: [N_wgd, edge_dim] - embeddings of edges with predicted WGD
            all_edges: [E, edge_dim] - all edge embeddings (candidate partners)
        Returns:
            scores: [N_wgd, E] - raw logits (NOT softmaxed, so loss can apply CE)
        """
        q = self.query_proj(query_edges)  # [N_wgd, D]
        k = self.key_proj(all_edges)       # [E, D]
        scores = torch.matmul(q, k.t()) / self.scale  # [N_wgd, E]
        return scores
