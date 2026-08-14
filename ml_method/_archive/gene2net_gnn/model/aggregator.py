"""Module 2: Cross-tree aggregator using multi-head attention pooling."""
import torch
import torch.nn as nn


class CrossTreeAggregator(nn.Module):
    """Aggregate variable-count embeddings for one species into fixed-size representation.

    Uses K learned query vectors. Each query attends over all input embeddings
    independently, producing K summary vectors that are concatenated.

    For a diploid species (1 copy per tree), all heads will learn similar summaries.
    For a tetraploid (2 copies per tree), different heads can specialize to
    capture different subgenome signals.
    """
    def __init__(self, embed_dim: int, n_heads: int = 8):
        super().__init__()
        self.embed_dim = embed_dim
        self.n_heads = n_heads
        self.queries = nn.Parameter(torch.randn(n_heads, embed_dim))
        self.key_proj = nn.Linear(embed_dim, embed_dim)
        self.scale = embed_dim ** 0.5

    def forward(self, embeddings: torch.Tensor) -> torch.Tensor:
        """
        Args:
            embeddings: [M, embed_dim] where M = total copies across all gene trees
        Returns:
            [n_heads * embed_dim] aggregated representation
        """
        if embeddings.shape[0] == 0:
            return torch.zeros(self.n_heads * self.embed_dim,
                             device=self.queries.device)

        keys = self.key_proj(embeddings)     # [M, D]
        # Attention: [K, D] × [D, M] → [K, M]
        attn = torch.matmul(self.queries, keys.t()) / self.scale
        attn = torch.softmax(attn, dim=-1)   # [K, M] normalize over M per head
        # Weighted sum: [K, M] × [M, D] → [K, D]
        output = torch.matmul(attn, embeddings)
        return output.reshape(-1)  # [K*D]
