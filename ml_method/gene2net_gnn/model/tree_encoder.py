"""Gene Tree Encoder: bottom-up then top-down message passing on tree structure."""
import torch
import torch.nn as nn
from torch_geometric.nn import MessagePassing


class TreeConv(MessagePassing):
    """One-directional message passing on a tree."""

    def __init__(self, in_dim: int, out_dim: int):
        super().__init__(aggr="add", flow="source_to_target")
        self.mlp = nn.Sequential(
            nn.Linear(in_dim + out_dim + 1, out_dim),  # source + target + branch_length
            nn.ReLU(),
            nn.Linear(out_dim, out_dim),
        )

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor, edge_attr=None) -> torch.Tensor:
        return self.propagate(edge_index, x=x, edge_attr=edge_attr)

    def message(self, x_j: torch.Tensor, x_i: torch.Tensor, edge_attr) -> torch.Tensor:
        """x_j = source features, x_i = target features, edge_attr = branch lengths."""
        if edge_attr is not None:
            bl = edge_attr.unsqueeze(-1)
        else:
            bl = torch.zeros(x_j.shape[0], 1, device=x_j.device)
        return self.mlp(torch.cat([x_j, x_i, bl], dim=-1))


class GeneTreeEncoder(nn.Module):
    """Encode a gene tree into per-leaf embeddings.

    Uses a bottom-up pass (child→parent) followed by a top-down pass
    (parent→child).  The bottom-up and top-down representations for each
    leaf node are concatenated and projected to *embed_dim*.

    Args:
        n_species: Number of distinct species (determines embedding table size).
        embed_dim: Hidden and output embedding dimensionality.
        n_rounds: Number of message-passing rounds in each direction.
    """

    def __init__(self, n_species: int, embed_dim: int = 64, n_rounds: int = 2):
        super().__init__()
        self.embed_dim = embed_dim
        self.species_embedding = nn.Embedding(n_species, embed_dim)
        self.internal_embedding = nn.Parameter(torch.randn(embed_dim))

        # Bottom-up layers (child → parent messages)
        self.bottom_up_layers = nn.ModuleList(
            [TreeConv(embed_dim, embed_dim) for _ in range(n_rounds)]
        )
        # Top-down layers (parent → child messages)
        self.top_down_layers = nn.ModuleList(
            [TreeConv(embed_dim, embed_dim) for _ in range(n_rounds)]
        )

        # Final projection from concatenated [bu || td] → embed_dim
        self.output_proj = nn.Linear(2 * embed_dim, embed_dim)

    def forward(
        self,
        edge_index: torch.Tensor,
        species_ids: torch.Tensor,
        branch_lengths: torch.Tensor,
    ) -> torch.Tensor:
        """Forward pass.

        Args:
            edge_index: [2, 2E] undirected edges produced by
                ``tree_to_edge_index``.  For every parent-child pair the
                function adds [parent, child] then [child, parent], so:
                  - even column indices (0, 2, …) → parent→child edges
                  - odd column indices  (1, 3, …) → child→parent edges
            species_ids: [N] species index per node; -1 for internal nodes.
            branch_lengths: [2E] branch length per directed edge (same ordering
                as edge_index columns).

        Returns:
            leaf_embeddings: [n_leaves, embed_dim] – one embedding per leaf,
                in the same order as leaf nodes appear in ``species_ids``.
        """
        n_nodes = species_ids.shape[0]

        # --- Node initialisation -------------------------------------------
        # All nodes start with the shared internal-node parameter; leaf nodes
        # are overwritten with their species embedding.
        x = self.internal_embedding.unsqueeze(0).expand(n_nodes, -1).clone()
        leaf_mask = species_ids >= 0
        x[leaf_mask] = self.species_embedding(species_ids[leaf_mask])

        # --- Split directed edges -------------------------------------------
        # even columns: parent→child  |  odd columns: child→parent
        parent_to_child = edge_index[:, 0::2]   # [2, E]
        child_to_parent = edge_index[:, 1::2]   # [2, E]
        bl_p2c = branch_lengths[0::2]
        bl_c2p = branch_lengths[1::2]

        # --- Bottom-up pass (child→parent) ----------------------------------
        x_bu = x.clone()
        for layer in self.bottom_up_layers:
            x_bu = x_bu + layer(x_bu, child_to_parent, bl_c2p)

        # --- Top-down pass (parent→child) -----------------------------------
        x_td = x.clone()
        for layer in self.top_down_layers:
            x_td = x_td + layer(x_td, parent_to_child, bl_p2c)

        # --- Output: leaves only --------------------------------------------
        combined = torch.cat([x_bu[leaf_mask], x_td[leaf_mask]], dim=-1)
        leaf_embeddings = self.output_proj(combined)

        return leaf_embeddings
