"""GIN-based gene tree encoder for WGD detection.

Encodes gene trees using Graph Isomorphism Network (GIN), pools per species
and per clade, then aggregates across gene trees to produce per-edge
representations for the species tree.

GIN is chosen over GCN because it's provably more expressive at distinguishing
tree topologies — critical for detecting WHERE paralogs sit in gene trees.
"""
import torch
import torch.nn as nn
from torch_geometric.nn import GINConv
from torch_geometric.utils import scatter
from typing import Dict, List, Set


class GeneTreeGIN(nn.Module):
    """Encode gene trees with GIN and aggregate per species tree edge.

    Architecture:
        1. Species embedding: each species gets a learnable embedding
        2. GIN layers: 2-layer GIN with residual connections
        3. Species pooling: scatter-based pooling per (gene tree, species)
        4. Clade pooling: vectorized gather + mean per species tree edge
        5. Gene tree aggregation: mean + std across gene trees per edge

    Args:
        n_species: Number of distinct species.
        embed_dim: Embedding dimension throughout.
        n_gin_layers: Number of GIN layers.
        dropout: Dropout rate.
    """

    def __init__(
        self,
        n_species: int,
        embed_dim: int = 64,
        n_gin_layers: int = 2,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.n_species = n_species
        self.embed_dim = embed_dim

        # Species embedding for gene tree leaves
        # +1 for internal nodes (index n_species)
        self.species_emb = nn.Embedding(n_species + 1, embed_dim)
        self.internal_idx = n_species

        # GIN layers
        self.gin_layers = nn.ModuleList()
        self.gin_norms = nn.ModuleList()
        for _ in range(n_gin_layers):
            mlp = nn.Sequential(
                nn.Linear(embed_dim, embed_dim),
                nn.ReLU(),
                nn.Linear(embed_dim, embed_dim),
            )
            self.gin_layers.append(GINConv(mlp, train_eps=True))
            self.gin_norms.append(nn.LayerNorm(embed_dim))

        self.dropout = nn.Dropout(dropout)

    def forward(
        self,
        gene_tree_edge_indices: List[torch.Tensor],
        gene_tree_species_ids: List[torch.Tensor],
        gene_tree_leaf_masks: List[torch.Tensor],
        edge_clades: Dict[int, Set[int]],
        n_edges: int,
    ) -> torch.Tensor:
        """Encode gene trees and aggregate per species tree edge.

        Args:
            gene_tree_edge_indices: List of [2, 2E_gt] per gene tree.
            gene_tree_species_ids: List of [N_gt] with species index per node
                (-1 for internal nodes).
            gene_tree_leaf_masks: List of [N_gt] boolean per gene tree.
            edge_clades: Dict mapping species tree edge index to set of
                species indices in the clade.
            n_edges: Number of undirected edges in species tree.

        Returns:
            edge_gene_features: [n_edges, 2 * embed_dim] aggregated gene tree
                features per species tree edge (mean + std).
        """
        device = self.species_emb.weight.device
        n_gt = len(gene_tree_edge_indices)

        if n_gt == 0:
            return torch.zeros(n_edges, 2 * self.embed_dim, device=device)

        # --- Step 1: Batch all gene trees via manual concatenation ---
        # (Avoids Data/Batch overhead from PyG)
        all_emb_ids = []
        all_sp_ids = []
        all_leaf_masks = []
        batch_ids = []
        all_ei = []
        offset = 0

        for i in range(n_gt):
            sp = gene_tree_species_ids[i]
            n_nodes = sp.shape[0]

            emb_ids = sp.clone()
            emb_ids[emb_ids < 0] = self.internal_idx
            emb_ids = emb_ids.clamp(0, self.internal_idx)

            all_emb_ids.append(emb_ids)
            all_sp_ids.append(sp)
            all_leaf_masks.append(gene_tree_leaf_masks[i])
            batch_ids.append(torch.full((n_nodes,), i, dtype=torch.long, device=device))
            all_ei.append(gene_tree_edge_indices[i] + offset)
            offset += n_nodes

        edge_index = torch.cat(all_ei, dim=1)
        emb_ids = torch.cat(all_emb_ids)
        sp_ids = torch.cat(all_sp_ids)
        leaf_mask = torch.cat(all_leaf_masks)
        batch_vec = torch.cat(batch_ids)

        # GIN message passing
        x = self.species_emb(emb_ids)
        for gin, norm in zip(self.gin_layers, self.gin_norms):
            x = x + gin(x, edge_index)  # residual
            x = norm(x)
            x = self.dropout(x)

        # --- Step 2: Vectorized species pooling with scatter ---
        # Instead of Python loops over nodes, use scatter_mean with composite key
        valid_leaf = leaf_mask & (sp_ids >= 0)
        leaf_idx = valid_leaf.nonzero(as_tuple=True)[0]

        if leaf_idx.shape[0] == 0:
            return torch.zeros(n_edges, 2 * self.embed_dim, device=device)

        leaf_x = x[leaf_idx]                # [N_leaves, embed_dim]
        leaf_gt = batch_vec[leaf_idx]        # [N_leaves] gene tree index
        leaf_sp = sp_ids[leaf_idx]           # [N_leaves] species index

        # Composite key for (gene_tree, species) pairs
        composite = leaf_gt * self.n_species + leaf_sp
        pool_size = n_gt * self.n_species

        # Scatter mean: pool all leaves of same species in same gene tree
        species_pool = scatter(leaf_x, composite, dim=0,
                               dim_size=pool_size, reduce='mean')
        # Track which (gt, sp) pairs have valid data
        sp_count = scatter(torch.ones(leaf_idx.shape[0], device=device),
                           composite, dim=0, dim_size=pool_size, reduce='sum')

        # Reshape to [n_gt, n_species, embed_dim] for efficient indexing
        sp_pool_2d = species_pool.view(n_gt, self.n_species, self.embed_dim)
        sp_valid_2d = sp_count.view(n_gt, self.n_species) > 0

        # --- Step 3: Clade pooling per edge, aggregated across gene trees ---
        # For each edge, gather clade species embeddings from all gene trees,
        # then compute mean + std across gene trees
        edge_features = torch.zeros(n_edges, 2 * self.embed_dim, device=device)

        for edge_idx in range(n_edges):
            clade_sp = edge_clades.get(edge_idx, set())
            if not clade_sp:
                continue

            sp_list = [sp for sp in clade_sp if sp < self.n_species]
            if not sp_list:
                continue

            sp_tensor = torch.tensor(sp_list, dtype=torch.long, device=device)

            # Gather: [n_gt, n_clade_sp, embed_dim]
            clade_embs = sp_pool_2d[:, sp_tensor]
            # Validity: [n_gt, n_clade_sp]
            clade_valid = sp_valid_2d[:, sp_tensor]

            # Number of valid species per gene tree
            n_valid_per_gt = clade_valid.sum(dim=1)  # [n_gt]
            gt_has_data = n_valid_per_gt > 0

            if gt_has_data.sum() == 0:
                continue

            # Mean pool over clade species per gene tree (masking invalid)
            masked_embs = clade_embs * clade_valid.unsqueeze(-1).float()
            gt_means = masked_embs.sum(dim=1) / n_valid_per_gt.unsqueeze(-1).clamp(min=1)
            # gt_means: [n_gt, embed_dim]

            # Keep only gene trees that had data
            valid_means = gt_means[gt_has_data]  # [n_valid_gt, embed_dim]

            # Aggregate across gene trees: mean + std
            mean_emb = valid_means.mean(dim=0)
            if valid_means.shape[0] > 1:
                std_emb = valid_means.std(dim=0)
            else:
                std_emb = torch.zeros(self.embed_dim, device=device)

            edge_features[edge_idx] = torch.cat([mean_emb, std_emb])

        return edge_features
