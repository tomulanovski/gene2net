"""GIN-based gene tree encoder for WGD detection.

Encodes gene trees using Graph Isomorphism Network (GIN), pools per species
and per clade, then aggregates across gene trees to produce per-edge
representations for the species tree.

GIN is chosen over GCN because it's provably more expressive at distinguishing
tree topologies — critical for detecting WHERE paralogs sit in gene trees.
"""
import torch
import torch.nn as nn
from torch_geometric.nn import GINConv, global_mean_pool
from torch_geometric.data import Data, Batch
from typing import Dict, List, Set


class GeneTreeGIN(nn.Module):
    """Encode gene trees with GIN and aggregate per species tree edge.

    Architecture:
        1. Species embedding: each species gets a learnable embedding
        2. GIN layers: 2-layer GIN with residual connections
        3. Species pooling: pool leaf embeddings per species per gene tree
        4. Clade pooling: pool species embeddings per clade (species tree edge)
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
        n_gene_trees = len(gene_tree_edge_indices)

        if n_gene_trees == 0:
            return torch.zeros(n_edges, 2 * self.embed_dim, device=device)

        # --- Step 1: Batch all gene trees and run GIN ---
        data_list = []
        tree_offsets = []  # track node offsets per tree in the batch
        offset = 0

        for gt_idx in range(n_gene_trees):
            ei = gene_tree_edge_indices[gt_idx].to(device)
            sp_ids = gene_tree_species_ids[gt_idx].to(device)
            n_nodes = sp_ids.shape[0]

            # Map species_ids to embedding indices
            # -1 (internal) → self.internal_idx
            emb_ids = sp_ids.clone()
            emb_ids[emb_ids < 0] = self.internal_idx
            # Clamp to valid range
            emb_ids = emb_ids.clamp(0, self.internal_idx)

            data_list.append(Data(
                edge_index=ei,
                emb_ids=emb_ids,
                sp_ids=sp_ids,
                leaf_mask=gene_tree_leaf_masks[gt_idx].to(device),
                num_nodes=n_nodes,
            ))
            tree_offsets.append(offset)
            offset += n_nodes

        batch = Batch.from_data_list(data_list)

        # Initial node features from species embeddings
        x = self.species_emb(batch.emb_ids)

        # GIN message passing
        for gin, norm in zip(self.gin_layers, self.gin_norms):
            x = x + gin(x, batch.edge_index)  # residual
            x = norm(x)
            x = self.dropout(x)

        # --- Step 2: Pool per species per gene tree ---
        # For each gene tree, for each species, mean-pool leaf embeddings
        # Result: per_species_emb[gt_idx][sp_idx] = embedding

        per_gt_species_emb = []
        for gt_idx in range(n_gene_trees):
            start = tree_offsets[gt_idx]
            end = tree_offsets[gt_idx + 1] if gt_idx + 1 < n_gene_trees else x.shape[0]
            gt_x = x[start:end]
            gt_sp_ids = batch.sp_ids[start:end]
            gt_leaf_mask = batch.leaf_mask[start:end]

            # Pool leaves by species
            species_embs = {}
            for node_idx in range(gt_x.shape[0]):
                if not gt_leaf_mask[node_idx]:
                    continue
                sp = gt_sp_ids[node_idx].item()
                if sp < 0:
                    continue
                if sp not in species_embs:
                    species_embs[sp] = []
                species_embs[sp].append(gt_x[node_idx])

            # Mean pool per species
            pooled = {}
            for sp, embs in species_embs.items():
                pooled[sp] = torch.stack(embs).mean(dim=0)

            per_gt_species_emb.append(pooled)

        # --- Step 3: Pool per clade and aggregate across gene trees ---
        # For each species tree edge, collect clade embeddings from each
        # gene tree, then compute mean + std across gene trees

        edge_features = torch.zeros(n_edges, 2 * self.embed_dim, device=device)

        for edge_idx in range(n_edges):
            clade_sp = edge_clades.get(edge_idx, set())
            if not clade_sp:
                continue

            gt_clade_embs = []
            for gt_idx in range(n_gene_trees):
                sp_embs = per_gt_species_emb[gt_idx]
                # Collect embeddings for species in this clade
                present = [sp_embs[sp] for sp in clade_sp if sp in sp_embs]
                if present:
                    # Mean pool across clade species
                    gt_clade_embs.append(torch.stack(present).mean(dim=0))

            if gt_clade_embs:
                stacked = torch.stack(gt_clade_embs)  # [n_valid_gt, embed_dim]
                mean_emb = stacked.mean(dim=0)
                std_emb = stacked.std(dim=0) if stacked.shape[0] > 1 else torch.zeros_like(mean_emb)
                edge_features[edge_idx] = torch.cat([mean_emb, std_emb])

        return edge_features
