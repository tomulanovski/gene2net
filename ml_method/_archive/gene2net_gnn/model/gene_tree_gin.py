"""Gene tree feature extraction for WGD detection.

Computes per-edge contrast features from gene tree copy counts.
The key WGD signal is the CONTRAST between copy counts below an edge
vs outside — WGD edges show ~2x copies below vs ~1x outside.

This replaces the GIN encoder which learned topology but couldn't
distinguish WGD edges from their ancestors (both show elevated copies).

Also retains a GIN encoder for learned topology features that may
capture patterns beyond simple copy counts.
"""
import torch
import torch.nn as nn
from torch_geometric.nn import GINConv
from torch_geometric.utils import scatter
from typing import Dict, List, Set

# 10 contrast features + 6 copy count stats from GIN
N_CONTRAST_FEATURES = 10


class GeneTreeGIN(nn.Module):
    """Encode gene trees with GIN + compute contrast features per edge.

    Two complementary feature sources:
        1. GIN learned embeddings: topology-aware features (mean+std per edge)
        2. Contrast features: below-vs-outside copy count comparison

    The contrast features capture the primary WGD signal that the GIN
    pooling erases: the ratio of copy counts in the child clade vs the
    rest of the tree.

    Output dimension: 2 * embed_dim + N_CONTRAST_FEATURES
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
        self.output_dim = 2 * embed_dim + N_CONTRAST_FEATURES

        # Species embedding for gene tree leaves
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

    def _compute_copy_counts(
        self,
        gene_tree_species_ids: List[torch.Tensor],
        gene_tree_leaf_masks: List[torch.Tensor],
        device: torch.device,
    ):
        """Compute copy count per (gene_tree, species) via scatter.

        Returns:
            sp_count_2d: [n_gt, n_species] copy counts.
            sp_valid_2d: [n_gt, n_species] boolean validity mask.
        """
        n_gt = len(gene_tree_species_ids)
        all_sp = []
        all_gt = []

        for i in range(n_gt):
            sp = gene_tree_species_ids[i]
            mask = gene_tree_leaf_masks[i]
            valid = mask & (sp >= 0)
            valid_sp = sp[valid]
            all_sp.append(valid_sp)
            all_gt.append(torch.full_like(valid_sp, i))

        all_sp_cat = torch.cat(all_sp)
        all_gt_cat = torch.cat(all_gt)
        composite = all_gt_cat * self.n_species + all_sp_cat
        pool_size = n_gt * self.n_species

        sp_count = scatter(torch.ones(all_sp_cat.shape[0], device=device),
                           composite, dim=0, dim_size=pool_size, reduce='sum')
        sp_count_2d = sp_count.view(n_gt, self.n_species)
        sp_valid_2d = sp_count_2d > 0

        return sp_count_2d, sp_valid_2d

    def _compute_contrast_features(
        self,
        sp_count_2d: torch.Tensor,
        sp_valid_2d: torch.Tensor,
        edge_clades: Dict[int, Set[int]],
        n_edges: int,
        device: torch.device,
    ) -> torch.Tensor:
        """Compute below-vs-outside contrast features per edge.

        For each edge, for each gene tree:
        - avg_copies_below: mean copy count for species in clade
        - avg_copies_outside: mean copy count for species outside clade
        - copy_ratio: below / outside (WGD → ~2.0, non-WGD → ~1.0)
        - frac_dup_below: fraction of clade species with >1 copy
        - frac_dup_outside: fraction of outside species with >1 copy

        Then aggregate across gene trees: mean + std of each.

        Returns:
            contrast_features: [n_edges, 10]
        """
        n_gt = sp_count_2d.shape[0]
        all_species = set(range(self.n_species))
        contrast = torch.zeros(n_edges, N_CONTRAST_FEATURES, device=device)

        for edge_idx in range(n_edges):
            clade_sp = edge_clades.get(edge_idx, set())
            if not clade_sp:
                continue

            below = [sp for sp in clade_sp if sp < self.n_species]
            outside = [sp for sp in (all_species - clade_sp) if sp < self.n_species]

            if not below or not outside:
                continue

            below_t = torch.tensor(below, dtype=torch.long, device=device)
            outside_t = torch.tensor(outside, dtype=torch.long, device=device)

            # Copy counts: [n_gt, n_below] and [n_gt, n_outside]
            counts_below = sp_count_2d[:, below_t]
            counts_outside = sp_count_2d[:, outside_t]
            valid_below = sp_valid_2d[:, below_t]
            valid_outside = sp_valid_2d[:, outside_t]

            n_valid_below = valid_below.sum(dim=1).float()
            n_valid_outside = valid_outside.sum(dim=1).float()
            gt_has_data = (n_valid_below > 0) & (n_valid_outside > 0)

            if gt_has_data.sum() == 0:
                continue

            # Per gene tree stats (masking invalid species)
            avg_below = (counts_below * valid_below.float()).sum(dim=1) / n_valid_below.clamp(min=1)
            avg_outside = (counts_outside * valid_outside.float()).sum(dim=1) / n_valid_outside.clamp(min=1)
            copy_ratio = avg_below / avg_outside.clamp(min=0.1)
            frac_dup_below = ((counts_below > 1) & valid_below).sum(dim=1).float() / n_valid_below.clamp(min=1)
            frac_dup_outside = ((counts_outside > 1) & valid_outside).sum(dim=1).float() / n_valid_outside.clamp(min=1)

            # Filter to gene trees with data on both sides
            ab = avg_below[gt_has_data]
            cr = copy_ratio[gt_has_data]
            fdb = frac_dup_below[gt_has_data]
            fdo = frac_dup_outside[gt_has_data]
            dup_contrast = fdb - fdo

            n = gt_has_data.sum().float()
            z = torch.tensor(0.0, device=device)

            contrast[edge_idx] = torch.stack([
                ab.mean(),
                ab.std() if n > 1 else z,
                cr.mean(),
                cr.std() if n > 1 else z,
                fdb.mean(),
                fdb.std() if n > 1 else z,
                fdo.mean(),
                fdo.std() if n > 1 else z,
                dup_contrast.mean(),
                dup_contrast.std() if n > 1 else z,
            ])

        return contrast

    def forward(
        self,
        gene_tree_edge_indices: List[torch.Tensor],
        gene_tree_species_ids: List[torch.Tensor],
        gene_tree_leaf_masks: List[torch.Tensor],
        edge_clades: Dict[int, Set[int]],
        n_edges: int,
    ) -> torch.Tensor:
        """Encode gene trees and compute contrast features per edge.

        Returns:
            edge_features: [n_edges, 2*embed_dim + 10] with GIN embeddings
                and contrast features.
        """
        device = self.species_emb.weight.device
        n_gt = len(gene_tree_edge_indices)

        if n_gt == 0:
            return torch.zeros(n_edges, self.output_dim, device=device)

        # --- Copy counts (used by both GIN pooling and contrast features) ---
        sp_count_2d, sp_valid_2d = self._compute_copy_counts(
            gene_tree_species_ids, gene_tree_leaf_masks, device
        )

        # --- Contrast features (the key WGD signal) ---
        contrast_feats = self._compute_contrast_features(
            sp_count_2d, sp_valid_2d, edge_clades, n_edges, device
        )

        # --- GIN encoding ---
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
            x = x + gin(x, edge_index)
            x = norm(x)
            x = self.dropout(x)

        # Vectorized species pooling with scatter
        valid_leaf = leaf_mask & (sp_ids >= 0)
        leaf_idx = valid_leaf.nonzero(as_tuple=True)[0]

        gin_features = torch.zeros(n_edges, 2 * self.embed_dim, device=device)

        if leaf_idx.shape[0] > 0:
            leaf_x = x[leaf_idx]
            leaf_gt = batch_vec[leaf_idx]
            leaf_sp = sp_ids[leaf_idx]

            composite = leaf_gt * self.n_species + leaf_sp
            pool_size = n_gt * self.n_species
            species_pool = scatter(leaf_x, composite, dim=0,
                                    dim_size=pool_size, reduce='mean')
            sp_pool_2d = species_pool.view(n_gt, self.n_species, self.embed_dim)

            for edge_idx in range(n_edges):
                clade_sp = edge_clades.get(edge_idx, set())
                if not clade_sp:
                    continue

                sp_list = [sp for sp in clade_sp if sp < self.n_species]
                if not sp_list:
                    continue

                sp_tensor = torch.tensor(sp_list, dtype=torch.long, device=device)
                clade_embs = sp_pool_2d[:, sp_tensor]
                clade_valid = sp_valid_2d[:, sp_tensor]

                n_valid_per_gt = clade_valid.sum(dim=1)
                gt_has_data = n_valid_per_gt > 0

                if gt_has_data.sum() == 0:
                    continue

                masked_embs = clade_embs * clade_valid.unsqueeze(-1).float()
                gt_means = masked_embs.sum(dim=1) / n_valid_per_gt.unsqueeze(-1).clamp(min=1)
                valid_means = gt_means[gt_has_data]

                mean_emb = valid_means.mean(dim=0)
                if valid_means.shape[0] > 1:
                    std_emb = valid_means.std(dim=0)
                else:
                    std_emb = torch.zeros(self.embed_dim, device=device)

                gin_features[edge_idx] = torch.cat([mean_emb, std_emb])

        # Concatenate GIN embeddings + contrast features
        return torch.cat([gin_features, contrast_feats], dim=-1)
