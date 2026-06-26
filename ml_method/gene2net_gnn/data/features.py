"""Hand-crafted features for gene tree discordance and copy number patterns."""
from collections import Counter
from itertools import combinations
from typing import Dict, List, Optional, Set, Tuple

from ete3 import Tree


def compute_copy_count_features(gene_trees: List[Tree], species: str) -> dict:
    """Compute copy count statistics for a species across all gene trees.

    Parameters
    ----------
    gene_trees : list of Tree
        Gene trees to analyse.
    species : str
        Target species name.

    Returns
    -------
    dict with keys:
        mean_copies, var_copies, mode_copies, p_absent, p_1_copy,
        p_2_copies, p_3plus_copies, max_copies
    """
    counts = []
    for tree in gene_trees:
        n = sum(1 for leaf in tree.get_leaves() if leaf.name == species)
        counts.append(n)

    n_trees = len(gene_trees)
    if n_trees == 0:
        return {
            "mean_copies": 0.0,
            "var_copies": 0.0,
            "mode_copies": 0,
            "p_absent": 1.0,
            "p_1_copy": 0.0,
            "p_2_copies": 0.0,
            "p_3plus_copies": 0.0,
            "max_copies": 0,
        }

    mean_copies = sum(counts) / n_trees
    var_copies = sum((c - mean_copies) ** 2 for c in counts) / n_trees

    counter = Counter(counts)
    mode_copies = counter.most_common(1)[0][0]

    p_absent = sum(1 for c in counts if c == 0) / n_trees
    p_1_copy = sum(1 for c in counts if c == 1) / n_trees
    p_2_copies = sum(1 for c in counts if c == 2) / n_trees
    p_3plus_copies = sum(1 for c in counts if c >= 3) / n_trees
    max_copies = max(counts) if counts else 0

    return {
        "mean_copies": mean_copies,
        "var_copies": var_copies,
        "mode_copies": mode_copies,
        "p_absent": p_absent,
        "p_1_copy": p_1_copy,
        "p_2_copies": p_2_copies,
        "p_3plus_copies": p_3plus_copies,
        "max_copies": max_copies,
    }


def compute_clustering_profile(
    gene_trees: List[Tree],
    species: str,
    all_species: Set[str],
) -> Dict[str, float]:
    """Compute co-clustering frequencies between a target species and all others.

    For each other species S, compute the fraction of gene trees where at least
    one leaf of 'species' has a sibling leaf that is species S (i.e., they share
    a direct common parent as leaves).

    Parameters
    ----------
    gene_trees : list of Tree
    species : str
        Target species.
    all_species : set of str
        All species to consider.

    Returns
    -------
    dict mapping each other species name to its co-clustering frequency (float in [0,1]).
    """
    other_species = all_species - {species}
    co_cluster_counts: Dict[str, int] = {s: 0 for s in other_species}

    for tree in gene_trees:
        # For this tree, find which other species co-cluster (are sisters) with any copy of target
        tree_co_clusters: Set[str] = set()
        for leaf in tree.get_leaves():
            if leaf.name != species:
                continue
            parent = leaf.up
            if parent is None:
                continue
            for sibling in parent.get_children():
                if sibling is leaf:
                    continue
                if sibling.is_leaf() and sibling.name in other_species:
                    tree_co_clusters.add(sibling.name)

        for s in tree_co_clusters:
            co_cluster_counts[s] += 1

    n_trees = len(gene_trees)
    if n_trees == 0:
        return {s: 0.0 for s in other_species}

    return {s: co_cluster_counts[s] / n_trees for s in other_species}


def compute_clustering_summary(
    gene_trees: List[Tree],
    species: str,
    all_species: Set[str],
) -> List[float]:
    """Compute fixed-size summary statistics of the clustering profile.

    Returns 5 values: [mean, std, max, min, median] of co-clustering
    frequencies. This is independent of the number of species, unlike
    the full clustering profile.
    """
    profile = compute_clustering_profile(gene_trees, species, all_species)
    values = list(profile.values())
    if not values:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / n
    std = var ** 0.5
    sorted_vals = sorted(values)
    median = sorted_vals[n // 2] if n % 2 == 1 else (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2
    return [mean, std, max(values), min(values), median]


def compute_species_tree_edge_features(
    species_tree: Tree,
    gene_trees: List[Tree],
) -> Dict[int, dict]:
    """Compute edge-level features for a species tree using gene tree concordance.

    For each edge (parent -> child) in the species tree, compute:
        - concordance_factor: fraction of gene trees supporting this edge's bipartition
        - branch_length: the edge's branch length
        - clade_size: number of leaves below this edge
        - depth: number of edges from root to this edge

    Bipartition concordance is computed at the species level: group all copies of
    a species in a gene tree to the same side as the species appears in the
    species tree bipartition.

    Parameters
    ----------
    species_tree : Tree
        The species tree (ete3 Tree with branch lengths).
    gene_trees : list of Tree

    Returns
    -------
    dict mapping edge_index (int) to a dict of features.
    """
    all_species = frozenset(species_tree.get_leaf_names())
    n_trees = len(gene_trees)

    # Build list of (edge_index, node) for all non-root nodes
    edge_nodes = []
    idx = 0
    for node in species_tree.traverse("preorder"):
        if node.is_root():
            continue
        edge_nodes.append((idx, node))
        idx += 1

    # Pre-compute species-level bipartitions for the species tree edges
    # bipartition: (below_set, above_set)
    species_bipartitions = []
    for edge_idx, node in edge_nodes:
        below = frozenset(node.get_leaf_names())
        above = all_species - below
        species_bipartitions.append((edge_idx, below, above))

    # For each gene tree, compute its species-level bipartitions
    def gene_tree_species_bipartitions(gtree: Tree) -> Set[frozenset]:
        """Return the set of bipartitions (as frozensets of below-species) for a gene tree."""
        gt_all_species = frozenset(leaf.name for leaf in gtree.get_leaves())
        biparts: Set[frozenset] = set()
        for node in gtree.traverse("preorder"):
            if node.is_root():
                continue
            below_sp = frozenset(leaf.name for leaf in node.get_leaves())
            # Only include if it's a proper bipartition (not the full set)
            if below_sp and below_sp != gt_all_species:
                biparts.add(below_sp)
        return biparts

    # Count concordance for each edge
    concordance_counts = [0] * len(edge_nodes)
    for gtree in gene_trees:
        gt_biparts = gene_tree_species_bipartitions(gtree)
        gt_species = frozenset(leaf.name for leaf in gtree.get_leaves())
        for i, (edge_idx, below, above) in enumerate(species_bipartitions):
            # Restrict bipartition to the species present in this gene tree
            restricted_below = below & gt_species
            restricted_above = above & gt_species
            if not restricted_below or not restricted_above:
                # Uninformative for this gene tree
                continue
            # Check if this restricted bipartition is supported
            if restricted_below in gt_biparts or restricted_above in gt_biparts:
                concordance_counts[i] += 1

    # Compute depth for each edge (number of edges from root)
    def compute_depth(node: Tree) -> int:
        depth = 0
        current = node
        while not current.is_root():
            depth += 1
            current = current.up
        return depth

    # Build output
    result: Dict[int, dict] = {}
    for i, (edge_idx, node) in enumerate(edge_nodes):
        n_informative = sum(
            1 for gtree in gene_trees
            if (species_bipartitions[i][1] & frozenset(leaf.name for leaf in gtree.get_leaves()))
            and (species_bipartitions[i][2] & frozenset(leaf.name for leaf in gtree.get_leaves()))
        )
        cf = concordance_counts[i] / n_informative if n_informative > 0 else 0.0

        result[edge_idx] = {
            "concordance_factor": cf,
            "branch_length": node.dist,
            "clade_size": len(node.get_leaves()),
            "depth": compute_depth(node),
        }

    return result


# ---------------------------------------------------------------------------
# WGD-detection features (per species-tree edge / clade).
#
# These target the signatures that copy-count and co-clustering summaries blur:
#   - duplication synchrony: do clade members get duplicated *together*? (WGD
#     duplicates the whole clade at once, so copy spikes are correlated;
#     small-scale duplication hits species independently)
#   - mirrored-sister fraction: the ((A,B),(A,B)) autopolyploidy signature
#   - copy-pair divergence: Ks-analogue from gene tree branch lengths
#     (substitutions per site). WGD copies share one "birthday" -> consistent,
#     synchronized divergence; scattered gene duplicates -> variable divergence.
# ---------------------------------------------------------------------------


def _index_gene_tree(tree: Tree):
    """Precompute per-tree structures reused by the detection features.

    Returns
    -------
    leaf_nodes : dict mapping species name -> list of leaf nodes (the copies).
    root_dist : dict mapping node -> cumulative branch-length distance to root.
    scale : float, mean root-to-leaf distance (used to normalise divergences so
        the per-sample substitution rate cancels out). Always > 0.
    dup_sister_sets : list of frozenset, species sets of nodes whose two
        children have identical species sets (the mirrored-sister pattern).
    """
    leaf_nodes: Dict[str, List] = {}
    for leaf in tree.get_leaves():
        leaf_nodes.setdefault(leaf.name, []).append(leaf)

    root_dist: Dict[object, float] = {}
    for node in tree.traverse("preorder"):
        if node.is_root():
            root_dist[node] = 0.0
        else:
            root_dist[node] = root_dist[node.up] + (node.dist or 0.0)

    leaf_root_dists = [root_dist[l] for l in tree.get_leaves()]
    scale = (sum(leaf_root_dists) / len(leaf_root_dists)) if leaf_root_dists else 1.0
    if scale <= 0:
        scale = 1.0

    dup_sister_sets: List[frozenset] = []
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        children = node.get_children()
        if len(children) != 2:
            continue
        s0 = frozenset(l.name for l in children[0].get_leaves())
        s1 = frozenset(l.name for l in children[1].get_leaves())
        if s0 == s1 and len(s0) >= 1:
            dup_sister_sets.append(s0)

    return leaf_nodes, root_dist, scale, dup_sister_sets


def _copy_pair_distance(copies: List, root_dist: Dict[object, float]) -> Optional[float]:
    """Smallest path distance between any two copies of a species in one tree.

    The minimum pair is taken because, with mixed WGD + gene duplications, the
    most recent split is the best single proxy for the WGD copy pair. Distance
    uses cumulative root distances: d(a,b) = root[a] + root[b] - 2*root[lca].
    Returns None if fewer than two copies.
    """
    if len(copies) < 2:
        return None
    best = None
    for a, b in combinations(copies, 2):
        lca = a.get_common_ancestor(b) if a is not b else a
        d = root_dist[a] + root_dist[b] - 2 * root_dist[lca]
        if best is None or d < best:
            best = d
    return best


def _mean_pairwise_correlation(indicator_matrix: List[List[int]]) -> float:
    """Mean pairwise Pearson correlation between columns (species) of a 0/1 matrix.

    Rows are gene trees, columns are clade species; entry is 1 if that species
    has >= 2 copies in that gene tree. Columns with zero variance are skipped.
    Returns 0.0 when fewer than two informative columns exist.
    """
    n_rows = len(indicator_matrix)
    if n_rows == 0:
        return 0.0
    n_cols = len(indicator_matrix[0])
    if n_cols < 2:
        return 0.0

    # Column means and variances
    means = [0.0] * n_cols
    for row in indicator_matrix:
        for j, v in enumerate(row):
            means[j] += v
    means = [m / n_rows for m in means]

    variances = [0.0] * n_cols
    for row in indicator_matrix:
        for j, v in enumerate(row):
            variances[j] += (v - means[j]) ** 2

    informative = [j for j in range(n_cols) if variances[j] > 1e-12]
    if len(informative) < 2:
        return 0.0

    corr_sum = 0.0
    n_pairs = 0
    for j, k in combinations(informative, 2):
        cov = 0.0
        for row in indicator_matrix:
            cov += (row[j] - means[j]) * (row[k] - means[k])
        denom = (variances[j] * variances[k]) ** 0.5
        if denom > 1e-12:
            corr_sum += cov / denom
            n_pairs += 1

    return corr_sum / n_pairs if n_pairs > 0 else 0.0


def compute_species_tree_edge_detection_features(
    species_tree: Tree,
    gene_trees: List[Tree],
) -> Dict[int, dict]:
    """Compute per-edge WGD-detection features keyed by the same edge index
    as ``compute_species_tree_edge_features`` (preorder, non-root nodes).

    For each edge with clade S (the species below it), returns:
        - dup_synchrony: mean pairwise correlation of "duplicated" indicators
          across clade species (high -> whole clade duplicated together -> WGD).
        - mirrored_sister_frac: fraction of gene trees where some subset of S
          appears as two identical sister subtrees (autopolyploidy signature).
        - copy_pair_div_mean: mean normalised copy-pair divergence across clade
          species and gene trees (Ks-analogue; copies' typical age).
        - copy_pair_div_cv: coefficient of variation of those divergences
          (low -> synchronised ages -> WGD; high -> scattered -> gene dup).
        - frac_clade_duplicated: mean fraction of clade species duplicated per
          gene tree (abundance of duplication within the clade).
    """
    # Edge enumeration must match compute_species_tree_edge_features.
    edge_nodes: List[Tuple[int, object]] = []
    idx = 0
    for node in species_tree.traverse("preorder"):
        if node.is_root():
            continue
        edge_nodes.append((idx, node))
        idx += 1

    # Precompute per gene tree once.
    indexed = [_index_gene_tree(gt) for gt in gene_trees]
    n_trees = len(gene_trees)

    result: Dict[int, dict] = {}
    for edge_idx, node in edge_nodes:
        clade = sorted(node.get_leaf_names())
        clade_set = frozenset(clade)

        # Per-tree duplication indicators (1 if species has >= 2 copies).
        indicator_matrix: List[List[int]] = []
        frac_dup_per_tree: List[float] = []
        # Copy-pair divergences gathered across clade species and trees.
        divergences: List[float] = []
        mirrored_count = 0

        for (leaf_nodes, root_dist, scale, dup_sets) in indexed:
            # Fixed columns = full clade, so rows stay aligned even when gene
            # loss makes some species absent in a tree (absent -> not duplicated).
            row = []
            n_present = 0
            n_dup = 0
            for sp in clade:
                copies = leaf_nodes.get(sp, [])
                if copies:
                    n_present += 1
                is_dup = 1 if len(copies) >= 2 else 0
                row.append(is_dup)
                if is_dup:
                    n_dup += 1
                    d = _copy_pair_distance(copies, root_dist)
                    if d is not None:
                        divergences.append(d / scale)
            if n_present:
                indicator_matrix.append(row)
                frac_dup_per_tree.append(n_dup / n_present)

            # Mirrored-sister: any duplicated sister set contained in this clade.
            if any(m and m <= clade_set for m in dup_sets):
                mirrored_count += 1

        dup_synchrony = _mean_pairwise_correlation(indicator_matrix)

        if divergences:
            n = len(divergences)
            mean_div = sum(divergences) / n
            var_div = sum((d - mean_div) ** 2 for d in divergences) / n
            std_div = var_div ** 0.5
            cv_div = (std_div / mean_div) if mean_div > 1e-12 else 0.0
        else:
            mean_div = 0.0
            cv_div = 0.0

        frac_clade_dup = (
            sum(frac_dup_per_tree) / len(frac_dup_per_tree)
            if frac_dup_per_tree else 0.0
        )
        mirrored_frac = mirrored_count / n_trees if n_trees > 0 else 0.0

        result[edge_idx] = {
            "dup_synchrony": dup_synchrony,
            "mirrored_sister_frac": mirrored_frac,
            "copy_pair_div_mean": mean_div,
            "copy_pair_div_cv": cv_div,
            "frac_clade_duplicated": frac_clade_dup,
        }

    return result


# ---------------------------------------------------------------------------
# Pairwise partner feature (for the allopolyploidy partner head).
#
# Allo signal is between two edges: clade i's duplicated copy lands next to the
# partner clade j, so in the gene trees species of clade i co-cluster with
# species of clade j. We measure species co-clustering, then aggregate it to an
# E x E pairwise feature. Computed from the stored gene-tree tensors (topology
# only) so no repackaging is needed. Imports torch lazily to keep the rest of
# this module dependency-light.
# ---------------------------------------------------------------------------


def species_coclustering_matrix(gene_tree_edge_indices, gene_tree_species_ids,
                                gene_tree_leaf_masks, n_species):
    """[n_species, n_species]: fraction of gene trees where two species appear
    as direct sibling leaves (same parent, both leaves)."""
    import torch
    C = torch.zeros(n_species, n_species)
    n_trees = len(gene_tree_edge_indices)
    if n_trees == 0:
        return C

    for ei, sp, lm in zip(gene_tree_edge_indices, gene_tree_species_ids, gene_tree_leaf_masks):
        # Group leaf-children species by parent (even-indexed parent->child edges).
        parent_leaf_species = {}
        for k in range(0, ei.shape[1], 2):
            p = int(ei[0, k]); c = int(ei[1, k])
            if bool(lm[c]) and int(sp[c]) >= 0:
                parent_leaf_species.setdefault(p, set()).add(int(sp[c]))
        seen = set()
        for specs in parent_leaf_species.values():
            specs = sorted(specs)
            for a_i in range(len(specs)):
                for b_i in range(a_i + 1, len(specs)):
                    seen.add((specs[a_i], specs[b_i]))
        for a, b in seen:
            C[a, b] += 1.0
            C[b, a] += 1.0

    return C / n_trees


def edge_clades_species(edge_index_pre, is_leaf, node_species):
    """Species-index set below each preorder edge's child node.

    edge_index_pre: reordered (preorder) edge_index.
    node_species: list/tensor length N; species index for leaf nodes else -1.
    """
    children = {}
    for k in range(0, edge_index_pre.shape[1], 2):
        p = int(edge_index_pre[0, k]); c = int(edge_index_pre[1, k])
        children.setdefault(p, []).append(c)

    memo = {}

    def below(n):
        if n in memo:
            return memo[n]
        if bool(is_leaf[n]):
            sidx = int(node_species[n])
            s = {sidx} if sidx >= 0 else set()
        else:
            s = set()
            for c in children.get(n, []):
                s |= below(c)
        memo[n] = s
        return s

    clades = []
    for k in range(0, edge_index_pre.shape[1], 2):
        clades.append(below(int(edge_index_pre[1, k])))
    return clades


def pairwise_partner_features(coclust, edge_clades):
    """[E, E, 2]: mean and max species co-clustering between clade i and clade j."""
    import torch
    E = len(edge_clades)
    feat = torch.zeros(E, E, 2)
    idx_lists = [torch.tensor(sorted(c), dtype=torch.long) for c in edge_clades]
    for i in range(E):
        ci = idx_lists[i]
        if ci.numel() == 0:
            continue
        for j in range(E):
            cj = idx_lists[j]
            if cj.numel() == 0:
                continue
            sub = coclust[ci][:, cj]
            feat[i, j, 0] = sub.mean()
            feat[i, j, 1] = sub.max()
    return feat


# ---------------------------------------------------------------------------
# Copy-aware cluster-support partner feature.
#
# The thin pairwise-sister feature above is symmetric and direct-sister only.
# This richer feature is copy-aware: for a species duplicated in a gene tree, one
# copy stays "home" (its local neighborhood is dominated by its own species-tree
# clade) and the other lands "away" near its allopolyploid partner. For each away
# copy we measure what fraction of its local neighborhood is each partner clade's
# species, accumulated over gene trees. This is the sharp allo signal the partner
# head needs, and it is computed from the same stored gene-tree tensors so no
# repackaging is required. Imports torch lazily like the functions above.
# ---------------------------------------------------------------------------


def gene_tree_copy_neighborhoods(edge_index, species_ids, leaf_mask, k):
    """Per-copy local neighborhoods for the duplicated species of one gene tree.

    A species is "duplicated" if it has >= 2 leaf copies in this gene tree. For
    each copy of each duplicated species, the neighborhood is the smallest
    enclosing clade grown upward while its total leaf count stays <= ``k`` (the
    copy's parent subtree at minimum). The focal copy itself is excluded.

    Parameters
    ----------
    edge_index : LongTensor [2, 2E]
        Undirected edges, parent->child at even indices (as stored on disk).
    species_ids : LongTensor [N]
        Species index per node; -1 for internal nodes.
    leaf_mask : BoolTensor [N]
    k : int
        Maximum neighborhood leaf count.

    Returns
    -------
    list of (species_id:int, neighborhood:Counter[species_id -> count])
        One entry per copy of each duplicated species.
    """
    children: Dict[int, List[int]] = {}
    parent: Dict[int, int] = {}
    for t in range(0, edge_index.shape[1], 2):
        p = int(edge_index[0, t]); c = int(edge_index[1, t])
        children.setdefault(p, []).append(c)
        parent[c] = p

    memo: Dict[int, List[int]] = {}

    def leaves_under(n: int) -> List[int]:
        if n in memo:
            return memo[n]
        if bool(leaf_mask[n]):
            r = [n]
        else:
            r = []
            for c in children.get(n, []):
                r += leaves_under(c)
        memo[n] = r
        return r

    leaf_nodes = [n for n in range(len(leaf_mask)) if bool(leaf_mask[n])]
    sp_count = Counter(int(species_ids[n]) for n in leaf_nodes if int(species_ids[n]) >= 0)
    duplicated = {s for s, ct in sp_count.items() if ct >= 2}

    results: List[Tuple[int, Counter]] = []
    for ln in leaf_nodes:
        a = int(species_ids[ln])
        if a < 0 or a not in duplicated or ln not in parent:
            continue
        # Grow the neighborhood upward while staying within the size cap.
        cand = parent[ln]
        while cand in parent and len(leaves_under(parent[cand])) <= k:
            cand = parent[cand]
        neigh: Counter = Counter()
        for x in leaves_under(cand):
            if x == ln:
                continue
            s = int(species_ids[x])
            if s >= 0:
                neigh[s] += 1
        results.append((a, neigh))
    return results


def copy_aware_cluster_support(gene_tree_edge_indices, gene_tree_species_ids,
                               gene_tree_leaf_masks, edge_clades, n_species,
                               k=10, away_threshold=0.5):
    """[E, E, 2]: copy-aware cluster support between clade i and clade j.

    For each away copy of a species in clade i (a copy whose neighborhood is less
    than ``away_threshold`` clade-i species), accumulate the fraction of its
    neighborhood that is clade-j species. Channel 0 is that sum normalized by the
    gene-tree count (support intensity); channel 1 is the per-copy max (peak
    support).
    """
    import torch
    E = len(edge_clades)
    support_sum = torch.zeros(E, E)
    support_max = torch.zeros(E, E)
    n_trees = len(gene_tree_edge_indices)
    if n_trees == 0 or E == 0:
        return torch.stack([support_sum, support_max], dim=-1)

    membership = torch.zeros(E, n_species)
    for i, clade in enumerate(edge_clades):
        for s in clade:
            if 0 <= s < n_species:
                membership[i, s] = 1.0

    for ei, sp, lm in zip(gene_tree_edge_indices, gene_tree_species_ids, gene_tree_leaf_masks):
        for a, neigh in gene_tree_copy_neighborhoods(ei, sp, lm, k):
            size = sum(neigh.values())
            if size == 0:
                continue
            counts = torch.zeros(n_species)
            for s, c in neigh.items():
                if 0 <= s < n_species:
                    counts[s] = c
            fracs = (membership @ counts) / size            # [E]: frac of neighborhood in each clade
            for i in range(E):
                if membership[i, a] == 0:
                    continue                                 # species a not in clade i
                if fracs[i].item() >= away_threshold:
                    continue                                 # home copy w.r.t. clade i
                support_sum[i] += fracs
                support_max[i] = torch.maximum(support_max[i], fracs)

    support_sum /= n_trees
    return torch.stack([support_sum, support_max], dim=-1)
