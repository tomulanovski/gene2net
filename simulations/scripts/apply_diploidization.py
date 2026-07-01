#!/usr/bin/env python3
"""Apply post-WGD diploidization (fractionation) to simulated gene trees + alignments.

Background
----------
SimPhy is fed a MUL species tree (``tree_10_mil.nex``) in which each polyploid
subgenome is a separate, uniquely-named lineage (``X_1``, ``X_2``).  Every gene
therefore starts with BOTH homeolog copies present.  Real polyploids undergo
rapid diploidization: one homeolog of many genes is lost (fractionation).

This script models that as a post-processing pass.  For each WGD event a
per-event retention rate ``q_e`` is drawn from a configurable distribution; for
each gene, with probability ``1 - q_e`` one homeolog copy of a species is dropped
(unbiased).  The same removals are applied to the gene tree and its matched
alignment so downstream method-prep consumes consistent data.  Ground-truth
networks are left untouched.

Copies are mapped to WGD events by reading the real ``tree_10_mil.nex``: copy
suffixes are stripped to recover the MUL structure, which is folded (Holm 2006,
via ``reticulate_tree.py``) to identify events and their isomorphic copy groups.
"""
import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from ete3 import Tree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reticulate_tree import ReticulateTree


def species_tree_leaf(gene_tip):
    """Recover the species-tree leaf name from a SimPhy gene-tree tip.

    SimPhy appends ``_<locus>_<individual>`` to each species-tree leaf name, and
    species names are underscore-free, so stripping the last two ``_``-fields
    yields the species-tree leaf (which keeps the polyploid ``_N`` copy suffix).

        ``Lamiumconfertum1_1_0_0`` -> ``Lamiumconfertum1_1``
        ``Lamiumalbum2_0_0``       -> ``Lamiumalbum2``
    """
    return gene_tip.rsplit("_", 2)[0]


def split_copy(leaf):
    """Split a species-tree leaf into (base species name, copy index or None).

    Species names are underscore-free, so any ``_`` marks the polyploid copy
    suffix (1-based).  Diploids have no suffix.

        ``Lamiumhybridum3_1`` -> ("Lamiumhybridum3", 1)
        ``Lamiumalbum2``      -> ("Lamiumalbum2", None)
    """
    if "_" in leaf:
        base, idx = leaf.rsplit("_", 1)
        return base, int(idx)
    return leaf, None


def _iso_leaf_map(a, b, canonical):
    """Map uniq leaf names of subtree ``a`` to the corresponding leaves of the
    isomorphic subtree ``b`` (children paired by canonical form)."""
    if a.is_leaf():
        return {a.uniq: b.uniq}
    a_children = sorted(a.children, key=lambda n: canonical[n])
    b_children = sorted(b.children, key=lambda n: canonical[n])
    mapping = {}
    for ca, cb in zip(a_children, b_children):
        mapping.update(_iso_leaf_map(ca, cb, canonical))
    return mapping


def build_events(species_tree):
    """Fold a unique-named species tree into its WGD events.

    ``species_tree`` is a Newick string or ete3 ``Tree`` whose polyploid copies
    carry ``_N`` suffixes (as in ``tree_10_mil.nex``).  Copy suffixes are stripped
    to recover the MUL structure, which is folded with the Holm 2006 grouping
    (mirroring ``ReticulateTree.multree_to_dag``): at each height, top-down,
    subtrees with identical canonical form are an isomorphic group = one WGD
    event.  After an event the non-leader copies are set aside so nested events
    inside the surviving (leader) copy are still found.

    Because the non-leader copies are set aside, an inner event is detected only
    in the leader block; it is then **mirrored** onto the other blocks of its
    immediately-enclosing event (via the block isomorphism) so that fractionation
    applies whichever outer block survives.  Mirrored inner events share the
    ``group`` (hence the drawn retention rate) of their source event.

    Returns a list of events, each a dict::

        {"sides":         [frozenset(species-tree leaf names), ...],  # the copies
         "autopolyploidy": bool,
         "height":         int,
         "group":          int}   # events sharing a group share one q_e
    """
    t = species_tree if isinstance(species_tree, Tree) else Tree(species_tree, format=1)
    t = t.copy()
    for leaf in t.iter_leaves():
        leaf.add_feature("uniq", leaf.name)
        leaf.name = split_copy(leaf.name)[0]

    height_map = ReticulateTree._compute_height_map(t)
    canonical = ReticulateTree._compute_canonical_forms(t)

    removed = set()  # ids of nodes excluded by an earlier (outer) fold
    events = []
    group = 0
    for h in sorted(height_map.keys(), reverse=True):
        nodes = [n for n in height_map[h] if id(n) not in removed]
        groups = defaultdict(list)
        for n in nodes:
            groups[canonical[n]].append(n)
        for members in groups.values():
            if len(members) <= 1:
                continue
            # stable leader: sort by sorted uniq-leaf-set
            members.sort(key=lambda m: sorted(l.uniq for l in m.iter_leaves()))
            sides = [frozenset(l.uniq for l in m.iter_leaves()) for m in members]
            parents = [m.up for m in members]
            autopolyploidy = any(parents.count(p) > 1 for p in parents)
            events.append({
                "sides": sides, "autopolyploidy": autopolyploidy,
                "height": h, "group": group, "_nodes": members,
            })
            group += 1
            for iso in members[1:]:
                removed.update(id(d) for d in iso.traverse())

    # Mirror each inner event onto the other blocks of its immediate enclosing event.
    for inner in list(events):
        inner_leaves = frozenset().union(*inner["sides"])
        enclosing = [
            e for e in events
            if e is not inner and e["height"] > inner["height"]
            and inner_leaves <= e["sides"][0]  # inside e's leader block
        ]
        if not enclosing:
            continue
        outer = min(enclosing, key=lambda e: len(e["sides"][0]))  # immediate parent
        leader_node = outer["_nodes"][0]
        for other_node in outer["_nodes"][1:]:
            imap = _iso_leaf_map(leader_node, other_node, canonical)
            if not all(l in imap for side in inner["sides"] for l in side):
                continue
            events.append({
                "sides": [frozenset(imap[l] for l in side) for side in inner["sides"]],
                "autopolyploidy": inner["autopolyploidy"],
                "height": inner["height"], "group": inner["group"], "_nodes": None,
            })

    for e in events:
        e.pop("_nodes", None)
    return events


def draw_retention_rates(n_events, dist="beta", rng=None, **params):
    """Draw one retention rate ``q_e`` per WGD event.

    ``dist`` is one of:

    - ``"fixed"``   : ``value`` for every event.
    - ``"uniform"`` : ``q ~ Uniform(low, high)``.
    - ``"beta"``    : ``q ~ Beta(alpha, beta)``; alternatively give ``mean`` and
      ``concentration`` (``alpha = mean * concentration``,
      ``beta = (1 - mean) * concentration``).

    ``rng`` is a ``numpy`` Generator (pass a seeded one for reproducibility).
    """
    if rng is None:
        rng = np.random.default_rng()

    if dist == "fixed":
        return [float(params["value"])] * n_events
    if dist == "uniform":
        return [float(x) for x in rng.uniform(params["low"], params["high"], n_events)]
    if dist == "beta":
        if "mean" in params:
            conc = params["concentration"]
            alpha = params["mean"] * conc
            beta = (1.0 - params["mean"]) * conc
        else:
            alpha, beta = params["alpha"], params["beta"]
        if alpha <= 0 or beta <= 0:  # degenerate Beta -> point mass
            return [1.0 if beta <= 0 else 0.0] * n_events
        return [float(x) for x in rng.beta(alpha, beta, n_events)]
    raise ValueError(f"unknown distribution: {dist!r}")


def plan_removals(events, q_by_group, present_leaves, rng):
    """Decide which species-tree leaves to fractionate out of one gene.

    For each event (outer -> inner) and each species, the species' copies are its
    leaves on each side (a "block"; usually one leaf, but a whole sub-block for an
    outer event of a nested polyploid).  Among the blocks still present, one is
    kept and each of the others is dropped with probability ``1 - q_e``.  The kept
    block is chosen uniformly at random (unbiased fractionation), so a species
    always retains at least one copy per gene.  ``q_by_group`` is indexed by each
    event's ``group`` (mirrored nested events share their source's group/rate).

    Returns the set of species-tree leaf names to remove from this gene.
    """
    removed = set()
    for event in events:
        q = q_by_group[event["group"]]
        per_species = defaultdict(list)  # species -> [block per side]
        for side in event["sides"]:
            by_species = defaultdict(set)
            for leaf in side:
                by_species[split_copy(leaf)[0]].add(leaf)
            for species, block in by_species.items():
                per_species[species].append(block)

        for blocks in per_species.values():
            present = [
                {l for l in block if l in present_leaves and l not in removed}
                for block in blocks
            ]
            present = [b for b in present if b]
            if len(present) < 2:
                continue
            order = list(rng.permutation(len(present)))
            for j in order[1:]:  # order[0] is the kept block
                if rng.random() < 1.0 - q:
                    removed.update(present[j])
    return removed


def prune_gene_tree(newick, removed_leaves):
    """Drop gene-tree tips whose species-tree leaf is in ``removed_leaves``.

    Branch lengths of surviving lineages are preserved (collapsed internal nodes
    have their lengths merged into the kept child).
    """
    t = Tree(newick, format=1)
    keep = [
        l.name for l in t.iter_leaves()
        if species_tree_leaf(l.name) not in removed_leaves
    ]
    t.prune(keep, preserve_branch_length=True)
    return t.write(format=5)


def read_phylip(text):
    """Parse PHYLIP sequential text into an ordered {seq_id: sequence} dict."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    seqs = {}
    for line in lines[1:]:  # line 0 is "ntax nchar"
        seq_id, sequence = line.split()
        seqs[seq_id] = sequence
    return seqs


def write_phylip(seqs):
    """Serialize an ordered {seq_id: sequence} dict to PHYLIP sequential text."""
    nchar = len(next(iter(seqs.values()))) if seqs else 0
    width = max((len(sid) for sid in seqs), default=0) + 2
    lines = [f"{len(seqs)} {nchar}"]
    lines += [f"{sid.ljust(width)}{seq}" for sid, seq in seqs.items()]
    return "\n".join(lines) + "\n"


def prune_alignment(text, removed_leaves):
    """Drop sequences whose species-tree leaf is in ``removed_leaves``.

    The ``ntax`` count in the header is updated; ``nchar`` is unchanged.
    """
    seqs = read_phylip(text)
    kept = {
        sid: seq for sid, seq in seqs.items()
        if species_tree_leaf(sid) not in removed_leaves
    }
    return write_phylip(kept)


def read_species_tree(path):
    """Read a species tree as a Newick string from a NEXUS or Newick file."""
    text = Path(path).read_text()
    if text.lstrip().upper().startswith("#NEXUS"):
        for line in text.splitlines():
            line = line.strip()
            if line.lower().startswith("tree") and "=" in line:
                return line.split("=", 1)[1].strip()
        raise ValueError(f"no tree found in NEXUS file: {path}")
    return text.strip()


def diploidize_replicate(species_tree_path, in_dir, out_dir, *,
                         dist, dist_params, seed):
    """Diploidize every gene of one replicate.

    Reads the species tree, identifies WGD events, draws one retention rate per
    event, then fractionates each ``g_trees{NNNN}.trees`` gene tree and its
    matched ``alignments/alignment_{NNNN}.phy`` consistently, writing pruned
    copies under ``out_dir`` (same layout).  Returns a summary dict (also useful
    for the on-disk ``diploidization_summary.json``).
    """
    in_dir, out_dir = Path(in_dir), Path(out_dir)
    events = build_events(read_species_tree(species_tree_path))
    n_groups = 1 + max((e["group"] for e in events), default=-1)
    rng = np.random.default_rng(seed)
    q_by_group = draw_retention_rates(n_groups, dist=dist, rng=rng, **dist_params)

    (out_dir / "alignments").mkdir(parents=True, exist_ok=True)
    copy_numbers = defaultdict(lambda: defaultdict(int))  # species -> {count: n_genes}
    n_genes = 0

    for gene_path in sorted(in_dir.glob("g_trees*.trees")):
        if gene_path.suffix != ".trees":  # skip g_trees*.trees.log
            continue
        gene_num = gene_path.stem.replace("g_trees", "")
        aln_path = in_dir / "alignments" / f"alignment_{gene_num}.phy"
        if not aln_path.exists():
            print(f"WARNING: no alignment for {gene_path.name}, skipping")
            continue

        newick = gene_path.read_text().strip()
        present = {species_tree_leaf(l.name) for l in Tree(newick, format=1).iter_leaves()}
        removed = plan_removals(events, q_by_group, present, rng)

        (out_dir / gene_path.name).write_text(prune_gene_tree(newick, removed) + "\n")
        (out_dir / "alignments" / aln_path.name).write_text(
            prune_alignment(aln_path.read_text(), removed)
        )

        kept = present - removed
        per_species = defaultdict(int)
        for leaf in kept:
            per_species[split_copy(leaf)[0]] += 1
        for species, count in per_species.items():
            copy_numbers[species][count] += 1
        n_genes += 1

    if n_genes == 0:
        print(f"WARNING: no g_trees*.trees processed in {in_dir} "
              f"(multi-batch layouts are not supported)")

    # One summary row per event group (mirrored nested events share a group).
    event_rows = []
    for group_id in sorted({e["group"] for e in events}):
        e = next(e for e in events if e["group"] == group_id)
        event_rows.append({
            "species": sorted({split_copy(l)[0] for side in e["sides"] for l in side}),
            "n_copies": len(e["sides"]),
            "autopolyploidy": e["autopolyploidy"],
            "q": q_by_group[group_id],
        })

    summary = {
        "dist": dist,
        "dist_params": dist_params,
        "seed": seed,
        "n_genes": n_genes,
        "events": event_rows,
        "copy_numbers": {
            sp: {str(c): n for c, n in sorted(counts.items())}
            for sp, counts in sorted(copy_numbers.items())
        },
    }
    with open(out_dir / "diploidization_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


# Same 21 networks as the rest of the pipeline (rescale_and_keep_ultrametric.py etc.)
NETWORKS = [
    "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
    "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019",
    "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011",
    "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014",
    "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014",
]


def stable_seed(base_seed, network, replicate):
    """A reproducible per-(network, replicate) seed derived from ``base_seed``.

    Uses a stable hash (not Python's salted ``hash()``) so reruns reproduce the
    same retention-rate draws and removals.
    """
    key = f"{base_seed}|{network}|{replicate}".encode()
    return int.from_bytes(hashlib.sha256(key).digest()[:8], "big")


def diploidize_config(base_dir, config, out_config, networks, replicates,
                      *, dist, dist_params, seed):
    """Diploidize every (network, replicate) of a config into ``out_config``.

    Layout (per the SimPhy pipeline)::

        {base_dir}/{network}/tree_10_mil.nex
        {base_dir}/{network}/data/{config}/replicate_{r}/1/   ->   .../{out_config}/...

    Returns ``{network: {replicate: summary}}``.
    """
    base_dir = Path(base_dir)
    results = defaultdict(dict)
    for network in networks:
        species_tree = base_dir / network / "tree_10_mil.nex"
        if not species_tree.exists():
            print(f"WARNING: missing species tree {species_tree}, skipping {network}")
            continue
        for r in replicates:
            in_dir = base_dir / network / "data" / config / f"replicate_{r}" / "1"
            if not in_dir.is_dir():
                print(f"WARNING: missing {in_dir}, skipping {network} replicate {r}")
                continue
            out_dir = base_dir / network / "data" / out_config / f"replicate_{r}" / "1"
            print(f"Diploidizing {network} replicate {r} -> {out_config}")
            results[network][r] = diploidize_replicate(
                species_tree, in_dir, out_dir,
                dist=dist, dist_params=dist_params,
                seed=stable_seed(seed, network, r),
            )
    return results


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Apply post-WGD diploidization (fractionation) to simulated "
                    "gene trees + alignments.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Beta-distributed retention, mean 0.5, into a new config
  python apply_diploidization.py conf_ils_low_10M --mean 0.5 --concentration 4

  # Fixed retention rate of 0.3
  python apply_diploidization.py conf_ils_low_10M --dist fixed --value 0.3

  # Sweep a single network / replicate
  python apply_diploidization.py conf_ils_low_10M --networks Bendiksby_2011 --replicates 1
""",
    )
    parser.add_argument("config", help="input config name (e.g. conf_ils_low_10M)")
    parser.add_argument("--base-dir", default=None,
                        help="dir with network subdirs (default: simulations/simulations "
                             "relative to this script)")
    parser.add_argument("--out-config", default=None,
                        help="output config name (default: {config}_dip{tag})")
    parser.add_argument("--networks", nargs="+", default=None,
                        help="subset of networks (default: all 21)")
    parser.add_argument("--replicates", nargs="+", type=int, default=[1, 2, 3, 4, 5],
                        help="replicate numbers (default: 1..5)")
    parser.add_argument("--seed", type=int, default=42, help="base random seed")

    parser.add_argument("--dist", choices=["beta", "uniform", "fixed"], default="beta",
                        help="retention-rate distribution (default: beta)")
    parser.add_argument("--mean", type=float, help="beta: mean retention rate")
    parser.add_argument("--concentration", type=float, default=4.0,
                        help="beta: concentration (default: 4.0)")
    parser.add_argument("--alpha", type=float, help="beta: alpha (overrides --mean)")
    parser.add_argument("--beta", type=float, help="beta: beta (overrides --mean)")
    parser.add_argument("--low", type=float, help="uniform: lower bound")
    parser.add_argument("--high", type=float, help="uniform: upper bound")
    parser.add_argument("--value", type=float, help="fixed: retention rate")
    return parser.parse_args(argv)


def _dist_params_from_args(args):
    if args.dist == "fixed":
        if args.value is None:
            raise SystemExit("--value is required for --dist fixed")
        return {"value": args.value}, f"fix{int(round(args.value * 100)):03d}"
    if args.dist == "uniform":
        if args.low is None or args.high is None:
            raise SystemExit("--low and --high are required for --dist uniform")
        return ({"low": args.low, "high": args.high},
                f"uni{int(round(args.low * 100)):03d}-{int(round(args.high * 100)):03d}")
    # beta
    if args.alpha is not None and args.beta is not None:
        mean = args.alpha / (args.alpha + args.beta)
        return ({"alpha": args.alpha, "beta": args.beta}, f"dip{int(round(mean * 100)):03d}")
    if args.mean is None:
        raise SystemExit("--mean (with --concentration) or --alpha/--beta required for --dist beta")
    return ({"mean": args.mean, "concentration": args.concentration},
            f"dip{int(round(args.mean * 100)):03d}")


def main(argv=None):
    args = _parse_args(argv)
    base_dir = (Path(args.base_dir) if args.base_dir
                else Path(__file__).resolve().parent.parent / "simulations")
    if not base_dir.is_dir():
        raise SystemExit(f"base dir not found: {base_dir}")

    dist_params, tag = _dist_params_from_args(args)
    out_config = args.out_config or f"{args.config}_{tag}"
    networks = args.networks or NETWORKS

    print(f"Base dir:    {base_dir}")
    print(f"Config:      {args.config} -> {out_config}")
    print(f"Networks:    {len(networks)}   Replicates: {args.replicates}")
    print(f"Retention:   dist={args.dist} params={dist_params} seed={args.seed}")

    diploidize_config(base_dir, args.config, out_config, networks, args.replicates,
                      dist=args.dist, dist_params=dist_params, seed=args.seed)


if __name__ == "__main__":
    main()
