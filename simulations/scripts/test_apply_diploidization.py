"""Tests for apply_diploidization.py (post-WGD diploidization / fractionation).

Run:  pytest simulations/scripts/test_apply_diploidization.py -v
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from collections import Counter, defaultdict as _dd

import numpy as np
from ete3 import Tree

import apply_diploidization as ad

NETWORKS_DIR = Path(__file__).resolve().parent.parent / "networks"


def _uniquify(newick):
    """Add _1.._N suffixes to repeated leaves, as the pipeline does for SimPhy."""
    t = Tree(newick, format=1)
    counts = Counter(l.name for l in t.iter_leaves())
    seen = _dd(int)
    for l in t.iter_leaves():
        if counts[l.name] > 1:
            seen[l.name] += 1
            l.name = f"{l.name}_{seen[l.name]}"
    return t


# ───────────────────────── name parsing ─────────────────────────

class TestSpeciesTreeLeaf:
    """gene-tree tip -> species-tree leaf name (strip SimPhy _locus_individual)."""

    def test_polyploid_copy_tip(self):
        assert ad.species_tree_leaf("Lamiumconfertum1_1_0_0") == "Lamiumconfertum1_1"

    def test_diploid_tip_with_trailing_digit_in_name(self):
        # 'Lamiumalbum2' is a diploid species whose NAME ends in 2 (no copy suffix)
        assert ad.species_tree_leaf("Lamiumalbum2_0_0") == "Lamiumalbum2"

    def test_plain_diploid_tip(self):
        assert ad.species_tree_leaf("Galeopsisspeciosa_0_0") == "Galeopsisspeciosa"


class TestSplitCopy:
    """species-tree leaf -> (base species name, copy index or None)."""

    def test_polyploid_copy(self):
        assert ad.split_copy("Lamiumhybridum3_1") == ("Lamiumhybridum3", 1)

    def test_polyploid_copy_two(self):
        assert ad.split_copy("Lamiumhybridum3_2") == ("Lamiumhybridum3", 2)

    def test_diploid_no_suffix(self):
        assert ad.split_copy("Lamiumalbum2") == ("Lamiumalbum2", None)

    def test_plain_diploid(self):
        assert ad.split_copy("Galeopsisladanum") == ("Galeopsisladanum", None)


# ───────────────────────── event detection ─────────────────────────

def _sides_as_set(event):
    """Sides of an event as an order-independent set of frozensets."""
    return {frozenset(s) for s in event["sides"]}


class TestBuildEvents:
    """build_events folds a unique-named species tree into WGD events."""

    def test_diploid_only_no_events(self):
        events = ad.build_events("((A,X),(B,Y));")
        assert events == []

    def test_autopolyploid_duplicated_clade(self):
        # clade (A,B) duplicated as identical siblings -> one autopolyploidy event
        events = ad.build_events("(((A_1,B_1),(A_2,B_2)),C);")
        assert len(events) == 1
        assert _sides_as_set(events[0]) == {
            frozenset({"A_1", "B_1"}),
            frozenset({"A_2", "B_2"}),
        }
        assert events[0]["autopolyploidy"] is True

    def test_allopolyploid_copies_under_different_parents(self):
        # species A duplicated, copies grafted into two different clades -> allo event
        events = ad.build_events("((A_1,X),(A_2,Y));")
        assert len(events) == 1
        assert _sides_as_set(events[0]) == {frozenset({"A_1"}), frozenset({"A_2"})}
        assert events[0]["autopolyploidy"] is False

    def test_flat_triple_copy(self):
        # three identical siblings -> one event with three sides
        events = ad.build_events("((X_1,X_2,X_3),Y);")
        assert len(events) == 1
        assert _sides_as_set(events[0]) == {
            frozenset({"X_1"}),
            frozenset({"X_2"}),
            frozenset({"X_3"}),
        }

    def test_nested_autopolyploid_outer_then_inner(self):
        # two rounds of WGD on X: ((X,X),(X,X)) -> outer event + inner event
        events = ad.build_events("(((X_1,X_2),(X_3,X_4)),Y);")
        assert len(events) == 2
        # ordered outer (taller) first
        outer, inner = events
        assert _sides_as_set(outer) == {
            frozenset({"X_1", "X_2"}),
            frozenset({"X_3", "X_4"}),
        }
        assert _sides_as_set(inner) == {frozenset({"X_1"}), frozenset({"X_2"})}


class TestBuildEventsRealNetwork:
    """Validate build_events against a real MUL-tree with known ground truth."""

    def test_bendiksby_six_binary_allo_events(self):
        # mul_tree_final_stats.csv: 6 polyploids, Max_Copies=2, H_Strict=6, auto=0
        nwk = (NETWORKS_DIR / "Bendiksby_2011.tre").read_text().strip()
        events = ad.build_events(_uniquify(nwk))

        assert len(events) == 6
        species = set()
        for e in events:
            assert len(e["sides"]) == 2          # all 2-copy
            assert e["autopolyploidy"] is False  # all allopolyploid
            bases = {ad.split_copy(n)[0] for side in e["sides"] for n in side}
            assert len(bases) == 1               # each event duplicates one species
            species.add(bases.pop())

        assert species == {
            "Lamiumamplexicaulevarorientale2", "Lamiumconfertum1",
            "Lamiumgaleobdolonsubspargentatum", "Lamiumgaleobdolonsubspmontanum",
            "Lamiumhybridum1", "Lamiumhybridum3",
        }


# ───────────────────────── retention-rate draws ─────────────────────────

class TestDrawRetentionRates:
    def test_fixed_returns_constant(self):
        q = ad.draw_retention_rates(5, dist="fixed", value=0.4,
                                    rng=np.random.default_rng(0))
        assert q == [0.4] * 5

    def test_beta_deterministic_for_same_seed(self):
        a = ad.draw_retention_rates(10, dist="beta", alpha=2, beta=2,
                                    rng=np.random.default_rng(7))
        b = ad.draw_retention_rates(10, dist="beta", alpha=2, beta=2,
                                    rng=np.random.default_rng(7))
        assert a == b
        assert len(a) == 10
        assert all(0.0 <= x <= 1.0 for x in a)

    def test_beta_mean_parametrization(self):
        # mean/concentration -> alpha=mean*conc, beta=(1-mean)*conc
        q = ad.draw_retention_rates(20000, dist="beta", mean=0.5, concentration=4,
                                    rng=np.random.default_rng(1))
        assert abs(np.mean(q) - 0.5) < 0.02

    def test_beta_degenerate_mean_one_is_point_mass(self):
        q = ad.draw_retention_rates(5, dist="beta", mean=1.0, concentration=4,
                                    rng=np.random.default_rng(0))
        assert q == [1.0] * 5

    def test_beta_degenerate_mean_zero_is_point_mass(self):
        q = ad.draw_retention_rates(5, dist="beta", mean=0.0, concentration=4,
                                    rng=np.random.default_rng(0))
        assert q == [0.0] * 5

    def test_uniform_within_bounds(self):
        q = ad.draw_retention_rates(1000, dist="uniform", low=0.2, high=0.8,
                                    rng=np.random.default_rng(2))
        assert all(0.2 <= x <= 0.8 for x in q)


# ───────────────────────── per-gene removal planner ─────────────────────────

AUTO_TREE = "(((A_1,B_1),(A_2,B_2)),C);"      # one auto event, clade (A,B) duplicated
ALLO_TREE = "((A_1,X),(A_2,Y));"               # one allo event on species A
NESTED_TREE = "(((X_1,X_2),(X_3,X_4)),Y);"     # tetraploid X, two nested events


def _present(tree_newick):
    return {l.name for l in Tree(tree_newick, format=1).iter_leaves()}


class TestPlanRemovals:
    def test_q1_removes_nothing(self):
        events = ad.build_events(AUTO_TREE)
        present = _present(AUTO_TREE)
        for seed in range(20):
            removed = ad.plan_removals(events, [1.0] * len(events), present,
                                       np.random.default_rng(seed))
            assert removed == set()

    def test_q0_collapses_binary_event_to_single_copy(self):
        events = ad.build_events(ALLO_TREE)
        present = _present(ALLO_TREE)
        removed = ad.plan_removals(events, [0.0] * len(events), present,
                                   np.random.default_rng(0))
        # exactly one of A_1 / A_2 removed; the other survives
        assert removed in ({"A_1"}, {"A_2"})

    def test_q0_keeps_at_least_one_copy_per_species_nested(self):
        events = ad.build_events(NESTED_TREE)
        present = _present(NESTED_TREE)
        for seed in range(20):
            removed = ad.plan_removals(events, [0.0] * len(events), present,
                                       np.random.default_rng(seed))
            survivors = present - removed
            x_survivors = {l for l in survivors if ad.split_copy(l)[0] == "X"}
            assert len(x_survivors) >= 1

    def test_unbiased_drop_is_roughly_fifty_fifty(self):
        events = ad.build_events(ALLO_TREE)
        present = _present(ALLO_TREE)
        rng = np.random.default_rng(123)
        dropped_a1 = 0
        n = 4000
        for _ in range(n):
            removed = ad.plan_removals(events, [0.0] * len(events), present, rng)
            if removed == {"A_1"}:
                dropped_a1 += 1
        assert abs(dropped_a1 / n - 0.5) < 0.05


# ───────────────────────── pruning gene trees + alignments ─────────────────────────

GENE_TREE = "((A_1_0_0:1,A_2_0_0:1):1,(B_0_0:1,C_0_0:1):1);"
PHYLIP = "3 4\nA_1_0_0   ACGT\nA_2_0_0   ACGA\nB_0_0   ACGC\n"


class TestPruneGeneTree:
    def test_removes_matching_tips_keeps_others(self):
        out = ad.prune_gene_tree(GENE_TREE, {"A_1"})
        tips = {l.name for l in Tree(out, format=1).iter_leaves()}
        assert tips == {"A_2_0_0", "B_0_0", "C_0_0"}

    def test_remove_nothing_keeps_all_tips(self):
        out = ad.prune_gene_tree(GENE_TREE, set())
        tips = {l.name for l in Tree(out, format=1).iter_leaves()}
        assert tips == {"A_1_0_0", "A_2_0_0", "B_0_0", "C_0_0"}

    def test_output_keeps_branch_lengths(self):
        out = ad.prune_gene_tree(GENE_TREE, {"A_1"})
        assert all(l.dist > 0 for l in Tree(out, format=1).iter_leaves())


class TestPhylip:
    def test_read_phylip_preserves_order_and_seqs(self):
        d = ad.read_phylip(PHYLIP)
        assert list(d.keys()) == ["A_1_0_0", "A_2_0_0", "B_0_0"]
        assert d["A_2_0_0"] == "ACGA"

    def test_prune_alignment_drops_seqs_and_updates_header(self):
        out = ad.prune_alignment(PHYLIP, {"A_1"})
        d = ad.read_phylip(out)
        assert list(d.keys()) == ["A_2_0_0", "B_0_0"]
        assert out.splitlines()[0].split() == ["2", "4"]

    def test_prune_alignment_roundtrips_sequences(self):
        out = ad.prune_alignment(PHYLIP, {"A_1"})
        d = ad.read_phylip(out)
        assert d["B_0_0"] == "ACGC"


# ───────────────────────── per-replicate orchestration ─────────────────────────

SPECIES_NEXUS = (
    "#NEXUS\nbegin trees;\n"
    "tree 1=((A_1:1,A_2:1):1,B:2);\n"  # A is an autopolyploid (2 copies)
    "end;\n"
)


def _make_replicate(tmp_path):
    """Build a minimal species tree + one gene (tree + alignment) on disk."""
    sp = tmp_path / "tree_10_mil.nex"
    sp.write_text(SPECIES_NEXUS)
    in_dir = tmp_path / "in"
    (in_dir / "alignments").mkdir(parents=True)
    (in_dir / "g_trees0001.trees").write_text(
        "((A_1_0_0:1,A_2_0_0:1):1,B_0_0:2);\n"
    )
    (in_dir / "g_trees0001.trees.log").write_text("ignore me\n")
    (in_dir / "alignments" / "alignment_0001.phy").write_text(
        "3 4\nA_1_0_0   ACGT\nA_2_0_0   ACGA\nB_0_0   ACGC\n"
    )
    return sp, in_dir


class TestDiploidizeReplicate:
    def test_fixed_zero_collapses_polyploid(self, tmp_path):
        sp, in_dir = _make_replicate(tmp_path)
        out_dir = tmp_path / "out"
        ad.diploidize_replicate(sp, in_dir, out_dir,
                                dist="fixed", dist_params={"value": 0.0}, seed=0)

        tips = {l.name for l in Tree((out_dir / "g_trees0001.trees").read_text(),
                                     format=1).iter_leaves()}
        a_copies = {t for t in tips if ad.species_tree_leaf(t).startswith("A_")}
        assert len(a_copies) == 1                       # collapsed to single copy
        assert any(ad.species_tree_leaf(t) == "B" for t in tips)

        aln = ad.read_phylip((out_dir / "alignments" / "alignment_0001.phy").read_text())
        assert len(aln) == 2                            # one A copy dropped from alignment
        assert set(aln.keys()) == set(tips)             # tree and alignment stay in sync

    def test_fixed_one_keeps_everything(self, tmp_path):
        sp, in_dir = _make_replicate(tmp_path)
        out_dir = tmp_path / "out"
        ad.diploidize_replicate(sp, in_dir, out_dir,
                                dist="fixed", dist_params={"value": 1.0}, seed=0)
        tips = {l.name for l in Tree((out_dir / "g_trees0001.trees").read_text(),
                                     format=1).iter_leaves()}
        assert tips == {"A_1_0_0", "A_2_0_0", "B_0_0"}

    def test_log_files_are_not_processed(self, tmp_path):
        sp, in_dir = _make_replicate(tmp_path)
        out_dir = tmp_path / "out"
        ad.diploidize_replicate(sp, in_dir, out_dir,
                                dist="fixed", dist_params={"value": 1.0}, seed=0)
        assert not (out_dir / "g_trees0001.trees.log").exists()

    def test_summary_records_event_rates_and_copy_numbers(self, tmp_path):
        sp, in_dir = _make_replicate(tmp_path)
        out_dir = tmp_path / "out"
        summary = ad.diploidize_replicate(sp, in_dir, out_dir,
                                          dist="fixed", dist_params={"value": 0.0}, seed=0)
        assert summary["n_genes"] == 1
        assert len(summary["events"]) == 1
        assert summary["events"][0]["q"] == 0.0
        # after full fractionation A has exactly one copy in the single gene
        assert summary["copy_numbers"]["A"] == {"1": 1}


class TestStableSeed:
    def test_deterministic(self):
        assert ad.stable_seed(42, "Bendiksby_2011", 1) == ad.stable_seed(42, "Bendiksby_2011", 1)

    def test_varies_by_network_replicate_and_base(self):
        s = ad.stable_seed(42, "Bendiksby_2011", 1)
        assert s != ad.stable_seed(42, "Bendiksby_2011", 2)
        assert s != ad.stable_seed(42, "Ding_2023", 1)
        assert s != ad.stable_seed(7, "Bendiksby_2011", 1)


class TestDiploidizeConfig:
    def _make_layout(self, tmp_path, network, config, replicate):
        net = tmp_path / network
        (net).mkdir()
        (net / "tree_10_mil.nex").write_text(SPECIES_NEXUS)
        rep = net / "data" / config / f"replicate_{replicate}" / "1"
        (rep / "alignments").mkdir(parents=True)
        (rep / "g_trees0001.trees").write_text("((A_1_0_0:1,A_2_0_0:1):1,B_0_0:2);\n")
        (rep / "alignments" / "alignment_0001.phy").write_text(
            "3 4\nA_1_0_0   ACGT\nA_2_0_0   ACGA\nB_0_0   ACGC\n"
        )

    def test_writes_diploidized_config_layout(self, tmp_path):
        self._make_layout(tmp_path, "Bendiksby_2011", "conf_ils_low_10M", 1)
        ad.diploidize_config(
            base_dir=tmp_path, config="conf_ils_low_10M",
            out_config="conf_ils_low_10M_dip050",
            networks=["Bendiksby_2011"], replicates=[1],
            dist="fixed", dist_params={"value": 0.0}, seed=42,
        )
        out = (tmp_path / "Bendiksby_2011" / "data" / "conf_ils_low_10M_dip050"
               / "replicate_1" / "1")
        assert (out / "g_trees0001.trees").exists()
        assert (out / "alignments" / "alignment_0001.phy").exists()
        assert (out / "diploidization_summary.json").exists()


class TestRealizedRetention:
    """Across many genes the realized fraction keeping both copies tracks q."""

    def test_fraction_keeping_both_copies_approximates_q(self):
        events = ad.build_events(ALLO_TREE)   # one binary event on species A
        present = _present(ALLO_TREE)
        rng = np.random.default_rng(99)
        q = 0.6
        n, both = 8000, 0
        for _ in range(n):
            removed = ad.plan_removals(events, [q] * len(events), present, rng)
            if not removed:                   # both A copies retained for this gene
                both += 1
        assert abs(both / n - q) < 0.03

