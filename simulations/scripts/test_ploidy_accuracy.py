"""Tests for ploidy_accuracy.py."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import ploidy_accuracy as pa

STATS = (
    "Filename,Num_Species,Num_Polyploids,Polyploid_Names\n"
    "Net1.tre,10,2,\"Aspecies, Bspecies\"\n"
    "Net2.tre,5,1,Csingle\n"
    "Net3.tre,4,0,\n"
)

DIST = (
    "Species\tRepresentativeCopyNumber\tDistribution\n"
    "Aspecies\t2\t2:480, 1:424\n"      # polyploid, correctly called
    "Bspecies\t1\t1:900, 2:90\n"       # polyploid, misclassified as diploid
    "Zdiploid\t1\t1:990\n"            # true diploid, ignored
)


class TestParseTruePolyploids:
    def test_parses_quoted_and_unquoted_names(self):
        tp = pa.parse_true_polyploids(STATS)
        assert tp["Net1"] == {"Aspecies", "Bspecies"}
        assert tp["Net2"] == {"Csingle"}

    def test_zero_polyploids_is_empty_set(self):
        tp = pa.parse_true_polyploids(STATS)
        assert tp["Net3"] == set()

    def test_strips_tre_suffix_from_network(self):
        tp = pa.parse_true_polyploids(STATS)
        assert "Net1" in tp and "Net1.tre" not in tp


class TestParseDistribution:
    def test_species_to_copy_number(self):
        d = pa.parse_distribution(DIST)
        assert d == {"Aspecies": 2, "Bspecies": 1, "Zdiploid": 1}


class TestRecordsForDistribution:
    def test_only_true_polyploids_scored(self):
        dist = pa.parse_distribution(DIST)
        recs = pa.records_for_distribution("Net1", "1", {"Aspecies", "Bspecies"}, dist)
        by_sp = {r["species"]: r for r in recs}
        assert set(by_sp) == {"Aspecies", "Bspecies"}      # Zdiploid not scored
        assert by_sp["Aspecies"]["present"] and by_sp["Aspecies"]["correct"]
        assert by_sp["Bspecies"]["present"] and not by_sp["Bspecies"]["correct"]

    def test_true_polyploid_absent_from_distribution(self):
        dist = {"Aspecies": 2}  # Bspecies missing entirely
        recs = pa.records_for_distribution("Net1", "1", {"Aspecies", "Bspecies"}, dist)
        by_sp = {r["species"]: r for r in recs}
        assert by_sp["Bspecies"]["present"] is False


class TestAggregate:
    def test_recall_overall_and_per_network(self):
        records = [
            {"network": "N", "replicate": "1", "species": "A", "present": True, "correct": True},
            {"network": "N", "replicate": "1", "species": "B", "present": True, "correct": False},
            {"network": "N", "replicate": "2", "species": "A", "present": True, "correct": True},
            {"network": "N", "replicate": "2", "species": "B", "present": True, "correct": True},
        ]
        agg = pa.aggregate(records)
        assert agg["overall"]["n_true"] == 4
        assert agg["overall"]["n_correct"] == 3
        assert agg["overall"]["recall"] == 0.75

    def test_absent_excluded_from_denominator_but_counted(self):
        records = [
            {"network": "N", "replicate": "1", "species": "A", "present": True, "correct": True},
            {"network": "N", "replicate": "1", "species": "B", "present": False, "correct": False},
        ]
        agg = pa.aggregate(records)
        assert agg["overall"]["n_true"] == 1      # only present ones in denominator
        assert agg["overall"]["n_correct"] == 1
        assert agg["overall"]["n_absent"] == 1
        assert agg["overall"]["recall"] == 1.0


class TestQFromConfig:
    def test_fix_suffix_parsed(self):
        assert pa.q_from_config("conf_dup_loss_medium_10M_ne1M_fix025") == 0.25
        assert pa.q_from_config("conf_dup_loss_medium_10M_ne1M_fix000") == 0.0

    def test_no_suffix_is_baseline_one(self):
        assert pa.q_from_config("conf_dup_loss_medium_10M_ne1M") == 1.0
