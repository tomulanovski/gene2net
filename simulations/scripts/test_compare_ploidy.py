"""Tests for compare_ploidy.py."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import compare_ploidy as cp

DIST = (
    "Species\tRepresentativeCopyNumber\tDistribution\n"
    "Xpoly\t2\t2:480, 1:424, 3:56\n"
    "Ypoly\t2\t2:463, 1:434\n"
    "Zdip\t1\t1:900, 2:90\n"
)


class TestParseDistribution:
    def test_parses_species_to_representative_copy_number(self):
        d = cp.parse_distribution(DIST)
        assert d == {"Xpoly": 2, "Ypoly": 2, "Zdip": 1}


class TestAggregate:
    def test_polyploid_retention_metrics(self):
        records = [
            {"network": "N", "replicate": 1, "species": "X", "a": 2, "b": 1},  # lost
            {"network": "N", "replicate": 1, "species": "Y", "a": 2, "b": 2},  # kept
            {"network": "N", "replicate": 1, "species": "Z", "a": 1, "b": 1},  # diploid
        ]
        agg = cp.aggregate(records)
        net = agg["per_network"]["N"]
        assert net["n_polyploid_orig"] == 2      # X and Y called polyploid originally
        assert net["n_still_polyploid"] == 1     # only Y stays polyploid
        assert net["n_lost"] == 1                # X drops to diploid
        assert net["retained_frac"] == 0.5

    def test_overall_rolls_up(self):
        records = [
            {"network": "N1", "replicate": 1, "species": "X", "a": 2, "b": 2},
            {"network": "N2", "replicate": 1, "species": "Y", "a": 2, "b": 1},
        ]
        agg = cp.aggregate(records)
        assert agg["overall"]["n_polyploid_orig"] == 2
        assert agg["overall"]["n_still_polyploid"] == 1
        assert agg["overall"]["retained_frac"] == 0.5
