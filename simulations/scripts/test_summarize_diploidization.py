"""Tests for summarize_diploidization.py."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import summarize_diploidization as sd


# One replicate summary: polyploid A has 40 single-copy and 60 two-copy genes;
# C is a diploid (not in events) and must be ignored in the polyploid stats.
S1 = {
    "events": [{"species": ["A"], "q": 0.6}],
    "copy_numbers": {"A": {"1": 40, "2": 60}, "C": {"1": 100}},
    # of 100 eligible genes, 40 had a copy dropped -> realized retention 0.6
    "removal": {"A": {"eligible": 100, "reduced": 40, "copies_removed": 40,
                      "realized_retention": 0.6}},
}


class TestAggregate:
    def test_per_network_metrics(self):
        records = [
            {"network": "net1", "summary": S1},
            {"network": "net1", "summary": S1},  # 2 identical replicates
        ]
        agg = sd.aggregate(records)
        net = agg["per_network"]["net1"]
        assert net["n_replicates"] == 2
        assert net["n_polyploids"] == 1
        assert net["n_obs"] == 200            # 2 reps * (40 + 60)
        assert net["single_copy_frac"] == 0.4  # 80 / 200
        assert abs(net["mean_copies"] - 1.6) < 1e-9  # (80*1 + 120*2)/200
        assert net["q"] == 0.6
        # diploidization-only: 2 reps * 40 reduced / (2 reps * 100 eligible)
        assert net["realized_retention"] == 0.6

    def test_diploids_excluded(self):
        agg = sd.aggregate([{"network": "n", "summary": S1}])
        # only A counts; C (diploid) is not in the polyploid observations
        assert agg["per_network"]["n"]["n_obs"] == 100

    def test_overall_rolls_up_networks(self):
        records = [
            {"network": "net1", "summary": S1},
            {"network": "net2", "summary": S1},
        ]
        agg = sd.aggregate(records)
        assert agg["overall"]["n_networks"] == 2
        assert agg["overall"]["n_obs"] == 200
        assert abs(agg["overall"]["mean_copies"] - 1.6) < 1e-9
