import csv
from scripts.summarize_oracle import summarize


def _write(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def test_summarize_means(tmp_path):
    scores = tmp_path / "reconstruction_scores.csv"
    mapping = tmp_path / "mapping_scores.csv"
    _write(scores, [
        {"edit_distance_multree": "0.30", "ret_leaf_jaccard": "0.20", "ret_sisters_jaccard": "0.35"},
        {"edit_distance_multree": "0.40", "ret_leaf_jaccard": "0.10", "ret_sisters_jaccard": "0.25"},
    ], ["edit_distance_multree", "ret_leaf_jaccard", "ret_sisters_jaccard"])
    _write(mapping, [
        {"sample": "0001", "type": "allo", "source": "true", "target": "1.0", "A": "0.8", "B": "1.0"},
        {"sample": "0002", "type": "auto", "source": "true", "target": "1.0", "A": "1.0", "B": "1.0"},
    ], ["sample", "type", "source", "target", "A", "B"])
    out = summarize(str(scores), str(mapping))
    assert abs(out["edit_distance_multree"] - 0.35) < 1e-9
    assert abs(out["allo_A_snap"] - 0.8) < 1e-9  # only the allo row counts
