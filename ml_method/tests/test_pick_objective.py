from gene2net_gnn.training.trainer_reconstruct import pick_objective


def test_pick_objective_respects_guard():
    hist = [
        {"f1": 0.80, "allo_acc": 0.40},
        {"f1": 0.78, "allo_acc": 0.60},   # excluded: f1 below guard
        {"f1": 0.81, "allo_acc": 0.50},
    ]
    assert pick_objective(hist, 0.79) == 0.50


def test_pick_objective_none_qualify_returns_negative():
    hist = [{"f1": 0.70, "allo_acc": 0.99}]
    assert pick_objective(hist, 0.79) == -1.0
