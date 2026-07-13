from scripts.metadata_event_stats import event_stats


def test_counts_auto_allo_single_and_clade():
    events = [
        {"event_type": "auto", "target_clade": ["sp1", "sp2"]},
        {"event_type": "allo", "target_clade": ["sp3"], "partner_clade": ["sp9"]},
        {"event_type": "allo", "target_clade": ["sp4", "sp5"], "partner_clade": ["sp9"]},
    ]
    s = event_stats(events)
    assert s["n_auto"] == 1
    assert s["n_allo"] == 2
    assert s["n_allo_single"] == 1
    assert s["n_allo_clade"] == 1
    assert s["max_target_copies"] == 2
