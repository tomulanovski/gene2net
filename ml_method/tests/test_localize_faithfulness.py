from scripts.localize_faithfulness import classify_sample


def test_auto_only():
    events = [{"event_type": "auto", "target_clade": ["s1", "s2"]}]
    assert classify_sample(events) == {"has_auto", "auto_only"}


def test_single_and_clade_allo():
    events = [
        {"event_type": "auto", "target_clade": ["s1"]},
        {"event_type": "allo", "target_clade": ["s3"], "partner_clade": ["s9"]},
        {"event_type": "allo", "target_clade": ["s4", "s5"], "partner_clade": ["s9"]},
    ]
    tags = classify_sample(events)
    assert "has_auto" in tags
    assert "has_single_allo" in tags
    assert "has_clade_allo" in tags
    assert "auto_only" not in tags


def test_no_events():
    assert classify_sample([]) == {"no_events"}
