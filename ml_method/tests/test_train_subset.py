import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from train_reconstruct import select_train_subset


def test_fraction_takes_prefix():
    xs = list(range(100))
    assert select_train_subset(xs, fraction=0.25) == list(range(25))


def test_max_samples_takes_prefix():
    xs = list(range(100))
    assert select_train_subset(xs, max_samples=10) == list(range(10))


def test_none_returns_all():
    xs = list(range(5))
    assert select_train_subset(xs) == xs


def test_fraction_wins_over_max():
    xs = list(range(100))
    assert select_train_subset(xs, fraction=0.5, max_samples=3) == list(range(50))
