"""load_model must fail loudly when no checkpoint exists.

Silently returning a randomly initialized model (the previous behavior) makes
evaluation report meaningless near-chance numbers that look like a real result.
"""
import pytest

from scripts.reconstruct_mul_tree import load_model


def test_load_model_raises_when_no_checkpoint(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_model(str(tmp_path), {}, "cpu")
