"""Make the in-tree build importable without installing anything.

`source env.sh` already puts `build/python` on PYTHONPATH; this is the
fallback for a bare `pytest python/tests` in a shell that did not.
"""
import os
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))          # the LiPolGen repo
_BUILD = os.path.join(_ROOT, "build", "python")

if _BUILD not in sys.path and os.path.isdir(_BUILD):
    sys.path.insert(0, _BUILD)


@pytest.fixture(scope="session")
def repo_root():
    return _ROOT


@pytest.fixture(scope="session")
def reference_dir(repo_root):
    d = os.path.join(repo_root, "validation", "reference")
    if not os.path.isdir(d):
        pytest.skip("validation/reference is not populated")
    return d


@pytest.fixture(scope="session")
def polligen(repo_root):
    """`polligen` from the sibling PolarizedLithiumSim checkout, or skip."""
    evgen = os.path.join(os.path.dirname(repo_root), "PolarizedLithiumSim",
                         "evgen")
    fastsim = os.path.join(os.path.dirname(repo_root), "PolarizedLithiumSim",
                           "fastsim")
    for p in (evgen, fastsim):
        if os.path.isdir(p) and p not in sys.path:
            sys.path.insert(0, p)
    try:
        import polligen
    except Exception as exc:                       # pragma: no cover
        pytest.skip("polligen is not importable: %s" % exc)
    return polligen
