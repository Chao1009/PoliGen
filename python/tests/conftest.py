"""Make the in-tree build importable without installing anything.

`source env.sh` already puts `build/python` on PYTHONPATH; this is the
fallback for a bare `pytest python/tests` in a shell that did not.

Set `LIPOLGEN_TESTS_USE_INSTALLED=1` to opt out of this and exercise whatever
`import lipolgen` already resolves to instead -- e.g. a `pip install`ed
wheel/editable install on the current `sys.path` -- which is how this same
suite doubles as the packaging gate (docs/OPEN_ITEMS_SOLUTIONS.md item 12).
Default behaviour (nothing set) is unchanged.
"""
import os
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))          # the LiPolGen repo
_BUILD = os.path.join(_ROOT, "build", "python")

if not os.environ.get("LIPOLGEN_TESTS_USE_INSTALLED") and os.path.isdir(_BUILD):
    # FIRST, unconditionally.  `env.sh` already exports it, but a bare
    # `cd python && pytest tests` puts the CWD (the SOURCE package, which has
    # no compiled `_lipolgen`) ahead of it, and the source copy then shadows
    # the staged one.  Moving it to the front is what makes both invocations
    # import the same package.
    while _BUILD in sys.path:
        sys.path.remove(_BUILD)
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
