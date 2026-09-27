"""Make the in-tree build importable without installing anything.

`source env.sh` already puts `build/python` on PYTHONPATH; this is the
fallback for a bare `pytest python/tests` in a shell that did not.

Set `LIPOLGEN_TESTS_USE_INSTALLED=1` to opt out of this and exercise whatever
`import lipolgen` already resolves to instead -- e.g. a `pip install`ed
wheel/editable install on the current `sys.path` -- which is how this same
suite doubles as the packaging gate (docs/OPEN_ITEMS_SOLUTIONS.md item 12).
Default behaviour (nothing set) is unchanged.

It also owns the benchmark opt-in (validation/benchmarks/README.md, "The
opt-in variable"; BENCHMARK_PLAN.md sec. 3): a test marked
`@pytest.mark.bench_expensive` -- an expensive benchmark row: one that costs
more than 30 s, or one the plan counts as expensive by kind (a full
generation, a chain leg, a generator-vs-generator sample, whatever it
costs) -- is SKIPPED, with a reason that names the variable, unless
`LIPOLGEN_BENCH=1`:

    LIPOLGEN_BENCH=1 python -m pytest python/tests -q

The marker is registered here (`pytest_configure`), not in pyproject.toml:
adding a `[tool.pytest.ini_options]` table there would move pytest's rootdir
and config-file discovery for every invocation.
"""
import os
import sys

import pytest

#: the environment variable that opts the expensive benchmark rows in, and
#: the marker that makes a test one of them
BENCH_ENV = "LIPOLGEN_BENCH"
BENCH_MARKER = "bench_expensive"
BENCH_SKIP_REASON = ("bench_expensive: an expensive benchmark row (> 30 s, or "
                     "a generator-vs-generator sample); set %s=1 to run it"
                     % BENCH_ENV)

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


def pytest_configure(config):
    """Register the opt-in marker (no PytestUnknownMarkWarning, and
    `--strict-markers` accepts it)."""
    config.addinivalue_line(
        "markers",
        "%s: an expensive benchmark row (> 30 s, or a generator-vs-generator "
        "sample); runs only with %s=1 (validation/benchmarks/README.md, "
        "'The opt-in variable')"
        % (BENCH_MARKER, BENCH_ENV))


def pytest_collection_modifyitems(config, items):
    """Skip every `bench_expensive` test unless LIPOLGEN_BENCH is exactly
    "1"; the skip reason names the variable."""
    if os.environ.get(BENCH_ENV) == "1":
        return
    skip = pytest.mark.skip(reason=BENCH_SKIP_REASON)
    for item in items:
        if item.get_closest_marker(BENCH_MARKER) is not None:
            item.add_marker(skip)


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
