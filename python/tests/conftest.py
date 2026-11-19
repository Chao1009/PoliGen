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


#: D7 (2026-10-27): set to 1 to run the suite without the shipped sets anyway
ALLOW_MISSING_SETS_ENV = "LIPOLGEN_TESTS_ALLOW_MISSING_SETS"


def pytest_sessionstart(session):
    """D7: the SHIPPED structure-function defaults are CT18NLO and
    NNPDFpol11_100 through the LHAPDF tier, so nearly every run-level test
    goes through them, and a missing set would surface as hundreds of
    identical failures.  The session stops at once instead, with the
    pre-flight's own message (it names the set, the searched paths and the
    toy escape) and return code 3.  A build WITHOUT the LHAPDF tier stops the
    same way: since D7 the tier is no longer optional for the shipped
    default, so this suite needs it too -- a consequence of the decision,
    stated in docs/PACKAGING.md, not a skip (a suite that skipped would pass
    without ever running the shipped default).  `LIPOLGEN_TESTS_ALLOW_
    MISSING_SETS=1` proceeds anyway: every default-path test then fails
    with that message and the toy-only tests pass.  CI never sets it."""
    if os.environ.get(ALLOW_MISSING_SETS_ENV) == "1":
        return
    try:
        import lipolgen as _lg
    except Exception:            # the import failure is reported by the tests
        return
    _l = _lg._lipolgen
    if not _l.HAVE_LHAPDF:
        pytest.exit("the shipped structure-function defaults (D7, 2026-10-27: "
                    "--unpol-sf ct18nlo, --pol-sf nnpdfpol) need the LHAPDF "
                    "tier, and this build has none (lipolgen.HAVE_LHAPDF is "
                    "False).  Rebuild with -DLIPOLGEN_WITH_LHAPDF=ON, or set "
                    "%s=1 to run the toy-only tests (the default-path ones then "
                    "fail by name)" % ALLOW_MISSING_SETS_ENV, returncode=3)
    problems = [p for p in (_l.lhapdf_set_problem(s)
                            for s in ("CT18NLO", "NNPDFpol11_100")) if p]
    if problems:
        pytest.exit("\n".join(problems) + "\n(set %s=1 to run the suite "
                    "anyway)" % ALLOW_MISSING_SETS_ENV, returncode=3)
