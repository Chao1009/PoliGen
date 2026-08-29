"""Throughput measured from Python, printed with -s.

Not a pass/fail performance gate (the machine varies); the assertions only
catch a collapse -- e.g. a binding that forgot to release the GIL, or a
per-event Python callback sneaking into the loop.  Run

    pytest python/tests/test_throughput.py -s

to see the numbers.
"""

import time

import pytest

import lipolgen as lg

N = 200000


def _rate(pipeline, n, nthreads=1):
    t0 = time.perf_counter()
    cols = pipeline.generate(n, False, nthreads)
    dt = time.perf_counter() - t0
    return cols["x"].size / dt, dt


@pytest.mark.parametrize("channel,floor", [
    ("inclusive", 5e4),
    ("tagged-alpha", 3e4),
    ("coherent", 3e4),
])
def test_columnar_throughput(channel, floor, capsys):
    cfg = lg.make_config(isotope="6Li", channel=channel, config=1, events=N,
                         seed=1)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    _rate(p, 20000)                                # warm the caches
    rate, dt = _rate(p, N)
    with capsys.disabled():
        print("\n  %-14s %8.0f ev/s single core (%d events, %.3f s)"
              % (channel, rate, N, dt))
    assert rate > floor, (channel, rate)


def test_threading_helps_or_at_least_does_not_hurt(capsys):
    cfg = lg.make_config(isotope="6Li", channel="tagged-alpha", config=1,
                         events=N, seed=1)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    _rate(p, 20000)
    one, _ = _rate(p, N, 1)
    many, _ = _rate(p, N, 8)
    with capsys.disabled():
        print("  %-14s %8.0f ev/s x1   %8.0f ev/s x8   (speedup %.2f)"
              % ("tagged-alpha", one, many, many / one))
    # The GIL is released around the loop, so 8 threads must not be SLOWER
    # than one by more than the block-scheduling noise.
    assert many > 0.7 * one


def test_event_object_and_hepmc_costs(capsys):
    """The two costs a user pays on purpose: materializing Event records and
    writing HepMC3 (docs/USAGE.md 8 measures ~18 k ev/s for the latter)."""
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=20000, seed=1)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    t0 = time.perf_counter()
    p.generate(20000, True, 1)
    dt_events = time.perf_counter() - t0
    with capsys.disabled():
        print("  %-14s %8.0f ev/s with Event records kept"
              % ("inclusive", 20000 / dt_events))
    assert 20000 / dt_events > 1e4
