"""InclusiveSampler through the binding: schema, determinism, thread safety."""

import numpy as np
import pytest

import lipolgen as lg


@pytest.fixture(scope="module")
def sampler():
    kern = lg.default_inclusive_kernel(lg.li6())
    cfg = lg.default_configs("6Li")[1]
    return lg.InclusiveSampler(kern, cfg, lg.generator_scenario(lg.Scenario()))


@pytest.fixture(scope="module")
def category():
    return lg.tensor_thirds_plan(0.7, 0.6).categories[0]


def test_polligen_sample_category_schema(sampler, category):
    ev = sampler.sample_n(category, 5000, seed=7, run=1, bunch=0)
    for key in ("x", "q2", "y", "phi", "m", "cell"):
        assert key in ev, key
        assert isinstance(ev[key], np.ndarray)
        assert ev[key].size == 5000
    assert ev["cell"].dtype.kind == "i"
    assert ev["category"] == category.name
    assert ev["lam_e"] == category.lam_e
    assert ((ev["phi"] >= 0.0) & (ev["phi"] < 2.0 * np.pi)).all()
    assert (ev["y"] > 0).all() and (ev["y"] < 1).all()


def test_in_cell_placement_matches_the_cell_index(sampler, category):
    """The polligen invariant `exp(logx_lo[cell]) <= x <= exp(logx_hi[cell])`
    (evgen/tests/test_sampler.py) -- it is what makes `cell` usable for the
    Mode-W reweighting."""
    ev = sampler.sample_n(category, 4000, seed=11)
    cell = ev["cell"]
    assert (ev["x"] >= np.exp(sampler.logx_lo[cell]) * (1 - 1e-12)).all()
    assert (ev["x"] <= np.exp(sampler.logx_hi[cell]) * (1 + 1e-12)).all()
    assert (ev["q2"] >= np.exp(sampler.logq2_lo[cell]) * (1 - 1e-12)).all()
    assert (ev["q2"] <= np.exp(sampler.logq2_hi[cell]) * (1 + 1e-12)).all()


def test_determinism_same_stream(sampler, category):
    a = sampler.sample_n(category, 3000, seed=20260713, run=2, bunch=5)
    b = sampler.sample_n(category, 3000, seed=20260713, run=2, bunch=5)
    for k in ("x", "q2", "y", "phi", "m", "cell"):
        assert np.array_equal(a[k], b[k]), k


def test_determinism_is_thread_count_independent(sampler, category):
    one = sampler.sample_n(category, 8000, seed=5, nthreads=1)
    many = sampler.sample_n(category, 8000, seed=5, nthreads=8)
    for k in ("x", "q2", "y", "phi", "m", "cell"):
        assert np.array_equal(one[k], many[k]), k


def test_a_different_stream_gives_different_events(sampler, category):
    a = sampler.sample_n(category, 2000, seed=1, run=1, bunch=0)
    b = sampler.sample_n(category, 2000, seed=2, run=1, bunch=0)
    assert not np.array_equal(a["x"], b["x"])
    c = sampler.sample_n(category, 2000, seed=1, run=1, bunch=1)
    assert not np.array_equal(a["x"], c["x"])


def test_event0_offsets_the_counter(sampler, category):
    whole = sampler.sample_n(category, 400, seed=3)
    tail = sampler.sample_n(category, 300, seed=3, event0=100)
    assert np.array_equal(whole["x"][100:], tail["x"])


def test_weights_for_consumes_the_events_dict(sampler):
    plan = lg.tensor_thirds_plan(0.7, 0.6)
    cats = list(plan.categories)
    ev = sampler.sample_n(cats[0], 500, seed=9)
    w = sampler.weights_for(ev, cats)
    assert w.shape == (500, len(cats))
    assert np.isfinite(w).all() and (w > 0).all()


def test_cross_sections_are_luminosity_share_invariant(sampler, category):
    """A share moves counts, never cross sections (bookkeeping.hpp)."""
    scaled = lg.SpinCategory(category.name, category.j, category.populations,
                             category.lam_e, category.pe, category.theta_s,
                             category.phi_s, 0.25 * category.lumi_fraction)
    assert sampler.sigma_tot_pb(scaled) == sampler.sigma_tot_pb(category)


def test_rng_stream_is_reproducible():
    a = [lg.Rng(42, 1, 2, i).uniform() for i in range(5)]
    b = [lg.Rng(42, 1, 2, i).uniform() for i in range(5)]
    assert a == b
    assert lg.Rng(42, 1, 2, 0).uniform() != lg.Rng(43, 1, 2, 0).uniform()
