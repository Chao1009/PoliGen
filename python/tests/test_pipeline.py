"""Pipeline generation through the binding: schema, physics gates, exports."""

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import export

N_TAGGED = 100000


@pytest.fixture(scope="module")
def tagged_run():
    """1e5 6Li alpha-tagged events at 10 GeV e x 99.5 GeV/u, YR optics."""
    cfg = lg.make_config(isotope="6Li", channel="tagged-alpha", config=1,
                         events=N_TAGGED, seed=1,
                         optics="yr-high-acceptance")
    plan = lg.tensor_thirds_plan(0.7, 0.6)
    p = lg.Pipeline(cfg, plan)
    return p, p.generate(0)


@pytest.fixture(scope="module")
def inclusive_run():
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=20000, seed=2)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    return p, p.generate(0)


# ------------------------------------------------------------------ schema

def test_columnar_schema(tagged_run):
    p, cols = tagged_run
    assert cols["x"].size == N_TAGGED
    for key in ("x", "q2", "y", "phi", "m", "m_ion", "m_struck", "k",
                "cos_theta_k", "phi_k", "pT", "theta", "p_lab", "R", "xL",
                "kx", "ky", "kz", "phi_spec", "route", "weight", "lam_e",
                "e_prime", "eta_e", "kp"):
        assert key in cols, key
    assert cols["kp"].shape == (N_TAGGED, 4)
    assert cols["category"].shape == (N_TAGGED,)
    assert set(np.unique(cols["category"])) <= set(cols["category_names"])
    assert cols["meta"]["channel"] == "tagged-6Li-alpha"
    assert cols["meta"]["sigma_gen_mb"] == pytest.approx(
        cols["meta"]["sigma_pb"] * 1e-9)


def test_kinematics_are_inside_the_generator_window(tagged_run):
    _, c = tagged_run
    sc = lg.generator_scenario(lg.Scenario())
    assert (c["q2"] >= sc.q2_min).all()
    assert (c["y"] >= sc.y_min).all() and (c["y"] <= sc.y_max).all()
    assert (c["x"] > 0).all() and (c["x"] <= 1).all()
    assert ((c["phi"] >= 0) & (c["phi"] < 2 * np.pi)).all()
    assert np.isfinite(c["kp"]).all()


def test_spin_labels_follow_the_run_plan(tagged_run):
    p, c = tagged_run
    names = list(c["category_names"])
    assert names == [cat.name for cat in p.plan.categories]
    for cat in p.plan.categories:
        mask = c["category"] == cat.name
        assert mask.any()
        assert np.allclose(c["pzz"][mask], c["pzz"][mask][0])
        assert set(np.unique(c["m_ion"][mask])) <= {-1.0, 0.0, 1.0}
        assert (c["lam_e"][mask] == cat.lam_e).all()
    # the struck cluster of the 6Li alpha tag is the embedded deuteron, S_c = 1
    assert set(np.unique(c["m_struck"])) <= {-1.0, 0.0, 1.0}


# ---------------------------------------------------------------- physics

def test_tag_fraction_matches_the_documented_yr_number(tagged_run):
    """docs/USAGE.md 5: 0.0248 at the Yellow Report high-acceptance optics,
    10 GeV e x 99.5 GeV/u, 6Li alpha tag (Python's own 0.0247)."""
    _, c = tagged_run
    frac = float(np.mean(export.rp_accepted(c)))
    assert frac == pytest.approx(0.0248, rel=0.20), frac


def test_tag_fraction_rises_at_the_tagging_optics(tagged_run):
    """The route is recomputed from the record, so one sample can be priced
    at several envelopes (docs/USAGE.md 5)."""
    p, _ = tagged_run
    p_u = p.beam_config.ion_momentum_per_nucleon
    tag = lg.optics_for(lg.OpticsChoice.Tagging, "6Li", p_u)
    evs = p.generate_range(0, 4000)
    seen = np.array([lg.rp_tagged(e, tag, p.pot_config) for e in evs])
    assert 0.15 < seen.mean() < 0.45


def test_four_momentum_and_charge_close(tagged_run):
    p, _ = tagged_run
    for i in (0, 7, 1234, N_TAGGED - 1):
        ev = p.event(i)
        r = lg.momentum_residual(ev)
        scale = lg.momentum_scale(ev)
        assert abs(r.e) < 1e-8 * max(scale.e, 1.0)
        assert abs(r.pz) < 1e-8 * max(abs(scale.pz), 1.0)
        assert abs(r.px) < 1e-8 and abs(r.py) < 1e-8
        assert lg.charge_residual(ev) == 0.0


def test_determinism_and_thread_independence():
    cfg = lg.make_config(isotope="6Li", channel="tagged-alpha", config=1,
                         events=20000, seed=99)
    plan = lg.tensor_thirds_plan(0.7, 0.6)
    a = lg.Pipeline(cfg, plan).generate(0, nthreads=1)
    b = lg.Pipeline(cfg, plan).generate(0, nthreads=8)
    for k in ("x", "q2", "y", "phi", "m_ion", "m_struck", "k", "route"):
        assert np.array_equal(a[k], b[k]), k


def test_generate_n_is_a_prefix_of_the_run(tagged_run):
    p, c = tagged_run
    head = p.generate(500)
    assert np.array_equal(head["x"], c["x"][:500])


def test_generate_lumi_rebuilds_the_run():
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=1000, seed=4)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    out = p.generate_lumi(1e-3, poisson=False)
    expected = sum(p.plan.lumi_share_vector(1e-3)[k] * s
                   for k, s in enumerate(p.sigma_per_category_pb()))
    assert out["x"].size == pytest.approx(expected, rel=0.02)
    assert out["meta"]["n_events"] == out["x"].size


# ----------------------------------------------------------------- exports

def test_tagged_dict_is_the_polligen_schema(tagged_run):
    _, c = tagged_run
    d = export.tagged_dict(c)
    for k in export.TAGGED_KEYS:
        assert k in d, k
    assert d["route"].dtype.kind == "i"
    # polligen.tagged.rp_accepted reads exactly this
    assert export.rp_accepted(d).sum() > 0


def test_inclusive_dict_is_the_polligen_schema(inclusive_run):
    _, c = inclusive_run
    d = export.inclusive_dict(c)
    for k in ("x", "q2", "y", "phi", "m", "category", "lam_e"):
        assert k in d, k
    # `cell` is a sampler internal, absent on the pipeline path (documented)
    assert "cell" not in d


def test_columns_from_events_reproduces_the_fast_path(tagged_run):
    """The Event-record path and the C++ columnar path are the same numbers."""
    p, c = tagged_run
    evs = p.generate_range(0, 200)
    d = export.columns_from_events(evs, p.optics, p.pot_config)
    for k in ("x", "q2", "y", "phi", "m_ion", "m_struck", "k", "cos_theta_k",
              "phi_k", "pT", "theta", "p_lab", "R", "xL", "kx", "ky", "kz",
              "phi_spec"):
        assert np.allclose(d[k], c[k][:200], rtol=1e-12, atol=0,
                           equal_nan=True), k
    assert np.array_equal(d["route"], c["route"][:200])
    assert np.array_equal(d["category"], c["category"][:200])


def test_columns_npz_round_trip(tagged_run, tmp_path):
    _, c = tagged_run
    path = tmp_path / "cols.npz"
    export.write_columns_npz(c, str(path))
    with np.load(str(path), allow_pickle=False) as f:
        assert np.array_equal(f["x"], c["x"])
        assert np.array_equal(f["route"], c["route"])
        import json
        meta = json.loads(str(f["meta"]))
    assert meta["channel"] == "tagged-6Li-alpha"


def test_coherent_channel_runs():
    cfg = lg.make_config(isotope="6Li", channel="coherent", config=1,
                         events=5000, seed=6, optics="tagging",
                         coherent={"f0": 0.04, "slope_b": 50.0, "amp": 0.01})
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    c = p.generate(0)
    assert c["x"].size == 5000
    assert (c["t"] > 0).all() and (c["t"] <= cfg.coherent_t_max).all()
    assert np.mean(c["t"]) == pytest.approx(1.0 / cfg.coherent.slope_b,
                                            rel=0.10)
