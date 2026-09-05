"""The npz `meta` must distinguish runs that differ in a knob (D3, D4, D5).

`python/bindings.cpp` states the rule on the `b1_*` block -- *without these
keys three otherwise identical npz files are indistinguishable, which is
exactly the "never quote a single row" rule failing silently* -- and until
2026-09-04 it had been applied to `b1_*` and to `rc_*` and to nothing else.
Measured consequences, all three reproduced by the tests below:

  * three FSI runs of one seed had **identical** `meta` while their total
    rate differed by 48 % / 42 % and individual weights by a factor 68.5
    (tagged-6Li-alpha, 20 000 events, seed 1234, tensor-thirds at
    P_z = 0.7 / P_zz = 0.6 / P_e = 0.7, sigma_XN = 40 mb);
  * fifteen `--pom-set` runs had **identical** `meta` while the whole
    hadronic final state moved -- the kaon fraction by a factor 2.8 over the
    twelve DPDF fits (3-10, 12-15) about set 6, 20 000 coherent events per
    set at 6Li config 1, seed 4242;
  * `coherent_t_max` moves the entire |t| spectrum, the tag acceptance and
    every c_2 in the file, and was neither in the `meta` nor reachable from
    the CLI at all.

Every block is CONDITIONAL, following the `rc` precedent: a default run's key
set does not move, which is what keeps the reference gates where they are.
"""

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import cli

PLAN = lg.make_plan("tensor-thirds", j=1.0, pz=0.7, pzz=0.6, pe=0.7)

KIN = ("k", "cos_theta_k", "phi_k", "x", "q2")


def _run(**kw):
    cols = lg.Pipeline(lg.make_config(**kw), PLAN).generate(0, False, 1)
    return cols, dict(cols["meta"])


# ------------------------------------------------------------------ D3, FSI

@pytest.fixture(scope="module")
def fsi_runs():
    out = {}
    for v in ("off", "glauber-cluster", "glauber-nucleon"):
        out[v] = _run(isotope="6Li", channel="tagged-6Li-alpha", events=2000,
                      seed=1234, fsi=v)
    return out


def test_fsi_runs_are_one_event_stream_reweighted(fsi_runs):
    """The premise of the measurement: the weight moves nothing else."""
    base = fsi_runs["off"][0]
    for v in ("glauber-cluster", "glauber-nucleon"):
        for c in KIN:
            assert np.array_equal(np.asarray(base[c]),
                                  np.asarray(fsi_runs[v][0][c])), (v, c)


def test_the_two_fsi_variants_are_not_one(fsi_runs):
    """D3: the inventory said they were identical.  They differ per event."""
    wc = np.asarray(fsi_runs["glauber-cluster"][0]["weight"])
    wn = np.asarray(fsi_runs["glauber-nucleon"][0]["weight"])
    assert not np.allclose(wc, wn, rtol=1e-14)
    r = wn / wc
    assert np.mean(np.abs(r - 1.0) > 0.01) > 0.98
    assert r.max() > 5.0
    assert r.min() < 0.9


def test_meta_distinguishes_the_three_fsi_runs(fsi_runs):
    m_off = fsi_runs["off"][1]
    m_c = fsi_runs["glauber-cluster"][1]
    m_n = fsi_runs["glauber-nucleon"][1]
    assert m_off != m_c and m_c != m_n
    # `off` carries today's key set exactly: the block is conditional.
    assert not [k for k in m_off if k.startswith("fsi")]
    assert m_c["fsi"] == "glauber-cluster"
    assert m_n["fsi"] == "glauber-nucleon"
    assert m_c["fsi_sigma_mb"] == m_n["fsi_sigma_mb"] == 40.0
    # ... and the numbers that ARE the difference between the variants:
    # shadowed 131.0 mb against the unshadowed A sigma_XN = 160.0.
    assert m_c["fsi_sigma_cluster_mb"] == pytest.approx(131.045, abs=1e-2)
    assert m_n["fsi_sigma_cluster_mb"] == pytest.approx(160.006, abs=1e-2)
    assert m_c["fsi_survival"] == pytest.approx(0.520239, abs=1e-5)
    assert m_n["fsi_survival"] == pytest.approx(0.582899, abs=1e-5)
    assert m_c["fsi_formation_ramp"] is False
    assert "fsi_ramp_w_lo" not in m_c


def test_the_sigma_xn_band_is_recorded():
    """`fsi_sigma_mb` carries a documented 20-40 mb band; both ends must be
    distinguishable in the record, or the band cannot be quoted at all."""
    _, m20 = _run(isotope="6Li", channel="tagged-6Li-alpha", events=500,
                  seed=1234, fsi="glauber-cluster", fsi_sigma_mb=20.0)
    _, m40 = _run(isotope="6Li", channel="tagged-6Li-alpha", events=500,
                  seed=1234, fsi="glauber-cluster", fsi_sigma_mb=40.0)
    assert m20 != m40
    assert m20["fsi_sigma_mb"] == 20.0 and m40["fsi_sigma_mb"] == 40.0
    assert m20["fsi_survival"] > m40["fsi_survival"]


# -------------------------------------------------------- D5, the |t| ceiling

def test_coherent_t_max_is_recorded_and_reachable():
    _, m = _run(isotope="6Li", channel="coherent", events=500, seed=7)
    assert m["coherent_t_max"] == lg.COHERENT_T_MAX_DEFAULT == 0.2
    assert m["coherent_eps_b0"] == -0.08
    assert m["coherent_slope_b"] == 50.0
    # The DERIVED positivity edge travels with the file, at the worst case.
    assert m["coherent_t_positivity_edge_pzz_m2"] == pytest.approx(0.245,
                                                                   abs=1e-12)
    cols, m2 = _run(isotope="6Li", channel="coherent", events=500, seed=7,
                    coherent_t_max=0.15)
    assert m2["coherent_t_max"] == 0.15
    assert m != m2
    assert np.asarray(cols["t"]).max() <= 0.15


def test_the_coherent_block_is_conditional():
    _, m = _run(isotope="6Li", channel="inclusive", events=200, seed=7)
    assert not [k for k in m if k.startswith("coherent")]


def test_t_positivity_edge_is_derived_not_fixed():
    """D5: the ceiling's SECOND reason moves with eps_b0; the anchor range,
    which is the first one, does not."""
    sc = lg.CoherentScenario()
    assert sc.t_positivity_edge(-2.0) == pytest.approx(0.245, abs=1e-12)
    assert sc.t_positivity_edge(1.0) == pytest.approx(0.495, abs=1e-12)
    assert sc.positivity_margin(sc.t_positivity_edge(-2.0), -2.0) == \
        pytest.approx(0.0, abs=1e-12)
    sc.eps_b0 = -0.0070          # the MEASURED 6Li quadrupole row
    assert sc.t_positivity_edge(-2.0) == pytest.approx(2.800, abs=1e-9)
    # ... i.e. 9.3x outside the |t| <= 0.30 the deformation input covers.
    assert sc.t_positivity_edge(-2.0) > 9.0 * 0.30


def test_the_cli_exposes_coherent_t_max():
    opts = cli.resolve(["--channel", "coherent", "--coherent-t-max", "0.1"])
    assert opts["coherent_t_max"] == 0.1
    # ... and it is NOT a CoherentScenario field, so it must not land there
    assert "t_max" not in (opts["coherent"] or {})
    assert cli.resolve(["--channel", "coherent"])["coherent_t_max"] is None
