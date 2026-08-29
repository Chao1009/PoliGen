"""T2 (PYTHIA) output exported in the polligen `HFSSample` .npz schema.

The gate is the HFS truth identity carried THROUGH the export: the file's own
`p4` column must reproduce

    Sum (E - p_z)_hadrons = (k + p_N - k')^-

(pythia_bridge.hpp `hfs_sigma_empz_exact`), so a transposed or truncated
column would show up here rather than in a C++ unit test.
"""

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import export

pytestmark = pytest.mark.skipif(not lg._lipolgen.HAVE_PYTHIA8,
                                reason="built without PYTHIA 8")

N_EVENTS = 200


@pytest.fixture(scope="module")
def hadronized():
    """200 inclusive 6Li events through the T2 tier, records kept."""
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=N_EVENTS, seed=20260713)
    beams = lg.default_configs(cfg.isotope)[cfg.beam_config]
    bridge = lg.PythiaBridge(beams, lg.PythiaBridgeOptions())
    lg.set_pythia_hadronizer(cfg, bridge)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    cols = p.generate(0, events=True, nthreads=1)   # bridge is not re-entrant
    return p, bridge, cols


def test_the_bridge_hadronized_most_events(hadronized):
    _, bridge, _ = hadronized
    st = bridge.stats
    assert st.n_called == N_EVENTS
    assert st.n_ok > 0.9 * N_EVENTS, (st.n_ok, st.n_failed)


def test_hadrons_are_on_the_record(hadronized):
    _, _, cols = hadronized
    n_with_hadrons = 0
    for ev in cols["events"]:
        had = [q for q in ev.particles
               if q.role == lg.Role.Hadron and q.status == lg.Status.Final]
        if had:
            n_with_hadrons += 1
            # the T0 pseudo-particle is demoted, never double counted
            x = [q for q in ev.particles if q.role == lg.Role.HadronicX]
            assert all(q.status == lg.Status.Intermediate for q in x)
    assert n_with_hadrons > 0.9 * N_EVENTS


def test_empz_identity_through_the_export(hadronized):
    """Sum(E - p_z) over the EXPORTED particle arrays == the exact truth."""
    _, _, cols = hadronized
    s = export.hfs_sample(cols["events"])
    p4, offs = s["p4"], s["offsets"]
    empz = p4[:, 0] - p4[:, 3]
    n_checked = 0
    for i, ev in enumerate(cols["events"]):
        lo, hi = offs[i], offs[i + 1]
        if hi == lo:
            continue                     # PYTHIA vetoed this one
        k = ev.find(lg.Role.BeamElectron).p
        kp = ev.find(lg.Role.ScatteredElectron).p
        pn = ev.find(lg.Role.StruckNucleon)
        assert pn is not None, "the bridge appends the struck nucleon"
        truth = lg.hfs_sigma_empz_exact(k, pn.p, kp)
        assert float(empz[lo:hi].sum()) == pytest.approx(truth, rel=1e-6), i
        n_checked += 1
    assert n_checked > 0.9 * N_EVENTS


def test_transverse_balance_through_the_export(hadronized):
    """Sum p_T of the hadrons balances the scattered electron's."""
    _, _, cols = hadronized
    s = export.hfs_sample(cols["events"])
    p4, offs = s["p4"], s["offsets"]
    for i, ev in enumerate(cols["events"]):
        lo, hi = offs[i], offs[i + 1]
        if hi == lo:
            continue
        kp = ev.find(lg.Role.ScatteredElectron).p
        px, py = p4[lo:hi, 1].sum(), p4[lo:hi, 2].sum()
        assert px == pytest.approx(-kp.px, abs=1e-6 * max(1.0, abs(kp.px)))
        assert py == pytest.approx(-kp.py, abs=1e-6 * max(1.0, abs(kp.py)))


def test_hfs_schema_and_kp(hadronized):
    _, _, cols = hadronized
    s = export.hfs_sample(cols["events"], meta=dict(cols["meta"]))
    assert s["offsets"].size == N_EVENTS + 1
    assert s["offsets"][0] == 0 and s["offsets"][-1] == s["pid"].size
    assert s["p4"].shape == (s["pid"].size, 4)
    assert s["charge"].size == s["pid"].size
    for k in ("x", "q2", "y", "weight"):
        assert s[k].size == N_EVENTS
    assert s["kp"].shape == (N_EVENTS, 4)
    # the scattered electron is NOT in the particle list (polligen convention)
    for i, ev in enumerate(cols["events"]):
        kp = ev.find(lg.Role.ScatteredElectron).p
        assert np.allclose(s["kp"][i], [kp.e, kp.px, kp.py, kp.pz])
    assert s["meta"]["sigma_gen_mb"] > 0.0


def test_npz_round_trip_without_polligen(hadronized, tmp_path):
    _, _, cols = hadronized
    path = tmp_path / "hfs.npz"
    written = export.write_hfs_npz(cols["events"], str(path),
                                   meta=dict(cols["meta"]))
    back = export.load_hfs_npz(str(path))
    for k in ("offsets", "pid", "charge", "p4", "x", "q2", "y", "kp",
              "weight"):
        assert np.array_equal(back[k], written[k]), k
    assert back["e_energy"] == pytest.approx(10.0)
    assert back["p_per_nucleon"] == pytest.approx(99.5)
    assert back["meta"]["sigma_gen_mb"] == written["meta"]["sigma_gen_mb"]


def test_polligen_hfssample_reads_the_file(hadronized, tmp_path, polligen):
    """The real consumer: `polligen.hfs.HFSSample.load` + `hadronic_sums`."""
    from polligen.hfs import HFSSample, hadronic_sums

    _, _, cols = hadronized
    path = tmp_path / "hfs_polligen.npz"
    export.write_hfs_npz(cols["events"], str(path), meta=dict(cols["meta"]))
    s = HFSSample.load(str(path))
    assert s.n_events == N_EVENTS
    assert s.e_energy == pytest.approx(10.0)
    assert s.p_per_nucleon == pytest.approx(99.5)
    assert s.s == pytest.approx(4.0 * 10.0 * 99.5)
    assert s.event_index().size == s.pid.size

    sigma, ptx, pty = hadronic_sums(s.p4, s.offsets)
    # polligen's own truth relation, Sigma ~ 2 E_e y + m_N^2/(E_N + p_N)
    ok = np.diff(s.offsets) > 0
    m_n = 0.9383
    e_n = np.sqrt(s.p_per_nucleon ** 2 + m_n ** 2)
    truth = 2.0 * s.e_energy * s.y + m_n ** 2 / (e_n + s.p_per_nucleon)
    assert np.allclose(sigma[ok], truth[ok], rtol=1e-3)
    # and the pT balance against the recorded scattered electron
    assert np.allclose(ptx[ok], -s.kp[ok, 1], atol=1e-6)
    assert np.allclose(pty[ok], -s.kp[ok, 2], atol=1e-6)
