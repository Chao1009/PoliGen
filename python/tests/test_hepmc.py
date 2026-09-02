"""HepMC3 output written from Python is readable by pyhepmc."""

import numpy as np
import pytest

import lipolgen as lg

pytestmark = pytest.mark.skipif(not lg._lipolgen.HAVE_HEPMC3,
                                reason="built without HepMC3")

pyhepmc = pytest.importorskip("pyhepmc")


@pytest.fixture(scope="module")
def written(tmp_path_factory):
    cfg = lg.make_config(isotope="6Li", channel="tagged-alpha", config=1,
                         events=300, seed=13)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    path = tmp_path_factory.mktemp("hepmc") / "tagged.hepmc"
    n = p.write_hepmc(str(path))
    assert n == 300
    return p, str(path)


def test_pyhepmc_reads_it_back(written):
    p, path = written
    n = 0
    with pyhepmc.open(path) as f:
        for ev in f:
            n += 1
            assert len(ev.particles) > 3
    assert n == 300


def test_four_momentum_conserves_on_readback(written):
    """Status 1 plus the ONE status-3 pdg-92 pseudo-particle is the whole
    final state (docs/USAGE.md 6), so it balances the beams."""
    p, path = written
    with pyhepmc.open(path) as f:
        for i, ev in enumerate(f):
            if i >= 25:
                break
            beams = [q for q in ev.particles if q.status == 4]
            out = [q for q in ev.particles
                   if q.status == 1 or (q.status == 3 and q.pid == 92)]
            assert len(beams) == 2
            tot_in = np.zeros(4)
            for q in beams:
                tot_in += (q.momentum.e, q.momentum.px, q.momentum.py,
                           q.momentum.pz)
            tot_out = np.zeros(4)
            for q in out:
                tot_out += (q.momentum.e, q.momentum.px, q.momentum.py,
                            q.momentum.pz)
            assert np.allclose(tot_in, tot_out, rtol=1e-7,
                               atol=1e-6 * max(1.0, abs(tot_in[0])))


def test_pdg_codes_round_trip(written):
    """10-digit ion codes survive the write/read cycle; charge does not live
    in a HepMC3 record at all, so `charge_residual` guards it on the
    generator side (test_pipeline) and this only pins the identities."""
    _, path = written
    with pyhepmc.open(path) as f:
        for i, ev in enumerate(f):
            if i >= 25:
                break
            pids = [(q.pid, q.status) for q in ev.particles]
            assert (11, 4) in pids                    # beam electron
            assert (1000030060, 4) in pids            # beam 6Li
            assert (11, 1) in pids                    # scattered electron
            assert (1000020040, 1) in pids            # tagged alpha
            assert (1000010020, 3) in pids            # off-shell struck d
            assert (92, 3) in pids                    # the T0 hadronic system


def test_writer_context_manager(tmp_path):
    cfg = lg.make_config(events=50, seed=17)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    path = tmp_path / "inclusive.hepmc"
    with lg.HepMC3Writer(str(path)) as w:
        for ev in p.generate_range(0, 50):
            w.write(ev)
    with pyhepmc.open(str(path)) as f:
        assert sum(1 for _ in f) == 50


def test_spin_attributes_are_on_the_event(written):
    """The plans/04 #17 convention: the ion spin labels ride on the GenEvent."""
    _, path = written
    with pyhepmc.open(path) as f:
        ev = next(iter(f))
        names = set(ev.attributes)
    for want in ("spin_J", "spin_M", "P_e", "lam_e", "P_z", "P_zz"):
        assert want in names, (want, sorted(names))


def test_electron_generated_mass_is_the_pdg_mass(written):
    """hepmc_writer.cpp's electron generated-mass rule: the core builds
    electrons massless (standard DIS kinematics), and the writer records
    generated_mass = the PDG electron mass for any |pdg| == 11 particle with
    Particle::mass == 0.0, so downstream readers see an on-shell electron."""
    _, path = written
    with pyhepmc.open(path) as f:
        for i, ev in enumerate(f):
            if i >= 25:
                break
            electrons = [q for q in ev.particles if abs(q.pid) == 11]
            assert electrons
            for q in electrons:
                assert q.generated_mass == pytest.approx(0.51099895e-3, rel=1e-9)
