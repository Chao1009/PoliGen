"""The module imports, and the bound surface is the one docs/USAGE.md uses."""

import numpy as np
import pytest

import lipolgen as lg


def test_import_and_version():
    assert lg.__version__ == "0.1.0"
    assert lg._lipolgen.HAVE_HEPMC3 in (True, False)


def test_constants_have_one_definition():
    # constants.hpp is mirrored read-only; nothing here recomputes a number.
    assert lg.TENSOR_LL_SIGN == 1.0
    assert lg.ALPHA_EM == pytest.approx(1.0 / 137.036, rel=0, abs=0)
    assert lg.M_NUCLEON == 0.9383
    assert lg.GEV2_TO_PB == 0.3894e9


@pytest.mark.parametrize("name", [
    # beams / spin / sf / xsec / bookkeeping / sampler / pipeline / io
    "Ion", "BeamConfig", "default_configs", "nucleus_mass",
    "rho_from_populations", "moments_along_axis", "populations_maxent",
    "ToyF2", "ToyG1", "MillerB1", "Li6B1", "NuclearF2", "toy_b1",
    "InclusiveKernel", "EventSpinState", "Amplitudes", "SFTables",
    "azz", "a_cos2phi", "a_parallel", "err_azz",
    "SpinCategory", "RunPlan", "tensor_thirds_plan", "helicity_flip_plan",
    "transverse_tensor_plan", "tensor_flip_plan",
    "Scenario", "generator_scenario", "InclusiveSampler",
    "PipelineConfig", "Pipeline", "PipelineChannel", "OpticsChoice",
    "Event", "Particle", "Vec4", "Rng", "Optics", "Route",
    "route_of", "rp_tagged", "momentum_residual", "charge_residual",
])
def test_bound_names(name):
    assert hasattr(lg, name), name


def test_optional_tiers_follow_the_build():
    if lg._lipolgen.HAVE_HEPMC3:
        assert hasattr(lg, "HepMC3Writer")
    if lg._lipolgen.HAVE_PYTHIA8:
        assert hasattr(lg, "PythiaBridge")
        assert hasattr(lg, "hfs_summary")
    if lg._lipolgen.HAVE_LHAPDF:
        assert hasattr(lg, "LhapdfSF")
        assert hasattr(lg, "Epps21Ratio")


def test_beams_match_the_documented_menu():
    cfgs = lg.default_configs("6Li")
    assert [c.ion_momentum_per_nucleon for c in cfgs] == [40.8, 99.5, 137.5]
    cfgs7 = lg.default_configs("7Li")
    assert [c.ion_momentum_per_nucleon for c in cfgs7] == [40.8, 99.5, 117.9]


def test_spin_maxent_rational_anchors():
    # (P_z, P_zz) = (8/13, 4/13) at J=1 and (0.7, 0.4) at J=3/2 -- the t=3
    # spin-temperature anchors of the validation matrix.
    p1 = lg.populations_maxent(1.0, 8.0 / 13.0)
    m1 = lg.moments_along_axis(1.0, p1)
    assert m1.vector == pytest.approx(8.0 / 13.0, rel=1e-12)
    assert m1.tensor == pytest.approx(4.0 / 13.0, rel=1e-10)
    p32 = lg.populations_maxent(1.5, 0.7)
    m32 = lg.moments_along_axis(1.5, p32)
    assert m32.vector == pytest.approx(0.7, rel=1e-12)
    assert m32.tensor == pytest.approx(0.4, rel=1e-10)


def test_arrays_are_numpy_and_own_their_buffers():
    r = lg.Rng(1, 0, 0, 0)
    a = r.uniforms(10)
    assert isinstance(a, np.ndarray) and a.dtype == np.float64
    assert a.base is not None            # the capsule holding the vector
    assert ((a > 0) & (a < 1)).all()


def test_vec4_algebra():
    a = lg.Vec4(5.0, 1.0, 2.0, 3.0)
    b = lg.Vec4(1.0, 0.0, 0.0, 1.0)
    c = a - b
    assert (c.e, c.px, c.py, c.pz) == (4.0, 1.0, 2.0, 2.0)
    assert a.pt() == pytest.approx(np.hypot(1.0, 2.0))
    assert np.allclose(a.as_array(), [5.0, 1.0, 2.0, 3.0])


# --- knobs that landed after the first binding pass (2026-08-29) -------------
#
# Each of these is a physics DEFAULT that is a choice, or an escape hatch for a
# known-broken path.  They are pinned here so that the next time one moves, a
# Python test says so instead of a downstream script silently changing answers.

def test_polarized_emc_transfer_names_its_baseline():
    """`cbt/tmt_valence_scale(baseline)` -- the header's own numbers."""
    assert lg.EMC_BASELINE_DEFAULT == lg.EmcBaseline.Epps21
    assert lg.cbt_valence_scale(lg.EmcBaseline.LegacyTable) == 1.0
    assert lg.tmt_valence_scale(lg.EmcBaseline.LegacyTable) == pytest.approx(
        0.397, abs=5e-4)
    assert lg.cbt_valence_scale(lg.EmcBaseline.Epps21) == pytest.approx(
        0.5105, abs=5e-4)
    assert lg.tmt_valence_scale(lg.EmcBaseline.Epps21) == pytest.approx(
        0.2027, abs=5e-4)
    assert lg.cbt_valence_scale() == lg.cbt_valence_scale(
        lg.EMC_BASELINE_DEFAULT)
    assert lg.emc_valence_depletion(lg.EmcBaseline.Epps21) == pytest.approx(
        lg.EMC_VALENCE_DEPLETION_EPPS21, rel=1e-12)


def test_kernel_defaults_are_the_python_ones():
    """target_mass and g2_scale: xsec.py's own defaults since 2026-08-29."""
    opt = lg.InclusiveKernel.Options()
    assert opt.target_mass is True
    assert opt.g2_scale == 1.0
    kern = lg.InclusiveKernel(lg.li6(), opt)
    assert kern.target_mass is True and kern.g2_scale == 1.0
    # g2_scale is the twist-3 handle: it multiplies g2 and nothing else
    opt.g2_scale = 0.0
    zero = lg.InclusiveKernel(lg.li6(), opt)
    assert zero.tables(0.1, 10.0, True).g2 == 0.0
    assert zero.tables(0.1, 10.0, True).g1 == kern.tables(0.1, 10.0, True).g1


def test_amplitudes_refuse_a_mismatched_spin():
    """The rank-2 branch reads state.j, tensor_moments reads ion().spin."""
    kern = lg.InclusiveKernel(lg.li6())            # J = 1
    bad = lg.EventSpinState(lam_e=1, pe=0.7, j=1.5, m=1.5)
    with pytest.raises(RuntimeError):
        kern.amplitudes(kern.tables(0.1, 10.0, True), 0.1, 10.0, 4000.0, bad)


def test_coherent_xpom_model_is_invertible():
    xp = lg.CoherentXpomModel()
    assert xp.m_x_min == lg.COHERENT_MX_MIN_DEFAULT
    q2, w2 = 10.0, 4000.0
    x_pom = xp.x_pom_min(q2, w2)
    m_x2 = lg.CoherentXpomModel.m_x2_of(x_pom, q2, w2)
    assert m_x2 == pytest.approx(xp.m_x_min ** 2, rel=1e-10)
    assert lg.CoherentXpomModel.x_pom_of(m_x2, q2, w2) == pytest.approx(
        x_pom, rel=1e-12)
    rng = lg.Rng(1, 0, 0, 0)
    draws = [xp.draw(q2, w2, rng) for _ in range(200)]
    assert all(x_pom * (1 - 1e-12) <= d <= xp.x_pom_max for d in draws)


def test_a_hadronizer_on_the_coherent_channel_is_refused():
    """C4: PythiaBridge v0 has no coherent-diffractive target, so the config
    refuses the pair rather than losing four-momentum silently."""
    cfg = lg.make_config(channel="coherent", events=10)
    assert cfg.hadronize_coherent is False
    cfg.hadronizer = lambda ev, rng: None
    with pytest.raises(RuntimeError, match="COHERENT"):
        cfg.validate()
    cfg.hadronize_coherent = True          # the deliberate escape hatch
    cfg.validate()


def test_optics_luminosity_factor_is_reachable():
    p = lg.Pipeline(lg.make_config(channel="tagged-alpha", events=10),
                    lg.tensor_thirds_plan(0.7, 0.6))
    assert p.optics_lumi_factor() == pytest.approx(1.0)
    cfg = lg.make_config(channel="tagged-alpha", events=10, optics="tagging")
    tagging = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    assert 0.0 < tagging.optics_lumi_factor() <= 1.0
