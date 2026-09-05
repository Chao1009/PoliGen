"""The module imports, and the bound surface is the one docs/USAGE.md uses."""

import math
import os
import re

import numpy as np
import pytest

import lipolgen as lg


def test_import_and_version(repo_root):
    # Not pinned to a literal: the version comes from CMakeLists.txt's
    # project() VERSION (single source of truth, see LIPOLGEN_VERSION) and
    # would otherwise drift out of sync with a hardcoded expectation here.
    # "Single source" is a claim, so it is VERIFIED rather than assumed: the
    # compiled-in string is checked for shape AND against that source, which
    # is the check that runs in the in-tree flow (where importlib.metadata
    # knows nothing).  When the package IS installed (wheel or
    # `pip install -e .`), the METADATA version scikit-build-core derived from
    # the same source is cross-checked too.
    assert re.fullmatch(r"\d+\.\d+\.\d+", lg.__version__), lg.__version__
    with open(os.path.join(repo_root, "CMakeLists.txt")) as fh:
        m = re.search(r"project\(\s*LiPolGen\s+VERSION\s+([0-9.]+)", fh.read())
    assert m, "CMakeLists.txt has no project(LiPolGen VERSION ...)"
    assert lg.__version__ == m.group(1)
    try:
        import importlib.metadata as importlib_metadata
        installed_version = importlib_metadata.version("lipolgen")
    except Exception:
        installed_version = None
    if installed_version is not None:
        assert lg.__version__ == installed_version
    assert lg._lipolgen.HAVE_HEPMC3 in (True, False)


def test_constants_have_one_definition():
    # constants.hpp is mirrored read-only; nothing here recomputes a number.
    assert lg.TENSOR_LL_SIGN == -1.0
    assert lg.ALPHA_EM == pytest.approx(1.0 / 137.036, rel=0, abs=0)
    assert lg.M_NUCLEON == 0.9383
    assert lg.GEV2_TO_PB == 0.3894e9
    # the 6Li cluster wave function -- beams.hpp, one definition, and the
    # inclusive eff_pol slots are a third of it each (plans/04 #6)
    assert lg.P_D_LI6 == 0.0867
    assert lg.P_D_DEUTERON == 0.045
    assert lg.LI6_CLUSTER_POLARIZATION == (
        (1.0 - 1.5 * lg.P_D_LI6) * (1.0 - 1.5 * lg.P_D_DEUTERON))
    assert lg.LI6_CLUSTER_POLARIZATION == pytest.approx(0.81123, abs=5e-6)
    li6 = lg.li6()
    assert li6.Z * li6.eff_pol_p == pytest.approx(
        lg.LI6_CLUSTER_POLARIZATION, rel=1e-15)
    assert li6.N * li6.eff_pol_n == pytest.approx(
        lg.LI6_CLUSTER_POLARIZATION, rel=1e-15)
    assert lg.deuteron().eff_pol_p == 1.0 - 1.5 * lg.P_D_DEUTERON
    # the retired Cloet convention stays reachable and is 1.233x this one
    assert lg.LI6_NAIVE_ONE_THIRD == 1.0 / 3.0
    assert lg.LI6_NAIVE_ONE_THIRD / (lg.LI6_CLUSTER_POLARIZATION / 3.0) == \
        pytest.approx(1.233, abs=5e-4)


@pytest.mark.parametrize("name", [
    # beams / spin / sf / xsec / bookkeeping / sampler / pipeline / io
    "Ion", "BeamConfig", "default_configs", "nucleus_mass",
    "rho_from_populations", "moments_along_axis", "populations_maxent",
    "ToyF2", "ToyG1", "MillerB1", "Li6B1", "NuclearF2", "toy_b1",
    "InclusiveKernel", "EventSpinState", "Amplitudes", "SFTables",
    "azz", "a_cos2phi", "a_parallel", "a_parallel_exact", "err_azz",
    "depolarization_effective", "gamma_squared", "epsilon_gamma",
    "depolarization_d_gamma", "eta_gamma", "theta_q_cos_sin",
    "cosyn_tensor_sfs", "cosyn_unpolarized_sfs", "TensorHarmonics",
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
        # MSTW2008 LO rides the PYTHIA tier, not the LHAPDF one: the grid is
        # PYTHIA's own pdfdata/mstw2008lo.00.dat (see python/tests/
        # test_mstw_sf.py).
        assert hasattr(lg, "MstwSF")
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
        0.5322, abs=5e-4)
    assert lg.tmt_valence_scale(lg.EmcBaseline.Epps21) == pytest.approx(
        0.2113, abs=5e-4)
    assert lg.cbt_valence_scale() == lg.cbt_valence_scale(
        lg.EMC_BASELINE_DEFAULT)
    assert lg.emc_valence_depletion(lg.EmcBaseline.Epps21) == pytest.approx(
        lg.EMC_VALENCE_DEPLETION_EPPS21, rel=1e-12)


def test_kernel_defaults_are_the_python_ones():
    """target_mass and g2_scale: xsec.py's own defaults since 2026-08-29."""
    opt = lg.InclusiveKernel.Options()
    assert opt.target_mass is True
    assert opt.g2_scale == 1.0
    assert opt.tensor_gamma is False      # the exact tensor sector is OPT-IN
    kern = lg.InclusiveKernel(lg.li6(), opt)
    assert kern.target_mass is True and kern.g2_scale == 1.0
    assert kern.tensor_gamma is False
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


def test_a_hadronizer_on_the_coherent_channel_is_accepted():
    """C4, closed 2026-08-30.  The coherent channel names its own T2 target
    (`Role.Pomeron`, P_IP = P_ion - P_recoil), so a hadronizer on it is an
    ordinary configuration -- no refusal, no `hadronize_coherent` escape
    hatch, and no `PipelineConfig` attribute of that name any more."""
    cfg = lg.make_config(channel="coherent", events=10)
    assert not hasattr(cfg, "hadronize_coherent")
    cfg.hadronizer = lambda ev, rng: None
    cfg.validate()


def test_the_coherent_record_carries_a_pomeron():
    cfg = lg.make_config(channel="coherent", events=8)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    ev = p.event(0)
    ip = [q for q in ev.particles if q.role == lg.Role.Pomeron]
    assert len(ip) == 1
    assert ip[0].pdg == 990
    assert ip[0].charge == 0.0
    assert ip[0].status == lg.Status.Intermediate
    rec = [q for q in ev.particles if q.role == lg.Role.IntactRecoil][0]
    ion = ev.particles[1]
    for c in ("e", "px", "py", "pz"):
        assert getattr(ip[0].p, c) == pytest.approx(
            getattr(ion.p, c) - getattr(rec.p, c), abs=1e-9)


@pytest.mark.skipif(not lg.HAVE_PYTHIA8, reason="no PYTHIA 8 tier")
def test_the_coherent_t2_options_are_bound():
    o = lg.PythiaBridgeOptions()
    assert o.coherent_t2 == lg.CoherentT2.Pomeron
    assert o.pom_set == 6
    assert o.pom_rescale == 1.0
    o.coherent_t2 = lg.CoherentT2.Off
    assert o.coherent_t2 == lg.CoherentT2.Off


def test_optics_luminosity_factor_is_reachable():
    p = lg.Pipeline(lg.make_config(channel="tagged-alpha", events=10),
                    lg.tensor_thirds_plan(0.7, 0.6))
    assert p.optics_lumi_factor() == pytest.approx(1.0)
    cfg = lg.make_config(channel="tagged-alpha", events=10, optics="tagging")
    tagging = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    assert 0.0 < tagging.optics_lumi_factor() <= 1.0


# ------------------------------------------------------- the VMC cluster waves


def _have_vmc():
    return os.path.isfile(os.path.join(lg._lipolgen.data_dir(),
                                       "vmc", "momenta", "li6_ad1.momentum"))


@pytest.mark.skipif(not _have_vmc(), reason="data/vmc is not present")
def test_vmc_cluster_waves_are_reachable_and_carry_their_provenance():
    ch = lg._lipolgen.li6_alpha_channel(
        source=lg._lipolgen.ClusterWaveSource.VmcAV18)
    assert "VMC" in ch.label
    ls = sorted(w.l for w in ch.waves)
    assert ls == [0, 2]
    for w in ch.waves:
        assert w.vmc is not None
        assert w.vmc.l == w.l
        assert "li6_ad1.momentum" in w.vmc.provenance
        # zero outside the table: the ANL grid stops at 5 fm^-1
        assert w.vmc(2.0) == 0.0
    # `p_d` is ignored on the VMC path -- P_D is a property of the wave function
    d = [w for w in ch.waves if w.l == 2][0]
    assert d.prob == pytest.approx(lg._lipolgen.VMC_P_D_LI6, rel=1e-12)
    assert d.prob == pytest.approx(0.01935, abs=2e-5)
    # the S-D relative sign the momentum densities alone cannot give
    s = [w for w in ch.waves if w.l == 0][0]
    assert s.vmc(0.05) > 0 > d.vmc(0.05)
    # ... and the Hulthen default carries no table at all
    assert all(w.vmc is None for w in lg._lipolgen.li6_alpha_channel().waves)


@pytest.mark.skipif(not _have_vmc(), reason="data/vmc is not present")
def test_make_config_takes_cluster_wave_by_name():
    assert lg.make_config(channel="tagged-alpha", events=10).cluster_wave \
        == lg._lipolgen.ClusterWaveSource.Hulthen
    cfg = lg.make_config(channel="tagged-alpha", events=500,
                         cluster_wave="vmc")
    assert cfg.cluster_wave == lg._lipolgen.ClusterWaveSource.VmcAV18
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    assert "VMC" in p.tagged_model.channel.label
    cols = p.generate(0, 500, 1)
    assert cols["x"].size == 500


def test_a_parallel_splits_massless_from_exact():
    """`a_parallel(..., g2=None)` is the massless D(y) g1/F1; with a g2 it is
    `a_parallel_exact`, the E143 finite-gamma form -- the same split the
    Python's `asymmetries.a_parallel` makes."""
    g1, g2, f1, y, x, q2 = 0.03, -0.01, 0.9, 0.05, 0.1413, 3.127
    massless = lg.a_parallel(g1, f1, y, x, q2)
    assert massless == pytest.approx(
        lg.depolarization_d(y, x, q2) * g1 / f1, rel=0, abs=0)
    exact = lg.a_parallel(g1, f1, y, x, q2, g2=g2)
    assert exact == lg.a_parallel_exact(g1, g2, f1, y, x, q2)
    assert exact != massless
    # D_eff is the divisor that inverts the exact A_par without an
    # O(gamma^2) bias; with no rho it is the massless D, bit for bit
    assert lg.depolarization_effective(y, x, q2) == lg.depolarization_d(
        y, x, q2)
    deff = lg.depolarization_effective(y, x, q2, g2_over_g1=g2 / g1)
    assert exact == pytest.approx(deff * g1 / f1, rel=1e-14)


def test_tensor_gamma_is_off_by_default_and_reversible():
    """Switching `tensor_gamma` on must not move a massless number, and the
    two paths must agree as gamma^2 -> 0 (Q^2 -> infinity at fixed x)."""
    x, q2 = 0.1, 1.0e8
    s = q2 / (1.0e-3 * x)
    opt = lg.InclusiveKernel.Options()
    opt.b1_func = lambda xx, qq, f1: lg.toy_b1(xx, qq, f1)
    opt.target_mass = False
    massless = lg.InclusiveKernel(lg.li6(), opt)
    opt.tensor_gamma = True
    exact = lg.InclusiveKernel(lg.li6(), opt)
    st = lg.EventSpinState(lam_e=0, pe=0.0, j=1.0, m=1.0, theta_s=0.0)
    w0 = massless.amplitudes(massless.tables(x, q2), x, q2, s, st).w_avg
    w1 = exact.amplitudes(exact.tables(x, q2), x, q2, s, st).w_avg
    assert w1 == pytest.approx(w0, rel=1e-9)
    # b3/b4 default to zero and are built on every path, but only the exact
    # one reads them
    t = massless.tables(x, q2)
    assert t.b3 == 0.0 and t.b4 == 0.0


# ------------------------------------------------ the triton spectral function


def test_triton_sf_bindings_round_trip():
    """The `TritonSfChoice` knob, the `BreakupOptions.triton_sf` slot and the
    `CiofiSimulaTriton` object itself are all reachable from Python."""
    cfg = lg.PipelineConfig()
    assert cfg.triton_sf == lg.TritonSfChoice.Hulthen           # the default
    cfg.triton_sf = lg.TritonSfChoice.CiofiSimula
    assert cfg.triton_sf == lg.TritonSfChoice.CiofiSimula
    assert lg.make_config(channel="tagged-7Li-alpha", events=10).triton_sf \
        == lg.TritonSfChoice.Hulthen
    assert lg.make_config(channel="tagged-7Li-alpha", events=10,
                          triton_sf="ciofi-simula").triton_sf \
        == lg.TritonSfChoice.CiofiSimula
    with pytest.raises(ValueError):
        lg.make_config(triton_sf="faddeev")

    sf = lg.CiofiSimulaTriton(lg.CiofiSimulaOptions())
    # S_0 = 0.6525 is the 3He(e,e'p)d spectroscopic factor; the CS Eq. (28)
    # sum closes on 1 untuned (test_triton_sf.cpp measures all of this at
    # full precision; here the binding is what is under test).
    assert sf.s0 == pytest.approx(0.6525, rel=1e-3)
    assert sf.s0 + sf.s1 == pytest.approx(1.0, abs=5e-4)
    # no bound nn exists, so a struck proton's two-body weight is exactly 0,
    # and a struck neutron's is the n_0/(n_0 + n_1) ratio
    assert sf.p_two_body(0.1, 2212) == 0.0
    assert sf.p_two_body(0.1, 2112) == \
        sf.n0_cs(0.1) / (sf.n0_cs(0.1) + sf.n1_cs(0.1))
    assert lg.KAPPA_PN_SINGLET == pytest.approx(0.0083122, rel=1e-9)
    assert lg.triton_channel_name(lg.TritonChannel.NeutronPnCont) == "n+(pn)"

    bo = lg.BreakupOptions()
    assert bo.triton_sf is None          # null = the sequential model, bit
    bo.triton_sf = sf                    # for bit; non-null = spectral fn
    assert bo.triton_sf.p_two_body(0.2, 2212) == 0.0


def _li7_t1_events(triton_sf, n=200):
    cfg = lg.make_config(channel="tagged-7Li-alpha", events=n, seed=7,
                         optics="tagging", triton_sf=triton_sf)
    p = lg.Pipeline(cfg, lg.helicity_flip_plan(1.5, 0.7, 0.7))
    return [p.event(i) for i in range(n)]


def test_triton_sf_pipeline_conserves_with_both_options():
    """A small tagged-7Li T1 run per option: per-event whole-record closure,
    same-seed determinism, different records between the two models, and the
    n + (pn) third channel only on ciofi-simula."""
    ciofi = _li7_t1_events("ciofi-simula")
    hulthen = _li7_t1_events("hulthen")
    for evs in (ciofi, hulthen):
        for ev in evs:
            r = lg.momentum_residual(ev)
            assert max(abs(r.e), abs(r.px), abs(r.py), abs(r.pz)) < 1e-9
            assert lg.charge_residual(ev) == 0.0

    def struck_pz(evs):
        return [[q for q in ev.particles
                 if q.role == lg.Role.StruckNucleon][0].p.pz for ev in evs]

    # determinism: rebuilding the same config gives the identical records
    assert struck_pz(_li7_t1_events("ciofi-simula")) == struck_pz(ciofi)
    # a different triton model gives different records (same T0 kinematics,
    # different breakup draws)
    same = sum(a == b for a, b in zip(struck_pz(ciofi), struck_pz(hulthen)))
    assert same < len(ciofi) // 10

    def n_pn_events(evs):
        n = 0
        for ev in evs:
            pn = [q for q in ev.particles
                  if q.role == lg.Role.StruckNucleon][0]
            partners = [q.pdg for q in ev.particles
                        if q.role == lg.Role.PartnerSpectator]
            if pn.pdg == 2112 and len(partners) == 2:
                assert partners.count(2212) == 1   # one p + one n
                n += 1
        return n

    assert n_pn_events(ciofi) > 0        # the channel the sequential model
    assert n_pn_events(hulthen) == 0     # has no room for


# ------------------------------------------------------------ the FSI weight


def test_fsi_weight_reproduces_the_pinned_table_point():
    """GlauberFsiWeight on the production S+D 6Li alpha channel: the k = 0.10,
    theta = 0 row (the pure-S prototype rows are pinned in tests/test_fsi.cpp),
    the shadowed profile, and the wrong-spectator guard."""
    f = lg.GlauberFsiWeight(lg.li6_alpha_channel())
    kin = lg.FsiKinematics(k=0.10, cos_theta_k=1.0)
    assert f.weight(kin) == pytest.approx(0.617, abs=5e-3)
    # Glauber shadowing: sigma_Xalpha = 131.0 mb, NOT 4 x 40 = 160
    assert f.sigma_cluster_mb(40.0) == pytest.approx(131.0, abs=0.1)
    assert 0.0 < f.survival() < 1.0
    assert f.weight_normalised(kin) == pytest.approx(
        f.weight(kin) / f.survival(), rel=1e-12)
    # a weight built for the alpha must not distort another spectator
    wrong = lg.FsiKinematics(k=0.10, cos_theta_k=1.0,
                             spectator_z=1, spectator_a=2)
    assert f.weight(wrong) == 1.0


def test_fsi_pipeline_weights_are_a_weight_never_a_shift():
    """FSI on vs off: bit-identical kinematics, only Event.weight moves, and
    its mean estimates the logged survival probability."""
    unpol = lg.SpinCategory("unpol", 1.0, [1 / 3, 1 / 3, 1 / 3])
    plan = lg.RunPlan([unpol])
    on = lg.Pipeline(lg.make_config(channel="tagged-alpha", events=4000,
                                    fsi="glauber-cluster"), plan)
    off = lg.Pipeline(lg.make_config(channel="tagged-alpha", events=4000),
                      plan)
    assert off.fsi_weight is None
    fw = on.fsi_weight
    assert fw is not None
    assert fw.sigma_eff_mb(0.0) == 40.0
    con, coff = on.generate(0, False, 1), off.generate(0, False, 1)
    # the weight moved NOTHING: every drawn quantity is bit-identical
    for col in ("k", "cos_theta_k", "x", "q2", "alpha_s", "pt_s", "route"):
        assert (con[col] == coff[col]).all(), col
    assert (coff["weight"] == 1.0).all()
    w = con["weight"]
    assert (w >= 0.0).all()
    assert (w != 1.0).all()
    assert w.mean() == pytest.approx(fw.survival(), abs=0.02)


def test_fsi_band_discipline_and_the_channel_guard():
    """sigma_XN = 20 mb absorbs less than 40 mb, and the config refuses FSI
    off the tagged channels."""
    cfg20 = lg.make_config(channel="tagged-alpha", events=10, fsi=
                           "glauber-cluster", fsi_sigma_mb=20.0)
    assert cfg20.fsi_sigma_mb == 20.0
    with pytest.raises(RuntimeError):
        lg.make_config(channel="inclusive", events=10, fsi="glauber-cluster")
    with pytest.raises(RuntimeError):
        lg.make_config(channel="coherent", events=10, fsi="glauber-nucleon")


def test_fsi_weight_reaches_the_npz_export(tmp_path):
    unpol = lg.SpinCategory("unpol", 1.0, [1 / 3, 1 / 3, 1 / 3])
    p = lg.Pipeline(lg.make_config(channel="tagged-alpha", events=200,
                                   fsi="glauber-cluster"), lg.RunPlan([unpol]))
    cols = p.generate(0, False, 1)
    path = str(tmp_path / "fsi.npz")
    lg.write_columns_npz(cols, path)
    d = np.load(path, allow_pickle=True)
    assert (d["weight"] != 1.0).all()
    assert (d["weight"] == cols["weight"]).all()


@pytest.mark.skipif(not lg.HAVE_HEPMC3, reason="no HepMC3 writer")
def test_fsi_weight_reaches_hepmc_weights0(tmp_path):
    unpol = lg.SpinCategory("unpol", 1.0, [1 / 3, 1 / 3, 1 / 3])
    p = lg.Pipeline(lg.make_config(channel="tagged-alpha", events=5,
                                   fsi="glauber-cluster"), lg.RunPlan([unpol]))
    path = str(tmp_path / "fsi.hepmc")
    events = [p.event(i) for i in range(5)]
    with lg.HepMC3Writer(path) as w:
        for ev in events:
            w.write(ev)
    # Asciiv3: the per-event weight vector is the "W" record (the run-info
    # header carries a "W nominal" names line -- skip the non-numeric one)
    def _num(tok):
        try:
            float(tok)
            return True
        except ValueError:
            return False
    wlines = [ln.split() for ln in open(path)
              if (ln.startswith("W ") or ln.startswith("W\t"))
              and _num(ln.split()[1])]
    assert len(wlines) == 5
    for ev, ln in zip(events, wlines):
        assert ev.weight != 1.0
        assert float(ln[1]) == pytest.approx(ev.weight, rel=1e-9)


# ------------------------------- open items C4, C5.1 - C5.5 (2026-09-04)
#
# The Python mirrors of tests/test_coherent.cpp T10b,
# tests/test_cluster_config.cpp T23 and tests/test_tagged.cpp T24-T26.
# Numbers: docs/open_items/run_2026-09-03/phase_C_numbers.md sections C4, C5.


def test_c4_delta_b_is_defined_and_eps_b0_is_a_deuteron_number():
    sc = lg._lipolgen.CoherentScenario()
    # a_2(m) = -(Delta B_m / 2)|t|, exactly, at every m
    for t in (0.01, 0.05, 0.2, 0.3):
        for m in (0, 1, -1):
            assert sc.a2_m_state(t, m) == pytest.approx(
                -0.5 * sc.delta_b_m(m) * t, rel=1e-14)
    # eps_b0 = 2 Delta B_{+-1}/B = -Delta B_0/B (the old label was off by a sign)
    assert sc.eps_b0 == pytest.approx(2.0 * sc.delta_b_m(1) / sc.slope_b,
                                      rel=1e-14)
    assert sc.eps_b0 == pytest.approx(-sc.delta_b_m(0) / sc.slope_b, rel=1e-14)
    # slope_at_azimuth IS the definition, and its phi average is B
    assert sc.slope_at_azimuth(0.0, 0) == pytest.approx(
        sc.slope_b + sc.delta_b_m(0), rel=1e-14)
    assert sc.slope_at_azimuth(math.pi / 4.0, 1) == pytest.approx(
        sc.slope_b, abs=1e-12)
    # the shipped default implies a 6Li quadrupole 11.4x the measured one
    q_matter = lg.quadrupole_from_a2_slope(sc.a2_m_state(1.0, 1), 6)
    assert q_matter == pytest.approx(-1.8690781872, rel=1e-9)
    assert 0.5 * q_matter / lg.LI6_QUADRUPOLE_FM2 == pytest.approx(
        11.4246832958, rel=1e-9)
    # and the measured-6Li band it should be compared against
    s = lg.ClusterConfigSampler()
    band = s.quadrupole_band_fm2()
    eps = [-4.0 * lg.a2_from_quadrupole(2.0 * q, 6, 1.0, 1) / sc.slope_b
           for q in band]
    assert eps[0] == pytest.approx(-0.0070023, rel=1e-4)
    assert eps[1] == pytest.approx(-0.0171209, rel=1e-4)
    assert eps[2] == pytest.approx(-0.0526853, rel=1e-4)
    assert abs(sc.eps_b0) > abs(eps[2])
    # eps_b0 and slope_b are NOT independent: the physics is their product
    hi = lg._lipolgen.CoherentScenario()
    hi.slope_b = 60.0
    assert (lg.quadrupole_from_a2_slope(hi.a2_m_state(1.0, 1), 6) / q_matter
            == pytest.approx(1.2, rel=1e-12))


def test_c5_1_eta_is_a_converter_not_a_second_dial():
    s = lg.ClusterConfigSampler()
    q = s.quadrupole_for_eta(lg.LI6_ETA_DS_GK)
    assert q == pytest.approx(-0.1842160146606, rel=1e-9)
    o = lg._lipolgen.ClusterConfigOptions()
    o.quadrupole_target_fm2 = q
    d = lg.ClusterConfigSampler(o)
    assert d.asymptotic_ds_ratio() == pytest.approx(lg.LI6_ETA_DS_GK, rel=1e-12)
    # exact on a dialled sampler too, because eta is linear in the dial
    assert d.quadrupole_for_eta(lg.LI6_ETA_DS_GK) == pytest.approx(q, rel=1e-12)
    # the GK band spans a SIGN CHANGE in Q
    sig = math.hypot(lg.LI6_ETA_DS_GK_STAT, lg.LI6_ETA_DS_GK_SYST)
    assert s.quadrupole_for_eta(lg.LI6_ETA_DS_GK - sig) == pytest.approx(
        -0.4004752268, rel=1e-9)
    assert s.quadrupole_for_eta(lg.LI6_ETA_DS_GK + sig) > 0.0
    with pytest.raises(RuntimeError):
        s.quadrupole_for_eta(-0.10)


@pytest.mark.skipif(not _have_vmc(), reason="data/vmc is not present")
def test_c5_2_the_vmc_tables_carry_their_monte_carlo_band():
    L = lg._lipolgen
    c0 = L.li6_alpha_channel(source=L.ClusterWaveSource.VmcAV18)
    for w in c0.waves:
        assert w.vmc.has_errors
        assert len(w.vmc.dpsi) == len(w.vmc.psi)
        cor, quad = w.vmc.norm2_error()
        assert cor > quad > 0.0
    cor, quad = c0.waves[1].vmc.norm2_error()
    assert cor / c0.waves[1].vmc.norm2() == pytest.approx(0.016176, rel=1e-2)
    # n_sigma = 0 is bit for bit the cached table
    assert L.li6_alpha_channel(source=L.ClusterWaveSource.VmcAV18,
                               vmc_mc_sigma=0.0).waves[1].prob \
        == c0.waves[1].prob
    # the band: P_D moves ~1.1 % per sigma, the tensor dilution 0.02 %
    m0 = L.TaggedModel(c0)
    mp = L.TaggedModel(L.li6_alpha_channel(
        source=L.ClusterWaveSource.VmcAV18, vmc_mc_sigma=1.0))
    assert mp.channel.waves[1].prob / lg.VMC_P_D_LI6 == pytest.approx(
        1.011346, rel=1e-4)
    assert abs(mp.tensor_dilution() / m0.tensor_dilution() - 1.0) < 3e-4
    # a band that would be silently ignored is refused
    with pytest.raises(RuntimeError):
        L.li6_alpha_channel(source=L.ClusterWaveSource.Hulthen,
                            vmc_mc_sigma=1.0)
    cfg = lg.make_config(channel="tagged-alpha", events=10)
    cfg.cluster_vmc_mc_sigma = 1.0
    with pytest.raises(RuntimeError):
        cfg.validate()
    cfg.cluster_wave = L.ClusterWaveSource.VmcAV18
    cfg.validate()


@pytest.mark.skipif(not _have_vmc(), reason="data/vmc is not present")
def test_c5_4_the_deuteron_control_on_av18():
    L = lg._lipolgen
    assert L.deuteron_av18_p_d() == pytest.approx(0.0575998919874, rel=1e-9)
    assert L.deuteron_av18_p_d() / lg.P_D_DEUTERON == pytest.approx(1.28,
                                                                    rel=1e-4)
    h = L.TaggedModel(L.deuteron_channel())
    v = L.TaggedModel(L.deuteron_channel(source=L.ClusterWaveSource.VmcAV18))
    assert h.channel.waves[1].prob == lg.P_D_DEUTERON
    assert v.channel.waves[1].prob == L.deuteron_av18_p_d()
    assert v.vector_dilution() / h.vector_dilution() - 1.0 == pytest.approx(
        -0.020268, rel=1e-3)
    assert v.tensor_dilution() / h.tensor_dilution() - 1.0 == pytest.approx(
        -0.011822, rel=1e-3)
    # the relative S-D sign does NOT flip: psi_2 = +W for the real deuteron
    assert v.channel.waves[0].vmc(0.05) > 0.0
    assert v.channel.waves[1].vmc(0.05) > 0.0
    assert not v.channel.waves[1].vmc.has_errors


@pytest.mark.skipif(not _have_vmc(), reason="data/vmc is not present")
def test_c5_5_the_inclusive_tagged_drift_and_why_substitution_is_not_the_fix():
    L = lg._lipolgen
    # on the Hulthen default the two ARE one wave function
    hul = L.TaggedModel(L.li6_alpha_channel())
    assert hul.vector_dilution() == pytest.approx(1.0 - 1.5 * lg.P_D_LI6,
                                                  rel=2e-5)
    assert lg.LI6_CLUSTER_POLARIZATION == lg.li6_cluster_polarization(
        lg.P_D_LI6, lg.P_D_DEUTERON)
    # under VmcAV18 they drift: +11.61 % vector, +6.58 % rank-2
    vmc = L.TaggedModel(L.li6_alpha_channel(
        source=L.ClusterWaveSource.VmcAV18))
    assert vmc.vector_dilution() / hul.vector_dilution() - 1.0 \
        == pytest.approx(0.116131, rel=1e-4)
    assert vmc.tensor_dilution() / lg.LI6_B1_RANK2_TRANSFER - 1.0 \
        == pytest.approx(0.065762, rel=1e-4)
    assert lg.LI6_CLUSTER_POLARIZATION_VMC / lg.LI6_CLUSTER_POLARIZATION - 1.0 \
        == pytest.approx(0.116119, rel=1e-4)
    # substitution is NOT the fix: it moves AWAY from the ab initio anchor
    ab = lg.LI6_POLARIZATION_VMC_SIX_BODY
    assert lg.LI6_CLUSTER_POLARIZATION / ab - 1.0 == pytest.approx(-0.04336,
                                                                   rel=1e-3)
    assert lg.LI6_CLUSTER_POLARIZATION_VMC / ab - 1.0 == pytest.approx(
        +0.06772, rel=1e-3)
    both = lg.li6_cluster_polarization(lg.VMC_P_D_LI6, L.deuteron_av18_p_d())
    assert both == pytest.approx(0.8870761569, rel=1e-8)
    assert abs(lg.LI6_CLUSTER_POLARIZATION / ab - 1.0) \
        < abs(lg.LI6_CLUSTER_POLARIZATION_VMC / ab - 1.0)


@pytest.mark.skipif(not _have_vmc(), reason="data/vmc is not present")
def test_c5_5b_one_run_one_deuteron():
    """`--cluster-wave vmc` selects ONE deuteron everywhere it is read.

    Until 2026-09-04 the alpha-d relative motion came from the ANL VMC
    AV18+UX overlap while the EMBEDDED deuteron stayed on the scenario
    Hulthen P_D = 0.045 in both places a tagged-alpha run reads it -- the
    struck cluster's g1 and the T1 spin draw -- so every polarized
    tagged-alpha observable was 2.069 % high against the wave function the
    flag claims to select.
    """
    L = lg._lipolgen
    h = L.li6_alpha_channel()
    v = L.li6_alpha_channel(source=L.ClusterWaveSource.VmcAV18)
    # the default is beams.DEUTERON() itself, bit for bit
    assert h.dis_target.eff_pol_p == lg.deuteron().eff_pol_p
    assert h.dis_target.eff_pol_p == lg.DEUTERON_VECTOR_POLARIZATION
    # the VMC one is the same Ion with the AV18 D state, and nothing else
    assert (v.dis_target.name, v.dis_target.A, v.dis_target.Z,
            v.dis_target.spin) == (h.dis_target.name, h.dis_target.A,
                                   h.dis_target.Z, h.dis_target.spin)
    assert v.dis_target.eff_pol_p == lg.vector_dilution_of(
        L.deuteron_av18_p_d())
    assert v.dis_target.eff_pol_n == v.dis_target.eff_pol_p
    # the size that was wrong, and it is EXACT: g1A is linear in eff_pol
    was_high = h.dis_target.eff_pol_p / v.dis_target.eff_pol_p
    assert was_high - 1.0 == pytest.approx(0.0206872095, rel=1e-8)
    k_now = L.InclusiveKernel(v.dis_target)
    k_before = L.InclusiveKernel(lg.deuteron())
    for x, q2 in ((0.05, 2.0), (0.1, 5.0), (0.3, 10.0), (0.5, 20.0)):
        assert k_before.tables(x, q2).g1 / k_now.tables(x, q2).g1 \
            == pytest.approx(was_high, rel=1e-12)
    # the T1 breakup follows too.  `BreakupOptions.source` defaults to
    # Hulthen and `Pipeline` forwards `cluster_wave` into it; the model it
    # builds is `deuteron_channel(beta, p_d, source)`, and the gate is that
    # the dilution the SPIN DRAW implies is the one the RATE uses.  (The
    # `ClusterBreakup` object itself is checked in C++ T27.)
    assert L.BreakupOptions().source == L.ClusterWaveSource.Hulthen
    assert L.TaggedModel(L.deuteron_channel()).vector_dilution() \
        == pytest.approx(h.dis_target.eff_pol_p, rel=1e-5)
    assert L.TaggedModel(L.deuteron_channel(
        source=L.ClusterWaveSource.VmcAV18)).vector_dilution() \
        == pytest.approx(v.dis_target.eff_pol_p, rel=1e-5)
    assert L.deuteron_av18().eff_pol_p == v.dis_target.eff_pol_p
    # what one --cluster-wave vmc run now says 6Li's polarization is: C5.5's
    # THIRD row, 0.887076, not the mongrel 0.905427 it used to be
    whole = L.TaggedModel(v).vector_dilution() * v.dis_target.eff_pol_p
    assert whole == pytest.approx(
        lg.li6_cluster_polarization(lg.VMC_P_D_LI6, L.deuteron_av18_p_d()),
        rel=2e-5)
    assert whole / lg.LI6_CLUSTER_POLARIZATION_VMC - 1.0 == pytest.approx(
        -0.020270, rel=2e-3)
    # 7Li deliberately does not move: no AV18 A = 3 wave function exists here
    assert L.li7_alpha_channel(source=L.ClusterWaveSource.VmcAV18) \
        .dis_target.eff_pol_p == L.li7_alpha_channel().dis_target.eff_pol_p


def test_c5_5b_cluster_wave_is_refused_where_it_is_never_read():
    L = lg._lipolgen
    for channel, isotope, tagged in (("inclusive", "6Li", False),
                                     ("coherent", "6Li", False),
                                     ("tagged-alpha", "6Li", True),
                                     ("tagged-d-p", "d", True)):
        cfg = lg.make_config(isotope=isotope, channel=channel, events=10)
        cfg.validate()                      # Hulthen is legal everywhere
        cfg.cluster_wave = L.ClusterWaveSource.VmcAV18
        if tagged:
            cfg.validate()
        else:
            with pytest.raises(RuntimeError):
                cfg.validate()


def test_shifted_by_sigma_is_the_identity_without_errors():
    """`VmcRadial.shifted_by_sigma` is bound in `python/bindings.cpp` and was
    covered only indirectly, through `li6_alpha_channel(vmc_mc_sigma=...)`.
    A table built from Python carries no error column, so every shift of it
    must return the table unchanged -- that is what makes the band's zero row
    bit for bit today's numbers.
    """
    L = lg._lipolgen
    k = [0.0, 0.1, 0.2, 0.4]
    psi = [1.0, 0.6, 0.25, 0.05]
    v = L.VmcRadial(k, psi, 0, "hand-built, no MC error column")
    assert not v.has_errors and list(v.dpsi) == []
    for n in (0.0, 1.0, -3.0):
        s = v.shifted_by_sigma(n)
        assert list(s.psi) == psi
        assert list(s.k) == k
        assert s.l == v.l
        assert s.provenance == v.provenance
        assert s.norm2() == v.norm2()


@pytest.mark.skipif(not _have_vmc(), reason="data/vmc is not present")
def test_shifted_by_sigma_moves_psi_by_n_times_dpsi_fully_correlated():
    """On a real ANL table: psi -> psi + n*dpsi elementwise and in the SAME
    direction at every k (the fully correlated envelope, not a per-point
    error), dpsi itself is carried through unchanged so the band can be
    re-applied, the provenance records the shift, and n = 0 returns the table
    itself.  The `vmc_mc_sigma` run knob must agree with it point by point.
    """
    L = lg._lipolgen
    w = L.li6_alpha_channel(source=L.ClusterWaveSource.VmcAV18).waves[1].vmc
    assert w.has_errors
    psi, dpsi = list(w.psi), list(w.dpsi)
    assert any(d != 0.0 for d in dpsi)

    assert list(w.shifted_by_sigma(0.0).psi) == psi
    assert w.shifted_by_sigma(0.0).provenance == w.provenance

    for n in (1.0, -1.0, 2.5):
        s = w.shifted_by_sigma(n)
        assert list(s.k) == list(w.k)
        assert s.l == w.l
        assert list(s.dpsi) == dpsi           # the band survives the shift
        assert "sigma_MC" in s.provenance and "CORRELATED" in s.provenance
        assert len(s.psi) == len(psi)
        for a, b, d in zip(list(s.psi), psi, dpsi):
            assert a == pytest.approx(b + n * d, rel=0, abs=1e-15)

    # what `--cluster-vmc-mc-sigma` does IS this function, not a second path
    shifted = L.li6_alpha_channel(source=L.ClusterWaveSource.VmcAV18,
                                  vmc_mc_sigma=1.0).waves[1].vmc
    assert list(shifted.psi) == pytest.approx(list(w.shifted_by_sigma(1.0).psi),
                                              rel=1e-15, abs=0)
