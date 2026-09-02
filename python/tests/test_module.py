"""The module imports, and the bound surface is the one docs/USAGE.md uses."""

import os

import numpy as np
import pytest

import lipolgen as lg


def test_import_and_version():
    assert lg.__version__ == "0.1.0"
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
