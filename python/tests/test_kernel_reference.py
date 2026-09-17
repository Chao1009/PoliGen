"""The bound kernel reproduces the polligen reference tables at rtol 1e-12.

Same JSON, same tolerance and same kernel construction as the C++ doctest
suite (`tests/test_reference.cpp`, `build_kernel`): 6Li with the TOY b1 shape
(`toy_b1(..., mode="toy")`) and the 1e-3 toy Delta, target_mass off.  The
point of running it from Python too is that the BINDING does not silently
reorder or truncate an argument on the way in.
"""

import json
import os

import pytest

import lipolgen as lg

RTOL = 1e-12


@pytest.fixture(scope="module")
def xsec_ref(reference_dir):
    path = os.path.join(reference_dir, "xsec.json")
    if not os.path.exists(path):
        pytest.skip("validation/reference/xsec.json absent")
    with open(path) as f:
        return json.load(f)


def _kernel_6li(target_mass=False):
    """`tests/test_reference.cpp::build_kernel("6Li", ...)`.

    `target_mass` is set EXPLICITLY: which way the C++ default points is a
    physics decision that has already moved once (it became true on
    2026-08-29), and a reference block is dumped for one specific choice.
    """
    opt = lg.InclusiveKernel.Options()
    opt.b1_func = lambda x, q2, f1: lg.toy_b1(x, q2, f1, lg.B1Mode.Toy)
    opt.delta_func = lambda x, q2, f1: lg.toy_delta_gluon(x, q2, f1, 1e-3)
    opt.target_mass = target_mass
    return lg.InclusiveKernel(lg.li6(), opt)


def _kernel_6li_tensor_gamma():
    """`build_kernel("6Li", "kernel_tensor_gamma")`: the exact finite-gamma
    tensor sector with both higher-twist slots filled by the dump's SCENARIO
    shapes (b3 and b4 are unmeasured)."""
    opt = lg.InclusiveKernel.Options()
    opt.b1_func = lambda x, q2, f1: lg.toy_b1(x, q2, f1, lg.B1Mode.Toy)
    opt.delta_func = lambda x, q2, f1: lg.toy_delta_gluon(x, q2, f1, 1e-3)
    opt.b3_func = lambda x, q2, f1: 0.05 * f1
    opt.b4_func = lambda x, q2, f1: -0.02 * f1
    opt.tensor_gamma = True
    opt.target_mass = False
    return lg.InclusiveKernel(lg.li6(), opt)


def test_constants_match_the_reference(xsec_ref):
    c = xsec_ref["constants"]
    assert lg.ALPHA_EM == c["ALPHA_EM"]
    assert lg.GEV2_TO_PB == c["GEV2_TO_PB"]
    assert lg.M_NUCLEON == c["M_NUCLEON"]
    assert lg.TENSOR_LL_SIGN == c["TENSOR_LL_SIGN"]
    assert lg.LI6_CLUSTER_POLARIZATION == pytest.approx(
        c["LI6_CLUSTER_POLARIZATION"], rel=RTOL)
    for name, row in c["ions"].items():
        ion = lg.ion_by_name(name)
        assert ion.eff_pol_p == pytest.approx(row["eff_pol_p"], rel=RTOL)
        assert ion.eff_pol_n == pytest.approx(row["eff_pol_n"], rel=RTOL)
    for name, row in c["ions"].items():
        ion = lg.ion_by_name(name)
        assert (ion.A, ion.Z, ion.N) == (row["A"], row["Z"], row["N"])
        assert ion.spin == row["spin"]
        assert ion.mass_per_nucleon() == pytest.approx(
            row["mass_per_nucleon"], rel=RTOL)


def test_structure_function_tables(xsec_ref):
    kern = _kernel_6li()
    pts = xsec_ref["grid_points"]
    ref = xsec_ref["results"]["6Li"]["kernel_default"]["tables"]
    for i, pt in enumerate(pts):
        t = kern.tables(pt["x"], pt["q2"], True)
        for key in ("f1", "f2", "g1", "g2", "b1", "b2", "b3", "b4", "delta"):
            got = getattr(t, key)
            want = ref[key][i]
            if want == 0.0:
                assert got == 0.0, key
            else:
                assert got == pytest.approx(want, rel=RTOL), (key, i)


@pytest.mark.parametrize("block_name", ["kernel_default",
                                        "kernel_default_target_mass",
                                        "kernel_tensor_gamma"])
def test_amplitudes_every_state(xsec_ref, block_name):
    """(w_avg, a1, a2) for every (lam_e, axis, m) row of the reference."""
    kern = (_kernel_6li_tensor_gamma() if block_name == "kernel_tensor_gamma"
            else _kernel_6li(block_name.endswith("target_mass")))
    pts = xsec_ref["grid_points"]
    res = xsec_ref["results"]["6Li"]
    s = res["beam_config"]["s"]
    block = res[block_name]
    tabs = [kern.tables(p["x"], p["q2"], True) for p in pts]
    n_checked = 0
    for amp in block["amplitudes"]:
        st = lg.EventSpinState(lam_e=int(amp["lam_e"]), pe=amp["pe"], j=1.0,
                               m=amp["m"], theta_s=amp["theta_s"],
                               phi_s=amp["phi_s"])
        for i, pt in enumerate(pts):
            got = kern.amplitudes(tabs[i], pt["x"], pt["q2"], s, st, True)
            for key, value in (("w_avg", got.w_avg), ("a1", got.a1),
                               ("a2", got.a2)):
                want = amp[key][i]
                if want == 0.0:
                    assert abs(value) < 1e-300, (key, amp["m"], i)
                else:
                    assert value == pytest.approx(want, rel=RTOL), \
                        (key, amp["m"], amp["theta_s"], i)
            n_checked += 3
    assert n_checked == 3 * len(block["amplitudes"]) * len(pts)
    assert n_checked > 100          # the table really was walked


def test_dsigma_and_dsigma_unpol(xsec_ref):
    kern = _kernel_6li()
    pts = xsec_ref["grid_points"]
    res = xsec_ref["results"]["6Li"]
    s = res["beam_config"]["s"]
    block = res["kernel_default"]
    for i, pt in enumerate(pts):
        assert kern.dsigma_unpol(pt["x"], pt["q2"], s) == pytest.approx(
            block["dsigma_unpol"][i], rel=RTOL)
    for row in block["dsigma"]:
        st = lg.EventSpinState(lam_e=int(row["lam_e"]), pe=row["pe"], j=1.0,
                               m=row["m"], theta_s=row["theta_s"],
                               phi_s=row["phi_s"])
        for i, pt in enumerate(pts):
            for jph, phi in enumerate(row["phi_grid"]):
                got = kern.dsigma(pt["x"], pt["q2"], phi, s, st, True)
                assert got == pytest.approx(
                    row["dsigma_per_point"][i][jph], rel=RTOL)


def test_asymmetry_helpers_match(xsec_ref):
    """asymmetries.hpp through the binding, against the same table."""
    res = xsec_ref["results"]["6Li"]
    s = res["beam_config"]["s"]
    ref = res["kernel_default"]["asymmetries"]
    pts = xsec_ref["grid_points"]
    kern = _kernel_6li()
    for i, pt in enumerate(pts):
        x, q2 = pt["x"], pt["q2"]
        y = q2 / (s * x)
        assert y == pytest.approx(ref["y"][i], rel=RTOL)
        assert lg.r_sigma_lt(x, q2) == pytest.approx(ref["R"][i], rel=RTOL)
        assert lg.depolarization_d(y, x, q2) == pytest.approx(
            ref["depolarization_d"][i], rel=RTOL)
        t = kern.tables(x, q2, True)
        assert lg.a_parallel(t.g1, t.f1, y, x, q2) == pytest.approx(
            ref["a_parallel"][i], rel=RTOL)
        if ref["azz"][i] is not None:
            assert lg.azz(t.b1, t.f1, t.f2, x, y, t.b2) == pytest.approx(
                ref["azz"][i], rel=RTOL)
        if ref["a_cos2phi"][i] is not None:
            assert lg.a_cos2phi(t.delta, t.f1, t.f2, x, y) == pytest.approx(
                ref["a_cos2phi"][i], rel=RTOL)


def test_tensor_gamma_block(xsec_ref):
    """The exact finite-gamma tensor sector through the binding: Cosyn
    Eqs. (16), (17a)-(17e), (24) and the (h0, h1, h2) harmonics, at the same
    rtol 1e-12 as everything else.

    `tensor_gamma` is OFF in a default kernel -- this block is the only one
    of the file's four built with it on.
    """
    res = xsec_ref["results"]["6Li"]
    block = res.get("kernel_tensor_gamma")
    if block is None or "tensor_gamma" not in block:
        pytest.skip("xsec.json carries no tensor_gamma block")
    s = res["beam_config"]["s"]
    pts = xsec_ref["grid_points"]
    assert lg.InclusiveKernel.Options().tensor_gamma is False
    kern = _kernel_6li_tensor_gamma()
    assert kern.tensor_gamma is True
    tg = block["tensor_gamma"]
    tabs = [kern.tables(p["x"], p["q2"], True) for p in pts]
    for i, pt in enumerate(pts):
        x, q2 = pt["x"], pt["q2"]
        y = q2 / (s * x)
        g2 = lg.gamma_squared(x, q2)
        cq, sq = lg.theta_q_cos_sin(y, g2)
        assert cq == pytest.approx(tg["cos_theta_q"][i], rel=RTOL)
        assert sq == pytest.approx(tg["sin_theta_q"][i], rel=RTOL)
        t = tabs[i]
        f = lg.cosyn_tensor_sfs(t.b1, t.b2, t.b3, t.b4, x, g2)
        for key, got in zip(("F_TLL_T", "F_TLL_L", "F_TLT", "F_TTT"), f):
            assert got == pytest.approx(tg[key][i], rel=RTOL), (key, i)
        fu = lg.cosyn_unpolarized_sfs(t.f1, t.f2, x, g2)
        for key, got in zip(("F_UU_T", "F_UU_L"), fu):
            assert got == pytest.approx(tg[key][i], rel=RTOL), (key, i)
    n = 0
    for row in tg["harmonics"]:
        st = lg.EventSpinState(lam_e=0, pe=0.0, j=1.0, m=row["m"],
                               theta_s=row["theta_s"], phi_s=row["phi_s"])
        for i, pt in enumerate(pts):
            y = pt["q2"] / (s * pt["x"])
            h = kern.tensor_harmonics_gamma(tabs[i], pt["x"], pt["q2"], y, st)
            for key, got in (("h0", h.h0), ("h1", h.h1), ("h2", h.h2)):
                want = row[key][i]
                if want == 0.0:
                    assert abs(got) < 1e-300, (key, i)
                else:
                    assert got == pytest.approx(want, rel=RTOL), (key, i)
                n += 1
    assert n > 100


def test_gamma_zero_identity_is_exact_for_any_b2():
    """Cosyn Eqs. (17) collapse to the massless HJM b-sector at gamma = 0 --
    for ANY b2, with b3 and b4 cancelling identically.  That identity is what
    makes the switch reversible against a massless reference at all."""
    x, y, r = 0.1413, 0.02713, 0.3
    f1, f2 = 1.0, 2.0 * x * (1.0 + r)
    b1 = 0.02
    for b2s, b3, b4 in ((2.0, 0.0, 0.0), (1.3, 0.0, 0.0), (2.0, 0.011, -0.007)):
        b2 = b2s * x * b1
        eps = lg.epsilon_gamma(y, 0.0)
        f_t, f_l, _f_lt, f_tt = lg.cosyn_tensor_sfs(b1, b2, b3, b4, x, 0.0)
        fu_t, fu_l = lg.cosyn_unpolarized_sfs(f1, f2, x, 0.0)
        got = -(f_t + eps * f_l) / (fu_t + eps * fu_l)
        kern = b1 + (1.0 - y) / (x * y * y) * b2
        dphi = f1 + (1.0 - y) / (x * y * y) * f2
        assert got == pytest.approx(kern / dphi, rel=1e-12)
        assert f_tt == 0.0


def test_reference_manifest_lists_exactly_the_tables_on_disk(reference_dir,
                                                             repo_root):
    """`_manifest.json` is the directory's own index, and nothing read it.

    A key walk over the eight reference blobs (2026-09-16) found every numeric
    block asserted somewhere -- except `xsec.json`'s `toy_formulas` (they are
    documentation strings) and `_manifest.json` itself, which `grep -rn
    _manifest tests python/tests` matched ZERO times: only the three dump /
    repin scripts ever wrote it.  So a table dropped from the directory, or
    added to it and left out of the index, went unnoticed by every gate.
    """
    import json
    import os

    with open(os.path.join(reference_dir, "_manifest.json")) as f:
        manifest = json.load(f)

    on_disk = {p for p in os.listdir(reference_dir)
               if p.endswith(".json") and p != "_manifest.json"}
    listed = set(manifest["files"])
    assert listed == on_disk, (
        "manifest lists %s; the directory holds %s"
        % (sorted(listed - on_disk), sorted(on_disk - listed)))

    # Every blob names the script that made it, and that script exists.
    generators = manifest["generators"]
    assert set(generators) == on_disk, sorted(set(generators) ^ on_disk)
    for blob, script in generators.items():
        # The paths are recorded repo-relative with the repo name in front.
        rel = script.split("LiPolGen/", 1)[-1]
        assert os.path.isfile(os.path.join(repo_root, rel)), (blob, script)
