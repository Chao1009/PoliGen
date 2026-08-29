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


def test_constants_match_the_reference(xsec_ref):
    c = xsec_ref["constants"]
    assert lg.ALPHA_EM == c["ALPHA_EM"]
    assert lg.GEV2_TO_PB == c["GEV2_TO_PB"]
    assert lg.M_NUCLEON == c["M_NUCLEON"]
    assert lg.TENSOR_LL_SIGN == c["TENSOR_LL_SIGN"]
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
        for key in ("f1", "f2", "g1", "g2", "b1", "b2", "delta"):
            got = getattr(t, key)
            want = ref[key][i]
            if want == 0.0:
                assert got == 0.0, key
            else:
                assert got == pytest.approx(want, rel=RTOL), (key, i)


@pytest.mark.parametrize("block_name", ["kernel_default",
                                        "kernel_default_target_mass"])
def test_amplitudes_every_state(xsec_ref, block_name):
    """(w_avg, a1, a2) for every (lam_e, axis, m) row of the reference."""
    kern = _kernel_6li(block_name.endswith("target_mass"))
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
