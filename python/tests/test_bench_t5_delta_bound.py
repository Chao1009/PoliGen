# SPDX-License-Identifier: GPL-3.0-or-later
"""validation/benchmarks/t5_delta_moment_bound.py exercised as tests (run
2026-10-26 phase C2; ~3 s per harness run, so not behind LIPOLGEN_BENCH).

The harness applies the n = 2 moment of the gluonic Soffer-type bound
(Detmold-Shanahan PRD 94 (2016) 014507 Eqs. (13), (32)) to the shipped
toy_delta_gluon at the pipeline's 1e-2 scale and to the letter's moment_A /
moment_B, under assumptions its header declares, plus the spin-1
helicity-state positivity -2 <= A_T <= 1 (Cosyn et al. EPJ A 61 (2025) 83
Eq. (2)) on the shipped kernel's own per-state rates.  Its second pass adds
the variation V6 (Delta in DS16's normalisation = 2 x the pre-D4 kernel's, at equal
cos 2phi amplitude), which the tests below re-derive from DS16's printed
Eqs. (3)-(5) with numpy alone; since D4 (2026-10-27) the kernel's is DS16's (DS16_OVER_TREE = 1, V6 = PRE_D4_OVER_JM = 2).

These tests pin what the harness MEASURES and that its verdict follows from
the numbers; none moves or asserts a new value for any tree default.  The
letter's models exceed the bound under V6 (and at dilution 1 under every
reading but the loose edge); that is RECORDED and pinned as recorded, never
"fixed".  Nothing here needs the network.  The numbers need PYTHIA's python
module, its LHAPDF6 plugin and the two LHAPDF sets: without the module or the
plugin the row is BLOCKED; without a set the harness raises LHAPDF's own
"Info file not found for PDF set" message, exits 2 and the report records
blocked (environment), as row 4 (README convention 6; since phase F of run
2026-10-26).  Either way the tests that read numbers skip, naming the reason.
"""
import importlib.util
import math
import os
import shutil
import subprocess
import sys

import pytest

np = pytest.importorskip("numpy")
lg = pytest.importorskip("lipolgen")

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")
_NAME = "t5_delta_moment_bound"
_DATA = "delta_moment_bound_sources.json"


def _load(name=_NAME):
    path = os.path.join(_BENCH, name + ".py")
    spec = importlib.util.spec_from_file_location("bench_" + name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _flip(path, old, new, tmp_path):
    """A copy of `path` with the one occurrence of `old` replaced."""
    raw = open(path, "rb").read()
    assert raw.count(old) == 1, old
    bad = tmp_path / os.path.basename(path)
    bad.write_bytes(raw.replace(old, new))
    return str(bad)


def _run_harness(script, cwd=None, args=()):
    """The harness as a script, with the imported lipolgen on its path."""
    env = dict(os.environ)
    pkg_parent = os.path.dirname(os.path.dirname(os.path.abspath(lg.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [pkg_parent] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    return subprocess.run([sys.executable, script] + list(args), capture_output=True,
                          text=True, env=env, cwd=cwd, timeout=300)


@pytest.fixture(scope="module")
def mod():
    return _load()


_MISSING_SET = "Info file not found for PDF set"


def _run_or_skip(mod, verbose=False):
    """run(), skipping -- naming the reason -- on a BLOCKED row or on a set
    that is not installed (the harness's exit-2 path); any other error fails."""
    try:
        r = mod.run(verbose=verbose)
    except RuntimeError as exc:
        if _MISSING_SET not in str(exc):
            raise
        pytest.skip("environment: %s" % str(exc).splitlines()[0])
    if r["status"] == "blocked":
        pytest.skip(r["reason"])
    return r


@pytest.fixture(scope="module")
def rep(mod):
    return _run_or_skip(mod)


def _sub(rep, name):
    hits = [s for s in rep["subrows"] if s["name"] == name]
    assert len(hits) == 1, name
    return hits[0]


# ------------------------------------------------------------ the vendored file

def test_vendored_file_guards(mod, tmp_path):
    doc = mod.read_vendored()
    body = doc["body"]
    q2 = body["grid"]["q2_gev2"]
    assert len(q2) == 3 and q2[0] == 2.0 and q2[-1] == 10.0
    assert [s["key"] for s in body["sources"]] == ["DS16", "SS90", "LOI18", "CvDM17",
                                                   "CRSZ25", "KS20"]
    assert [s["name"] for s in body["pdf_sets"]] == ["CT18NLO", "NNPDFpol11_100"]
    assert body["pdf_sets"][0]["lhapdf_set_index"] == 14400
    assert body["pdf_sets"][-1]["x_min"] == 1e-05
    pos = body["grid"]["positivity"]
    assert len(pos["x"]) == 7 and pos["x"][0] == 0.001 and pos["x"][-1] == 0.9
    assert body["letter_models"]["moment_A"]["A_as_printed"] == -0.2924
    assert body["letter_models"]["moment_A"]["q2_ref_gev2"] == 4.474
    prov = doc["_provenance"]
    assert "NOT vendored" in prov["documents_outside_repository"]
    assert "PER NUCLEON" in prov["publishing_convention"]
    assert "second pass" in prov["compiled"]
    # one flipped byte in the body is refused ...
    bad = _flip(mod.DATA, b'"q2_ref_gev2": 4.474', b'"q2_ref_gev2": 4.475', tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_reference(bad)
    # ... and one in the provenance (every byte is digested)
    bad = _flip(mod.DATA, b"Facts quoted with attribution", b"Facts quoted with attributiom",
                tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.read_vendored(bad)


def test_vendored_equations_are_the_ones_the_header_applies(mod):
    src = {s["key"]: s for s in mod.load_reference()["sources"]}
    ds = src["DS16"]["equations"]
    assert ds["(13)"].startswith("Delta(x, Q^2) = -(alpha_s(Q^2)/(2 pi)) Tr Q^2 x^2 "
                                 "int_x^1 (dy/y^3) deltaG(y, Q^2)")
    assert ds["(32)"] == "|deltaG(x)| <= (1/2) (f1(x) + (1/2) f1LL(x) + g1(x))"
    assert ds["(9)"] == "C_n(Q^2, mu^2) = -(alpha_s(Q^2)/(2 pi)) Tr Q^2 * 2/(n+2)"
    # the definition the second pass derives V6 from
    assert "P(hH, h'H')_{mu nu, alpha beta} = eps*_mu(h') E*_alpha(H') E_beta(H) " \
           "eps_nu(h)" in ds["(3)"]
    assert "A_{+-,-+} = A_{-+,+-}" in ds["(4)"]
    assert ds["(5)"].startswith("W^{Delta=2}_{mu nu} = (1/2) {")
    assert src["CvDM17"]["equations"]["(36)"] == "|h1TT| <= (1/2) (f1 + f1LL/2 + g1)"
    assert "-0.012 alpha_s(Q^2)" in src["SS90"]["equations"]["(3.5)"]
    assert "a(x, Q^2) = -2 [g_xhat" in src["SS90"]["equations"]["(2.11)"]
    assert "circularly polarized" in src["SS90"]["notes"]
    assert "diag(2/3, -1/3, -1/3)" in src["LOI18"]["equations"]["(7)"]
    assert "- (x (1-y)/2) Delta(x, Q^2) cos(2 phi)" in src["LOI18"]["equations"]["(6)"]
    assert "[-2 <= A_T <= 1]" in src["CRSZ25"]["equations"]["(2)"]
    ks = src["KS20"]["equations"]["(68)"]
    assert "deltaG(x) = -g_xhat/xhat(x) + g_yhat/xhat(x) [28, 48]" in ks
    assert "h1TT,g(x) = -g_xhat/xhat(x) + g_yhat/xhat(x) [39, 41, 49]" in ks
    assert "E_x the linear polarization" in src["KS20"]["equations"]["(22)"]
    for s in src.values():
        assert s["file_sha256"] and len(s["file_sha256"]) == 64


# ------------------------------------------------------------ the algebra

def _trap(f, a, b, n=20001):
    h = (b - a) / (n - 1)
    xs = [a + i * h for i in range(n)]
    ys = [f(x) for x in xs]
    return h * (sum(ys) - 0.5 * (ys[0] + ys[-1]))


def test_n2_moment_of_the_relation(mod):
    """DS16 Eq. (13)'s n = 2 moment is -(alpha_s/(8 pi)) Tr Q^2 int y deltaG:
    checked numerically on deltaG = y^0.5 (1-y)^3 (alpha_s Tr Q^2/(2 pi) = 1)."""
    dg = lambda y: y ** 0.5 * (1.0 - y) ** 3              # noqa: E731
    inner = {}

    def delta(x):
        if x not in inner:
            inner[x] = _trap(lambda y: dg(y) / y ** 3, x, 1.0, 4001)
        return -x * x * inner[x]

    lhs = _trap(lambda x: x * delta(x), 1e-3, 1.0, 801)
    rhs = -0.25 * _trap(lambda y: y * dg(y), 0.0, 1.0)
    # the x < 1e-3 part of the lhs is O(1e-6); trapezoid on both sides
    assert lhs == pytest.approx(rhs, rel=2e-3)


def test_rule3_relation_is_the_same_per_nucleon():
    """Header A8 (rule 3): with x = A x_A and per-nucleon densities D(x) =
    D_A(x_A)/A^2, Eq. (13) in nuclear variables is Eq. (13) per nucleon."""
    a_nuc = 6.0
    dg_n = lambda y: y ** 0.7 * (1.0 - y) ** 4 if y < 1.0 else 0.0   # noqa: E731
    dg_a = lambda ya: a_nuc ** 2 * dg_n(a_nuc * ya)                   # noqa: E731
    for x in (0.02, 0.1, 0.3):
        per_nucleon = -x * x * _trap(lambda y: dg_n(y) / y ** 3, x, 1.0, 20001)
        xa = x / a_nuc
        nuclear = -xa * xa * _trap(lambda ya: dg_a(ya) / ya ** 3, xa, 1.0 / a_nuc, 20001)
        assert nuclear / a_nuc ** 2 == pytest.approx(per_nucleon, rel=1e-9)


# --- the second pass's two findings, re-derived from the printed equations

_G = np.diag([1.0, -1.0, -1.0, -1.0])


def _dot(a, b):
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def _helicity(h):
    """Transverse helicity vectors along z (upper index)."""
    return (-np.array([0, 1, 1j, 0]) if h == 1 else np.array([0, 1, -1j, 0])) / np.sqrt(2)


_XHAT = np.array([0, 1.0, 0, 0], dtype=complex)
_YHAT = np.array([0, 0, 1.0, 0], dtype=complex)
_ZHAT = np.array([0, 0, 0, 1.0], dtype=complex)
#: m_x = +1 along xhat (spin along xhat: circular in the y-z plane); E = xhat is m_x = 0
_MX_PLUS = -(_YHAT + 1j * _ZHAT) / np.sqrt(2)


def test_s2_ds16_projection_operators_give_eq5_and_the_f1_normalisation():
    """DS16 Eq. (3)'s P(hH,h'H') = eps*_mu(h') E*_alpha(H') E_beta(H) eps_nu(h)
    with Eq. (4)'s Delta = A_{+-,-+} = A_{-+,+-} = 1 give Eq. (5)'s transverse
    block (W_xx = -W_yy = 1/2 for E = xhat), and with A_{hH,hH} = 1 the
    unpolarized W_ij = delta_ij, i.e. -F1 g_ij with F1 = (1/3) sum_H A_{+H,+H}."""
    target = {1: _helicity(1), -1: _helicity(-1), 0: _ZHAT}

    def w_of(e, amps):
        out = np.zeros((4, 4), complex)
        for (h, hh, hp, hhp), a in amps.items():
            left = np.conj(_G @ _helicity(hp)) * _dot(np.conj(e), np.conj(target[hhp]))
            right = (_G @ _helicity(h)) * _dot(e, target[hh])
            out += a * np.outer(left, right)
        return out

    w = w_of(_XHAT, {(1, -1, -1, 1): 1.0, (-1, 1, 1, -1): 1.0})
    assert (w[1, 1].real, w[2, 2].real, abs(w[1, 2])) == pytest.approx((0.5, -0.5, 0.0),
                                                                       abs=1e-15)
    w = w_of(_MX_PLUS, {(1, -1, -1, 1): 1.0, (-1, 1, 1, -1): 1.0})
    assert (w[1, 1].real, w[2, 2].real) == pytest.approx((-0.25, 0.25), abs=1e-15)
    unpol = {(h, hh, h, hh): 1.0 for h in (1, -1) for hh in (1, 0, -1)}
    wu = sum(w_of(e, unpol) for e in (_XHAT, _YHAT, _ZHAT)) / 3.0
    assert (wu[1, 1].real, wu[2, 2].real, abs(wu[1, 2])) == pytest.approx((1.0, 1.0, 0.0),
                                                                          abs=1e-15)


def _ds16_cos2phi_over_f1(x, y, q2, e, mass=0.938):
    """DS16 Eq. (5) (Delta = 1) and W = -F1 g (F1 = 1) contracted with the
    unpolarized lepton tensor, target at rest, q along +z, the target state
    `e`; returns the cos(2 phi) coefficient over the phi-average of the F1
    term (phi the lepton-plane azimuth about q, from xhat)."""
    nu_l = q2 / (2 * mass * x)
    beam = nu_l / y
    qz = math.sqrt(nu_l ** 2 + q2)
    p = np.array([mass, 0, 0, 0.0])
    q = np.array([nu_l, 0, 0, qz])
    kz = (beam * nu_l + q2 / 2) / qz
    kt = math.sqrt(beam * beam - kz * kz)
    nu = _dot(p, q)
    qq = _dot(q, q)
    kap = 1.0 - mass * mass * qq / nu ** 2
    pl, ql = _G @ p, _G @ q
    el = _G @ e
    qe = _dot(q, e)
    a = np.conj(el) - np.conj(qe) / (kap * nu) * (pl - mass * mass / nu * ql)
    b = el - qe / (kap * nu) * (pl - mass * mass / nu * ql)
    ptl = pl - nu / qq * ql
    proj = _G - np.outer(ql, ql) / qq + qq / (kap * nu ** 2) * np.outer(ptl, ptl)
    scal = _dot(np.conj(e), e) + mass * mass / (kap * nu ** 2) * np.conj(qe) * qe
    w_delta = 0.5 * (np.outer(a, b) + np.outer(b, a) - proj * scal)
    phis = np.linspace(0.0, 2.0 * math.pi, 361)[:-1]
    d, f = [], []
    for ph in phis:
        k = np.array([beam, kt * math.cos(ph), kt * math.sin(ph), kz])
        kp = k - q
        lep = 2.0 * (np.outer(k, kp) + np.outer(kp, k) - _G * _dot(k, kp))
        d.append(np.real(np.sum(lep * w_delta)))
        f.append(np.real(np.sum(lep * -_G)))
    return 2.0 * float(np.mean(np.array(d) * np.cos(2.0 * phis))) / float(np.mean(f))


def test_s2_the_kernel_coefficient_is_ds16_eq5(mod):
    """SECOND PASS S2, re-taken at D4 (2026-10-27): DS16 Eq. (5) gives -(1/2)(1-y)/y^2 per unit
    Delta for a pure m = +-1 state along an axis perpendicular to q (LOI18 Eq. (6)), and since D4
    so does the kernel, a_cos2phi = -(1/2)(1-y)/y^2 (F2 = 0 isolates the F1 term): DS16_OVER_TREE
    = 1 (the pre-D4 kernel's -(1-y)/y^2 gave 1/2 and the constant 2, now PRE_D4_OVER_JM, V6's).
    The Bjorken limit is approached as 1/Q^2; m_x = 0 carries -2 times m_x = +-1 on both sides."""
    assert (mod.DS16_OVER_TREE, mod.PRE_D4_OVER_JM) == (1.0, 2.0)
    for x, y in ((0.1, 0.5), (0.3, 0.3), (0.01, 0.7)):
        tree = lg.a_cos2phi(1.0, 1.0, 0.0, x, y)
        assert tree == pytest.approx(-0.5 * (1.0 - y) / y ** 2, rel=1e-15)  # D4 (no 0.5 before)
        ratio = _ds16_cos2phi_over_f1(x, y, 1e5, _MX_PLUS) / tree
        assert ratio == pytest.approx(1.0 / mod.DS16_OVER_TREE, abs=2e-5)
        lin = _ds16_cos2phi_over_f1(x, y, 1e5, _XHAT)
        assert lin / _ds16_cos2phi_over_f1(x, y, 1e5, _MX_PLUS) == pytest.approx(-2.0,
                                                                                rel=1e-9)
    # at the row's own Q2 the finite-mass terms are a percent: 0.99214 (0.4961 before D4)
    assert _ds16_cos2phi_over_f1(0.1, 0.5, 5.0, _MX_PLUS) / lg.a_cos2phi(
        1.0, 1.0, 0.0, 0.1, 0.5) == pytest.approx(1.0, rel=0.02)


def test_s2_the_kernel_carries_the_a_cos2phi_coefficient():
    """The shipped kernel's a2 at m = +1, theta_S = 90 deg is the standalone
    a_cos2phi = -(1/2)(1-y)/y^2 Delta/D_phi since D4 (PHYSICS_CHANNELS.md sec. 2)."""
    k = lg.default_inclusive_kernel(lg.li6())
    for x, q2, y in ((0.1, 5.0, 0.5), (0.01, 2.0, 0.7), (0.3, 10.0, 0.3)):
        t = k.tables(x, q2)
        a = k.amplitudes(t, x, q2, q2 / (x * y), lg.EventSpinState(
            lam_e=0, pe=0.0, j=1.0, m=1.0, theta_s=math.pi / 2, phi_s=0.0))
        assert a.a2 == pytest.approx(lg.a_cos2phi(t.delta, t.f1, t.f2, x, y), rel=1e-14)


def test_s1_circular_state_carries_minus_half_the_linear_tensor_asymmetry():
    """SECOND PASS S1: for spin 1, T_ij = <(S_i S_j + S_j S_i)/2> - (2/3)
    delta_ij; T_xx - T_yy = +1/2 in the m_x = +1 (circular) state and -1 in
    the E_x (m_x = 0, linear) state, so (g_x - g_y)_circular = -(1/2)
    (g_x - g_y)_linear and SS90's a = -2 (g_x - g_y)_circular is the linear
    difference DS16's deltaG is (KS20 Eq. (68))."""
    eps = np.zeros((3, 3, 3))
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        eps[i, j, k], eps[i, k, j] = 1.0, -1.0
    s = [-1j * eps[i] for i in range(3)]        # (S_i)_{ab} = -i eps_{iab}

    def t_xx_minus_t_yy(v):
        v = v / np.linalg.norm(v)
        tens = [[np.real(np.conj(v) @ (0.5 * (s[i] @ s[j] + s[j] @ s[i])) @ v)
                 for j in range(3)] for i in range(3)]
        return tens[0][0] - tens[1][1]

    circular = -(np.array([0, 1.0, 0]) + 1j * np.array([0, 0, 1.0])) / np.sqrt(2)
    linear = np.array([1.0, 0, 0], dtype=complex)
    assert t_xx_minus_t_yy(circular) == pytest.approx(0.5, abs=1e-15)
    assert t_xx_minus_t_yy(linear) == pytest.approx(-1.0, abs=1e-15)


# ------------------------------------------------------------ the measurement

def test_row_passes_and_the_headline_is_pinned(rep):
    assert rep["status"] == "pass"
    toy = _sub(rep, "toy_1e-2[shipped]")
    assert toy["status"] == "pass" and toy["in_verdict"] is True
    assert (toy["status_central"], toy["status_v6"]) == ("pass", "pass")
    assert toy["q2"] == [2.0, 5.0, 10.0]
    assert toy["m2"] == pytest.approx([1.238347186e-4, 1.310790196e-4, 1.380460159e-4],
                                      rel=1e-7)
    assert toy["bound"] == pytest.approx([1.921562687e-3, 1.633121332e-3, 1.467657982e-3],
                                         rel=1e-7)
    assert toy["saturation"] == pytest.approx([0.06444479771, 0.08026287884, 0.09405870963],
                                              rel=1e-7)
    assert toy["saturation_v6"] == pytest.approx([0.1288895954, 0.1605257577,
                                                  0.1881174193], rel=1e-7)
    assert toy["saturation_tight"] == pytest.approx([0.2719227491, 0.3369110821,
                                                     0.3936582914], rel=1e-7)
    assert toy["saturation_loose"] == pytest.approx([0.02601242937, 0.03237320575,
                                                     0.03792158164], rel=1e-7)
    assert toy["f1ll_critical"] == pytest.approx([-1.923123096, -1.886331427,
                                                  -1.855643673], rel=1e-7)
    assert toy["f1ll_critical_v6"] == pytest.approx([-1.790650658, -1.721716574,
                                                     -1.66298279], rel=1e-7)
    assert all(e < 1e-15 for e in toy["m2_numerical_error"])
    assert all(t < 1e-12 for t in toy["m2_tail_x_lo"])


def test_the_shipped_delta_is_the_pipelines_toy_at_1e_2(rep):
    """The scale is READ from the kernel: the pipeline's 1e-2, not the
    function's 1e-3 default."""
    assert _sub(rep, "toy_1e-2[shipped]")["scale_measured"] == pytest.approx(1e-2, rel=1e-14)
    assert rep["checks"]["max_rel_dev_from_shipped_scale"] < 1e-12


def test_toy_moment_by_an_independent_quadrature(rep):
    """The toy's M2 again, with scipy's adaptive quad in ln x on the kernel's
    own Delta (not the harness's Gauss-Legendre rule)."""
    integrate = pytest.importorskip("scipy.integrate")
    k = lg.default_inclusive_kernel(lg.li6())
    toy = _sub(rep, "toy_1e-2[shipped]")
    for q2, m2 in zip(toy["q2"], toy["m2"]):
        def f(u, q2=q2):
            x = math.exp(u)
            return x * x * k.tables(x, q2).delta
        lo, _ = integrate.quad(f, math.log(1e-9), math.log(0.1), limit=200,
                               epsabs=0.0, epsrel=1e-10)
        hi, _ = integrate.quad(lambda x, q2=q2: x * k.tables(x, q2).delta, 0.1, 1.0,
                               limit=200, epsabs=0.0, epsrel=1e-10)
        assert lo + hi == pytest.approx(m2, rel=1e-8)


def test_bound_pieces_and_closed_forms(rep, mod):
    ion = lg.li6()
    for b in rep["bound_detail"]:
        assert b["p_eff"] == pytest.approx(lg.LI6_CLUSTER_POLARIZATION / 3.0, rel=1e-14)
        assert b["p_eff"] == pytest.approx(
            (ion.Z * ion.eff_pol_p + ion.N * ion.eff_pol_n) / ion.A, rel=1e-15)
        k = b["alpha_s"] / (8.0 * math.pi)
        assert b["B"] == pytest.approx(k * (2.0 / 3.0) * 0.5 * (b["xg"] + b["g1_moment"]),
                                       rel=1e-14)
        assert b["g1_moment"] == pytest.approx(b["p_eff"] * b["xdg"], rel=1e-14)
        assert b["B_dg0"] == pytest.approx(k * (2.0 / 3.0) * 0.5 * b["xg"], rel=1e-14)
        assert b["B_ss90"] == pytest.approx(2.0 * b["B"], rel=1e-14)
        assert b["B_nf4"] == pytest.approx(b["B"] * (10.0 / 9.0) / (2.0 / 3.0), rel=1e-14)
        assert b["B_ll_minus"] == pytest.approx(
            k * (2.0 / 3.0) * 0.5 * (0.5 * b["xg"] + b["g1_moment"]), rel=1e-14)
        assert b["B_peff_naive"] == pytest.approx(
            k * (2.0 / 3.0) * 0.5 * (b["xg"] + b["xdg"] / 3.0), rel=1e-14)
        # the edges (SECOND PASS S3): V4 in neither, V6 on the model side
        xg_lo, xg_hi = min(b["xg"], b["xg_epps21"]), max(b["xg"], b["xg_epps21"])
        assert b["B_tight"] == pytest.approx(k * (2.0 / 3.0) * 0.5 * 0.5 * xg_lo,
                                             rel=1e-14)
        assert b["B_loose"] == pytest.approx(
            k * (10.0 / 9.0) * 0.5 * (1.5 * xg_hi + max(b["g1_moment"], 0.0)), rel=1e-14)
        assert b["B_tight"] < b["B_dg0"] < b["B"] < b["B_loose"]
        assert b["xg_err"] < 1e-8 and b["xdg_err"] < 1e-8
        # |g1^g| <= f1^g / 2 fails only at x > 0.94, carrying < 3e-7 of <x>_g
        assert 0.94 < b["x_g1_above_half_f1"] < b["x_g1_above_f1"] < 0.98
        assert b["xg_share_above_x_half"] < 3e-7
    a_s = [b["alpha_s"] for b in rep["bound_detail"]]
    xg = [b["xg"] for b in rep["bound_detail"]]
    xdg = [b["xdg"] for b in rep["bound_detail"]]
    assert a_s == pytest.approx([0.354535328227, 0.283892922592, 0.247419816763], rel=1e-9)
    assert xg == pytest.approx([0.397601929906, 0.422961376441, 0.43670378769], rel=1e-8)
    assert xdg == pytest.approx([0.0408730005353, 0.0398438516164, 0.0390052607519], rel=1e-7)
    assert [b["xg_epps21"] for b in rep["bound_detail"]] == pytest.approx(
        [0.387398975277, 0.413318103314, 0.427455686421], rel=1e-7)
    assert [b["B_loose"] for b in rep["bound_detail"]] == pytest.approx(
        [4.76059797571e-3, 4.04899720355e-3, 3.64030216937e-3], rel=1e-7)
    assert [b["B_peff_naive"] / b["B"] - 1.0 for b in rep["bound_detail"]] == pytest.approx(
        [6.293e-3, 5.780e-3, 5.488e-3], abs=5e-6)


def test_bound_side_subrow_carries_the_pieces(rep):
    s = _sub(rep, "bound_side[central]")
    assert s["status"] == "recorded" and s["in_verdict"] is False
    assert s["tr_q2"] == pytest.approx(2.0 / 3.0, rel=1e-15)
    assert (s["ds16_over_tree"], s["pre_d4_over_jm"]) == (1.0, 2.0)  # D4, 2026-10-27 (2.0, -)
    for k in ("alpha_s", "xg", "xdg", "B", "B_tight", "B_loose", "B_peff_naive"):
        assert s[k] == [b[k] for b in rep["bound_detail"]]
    assert s["B"] == pytest.approx(_sub(rep, "toy_1e-2[shipped]")["bound"], rel=1e-15)


def test_the_plugin_reads_the_trees_grids(rep):
    g = rep["checks"]["same_grid"]
    assert g["points"] == 6
    assert g["worst_rel"]["f2p"] <= 1e-12 and g["worst_rel"]["g1p"] <= 1e-12


def test_letter_models_recorded_with_their_own_status(rep):
    """moment_A / moment_B at the letter's dilution 1/3 use ~70 % of the room
    under the central assumptions (A7, the kernel's own Delta since D4) and ~140 % under V6
    (the pre-D4 kernel / the letter's printed amplitudes, x2); undiluted they exceed it under
    both -- recorded, not in the verdict."""
    a = _sub(rep, "moment_A[letter]")
    b = _sub(rep, "moment_B[letter]")
    for s in (a, b):
        assert (s["status"], s["in_verdict"]) == ("recorded", False)
        assert (s["status_central"], s["status_v6"]) == ("pass", "fail")
    assert a["a_solved"] == pytest.approx(-0.29239892903, rel=1e-9)
    assert a["a_reproduced"] is True and a["dilution"] == "1/3"
    assert a["saturation"] == pytest.approx([0.7146852971, 0.6987238075, 0.7020124259],
                                            rel=1e-7)
    assert b["saturation"] == pytest.approx([0.7380145976, 0.6953382265, 0.6743255442],
                                            rel=1e-7)
    assert a["saturation_v6"] == pytest.approx([1.429370594, 1.397447615, 1.404024852],
                                               rel=1e-7)
    assert b["saturation_v6"] == pytest.approx([1.476029195, 1.390676453, 1.348651088],
                                               rel=1e-7)
    assert a["saturation_tight"] == pytest.approx([3.015591601, 2.932959762, 2.938090616],
                                                  rel=1e-7)
    assert a["saturation_loose"] == pytest.approx([0.2884748106, 0.2818230533,
                                                   0.2830298399], rel=1e-7)
    assert a["f1ll_critical"] == pytest.approx([-0.5864916293, -0.6179012863,
                                                -0.6103693059], rel=1e-7)
    # under V6 the scenario needs the 6Li gluons tensor-polarized at
    # f1LL^g >~ +0.8 f1^g, short of CvDM17 (35)'s f1LL <= f1
    assert a["f1ll_critical_v6"] == pytest.approx([0.8826122761, 0.8151437073,
                                                   0.827565945], rel=1e-7)
    assert b["f1ll_critical_v6"] == pytest.approx([0.978523488, 0.801256418,
                                                   0.714143613], rel=1e-7)
    assert a["saturation_model_lo_alpha_s"] == pytest.approx(
        [0.8168452313, 0.8002764808, 0.80261216], rel=1e-7)
    # moment_B's moment is dilution * C_BAG * alpha_s in closed form
    for m2, bd in zip(b["m2"], rep["bound_detail"]):
        assert m2 == pytest.approx(lg.C_BAG * bd["alpha_s"] / 3.0, rel=1e-12)
    for name in ("moment_A[dilution 1]", "moment_B[dilution 1]"):
        s = _sub(rep, name)
        assert (s["status"], s["in_verdict"]) == ("recorded", False)
        assert (s["status_central"], s["status_v6"]) == ("fail", "fail")
        assert all(v > 2.0 for v in s["saturation"])           # excluded centrally
        assert all(v > 1.0 for v in s["saturation_ss90"])      # even with V4's x2
        assert all(0.8 < v < 0.9 for v in s["saturation_loose"])   # allowed at the loose edge
        assert all(v > 2.0 for v in s["f1ll_critical"])        # above f1LL <= f1 (CvDM17 (35))
    assert _sub(rep, "moment_A[dilution 1]")["saturation"] == pytest.approx(
        [2.144055891, 2.096171422, 2.106037278], rel=1e-7)


def test_helicity_positivity_subrow(rep):
    p = _sub(rep, "helicity_positivity[shipped kernel]")
    assert p["status"] == "pass" and p["in_verdict"] is True
    assert (p["n_points"], p["n_rates"]) == (567, 1701)
    assert p["w_min"] == pytest.approx(0.9493997990, rel=1e-9)  # D4, 2026-10-27 (0.9493997515)
    assert p["a_t_min"] == pytest.approx(-0.1012003069, rel=1e-8)  # theta_S = 0: no Delta term
    assert p["a_t_max"] == pytest.approx(0.05060020096, rel=1e-8)  # D4 (0.05060024846 before)
    assert p["w_min_at"] == "x=0.9 Q2=2 y=0.95 theta_S=90 phi'=90 m=+0"
    assert -2.0 <= p["a_t_min"] <= p["a_t_max"] <= 1.0


def test_verdict_follows_the_numbers(mod, rep):
    toy = _sub(rep, "toy_1e-2[shipped]")
    for s, s6, m, b in zip(toy["saturation"], toy["saturation_v6"], toy["m2"],
                           toy["bound"]):
        assert s == pytest.approx(abs(m) / b, rel=1e-15)
        assert s6 == mod.PRE_D4_OVER_JM * s      # V6 since D4 (2026-10-27)
    bound = rep["bound_detail"]
    # synthetic models at and just beyond the bound: the -1e-12 margin rule
    at = [dict(m2=b["B"], err=0.0, tail=0.0) for b in bound]
    row = mod.saturation_row("at", at, bound, True)
    assert (row["status"], row["status_v6"]) == ("pass", "fail")
    over = [dict(m2=-b["B"] * (1.0 + 1e-9), err=0.0, tail=0.0) for b in bound]
    assert mod.saturation_row("over", over, bound, True)["status"] == "fail"
    one_over = at[:2] + over[2:]
    assert mod.saturation_row("one", one_over, bound, True)["status"] == "fail"
    half = [dict(m2=0.5 * b["B"], err=0.0, tail=0.0) for b in bound]
    assert mod.saturation_row("half", half, bound, False)["status_v6"] == "pass"
    assert mod.saturation_row("half", half, bound, False)["status"] == "recorded"
    # f1ll_critical is where the saturation reaches 1, f1ll_critical_v6 where
    # twice it does (PRE_D4_OVER_JM, V6 since D4)
    for m, b, ll, ll6 in zip(toy["m2"], bound, toy["f1ll_critical"],
                             toy["f1ll_critical_v6"]):
        k = b["alpha_s"] / (8.0 * math.pi) * (2.0 / 3.0) * 0.5
        assert k * (b["xg"] * (1.0 + 0.5 * ll) + b["g1_moment"]) == pytest.approx(
            abs(m), rel=1e-12)
        assert k * (b["xg"] * (1.0 + 0.5 * ll6) + b["g1_moment"]) == pytest.approx(
            2.0 * abs(m), rel=1e-12)


def test_a_set_of_another_version_is_a_blocked_row(mod, monkeypatch, capsys):
    """A set present in another version (or without its member file) is a
    BLOCKED row with an environment: reason, exit 0 -- the row ran and says
    why it cannot measure."""
    def other_version(entry):
        raise mod.EnvironmentBlocked("environment: LHAPDF set %s at /x is not the "
                                     "version this row was measured with (test)"
                                     % entry["name"])
    monkeypatch.setattr(mod, "check_set", other_version)
    r = mod.run(verbose=True)
    assert r["status"] == "blocked" and r["reason"].startswith("environment: LHAPDF set "
                                                                "CT18NLO")
    out = capsys.readouterr().out.splitlines()
    assert [l for l in out if l.startswith("REPORT |")][0].endswith("| BLOCKED")


def test_a_missing_set_is_blocked_by_the_environment(mod, monkeypatch, tmp_path):
    """A set that is not installed where LHAPDF searches: the row's only
    configuration needs it, so run() raises LHAPDF's own message, main()
    exits 2 and make_report's classifier records 'blocked (environment)' with
    the convention's reason -- README convention 6, as row 4 (phase F of run
    2026-10-26; a BLOCKED row with exit 0 before).  A trailing '::' keeps
    LHAPDF out of its install datadir, as LHAPDF::paths() does."""
    if not getattr(lg, "HAVE_LHAPDF", False):
        pytest.skip("environment: lipolgen built without the LHAPDF tier")
    monkeypatch.setenv("LHAPDF_DATA_PATH", str(tmp_path) + "::")
    monkeypatch.delenv("LHAPATH", raising=False)
    assert mod.lhapdf_dirs() == [str(tmp_path)]
    with pytest.raises(RuntimeError, match=r"Info file not found for PDF set 'CT18NLO'"
                       ) as exc:
        mod.run(verbose=False)
    rp = _load("_report")
    assert rp.classify_exception(exc.value) == (
        rp.ENV_BLOCKED, "environment: LHAPDF set CT18NLO not installed")
    assert mod.main([]) == 2


def test_lhapdf_dirs_follow_lhapdf_search_order(mod, monkeypatch, tmp_path):
    """LHAPDF::paths(): $LHAPDF_DATA_PATH (else $LHAPATH) split on ':', empty
    entries dropped, then the install datadir unless the variable ends in
    '::'.  Until phase F the datadir was appended always."""
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    monkeypatch.delenv("LHAPATH", raising=False)
    monkeypatch.setenv("LHAPDF_DATA_PATH", "%s:%s::" % (a, b))
    assert mod.lhapdf_dirs() == [a, b]
    monkeypatch.setenv("LHAPDF_DATA_PATH", "%s:%s" % (a, b))
    d = mod.lhapdf_dirs()
    assert d[:2] == [a, b] and len(d) <= 3          # + lhapdf-config's datadir, if any
    monkeypatch.delenv("LHAPDF_DATA_PATH")
    monkeypatch.setenv("LHAPATH", a + "::")
    assert mod.lhapdf_dirs() == [a]


def test_report_line(mod, capsys):
    r = _run_or_skip(mod, verbose=True)
    out = capsys.readouterr().out.splitlines()
    rows = [l for l in out if l.startswith("REPORT |")]
    subs = [l for l in out if l.startswith("SUBROW |")]
    assert len(rows) == 1 and len(subs) == 7
    assert subs[-1].startswith("SUBROW | t5_delta_moment_bound[bound_side[central]] | "
                               "B = 1.921563e-03 / 1.633121e-03 / 1.467658e-03")
    assert "PASS centrally (A7" in subs[2] and "FAIL under V6 (the pre-D4 kernel" in subs[2]  # D4
    assert rows[0].startswith("REPORT | t5_delta_moment_bound | shipped toy_delta_gluon (scale 0.01) "
                              "per nucleon, 6Li: M2 = int x Delta dx = 1.2383e-04 / 1.3108e-04 / "
                              "1.3805e-04 at Q2 = 2 / 5 / 10 GeV2, |M2|/B = 0.0644 / 0.0803 / 0.0941 "
                              "(Jaffe-Manohar Delta, the kernel's since D4; V6, the pre-D4 kernel / "
                              "the letter as printed, x2: 0.1289 / 0.1605 / 0.1881);")
    assert rows[0].endswith("| PASS")
    assert "np.float64" not in "\n".join(out)
    assert r["runtime_s"] < 30.0


def test_exit_codes(tmp_path):
    ok = _run_harness(os.path.join(_BENCH, _NAME + ".py"))
    if _MISSING_SET in ok.stderr:      # a set not installed: exit 2, no row
        assert ok.returncode == 2 and "REPORT |" not in ok.stdout
    else:
        assert ok.returncode == 0, ok.stderr
        assert "REPORT | t5_delta_moment_bound |" in ok.stdout
    d = tmp_path / "bench"
    (d / "data").mkdir(parents=True)
    shutil.copy(os.path.join(_BENCH, _NAME + ".py"), d / (_NAME + ".py"))
    raw = open(os.path.join(_BENCH, "data", _DATA), "rb").read()
    old = b'"A_as_printed": -0.2924'
    assert raw.count(old) == 1
    (d / "data" / _DATA).write_bytes(raw.replace(old, b'"A_as_printed": -0.2925'))
    bad = _run_harness(str(d / (_NAME + ".py")))
    assert bad.returncode == 2 and "sha256" in bad.stderr
    assert "REPORT |" not in bad.stdout


def test_plot_writes_a_png(mod, rep, tmp_path):
    pytest.importorskip("matplotlib")
    path = mod.plot(rep, out_dir=str(tmp_path))
    assert path == str(tmp_path / (_NAME + ".png"))
    raw = open(path, "rb").read()
    assert len(raw) > 10000 and raw[:8] == b"\x89PNG\r\n\x1a\n"
