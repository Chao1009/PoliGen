# SPDX-License-Identifier: GPL-3.0-or-later
"""The two STATIC-OBSERVABLE harnesses under validation/benchmarks/,
exercised as tests (scout items 5 and 6; both cheap, < 2 s together, so
neither sits behind LIPOLGEN_BENCH):

  t2_deuteron_static.py   eta_d and Q_d vs the tree's three deuteron waves
  t3_li6_rms_devries.py   6Li charge rms of HoSpin1FF vs de Vries Table I

A harness that DISAGREES with the tree is not a failing test: the tests pin
what each harness MEASURES and that its verdict follows from the numbers.
No test here moves, or asserts a new value for, any tree default -- the
Hulthen deuteron's eta is a recorded FAIL, pinned as one.
"""
import importlib.util
import math
import os
import shutil
import subprocess
import sys

import pytest

lg = pytest.importorskip("lipolgen")

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")


def _load(name):
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


def _run_harness(script, cwd=None):
    """The harness as a script, with the imported lipolgen on its path."""
    env = dict(os.environ)
    pkg_parent = os.path.dirname(os.path.dirname(os.path.abspath(lg.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [pkg_parent] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    return subprocess.run([sys.executable, script], capture_output=True,
                          text=True, env=env, cwd=cwd, timeout=300)


def _exit_code_on_altered_data(tmp_path, name, data_name, old, new):
    """Copy the harness next to a data/ holding a one-byte-altered copy of
    its vendored file: the harness must refuse it and exit 2."""
    d = tmp_path / "bench"
    (d / "data").mkdir(parents=True)
    shutil.copy(os.path.join(_BENCH, name + ".py"), d / (name + ".py"))
    src = os.path.join(_BENCH, "data", data_name)
    raw = open(src, "rb").read()
    assert raw.count(old) == 1
    (d / "data" / data_name).write_bytes(raw.replace(old, new))
    return _run_harness(str(d / (name + ".py")))


# ============================================================ t2_deuteron_static

@pytest.fixture(scope="module")
def dmod():
    return _load("t2_deuteron_static")


@pytest.fixture(scope="module")
def ds(dmod):
    fdeut = lg.data_path(lg.VMC_DEUTERON_WAVE)
    if not os.path.isfile(fdeut):
        pytest.skip("environment: the AV18 deuteron table %s is not at %s "
                    "(LIPOLGEN_DATA_DIR)" % (lg.VMC_DEUTERON_WAVE, fdeut))
    return dmod.run(verbose=False)


def test_deuteron_static_vendored_file_is_the_transcription(dmod, tmp_path):
    doc = dmod.read_vendored()
    pts = doc["body"]["points"]
    assert len(pts) == 2
    assert (pts[0]["quantity"], pts[0]["value"], pts[0]["sigma"]) == ("eta_d", 0.0256, 0.0004)
    assert (pts[-1]["quantity"], pts[-1]["value"], pts[-1]["sigma"]) == ("Q_d", 0.2859, 0.0003)
    prov = doc["_provenance"]
    assert prov["transcription"].startswith("SECOND-HAND")
    assert "NEITHER SOURCE PAPER WAS RE-READ" in prov["transcription"]
    assert prov["fetched"].startswith("2026-09-26")
    ref = dmod.load_reference()
    assert ref["eta_d"]["as_printed"] == "0.0256(4)"
    # one flipped byte in the body is refused ...
    bad = _flip(dmod.DATA, b'"value": 0.0256,', b'"value": 0.0257,', tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        dmod.load_reference(bad)
    # ... and so is one in the provenance (the whole file is digested)
    bad = _flip(dmod.DATA, b"SECOND-HAND", b"SECOND-HANd", tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        dmod.read_vendored(bad)


def test_deuteron_static_rule3_per_deuteron_normalisation(dmod, ds):
    """BENCHMARK_PLAN.md sec. 8 rule 3, ASSERTED: Q_d is the per-deuteron
    one-body operator on a wave function normalised to 1 -- invariant under
    u, w -> c u, c w; zero for a pure S wave; and the SAME operator as the
    tree's own `alpha_d_quadrupole_fm2` (at converged numerics)."""
    r, h = dmod.r_grid(40.0, 4000)
    u = r * dmod.np.exp(-0.3 * r)
    w = 0.05 * r ** 3 * dmod.np.exp(-0.35 * r)
    q1 = dmod.quadrupole_ia(r, u, w, h)[0]
    q3 = dmod.quadrupole_ia(r, 3.0 * u, 3.0 * w, h)[0]
    assert q3 == pytest.approx(q1, rel=1e-13)
    assert dmod.quadrupole_ia(r, u, 0.0 * w, h)[0] == 0.0
    m = ds["models"]
    assert abs(m["cdbonn"]["q_d"] - m["cdbonn"]["cross_checks"]["q_tree_converged"]) < 5e-7
    # every model is divided by its own norm, measured and pinned
    assert m["cdbonn"]["norm"] == pytest.approx(0.99999983, abs=1e-7)
    assert m["av18"]["norm"] == pytest.approx(0.99999813, abs=1e-7)
    assert m["hulthen"]["norm"] == pytest.approx(1.000144, abs=1e-6)
    # eta is a ratio: the Hulthen value does not depend on the overall scale
    assert m["hulthen"]["eta"] == pytest.approx(
        -m["hulthen"]["n2"] * dmod.hulthen_forms(
            lg.deuteron_channel().base.kappa(), 0.30)[2] / m["hulthen"]["n0"],
        rel=1e-14)


def test_deuteron_static_eta_values_and_cross_checks(ds):
    m = ds["models"]
    c, a, h = m["cdbonn"], m["av18"], m["hulthen"]
    # CD-Bonn: the tree's own accessor, as tests/test_cluster.cpp pins it
    assert c["eta"] == pytest.approx(0.0255713786530431, rel=1e-13)
    assert c["eta_check_r60"] == pytest.approx(c["eta"], rel=1e-12)
    # AV18: from the file's r block, against the file's own header
    assert a["eta"] == pytest.approx(0.02504467, abs=5e-8)
    assert abs(a["eta"] - a["cross_checks"]["header_eta"]) <= 5e-7
    assert a["eta_numerical_error"] < 1e-6
    assert a["a_s"] == pytest.approx(a["cross_checks"]["header_as"], abs=2e-6)
    assert a["kappa_fm"] == pytest.approx(0.2316065, abs=1e-6)
    # the k block the generator reads, through the tree's spline, at 20 fm
    assert a["cross_checks"]["eta_kblock_r"] == pytest.approx(0.024984, abs=2e-6)
    assert abs(a["cross_checks"]["eta_kblock_r"] / a["eta"] - 1) < 5e-3
    # Hulthen: pole residues, and the large-r ratio of the closed forms
    assert h["eta"] == pytest.approx(0.03247843, rel=1e-6)
    assert h["eta_check_r80"] == pytest.approx(h["eta"], rel=1e-8)
    assert h["beta_gev"] == 0.30 and h["p_d_nominal"] == lg.P_D_DEUTERON


def test_deuteron_static_qd_values_and_in_tree_cross_checks(ds):
    m = ds["models"]
    c, a, h = m["cdbonn"], m["av18"], m["hulthen"]
    assert c["q_d"] == pytest.approx(0.2704934, abs=2e-7)
    # tests/test_b1_nuclear.cpp item 5 prints/pins these two
    assert c["cross_checks"]["q_tree_default"] == pytest.approx(0.270178, abs=2e-6)
    assert a["cross_checks"]["q_tree_kblock"] == pytest.approx(0.269362, abs=2e-6)
    assert a["q_d"] == pytest.approx(0.2696703, abs=2e-7)
    assert a["q_d"] == pytest.approx(a["cross_checks"]["header_qm"], abs=5e-6)
    assert a["p_d"] == pytest.approx(a["cross_checks"]["header_dstate"], abs=2e-6)
    assert h["q_d"] == pytest.approx(0.290349, abs=1e-6)
    assert h["cross_checks"]["q_tree_tables"] == pytest.approx(0.289791, abs=2e-6)
    assert abs(h["q_d"] / h["cross_checks"]["q_tree_tables"] - 1) < 3e-3
    assert max(x["q_d_quad_error"] for x in m.values()) < 1e-8
    # the recorded distances (impulse approximation: about -5 %)
    assert c["q_d_rel_distance"] == pytest.approx(-0.05389, abs=1e-4)
    assert a["q_d_rel_distance"] == pytest.approx(-0.05677, abs=1e-4)
    assert h["q_d_rel_distance"] == pytest.approx(+0.01556, abs=1e-4)


def test_deuteron_static_verdict_follows_the_numbers(dmod, ds):
    ref = ds["reference_points"]
    e = ref["eta_d"]
    for name, m in ds["models"].items():
        assert m["eta_pull"] == pytest.approx((m["eta"] - e["value"]) / e["sigma"],
                                              rel=1e-14)
        assert m["eta_status"] == ("pass" if abs(m["eta_pull"]) <= 2.0 else "fail")
    assert ds["status"] == ("pass" if all(m["eta_status"] == "pass"
                                          for m in ds["models"].values()) else "fail")
    # measured: CD-Bonn and AV18 inside, the DEFAULT Hulthen pair far outside
    assert ds["models"]["cdbonn"]["eta_pull"] == pytest.approx(-0.0715, abs=1e-3)
    assert ds["models"]["av18"]["eta_pull"] == pytest.approx(-1.388, abs=1e-3)
    assert ds["models"]["hulthen"]["eta_pull"] == pytest.approx(17.196, abs=1e-3)
    assert ds["status"] == "fail" and ds["failing"] == ["hulthen"]
    # Q_d never enters the verdict: a reference Q_d 100x off changes nothing
    fake = {"eta_d": dict(e), "Q_d": dict(ref["Q_d"], value=28.59)}
    models = {k: dict(v) for k, v in ds["models"].items()}
    assert dmod.verdict(models, fake) == ds["status"]
    # and with every eta inside, the row would pass
    for v in models.values():
        v["eta"] = e["value"] + 1.9 * e["sigma"]
    assert dmod.verdict(models, fake) == "pass"
    # the beta band is recorded, never in the verdict
    band = ds["hulthen_beta_band"]
    assert band[0.20]["eta"] == pytest.approx(0.05940, abs=5e-5)
    assert band[0.40]["eta"] == pytest.approx(0.02143, abs=5e-5)


def test_deuteron_static_hulthen_closed_form_is_the_transform(dmod):
    """The partial-fraction r-space forms against a direct numerical
    bare-j_L transform of the momentum-space pair (QAWF for the tail)."""
    pytest.importorskip("scipy", reason="scipy (the numerical Fourier-Bessel "
                                        "check of the Hulthen closed forms)")
    from scipy import integrate
    from scipy.special import spherical_jn
    ch = lg.deuteron_channel()
    kap = ch.base.kappa()
    # beta read from the channel, not a literal: the closed forms take one
    # beta for both waves, and the harness refuses a pair that differs
    beta = ch.waves[1].beta
    assert ch.waves[0].beta == beta == 0.30
    u_f, w_f, _ = dmod.hulthen_forms(kap, beta)
    np = dmod.np

    def transform(l, psi, r, a=2.0):
        # r sqrt(2/pi) INT k^2 j_l(kr) psi(k) dk; k^2 j_l(kr) r in sin/cos parts
        reg = integrate.quad(lambda k: k * k * spherical_jn(l, k * r) * r * psi(k),
                             0, a, limit=2000, epsabs=1e-16)[0]
        if l == 0:
            s = integrate.quad(lambda k: k * psi(k), a, np.inf, weight="sin", wvar=r)[0]
            c = 0.0
        else:
            s = integrate.quad(lambda k: (3 / (k * r * r) - k) * psi(k), a, np.inf,
                               weight="sin", wvar=r)[0]
            c = integrate.quad(lambda k: -(3 / r) * psi(k), a, np.inf,
                               weight="cos", wvar=r)[0]
        return math.sqrt(2 / math.pi) * (reg + s + c)

    psi0 = lambda k: 1 / (k * k + kap * kap) - 1 / (k * k + beta * beta)  # noqa: E731
    psi2 = lambda k: k * k / ((k * k + kap * kap) * (k * k + beta * beta) ** 2)  # noqa: E731
    # the pair transformed here IS the channel's Wave.radial, both waves
    for k in (0.05, 0.2, 0.7):
        assert ch.waves[0].radial(k, kap) == pytest.approx(psi0(k), rel=1e-13)
        assert ch.waves[1].radial(k, kap) == pytest.approx(psi2(k), rel=1e-13)
    for r_fm in (0.1, 1.0, 3.0, 10.0):
        r = r_fm / lg.HBARC_GEV_FM
        assert u_f(np.array([r]))[0] == pytest.approx(transform(0, psi0, r), rel=1e-8)
        assert w_f(np.array([r]))[0] == pytest.approx(transform(2, psi2, r), rel=1e-8)


def test_deuteron_static_hulthen_refuses_a_d_wave_with_its_own_beta(dmod,
                                                                     monkeypatch):
    """hulthen()'s closed forms and eta = -N_2 A/N_0 take the S wave's beta
    for the D wave too: a channel whose D wave carries another beta is
    refused, not silently computed with the wrong one."""
    real = lg.deuteron_channel()

    class _Wave:
        def __init__(self, w, beta):
            self._w, self.beta = w, beta

        def __getattr__(self, name):
            return getattr(self._w, name)

    class _Channel:
        base = real.base
        waves = [real.waves[0], _Wave(real.waves[1], 0.35)]

    monkeypatch.setattr(lg, "deuteron_channel", lambda **kw: _Channel())
    with pytest.raises(RuntimeError, match="one beta for both"):
        dmod.hulthen()


def test_deuteron_static_report_line(dmod, capsys):
    rep = dmod.run(verbose=True)
    out = capsys.readouterr().out.splitlines()
    rows = [l for l in out if l.startswith("REPORT |")]
    subs = [l for l in out if l.startswith("SUBROW |")]
    assert len(rows) == 1 and len(subs) == 3
    assert rows[0].startswith("REPORT | t2_deuteron_static | eta: cdbonn 0.025571 "
                              "(-0.07 sigma), av18 0.025045 (-1.39 sigma), hulthen "
                              "0.032478 (+17.20 sigma)")
    assert rows[0].endswith("| FAIL (recorded: hulthen +17.2 sigma outside; "
                            "nothing moved)")
    assert "np.float64" not in "\n".join(out)
    assert rep["runtime_s"] < 30.0


def test_deuteron_static_exit_codes(tmp_path):
    ok = _run_harness(os.path.join(_BENCH, "t2_deuteron_static.py"))
    assert ok.returncode == 0, ok.stderr
    assert "REPORT | t2_deuteron_static |" in ok.stdout
    bad = _exit_code_on_altered_data(tmp_path, "t2_deuteron_static",
                                     "deuteron_static_measured.json",
                                     b'"value": 0.2859,', b'"value": 0.2869,')
    assert bad.returncode == 2 and "sha256" in bad.stderr
    assert "REPORT |" not in bad.stdout


# ============================================================ t3_li6_rms_devries

@pytest.fixture(scope="module")
def rmod():
    return _load("t3_li6_rms_devries")


@pytest.fixture(scope="module")
def rms(rmod):
    return rmod.run(verbose=False)


def test_devries_vendored_file_is_the_survey_transcription(rmod, tmp_path):
    doc = rmod.read_vendored()
    pts = doc["body"]["points"]
    assert len(pts) == 4
    assert (pts[0]["rms_fm"], pts[0]["sigma_fm"], pts[0]["source_key"]) == (2.54, 0.05, "Su67")
    assert [p["rms_fm"] for p in pts] == [2.54, 2.56, 2.57, 2.589]
    assert (pts[-1]["rms_fm"], pts[-1]["sigma_fm"]) == (2.589, 0.039)
    assert pts[-1]["compilation"].startswith("Angeli")
    prov = doc["_provenance"]
    assert "SURVEY TRANSCRIPTION" in prov["transcription"]
    assert "NOT RE-READ" in prov["transcription"]
    dv, am = rmod.load_reference()
    assert len(dv) == 3 and am["as_printed"] == "2.589(39)"
    bad = _flip(rmod.DATA, b'"rms_fm": 2.56,', b'"rms_fm": 2.66,', tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        rmod.load_reference(bad)


def test_devries_rule3_per_nucleus_and_isoscalar_folding(rmod, rms):
    """Rule 3: F_c(0) = Z = 3 and the slope is of F/F(0), so Z cancels.
    Rule 2's analogue: the fold is ISOSCALAR, r^2 = r^2_pt + r^2_p + r^2_n."""
    t = rms["tree"]
    assert t["fc0"] == 3.0 == t["z"]
    ff = lg.HoSpin1FF.for_ion(lg.li6())
    hc2 = lg.HBARC_GEV_FM ** 2
    r2_twice, _, _ = rmod.richardson_r2(lambda s: 2.0 * ff.fc(s * hc2))
    assert r2_twice == pytest.approx(t["r2"], rel=1e-13)
    assert abs(t["decomposition_residual"]) < 1e-9
    assert t["r2_n"] < 0.0
    assert t["r_proton_only_fold"] == pytest.approx(2.595567, abs=1e-6)
    # the nucleon slopes are rc.hpp's dipole and Galster (MU_NEUTRON =
    # -1.9130427 is not bound; typed here as a pin of rc.hpp's value)
    assert t["r2_p"] == pytest.approx(12.0 / 0.71 * hc2, rel=1e-8)
    assert t["r2_n"] == pytest.approx(
        -6.0 * 1.9130427 / (4.0 * lg.PROTON_MASS ** 2) * hc2, rel=1e-8)
    # the point radius IS the named constant, up to (a, alpha)'s rounding
    assert t["r2_point"] == pytest.approx(lg.LI6_R2_POINT_FM2, abs=1e-4)


def test_devries_slope_numerics(rms):
    t = rms["tree"]
    assert t["r2_error"] < 1e-8
    assert abs(t["closed_form_residual"]) < 1e-10
    assert t["r_one_sided"] == pytest.approx(t["r"], abs=1e-6)
    d = t["richardson_diag"]
    assert all(abs(d[i + 1] - d[-1]) < abs(d[i] - d[-1]) for i in range(len(d) - 1))


def test_devries_radius_pinned(rms):
    t = rms["tree"]
    assert t["r"] == pytest.approx(2.5710008, abs=2e-7)
    assert math.sqrt(t["r2_point"]) == pytest.approx(2.465531, abs=1e-6)
    assert t["r_vmcft_edge"] == pytest.approx(2.570992, abs=1e-5)


def test_devries_verdict_follows_the_numbers(rmod, rms):
    r = rms["tree"]["r"]
    for p in rms["de_vries"]:
        assert p["pull"] == pytest.approx((r - p["rms_fm"]) / p["sigma_fm"], rel=1e-14)
        assert p["within"] == (abs(p["pull"]) <= 2.0)
    assert rms["status"] == ("pass" if all(p["within"] for p in rms["de_vries"])
                             else "fail")
    assert rms["status"] == "pass"
    assert [round(p["pull"], 2) for p in rms["de_vries"]] == [0.62, 0.22, 0.01]
    assert rms["window"] == pytest.approx((2.46, 2.64))
    # the verdict is the three de Vries rows only; synthetic radii either side
    dv, _ = rmod.load_reference()
    assert rmod.verdict(2.47, dv)[1] == "pass"
    assert rmod.verdict(2.45, dv)[1] == "fail"      # Li71a at -2.2 sigma
    assert rmod.verdict(2.65, dv)[1] == "fail"      # Su67 at +2.2 sigma
    # Angeli: a recorded closure, decomposed exactly
    c = rms["closure"]
    assert c["pull"] == pytest.approx(-0.4615, abs=1e-3)
    assert c["dr2"] == pytest.approx(-0.09288, abs=1e-5)
    assert abs(c["parts_residual"]) < 1e-9
    assert c["parts"]["darwin_foldy"] == -0.033


def test_devries_report_line(rmod, capsys):
    rep = rmod.run(verbose=True)
    out = capsys.readouterr().out.splitlines()
    rows = [l for l in out if l.startswith("REPORT |")]
    assert len(rows) == 1
    assert rows[0].startswith("REPORT | t3_li6_rms_devries | HoSpin1FF.for_ion(li6()) "
                              "C0 slope, nucleon-folded (G_E^p + G_E^n): r_ch = "
                              "2.57100 fm")
    assert rows[0].endswith("| PASS")
    assert "np.float64" not in "\n".join(out)
    assert rep["runtime_s"] < 30.0


def test_devries_exit_codes(tmp_path):
    ok = _run_harness(os.path.join(_BENCH, "t3_li6_rms_devries.py"))
    assert ok.returncode == 0, ok.stderr
    assert "REPORT | t3_li6_rms_devries |" in ok.stdout
    bad = _exit_code_on_altered_data(tmp_path, "t3_li6_rms_devries",
                                     "devries1987_li6_rms.json",
                                     b'"rms_fm": 2.54,', b'"rms_fm": 2.64,')
    assert bad.returncode == 2 and "sha256" in bad.stderr
    assert "REPORT |" not in bad.stdout
