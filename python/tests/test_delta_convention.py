# SPDX-License-Identifier: GPL-3.0-or-later
"""D4 (2026-10-27): the kernel's cos 2phi coefficient is Jaffe-Manohar's.

The author's decision D4 (docs/open_items/run_2026-10-27/DECISIONS.md;
docs/CONVENTIONS.md, "Delta normalisation") gives the double-helicity-flip
term the coefficient every source that derives it prints -- Jaffe-Manohar
PLB 223 (1989) 218 Eq. (16) (preprint CTP#1712 numbering), Sather-Schmidt
PRD 42 (1990) 1424 Eq. (A7), Detmold-Shanahan PRD 94 (2016) 014507 Eq. (5)
contracted, Ma-Wang-Zhang arXiv:1306.6693 Eq. (23):

    a2 = -(1/2) (1-y)/y^2 c_eff(m) sin^2(theta_S) Delta/D_phi,

half the JLab LOI's Eq. (9), which the tree carried until that date.  Every
other test that reads a2 compares the kernel with `a_cos2phi`, and the two
move together, so these pin the coefficient itself: the constant, the closed
form written out with a LITERAL 1/2, four absolute amplitudes, the m pattern,
and the provenance trace (`meta["delta_convention"]` of `Pipeline.generate()`,
so of every .npz the CLI writes, and the CLI run banner; a HepMC3 file carries
none), without which such a file made after D4 could not be told from one made
before.  The C++ twin is tests/test_xsec.cpp, "D4: the cos 2phi coefficient is
Jaffe-Manohar's".
"""
import math

import pytest

lg = pytest.importorskip("lipolgen")

#: the 6Li sweet-spot cells (x, Q2) of the letter's money plot, at
#: default_configs('6Li')[1] (s = 3980 GeV2), m = +1, theta_S = 90 deg,
#: unpolarized electrons, on default_inclusive_kernel(li6()) -- the shipped
#: toy Delta at 1e-2 on ToyF2 (a backend change moves these; re-measure)
CELLS = ((0.0112, 1.14), (0.0282, 1.14), (0.0708, 3.13), (0.141, 14.3))
#: a2 there, measured on the D4 build 2026-10-27: exactly 1/2 of the pre-D4
#: -1.055948214152054e-03 / -1.299961141425631e-03 / -1.440230754221272e-03 /
#: -1.326442680639959e-03 (measured on ca9b23a, plan_D4 sec. 2.1)
A2_D4 = (-5.279741070760268e-04, -6.499805707128153e-04,
         -7.201153771106358e-04, -6.632213403199796e-04)


def _state(m, theta_s=math.pi / 2):
    return lg.EventSpinState(lam_e=0, pe=0.0, j=1.0, m=m, theta_s=theta_s,
                             phi_s=0.0)


def test_the_constant_and_its_label():
    assert lg.DELTA_COS2PHI_NORM == 0.5
    assert lg.DELTA_CONVENTION == ("jaffe-manohar-1989; cos2phi coefficient "
                                   "-(1/2)(1-y)/y^2 (D4, 2026-10-27)")


def test_a_cos2phi_is_the_jaffe_manohar_closed_form():
    """m = +-1 along an axis at 90 deg: -(1/2)(1-y)/y^2 Delta/D_phi, the 1/2
    a literal (the LOI's Eq. (9), -(1-y)/y^2, was the tree's before D4)."""
    for x, y in ((0.01, 0.7), (0.1, 0.5), (0.3, 0.3)):
        for f1, f2, delta in ((1.0, 0.0, 1.0), (0.8, 0.5, -3e-3)):
            dphi = f1 + (1.0 - y) / (x * y * y) * f2
            assert lg.a_cos2phi(delta, f1, f2, x, y) == pytest.approx(
                -0.5 * (1.0 - y) / y ** 2 * delta / dphi, rel=1e-14)


def test_the_kernel_carries_the_closed_form_at_every_m_and_angle():
    k = lg.default_inclusive_kernel(lg.li6())
    s = lg.default_configs("6Li")[1].s_per_nucleon()
    for x, q2 in CELLS:
        t = k.tables(x, q2)
        y = q2 / (s * x)
        dphi = t.f1 + (1.0 - y) / (x * y * y) * t.f2
        for m in (1.0, 0.0, -1.0):
            c_eff = 3.0 * m * m - 2.0
            for th in (0.7, math.pi / 2):
                want = (-0.5 * (1.0 - y) / y ** 2 * c_eff * math.sin(th) ** 2
                        * t.delta / dphi)
                got = k.amplitudes(t, x, q2, s, _state(m, th)).a2
                assert got == pytest.approx(want, rel=1e-12), (x, q2, m, th)


def test_the_sweet_spot_amplitudes():
    """The four absolute numbers, and the m pattern: m = 0 carries -2 times
    m = +1, m = -1 the same as m = +1 (c_eff = 3m^2 - 2)."""
    k = lg.default_inclusive_kernel(lg.li6())
    s = lg.default_configs("6Li")[1].s_per_nucleon()
    for (x, q2), want in zip(CELLS, A2_D4):
        t = k.tables(x, q2)
        a_plus = k.amplitudes(t, x, q2, s, _state(1.0)).a2
        assert a_plus == pytest.approx(want, rel=1e-12), (x, q2)
        assert a_plus == lg.a_cos2phi(t.delta, t.f1, t.f2, x, q2 / (s * x))
        assert k.amplitudes(t, x, q2, s, _state(0.0)).a2 == -2.0 * a_plus
        assert k.amplitudes(t, x, q2, s, _state(-1.0)).a2 == a_plus


def test_a_generated_file_says_which_convention_it_carries():
    """`meta["delta_convention"]` on every run (a constant of the build, so
    it moves no diff between two runs), and the same string in the banner."""
    plan = lg.make_plan("transverse-tensor", j=1.0, pzz=0.6)
    cfg = lg.make_config(isotope="6Li", events=200, seed=7)
    cols = lg.Pipeline(cfg, plan).generate(0, False, 1)
    assert cols["meta"]["delta_convention"] == lg.DELTA_CONVENTION
    assert "Jaffe-Manohar" in cols["meta"]["rank2_input"]


def test_the_cli_banner_prints_the_convention(tmp_path, capsys):
    from lipolgen import cli
    rc = cli.main(["--isotope", "6Li", "--events", "200", "--seed", "7",
                   "--npz", str(tmp_path / "d4.npz")])
    assert rc == 0
    out = capsys.readouterr().out
    assert ("  Delta (cos 2phi) convention: %s" % lg.DELTA_CONVENTION) in out
