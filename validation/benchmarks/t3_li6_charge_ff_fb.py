#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 5 -- the UVa 6Li Fourier-Bessel charge
density vs the tree's C0 form factor (tier T3).

WHAT IS COMPARED.  The seven FB coefficients and R = 6.0 fm of the 6Li row of
the UVa Nuclear Charge Density archive's `FB_data.dat` (vendored, one row,
byte for byte, in `data/uva_ncd_fb_data_li6.dat` with the file's sha256 and
the row's) give rho(r) = SUM a_n j0(n pi r/R), r <= R.  From it this harness
computes, in closed form (each term's Fourier-Bessel transform is analytic):

  * F_C(q) = 4 pi INT rho j0(q r) r^2 dr, normalised to F_C(0) = Z = 3 (the
    truncation integrates to 2.99117, 0.29 % low; normalising by the row's
    own F_C(0) removes that and nothing else);
  * its FIRST ZERO q_FB (bisection to 1e-12 fm^-1);
  * <r^2>^1/2 = [INT rho r^4 / INT rho r^2]^1/2.

and compares with the tree's `HoSpin1FF.for_ion(li6())` C0 sector, SHIPPED
defaults, nothing moved:

  * the first zero of `fc(t)` (t = (q hbar c)^2).  `fc` is Z F_point(q)
    [G_E^p + G_E^n], and the nucleon factor has no zero, so this IS the
    shipped q0 = 3.0998 fm^-1 that tests/test_rc.cpp T11 gates in
    [2.9, 3.3] fm^-1 -- that window is NOT moved by this harness;
  * the minimum of the tree's own longitudinal |F_L|^2 = F_C0^2 + F_C2^2,
    F_C2 = (2 sqrt 2/3) eta_A F_q (rc.hpp's multipole conversion), which is
    the quantity Li et al. 1971 actually saw a minimum in;
  * the charge rms of the FOLDED `fc`, from its q -> 0 slope.

THE MUST-NOT-CONTRADICT BAND.  Two data-derived numbers bound where the 6Li
C0 zero / |F_L|^2 minimum is:
  q_FB  = 2.6944 fm^-1   (this file; 06_critic.md sec. 4.5 printed 2.6950)
  q_Li  = sqrt(8) = 2.8284 fm^-1  (Li, Sick, Whitney, Yearian, NPA 162
          (1971) 583: diffraction minimum at q^2 = 8 fm^-2, as quoted by
          03_data_nuclear.md sec. 2 -- the paper is NOT held; its precision
          is not known here)
The row PASSES if the tree's zero lies inside [q_FB, q_Li]; otherwise it is
a recorded FAIL carrying its distance from the band's upper edge.  A FAIL
here tunes nothing: the HO refit is open item Q1, and this harness only
prices it (docs/open_items/run_2026-09-23/phase_B2_nuclear.md).

CAVEATS, stated so the band is not over-read:
  (i)  q_Li is a minimum of a MEASURED cross section: Coulomb distortion
       (small at Z = 3, not computed here) and C2 fill-in enter it.  C2
       fill-in cannot rescue the tree's zero: `HoSpin1FF`'s F_q shares
       F_point with F_c, so its C2 vanishes at the SAME q0 and the tree's
       |F_L|^2 minimum sits on q0 (measured 2026-09-23: |q_min - q0| =
       1.7e-11 fm^-1, the golden-section tolerance; |F_L|^2 = 6.7e-30 there).
  (ii) the FB row's PROVENANCE IS UNVERIFIED (06 sec. 4.5: 6Li is in neither
       the ADNDT 14 nor the ADNDT 36 Table IV) -- verified content,
       unverified source.  Its zero is inside the fitted data's q range only
       if the underlying fit reached q ~ 2.7 fm^-1 (Li71 reached 3.66).
  (iii) the FB density is a CHARGE density (nucleon size folded in); the
       tree's zero is the zero of F_point, identical to that of the folded
       `fc` -- so the zeros are like for like, the radii are compared as
       charge radii only.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: no b1, no A-scaling assumption.
  2. isoscalar denominator -- NOT APPLICABLE: no structure function.  (The
     tree folds its point shape with the ISOSCALAR G_E^p + G_E^n, N = Z = 3;
     the FB density already contains the nucleons' charge distribution.)
  3. per nucleus / per nucleon -- both sides are PER NUCLEUS: the FB density
     is normalised to 4 pi INT rho r^2 dr = Z e; the tree's F_c(0) = Z = 3.
     Both are compared as F_C(q)/F_C(0), so the normalisation cancels in the
     zero and the rms.

Run:  source env.sh && python3 validation/benchmarks/t3_li6_charge_ff_fb.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import hashlib
import math
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "uva_ncd_fb_data_li6.dat")
MARKER = "# ---- BEGIN VERBATIM FB_data.dat LINES ----"
ROW_SHA256 = "b9655bdc97151a58614cd949c5d5a03c54abfee8ae4609c2053c65e482b24a00"

NAME = "t3_li6_charge_ff_fb"
Q_LI71 = math.sqrt(8.0)          # fm^-1, Li et al. 1971 |F_L|^2 minimum
T11_WINDOW = (2.9, 3.3)          # tests/test_rc.cpp T11 -- NOT moved
DE_VRIES_RMS = (2.54, 2.57)      # ADNDT 36 (1987) Table I, 6Li, three analyses
ANGELI_RMS = 2.589               # Angeli & Marinova ADNDT 99 (2013), rc.hpp


def load_row(path=DATA):
    raw = open(path, "rb").read()
    _, sep, body = raw.partition((MARKER + "\n").encode())
    if not sep:
        raise RuntimeError("%s: provenance marker missing" % path)
    lines = body.split(b"\n")
    row = [l for l in lines if l.strip().startswith(b"6Li ")]
    if len(row) != 1:
        raise RuntimeError("%s: expected exactly one 6Li row" % path)
    sha = hashlib.sha256(row[0]).hexdigest()
    if sha != ROW_SHA256:
        raise RuntimeError("%s: 6Li row sha256 %s != recorded %s"
                           % (path, sha, ROW_SHA256))
    f = row[0].split()
    # the archive's header says "Z A" but the values are A, Z (see the file)
    a_num, z_num = int(f[1]), int(f[2])
    coeffs = [float(v) for v in f[3:-1]]
    radius = float(f[-1])
    while coeffs and coeffs[-1] == 0.0:
        coeffs.pop()
    return dict(A=a_num, Z=z_num, a=coeffs, R=radius)


class FourierBessel:
    """rho(r) = SUM_n a_n j0(n pi r/R), r <= R; everything in closed form."""

    def __init__(self, a, R):
        self.a, self.R = list(a), float(R)
        self.k = [(n + 1) * math.pi / self.R for n in range(len(self.a))]

    def rho(self, r):
        if r > self.R:
            return 0.0
        return sum(an * (1.0 if r == 0 else math.sin(k * r) / (k * r))
                   for an, k in zip(self.a, self.k))

    def charge(self):
        """4 pi INT_0^R rho r^2 dr = 4 pi SUM a_n R^3 (-1)^(n+1)/(n pi)^2."""
        return 4 * math.pi * sum(
            an * self.R ** 3 * (-1) ** (n + 1) / (n * math.pi) ** 2
            for n, an in zip(range(1, len(self.a) + 1), self.a))

    def r4(self):
        """4 pi INT_0^R rho r^4 dr, termwise closed form."""
        tot = 0.0
        for n, an in zip(range(1, len(self.a) + 1), self.a):
            k = n * math.pi / self.R
            # INT_0^R sin(kr) r^3 dr / k, with sin(kR) = 0, cos(kR) = (-1)^n
            c = (-1) ** n
            integral = (-self.R ** 3 * c / k + 6.0 * self.R * c / k ** 3) / k
            tot += an * integral
        return 4 * math.pi * tot

    def ff(self, q):
        """4 pi INT rho j0(qr) r^2 dr (un-normalised), termwise closed form."""
        if q < 1e-9:
            return self.charge()
        tot = 0.0
        for an, k in zip(self.a, self.k):
            d = k - q
            s1 = self.R / 2.0 if abs(d) < 1e-12 else math.sin(d * self.R) / (2 * d)
            s2 = math.sin((k + q) * self.R) / (2 * (k + q))
            tot += an * (s1 - s2) / (k * q)
        return 4 * math.pi * tot

    def ff_norm(self, q):
        return self.ff(q) / self.charge()

    def rms(self):
        return math.sqrt(self.r4() / self.charge())


def bisect(f, lo, hi, tol=1e-12):
    flo = f(lo)
    if flo * f(hi) > 0:
        raise ValueError("no sign change in [%g, %g]" % (lo, hi))
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def first_zero(f, q_max=5.0, step=0.01):
    q = step
    fq = f(q)
    while q < q_max:
        q2 = q + step
        f2 = f(q2)
        if fq * f2 <= 0:
            return bisect(f, q, q2)
        q, fq = q2, f2
    return None


def golden_min(f, lo, hi, tol=1e-10):
    g = (math.sqrt(5) - 1) / 2
    c, d = hi - g * (hi - lo), lo + g * (hi - lo)
    while hi - lo > tol:
        if f(c) < f(d):
            hi = d
        else:
            lo = c
        c, d = hi - g * (hi - lo), lo + g * (hi - lo)
    return 0.5 * (lo + hi)


def tree_side():
    import lipolgen as lg  # noqa: WPS433
    hbarc = lg.HBARC_GEV_FM
    ion = lg.li6()
    ff = lg.HoSpin1FF.for_ion(ion)
    m_a = ion.mass()
    t_of = lambda q: (q * hbarc) ** 2
    fc = lambda q: ff.fc(t_of(q))
    q0 = first_zero(fc)

    def fl2(q):
        t = t_of(q)
        eta = t / (4.0 * m_a * m_a)
        c2 = 2.0 * math.sqrt(2.0) / 3.0 * eta * ff.fq(t)
        return ff.fc(t) ** 2 + c2 ** 2

    q_min_fl2 = golden_min(fl2, q0 - 0.3, q0 + 0.3)
    eps = 1e-3
    rms_folded = math.sqrt(-6.0 * (fc(eps) / fc(0.0) - 1.0) / eps ** 2)
    return dict(q0=q0, q_min_fl2=q_min_fl2, fl2_at_min=fl2(q_min_fl2),
                rms_folded=rms_folded, fc0=fc(0.0),
                provenance=ff.provenance())


def measure(path=DATA):
    row = load_row(path)
    fb = FourierBessel(row["a"], row["R"])
    q_fb = first_zero(fb.ff_norm)
    tree = tree_side()
    band = (q_fb, Q_LI71)
    inside = band[0] <= tree["q0"] <= band[1]
    dist_upper = tree["q0"] - band[1]
    return dict(
        row=row, fb_charge=fb.charge(), fb_rms=fb.rms(), q_fb=q_fb,
        fb_at_tree_q0=row["Z"] * fb.ff_norm(tree["q0"]),
        band=band, inside=inside, distance_above_band=dist_upper,
        rel_above_q_li=dist_upper / band[1],
        rel_above_q_fb=(tree["q0"] - band[0]) / band[0],
        t11_window=T11_WINDOW,
        band_meets_t11=not (band[1] < T11_WINDOW[0] or band[0] > T11_WINDOW[1]),
        gap_band_to_t11=T11_WINDOW[0] - band[1],
        de_vries_rms=DE_VRIES_RMS, angeli_rms=ANGELI_RMS,
        tree=tree,
        fb_table=[(q, row["Z"] * fb.ff_norm(q)) for q in
                  (0.5, 1.0, 1.5, 2.0, 2.5, q_fb, Q_LI71, tree["q0"], 3.5)],
    )


def run(path=DATA, verbose=True):
    m = measure(path)
    status = "pass" if m["inside"] else "fail"
    report = dict(
        name=NAME,
        generator="HoSpin1FF first C0 zero q0 = %.4f fm^-1 (shipped; T11 window "
                  "[%.1f, %.1f] not moved)" % (m["tree"]["q0"], *T11_WINDOW),
        reference="band [%.4f (UVa FB zero), %.4f (Li71 minimum)] fm^-1"
                  % m["band"],
        tolerance="q0 inside the band (must not contradict)",
        status=status,
        **{k: v for k, v in m.items()},
    )
    if verbose:
        t = m["tree"]
        print("FB row: A = %d, Z = %d, R = %.1f fm, %d coefficients; "
              "4 pi INT rho r^2 dr = %.6f" % (m["row"]["A"], m["row"]["Z"],
                                              m["row"]["R"], len(m["row"]["a"]),
                                              m["fb_charge"]))
        print("   q [fm^-1]   F_C^FB(q) (F_C(0) = Z)")
        for q, v in m["fb_table"]:
            print("   %8.4f    %+.6e" % (q, v))
        print("FB: first zero %.6f fm^-1, <r^2>^1/2 = %.4f fm "
              "(de Vries 2.54-2.57; Angeli %.3f)" % (m["q_fb"], m["fb_rms"],
                                                      ANGELI_RMS))
        print("tree HoSpin1FF: C0 zero %.6f fm^-1; |F_L|^2 minimum at %.6f "
              "(value %.2e: C2 shares the zero); folded charge rms %.4f fm"
              % (t["q0"], t["q_min_fl2"], t["fl2_at_min"], t["rms_folded"]))
        print("distance: q0 - q_Li71 = %+.4f fm^-1 (%+.1f %%), q0 - q_FB = %+.4f "
              "(%+.1f %%); band meets T11 window: %s (gap %.3f fm^-1)"
              % (m["distance_above_band"], 100 * m["rel_above_q_li"],
                 t["q0"] - m["q_fb"], 100 * m["rel_above_q_fb"],
                 m["band_meets_t11"], m["gap_band_to_t11"]))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, report["generator"], report["reference"], report["tolerance"],
            status.upper() + ("" if m["inside"] else
                              " (recorded, %+.4f fm^-1 above the band; open item Q1; "
                              "nothing moved)" % m["distance_above_band"])))
    return report


if __name__ == "__main__":
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    try:
        run()
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        sys.exit(2)
    sys.exit(0)
