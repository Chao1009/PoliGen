#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 4 -- NMC F2(6Li)/F2(D) vs the tree's EPPS21
baseline (tier T3).  The only DIS measurement on 6Li in existence.

WHAT IS COMPARED.  HEPData ins394050 Table 1 (24 points, CC0, vendored in
`data/hepdata_ins394050_table1.csv` with provenance and a sha256 of the
download) against `lipolgen.Epps21Ratio().f2_per_nucleon(x, Q2)` -- the grid
`EPPS21nlo_CT18Anlo_Li6` read through the tree's own class, default proton
baseline CT18NLO -- divided by the ISOSCALAR free nucleon
(F2p + F2n)/2 from `lipolgen.LhapdfSF("CT18NLO")`, at each point's OWN
(x, Q2).  Errors: stat (+) sys in quadrature; NMC's 0.004 overall
normalisation is NOT included (the table says so).  ZERO free parameters,
NO tuning: nothing in the tree is adjusted to this comparison.

WHICH POINTS.  Only points with Q2 >= max(q2_min) of the two LHAPDF grids
(EPPS21: 1.69 GeV2; CT18NLO: 1.677 GeV2) are compared -- 15 of 24.  The nine
points below (x <= 0.0085, Q2 = 0.034 .. 1.4 GeV2) are REPORTED as
below-grid and never enter a chi2: an LHAPDF extrapolation below the grid is
not the fit.  (The survey's 03_data_nuclear.md sec. 1.2 table printed the
x = 0.0085, Q2 = 1.4 point as on-grid; 1.4 < 1.69, so it is not.)

THE ROW.  The report row is the VALENCE subset x >= 0.30 (4 points), the
window `EMC_VALENCE_DEPLETION_EPPS21` averages over: chi2/ndf with the
tolerance p(chi2, ndf) >= 0.01.  The whole on-grid set (15 points) is
computed and carried in the report with the same tolerance, and the row
passes only if BOTH do.  The threshold p >= 0.01 is a conventional choice
made by this harness, AFTER the survey's 4.63/4 was known; it is stated here
so that it can be argued with, not tuned.

THE TRAP (03_data_nuclear.md sec. 1.3, printed once by `run()`): the EPPS21
LHAPDF grid is the ISOSCALAR-AVERAGED nucleon of 6Li (R_u = 0.578,
R_d = 2.829 at x = 0.65, Q2 = 5 are the n/p isospin average, not a nuclear
effect).  Dividing by the FREE PROTON gives a plausible-looking EMC curve
~9x too deep.  `run()` computes that WRONG number once, labelled WRONG, so
the size of the trap is on the record; it never enters the verdict.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: no b1 on either side, and no
     A-scaling assumption (no 2/6, no 0.5) is used anywhere in this row.
  2. isoscalar denominator -- APPLIED AND ASSERTED.  Generator side:
     F2^{6Li}/nucleon from the EPPS21 grid (the isoscalar-averaged bound
     nucleon) over (F2p + F2n)/2 of the FREE nucleon (CT18NLO; EPPS21's own
     baseline CT18ANLO is not installed here, which moves <1 - R> by 4.0 %
     relative over the valence window -- 03 sec. 1.2).  Reference side:
     NMC's ratio is itself corrected for "the slight non-isoscalarity of the
     targets" (4.5 % 7Li in the 6Li, 1.5 % H in the D; hep-ex/9504002).
     `isoscalar_denominator()` is the one place the denominator is built;
     `python/tests/test_bench_nuclear.py` asserts it IS (F2p + F2n)/2.
  3. per nucleus / per nucleon -- both sides are PER NUCLEON structure
     functions, so the ratio is 1 for no nuclear effect.  NOT like for like
     in ONE respect, stated rather than corrected (03 sec. 1.4): NMC divides
     by the DEUTERON, which carries its own few-percent EMC/binding effect at
     x >~ 0.5; the tree divides by the FREE isoscalar nucleon.  The true
     6Li-to-free-nucleon depletion is therefore LARGER than NMC's ratio
     shows, in the direction that would INCREASE 0.031.  No deuteron
     correction is applied (the tree has no model for one).

WHAT THIS DOES NOT VALIDATE (03 sec. 1.4): the polarized EMC effect (the
CBT/TMT curves this baseline scales), 7Li (no 7Li DIS data exists), or
shadowing point by point (NMC's x < 0.002 points are at Q2 = 0.03 .. 0.34,
below the grid).

Run:  source env.sh && python3 validation/benchmarks/t3_nmc_li6_over_d.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import csv
import hashlib
import math
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "hepdata_ins394050_table1.csv")
MARKER = "# ---- BEGIN VERBATIM HEPDATA CSV ----"
BODY_SHA256 = "3d8b92785ac82f3fccc3b440e66e8eb0852375e6bdc7fc31dfef1b14e4c1b376"

NAME = "t3_nmc_li6_over_d"
P_MIN = 0.01            # tolerance: p(chi2, ndf) >= 0.01
X_VALENCE = 0.30        # the valence subset (EMC_VALENCE_DEPLETION window)
NUCLEAR_SET = "EPPS21nlo_CT18Anlo_Li6"
PROTON_SET = "CT18NLO"


def load_table(path=DATA):
    """The 24 NMC points as dicts; refuses a body whose sha256 moved."""
    raw = open(path, "rb").read()
    head, sep, body = raw.partition((MARKER + "\n").encode())
    if not sep:
        raise RuntimeError("%s: provenance marker missing" % path)
    sha = hashlib.sha256(body).hexdigest()
    if sha != BODY_SHA256:
        raise RuntimeError("%s: vendored body sha256 %s != recorded %s"
                           % (path, sha, BODY_SHA256))
    lines = [l for l in body.decode().splitlines()
             if l.strip() and not l.startswith("#")]
    rows = list(csv.reader(lines))
    hdr = rows[0]
    assert hdr[:4] == ["X", "Q**2 [GEV**2]", "Y", "F2(Q=LI)/F2(Q=DEUT)"], hdr
    out = []
    for r in rows[1:]:
        x, q2, y, val, sp, sm, yp, ym = (float(v) for v in r)
        assert sp == -sm and yp == -ym  # symmetric, as the table says
        out.append(dict(x=x, q2=q2, y=y, ratio=val, stat=sp, sys=yp,
                        err=math.hypot(sp, yp)))
    return out


class Tree:
    """The tree side: Epps21Ratio over the ISOSCALAR free nucleon."""

    def __init__(self):
        import lipolgen as lg  # noqa: WPS433 (env.sh puts it on the path)
        self.lg = lg
        self.epps = lg.Epps21Ratio(NUCLEAR_SET, 0, PROTON_SET, 0)
        self.free = lg.LhapdfSF(PROTON_SET)
        self.q2_min = max(lg.LhapdfSF(NUCLEAR_SET).q2_min, self.free.q2_min)

    def isoscalar_denominator(self, x, q2):
        """RULE 2: (F2p + F2n)/2 of the free nucleon -- the ONLY denominator."""
        return 0.5 * (self.free.f2p(x, q2) + self.free.f2n(x, q2))

    def ratio(self, x, q2):
        return self.epps.f2_per_nucleon(x, q2) / self.isoscalar_denominator(x, q2)

    def ratio_WRONG_free_proton(self, x, q2):
        """The trap of 03 sec. 1.3.  Printed once; never in a verdict."""
        return self.epps.f2_per_nucleon(x, q2) / self.free.f2p(x, q2)


def chi2_sf(chi2, ndf):
    """Upper-tail probability of chi2 with ndf dof (scipy if present)."""
    try:
        from scipy.stats import chi2 as _c
        return float(_c.sf(chi2, ndf))
    except Exception:                               # pragma: no cover
        # regularised upper incomplete gamma by series, adequate here
        a, xx = ndf / 2.0, chi2 / 2.0
        term = s = 1.0 / a
        k = 1
        while term > 1e-15 * s:
            term *= xx / (a + k)
            s += term
            k += 1
        return 1.0 - math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * s


def measure(path=DATA):
    tree = Tree()
    pts = load_table(path)
    for p in pts:
        p["on_grid"] = p["q2"] >= tree.q2_min
        p["tree"] = tree.ratio(p["x"], p["q2"]) if p["on_grid"] else None
        p["tree_wrong"] = (tree.ratio_WRONG_free_proton(p["x"], p["q2"])
                           if p["on_grid"] else None)
        p["pull"] = ((p["ratio"] - p["tree"]) / p["err"]) if p["on_grid"] else None

    def chi2(sel, key="tree"):
        c = sum(((p["ratio"] - p[key]) / p["err"]) ** 2 for p in sel)
        return c, len(sel)

    grid = [p for p in pts if p["on_grid"]]
    val = [p for p in grid if p["x"] >= X_VALENCE]
    c_val, n_val = chi2(val)
    c_all, n_all = chi2(grid)
    c_wrong, _ = chi2(val, "tree_wrong")
    w = [1.0 / p["err"] ** 2 for p in val]
    wsum = sum(w)
    nmc_mean = sum(p["ratio"] * wi for p, wi in zip(val, w)) / wsum
    tree_mean = sum(p["tree"] * wi for p, wi in zip(val, w)) / wsum
    # the library constant's own definition, recomputed on the SHIPPED proton
    xs = [0.35 + 0.30 * i / 300.0 for i in range(301)]
    dep_iso = sum(1.0 - tree.ratio(x, 5.0) for x in xs) / len(xs)
    dep_wrong = sum(1.0 - tree.ratio_WRONG_free_proton(x, 5.0) for x in xs) / len(xs)
    return dict(
        points=pts, q2_min=tree.q2_min,
        chi2_valence=c_val, ndf_valence=n_val, p_valence=chi2_sf(c_val, n_val),
        chi2_ongrid=c_all, ndf_ongrid=n_all, p_ongrid=chi2_sf(c_all, n_all),
        chi2_valence_WRONG=c_wrong,
        nmc_valence_mean=nmc_mean, nmc_valence_mean_err=1.0 / math.sqrt(wsum),
        tree_valence_mean=tree_mean,
        depletion_iso_ct18nlo=dep_iso, depletion_WRONG=dep_wrong,
        depletion_constant=tree.lg.EMC_VALENCE_DEPLETION_EPPS21,
        r_u_065=tree.epps.ratio(2, 0.65, 5.0), r_d_065=tree.epps.ratio(1, 0.65, 5.0),
        n_below_grid=sum(1 for p in pts if not p["on_grid"]),
    )


def run(path=DATA, verbose=True):
    m = measure(path)
    ok = m["p_valence"] >= P_MIN and m["p_ongrid"] >= P_MIN
    status = "pass" if ok else "fail"
    report = dict(
        name=NAME,
        generator="chi2/ndf = %.4f/%d (x >= 0.30, p = %.3f); on-grid %.3f/%d (p = %.3f)"
                  % (m["chi2_valence"], m["ndf_valence"], m["p_valence"],
                     m["chi2_ongrid"], m["ndf_ongrid"], m["p_ongrid"]),
        reference="NMC F2(6Li)/F2(D), HEPData ins394050 Table 1 (%d of 24 on grid)"
                  % m["ndf_ongrid"],
        tolerance="p(chi2, ndf) >= %.2f on both sets" % P_MIN,
        status=status,
        **{k: v for k, v in m.items() if k != "points"},
        points=m["points"],
    )
    if verbose:
        print("   x        Q2     NMC F2(Li)/F2(D)     tree R_iso    pull")
        for p in m["points"]:
            if p["on_grid"]:
                print("%8.5f %7.3f   %.3f +- %.3f      %.4f      %+5.2f"
                      % (p["x"], p["q2"], p["ratio"], p["err"], p["tree"], p["pull"]))
            else:
                print("%8.5f %7.3f   %.3f +- %.3f      below grid (Q2 < %.3f)"
                      % (p["x"], p["q2"], p["ratio"], p["err"], m["q2_min"]))
        print("valence x >= 0.30: NMC weighted mean %.4f +- %.4f; tree at the same "
              "points %.4f" % (m["nmc_valence_mean"], m["nmc_valence_mean_err"],
                               m["tree_valence_mean"]))
        print("<1 - R_iso> over x in [0.35, 0.65], Q2 = 5, CT18NLO proton: %.6f "
              "(library constant, CT18ANLO proton: %.6f)"
              % (m["depletion_iso_ct18nlo"], m["depletion_constant"]))
        print("WRONG DENOMINATOR (free proton, the 03 sec. 1.3 trap) -- printed once, "
              "never used: <1 - R> = %.4f, chi2 = %.2f/%d on the valence set; "
              "R_u = %.3f, R_d = %.3f at x = 0.65, Q2 = 5"
              % (m["depletion_WRONG"], m["chi2_valence_WRONG"], m["ndf_valence"],
                 m["r_u_065"], m["r_d_065"]))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, report["generator"], report["reference"], report["tolerance"],
            status.upper()))
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
