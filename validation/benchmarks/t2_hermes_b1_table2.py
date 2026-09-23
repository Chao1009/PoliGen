#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 3 -- HERMES b1^d Table II as an assertion
(tier T2).

WHAT IS COMPARED.  The six (<x>, <Q^2>) rows of HERMES PRL 95 (2005) 242001
Table II (`data/hermes2005_b1d_table2.json`) against LiPolGen's deuteron b1,
per nucleon, at each row's own (<x>, <Q^2>), in five configurations and one
null baseline, with NO tuning of any of them:

  cdks_conv_mstw     DeuteronConvolutionB1 with unpol = MstwSF() -- the A = 2
                     gate's passing configuration (AV18 `fdeut.av18`, CDKS
                     Eq. (21) finite-|q| delta, r1998, CDKS Eq. (22) target
                     mass): the MstwSF path.  Needs the PYTHIA 8 tier.
  cdks_conv_default  DeuteronConvolutionB1() -- the class's shipped default
                     Options (unpol = ToyF2, otherwise as above).
  miller_shipped     MillerB1() -- the deuteron b1 the generator ships
                     (`--b1-model miller`): Miller's per-DEUTERON curve times
                     B1_MILLER_TABLE_TO_PER_NUCLEON = 0.5 (registry row 2,
                     status quo).
  miller_x1          the same curve WITHOUT the 0.5 (registry row 2's other
                     option), i.e. MillerB1() / B1_MILLER_TABLE_TO_PER_NUCLEON.
                     A PRICE, not a proposal.
  cdks_digitized     CdksB1() -- the digitized CDKS Fig. 4 column (Q^2 = 2.5,
                     per nucleon, B1_CDKS_TABLE_TO_PER_NUCLEON = 1), read at
                     every row's x regardless of its Q^2.
  null_b1_zero       b1 = 0 everywhere: NOT a LiPolGen model, the baseline a
                     chi2 must beat to say the data prefer a model to nothing.

A_zz^d is formed the way HERMES formed it, inverted: A_zz = -(2/3) b1/F1^d
with F1^d = (1 + gamma^2) F2^d/(2x(1 + R1990)), gamma^2 = 4 M^2 x^2/Q^2,
F2^d = (F2^p + F2^n)/2 per nucleon.  HERMES took F2^p from ALLM97 and
F2^n/F2^p from NMC 1992, neither of which is in the tree; this harness uses
MstwSF (ToyF2 if the PYTHIA tier is absent) and PRICES the difference against
the F1^d HERMES's own table implies, F1_impl = -(2/3) b1/A_zz (rounded
3-digit entries, so that ratio is itself good to ~1-5 %).  R1990 is
`data/whitlow1990_r1990.json`; the R1998 alternative is priced alongside.

TOLERANCE.  Per bin, the quadrature sum sqrt(stat^2 + syst^2) of the row
(the systematic includes a ~1e-3 A_zz normalisation that HERMES says is
correlated over bins; treated here as uncorrelated, which is what the plan
row asks).  Per configuration: chi2 = sum (model - data)^2/sigma^2 over the
six bins, ndf = 6 (no parameter is fitted).  STATUS: "pass" iff every bin is
within its tolerance; the chi2 is reported either way.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations, both sides named.  HERMES: PER NUCLEON (their F2^d =
     (F2^p + F2^n)/2).  DeuteronConvolutionB1 / CdksB1: PER NUCLEON, the 0.5
     dropped (CDKS Eq. (10)'s explicit 1/A; B1_CDKS_TABLE_TO_PER_NUCLEON = 1,
     CERTAIN).  MillerB1: Miller's table is PER DEUTERON and the shipped
     value keeps the 0.5 (B1_MILLER_TABLE_TO_PER_NUCLEON; LIKELY, NOT CERTAIN
     -- registry row 2).  A-scaling: none -- A = 2 on both sides; no 2/6
     (that is the 6Li inert-alpha assumption and is NOT applied here).  x is
     Bjorken x on the NUCLEON mass on both sides (runs to 2 for the deuteron;
     CDKS's convention, and HERMES's).
  2. isoscalar denominator: F2^d, F1^d are the isoscalar nucleon
     (F2^p + F2^n)/2 on both sides -- asserted in `measure()` by recomputing
     F1^d from f2p and f2n explicitly, never from a nuclear grid.
  3. structure functions per nucleon (no cross section enters).

Born level: HERMES's A_zz is RC-unfolded (RADGEN); the comparison is to a
Born structure function.  Q^2 = 0.51 at bin 1 is below MSTW2008's 1 GeV^2
grid floor (PYTHIA's reader extrapolates) and below R1990's SLAC x range.

Run:  source env.sh && python3 validation/benchmarks/t2_hermes_b1_table2.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import json
import math
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "hermes2005_b1d_table2.json")
R1990_DATA = os.path.join(HERE, "data", "whitlow1990_r1990.json")

NAME = "t2_hermes_b1_table2"
M_NUCLEON = 0.938272   # GeV; gamma^2 = 4 M^2 x^2 / Q^2 = Q^2/nu^2
CONFIGS = ("cdks_conv_mstw", "cdks_conv_default", "miller_shipped",
           "miller_x1", "cdks_digitized", "null_b1_zero")


def load_table(path=DATA):
    with open(path) as f:
        d = json.load(f)
    s = d["scale"]
    out = []
    for r in d["rows"]:
        x, q2, a, sa, ya, b, sb, yb = r
        out.append(dict(x=x, q2=q2, azz=a * s, azz_err=math.hypot(sa, ya) * s,
                        b1=b * s, b1_err=math.hypot(sb, yb) * s,
                        f1_implied=-(2.0 / 3.0) * b / a))
    return out


def r1990_factory(path=R1990_DATA):
    with open(path) as f:
        p = json.load(f)

    def r1990(x, q2):
        th = 1.0 + 12.0 * (q2 / (q2 + 1.0)) * (0.125 ** 2 / (0.125 ** 2 + x * x))
        return (p["b1"] / math.log(q2 / p["Lambda2_GeV2"]) * th + p["b2"] / q2
                + p["b3"] / (q2 * q2 + 0.09))
    return r1990


def have_mstw():
    import lipolgen._lipolgen as _l
    if not getattr(_l, "HAVE_PYTHIA8", False):
        return False, "PYTHIA 8 tier not built"
    d = _l.pythia8_pdfdata_dir()
    if not (d and os.path.isfile(os.path.join(d, "mstw2008lo.00.dat"))):
        return False, "PYTHIA 8 built without pdfdata/mstw2008lo.00.dat"
    return True, ""


def build_models():
    import lipolgen._lipolgen as _l
    models, blocked = {}, {}
    ok, why = have_mstw()
    if ok:
        o = _l.DeuteronConvolutionB1.Options()
        o.unpol = _l.MstwSF()
        models["cdks_conv_mstw"] = _l.DeuteronConvolutionB1(o)
    else:
        blocked["cdks_conv_mstw"] = why
    models["cdks_conv_default"] = _l.DeuteronConvolutionB1()
    mil = _l.MillerB1()
    models["miller_shipped"] = mil
    models["miller_x1"] = mil          # divided below, at evaluation
    models["cdks_digitized"] = _l.CdksB1()
    models["null_b1_zero"] = None
    f2 = _l.MstwSF() if ok else _l.ToyF2()
    return models, blocked, f2, ("MstwSF" if ok else "ToyF2")


def measure():
    import lipolgen._lipolgen as _l
    rows = load_table()
    r1990 = r1990_factory()
    models, blocked, f2, f2_name = build_models()
    miller_undo = 1.0 / _l.B1_MILLER_TABLE_TO_PER_NUCLEON
    per_bin = []
    for r in rows:
        x, q2 = r["x"], r["q2"]
        f2p, f2n = f2.f2p(x, q2), f2.f2n(x, q2)
        f2d = 0.5 * (f2p + f2n)                 # rule 2: the isoscalar nucleon
        g2 = 4.0 * M_NUCLEON ** 2 * x * x / q2
        r90, r98 = r1990(x, q2), _l.r1998(x, q2)
        f1d = (1.0 + g2) * f2d / (2.0 * x * (1.0 + r90))
        f1d_r98 = (1.0 + g2) * f2d / (2.0 * x * (1.0 + r98))
        b = {}
        for name, m in models.items():
            v = 0.0 if m is None else m.b1(x, q2, 0.0)  # f1 arg unused here
            b[name] = v * miller_undo if name == "miller_x1" else v
        per_bin.append(dict(
            x=x, q2=q2, data=r, f2d=f2d, f1d=f1d, f1d_r1998=f1d_r98,
            r1990=r90, r1998=r98, f1_ratio_to_implied=f1d / r["f1_implied"],
            b1=b,
            azz={k: -(2.0 / 3.0) * v / f1d for k, v in b.items()},
            azz_r1998={k: -(2.0 / 3.0) * v / f1d_r98 for k, v in b.items()},
        ))
    return per_bin, blocked, f2_name


def chi2_of(per_bin, name, obs, key=None):
    key = key or obs
    c, within = 0.0, 0
    for p in per_bin:
        d = p["data"][obs]
        e = p["data"][obs + "_err"]
        m = p[key][name]
        c += ((m - d) / e) ** 2
        within += abs(m - d) <= e
    return c, within


def run(verbose=True):
    per_bin, blocked, f2_name = measure()
    configs = {}
    for name in CONFIGS:
        if name in blocked:
            configs[name] = dict(status="blocked", reason=blocked[name])
            continue
        cb, wb = chi2_of(per_bin, name, "b1")
        ca, wa = chi2_of(per_bin, name, "azz")
        ca98, _ = chi2_of(per_bin, name, "azz", key="azz_r1998")
        configs[name] = dict(
            chi2_b1=cb, within_b1=wb, chi2_azz=ca, within_azz=wa,
            chi2_azz_r1998=ca98, ndf=len(per_bin),
            status="pass" if (wb == len(per_bin) and wa == len(per_bin)) else "fail")
    head = "cdks_conv_mstw" if "cdks_conv_mstw" not in blocked else "cdks_conv_default"
    h = configs[head]
    report = dict(
        name=NAME,
        generator="%s: chi2/ndf b1 = %.2f/%d, A_zz = %.2f/%d; %d/%d b1 bins "
                  "within tol" % (head, h["chi2_b1"], h["ndf"], h["chi2_azz"],
                                  h["ndf"], h["within_b1"], h["ndf"]),
        reference="HERMES PRL 95 (2005) 242001 Table II, 6 bins, per nucleon",
        tolerance="per bin sqrt(stat^2+syst^2)",
        status=h["status"], headline=head, configs=configs,
        per_bin=per_bin, f2_source=f2_name,
    )
    if verbose:
        print("F1^d rebuilt as HERMES did, with F2 = (F2p+F2n)/2 from %s and "
              "R1990, vs the F1^d implied by their table (-2/3 b1/A_zz):" % f2_name)
        for p in per_bin:
            print("  x=%.3f Q2=%.2f  F1_rebuilt=%.4f F1_implied=%.4f ratio=%.4f"
                  "  R1990=%.3f R1998=%.3f (1+R90)/(1+R98)=%.4f" % (
                      p["x"], p["q2"], p["f1d"], p["data"]["f1_implied"],
                      p["f1_ratio_to_implied"], p["r1990"], p["r1998"],
                      (1 + p["r1990"]) / (1 + p["r1998"])))
        print("b1^d per nucleon, units 1e-2 (HERMES +- quadrature error):")
        print("  x      HERMES          " + " ".join("%17s" % n for n in CONFIGS
                                                     if n not in blocked))
        for p in per_bin:
            d = p["data"]
            print("  %.3f  %+6.2f +- %5.2f  " % (p["x"], 100 * d["b1"], 100 * d["b1_err"])
                  + " ".join("%17.4f" % (100 * p["b1"][n]) for n in CONFIGS
                             if n not in blocked))
        print("A_zz^d, units 1e-2:")
        for p in per_bin:
            d = p["data"]
            print("  %.3f  %+6.2f +- %5.2f  " % (p["x"], 100 * d["azz"], 100 * d["azz_err"])
                  + " ".join("%17.4f" % (100 * p["azz"][n]) for n in CONFIGS
                             if n not in blocked))
        for n in CONFIGS:
            c = configs[n]
            if c["status"] == "blocked":
                print("REPORT | %s[%s] | BLOCKED: %s" % (NAME, n, c["reason"]))
                continue
            print("REPORT | %s[%s] | chi2/ndf b1 = %.3f/%d (%d/6 within), A_zz = "
                  "%.3f/%d (%d/6 within; with R1998 %.3f) | HERMES Table II | "
                  "tol per-bin quadrature | %s" % (
                      NAME, n, c["chi2_b1"], c["ndf"], c["within_b1"],
                      c["chi2_azz"], c["ndf"], c["within_azz"],
                      c["chi2_azz_r1998"], c["status"].upper()))
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
