#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T1 -- NMC F2 of the deuteron (survey row F-1) against the tree's
per-nucleon F2^d (tier T1, BENCHMARK_PLAN.md sec. 4 "then, by tier").

WHAT IS COMPARED.  NMC's F2^d PER NUCLEON -- M. Arneodo et al., Nucl. Phys. B
483 (1997) 3, hep-ph/9610231; HEPData ins424154 tables 17-32 (CC0), 158
points, 0.0045 <= x <= 0.5, 0.75 <= Q2 <= 65 GeV2 -- vendored VERBATIM from
the nnpdf-data 4.1.5 wheel (set NMC_NC_NOTFIXED_D, observable EM-F2-HEPDATA)
by `vendor_nnpdf_commondata.py` into `data/nnpdf_NMC_NC_NOTFIXED_D.json`,
whose `_provenance` block carries the wheel URL and sha256, the path inside
the wheel, the HEPData / INSPIRE record and the licences.  Against it: the
tree's per-nucleon F2 of the deuteron, (F2p + F2n)/2, at each point's OWN
(x, Q2), for three unpolarized backends (the `--unpol-sf` choices):

  toy      ToyF2() -- the SHIPPED default, i.e. exactly what the generator's
           `default_inclusive_kernel(deuteron()).tables(x, Q2).f2` returns.
           THE REPORT ROW'S VERDICT IS THIS CONFIGURATION'S.
  mstw     MstwSF() -- MSTW2008 LO read from PYTHIA 8's own pdfdata
           (`--unpol-sf mstw`).  Sub-row.  BLOCKED (environment) only if
           the PYTHIA 8 tier or its pdfdata grid is absent.
  ct18nlo  LhapdfSF("CT18NLO") (`--unpol-sf ct18nlo`).  Sub-row; BLOCKED
           (environment) when LHAPDF reports the set not installed --
           reported with that reason, never dropped, never faked.

ZERO free parameters, NO tuning: nothing in the tree is adjusted to this.

WHICH POINTS.  The chi2 is computed on ALL points and on THE DIS CUT,
declared here before any chi2 of this harness was computed: Q2 >= 1 GeV2 and
W2 = M^2 + Q2 (1 - x)/x >= 4 GeV2 (M = 0.938272 GeV), and in addition
Q2 >= the backend's own grid floor where it has one -- CT18NLO's
`LhapdfSF.q2_min()` (1.677); MSTW2008's grid starts at Q2 = 1 GeV2, so the
cut coincides with it; ToyF2 has no grid (it freezes its Q2 dependence below
1.1).  The W2 cut removes no point of this set (min W2 = 6.38 GeV2 above
Q2 = 1); 155 of the 158 points pass the cut for toy and mstw.

ERRORS AND WHICH CHI2 IS THE VERDICT.  PRIMARY (the verdict): the DIAGONAL
chi2, every uncertainty component of a point in quadrature (stat (+) sys;
the systematic treated as uncorrelated).  RECORDED beside it, never in the
verdict: the FULL-COVARIANCE chi2 with the components as NNPDF labels them,
C = diag(sum of UNCORR^2) + sum_k s_k s_k^T over the CORR components.  Here
the one 'sys' component is labelled MULT, CORR, so C makes it 100 %
correlated over all 158 points -- NNPDF's label, not a statement NMC makes
in the vendored file, which is why it is not the verdict.  Uncertainties are
the ABSOLUTE values given (a MULT component already multiplied by the data
central value; no t0 rescaling).  Whether NMC's overall normalisation
uncertainty is inside 'sys' is not stated (UNVERIFIED).

TOLERANCE.  p(chi2, ndf) >= 0.01 on the DIS-cut set, diagonal errors,
ndf = number of points (nothing is fitted).  The row's STATUS is the shipped
default's (toy); mstw and ct18nlo are sub-rows with their own status.  The
threshold follows row 4 (`t3_nmc_li6_over_d.py`) and was chosen AFTER the
plan scout's previews had been seen (scout, 2026-09-26: ToyF2 chi2 = 8308
and MSTW 417 over 155 points, diagonal, Q2 >= 1); it is stated here so that
it can be argued with, not tuned.  ToyF2 is labelled a TOY in sf.hpp
("anchored by eye ... factor-1.5 rate estimates ONLY"): its failure is a
RECORDED FAIL with its distance, not a defect of this harness, and moves
nothing.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations / A-scaling -- NOT APPLICABLE: no b1 on either side;
     A = 2 on both sides and no A-scaling assumption (no 0.5, no 2/6) is
     used anywhere in this row.
  2. isoscalar denominator -- APPLIED AND ASSERTED.  The tree side is the
     ISOSCALAR nucleon (F2p + F2n)/2 of the free nucleon, built in ONE place
     (`isoscalar_nucleon`); it is also, bit for bit, what the generator's
     kernel returns for the deuteron (`InclusiveKernel(deuteron).tables().f2`
     = NuclearF2 (Z F2p + N F2n)/A, no EMC hook) -- `run()` records the
     largest relative difference between the two, and
     python/tests/test_bench_t1_nucleon.py asserts both.  Reference side:
     the deuteron is isoscalar; no non-isoscalarity correction exists.
  3. per nucleus / per nucleon -- both sides are PER NUCLEON structure
     functions (NMC's values are the size of ONE nucleon's F2, 0.25-0.44 at
     small x, which excludes a per-deuteron reading).  NOT like for like in
     one respect, stated rather than corrected: NMC's F2d is the DEUTERON's
     (Fermi motion, binding, shadowing included; not corrected -- UNVERIFIED
     offline), the tree's is the free isoscalar nucleon; the tree has no
     deuteron-structure correction for F2.

WHAT THIS DOES NOT VALIDATE.  R = sigma_L/sigma_T (NMC's R points are not in
the package, and the R NMC extracted F2 with is not stated in the file),
F2 of the proton or neutron separately (only their sum enters), anything
nuclear beyond A = 2, and the polarized or tensor sector.

INTEGRITY.  The sha256 of every byte of the vendored file (`FILE_SHA256`)
and of its canonical body (`BODY_SHA256`, also recorded inside the file)
are re-checked on EVERY read: an altered table is refused with a
RuntimeError, so `__main__` exits 2, never 0 on moved numbers.

Run:  source env.sh && python3 validation/benchmarks/t1_nmc_f2d.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import hashlib
import json
import math
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "nnpdf_NMC_NC_NOTFIXED_D.json")
#: sha256 of every byte of the vendored file and of its canonical body, as
#: written by vendor_nnpdf_commondata.py from the pinned wheel (2026-09-26).
FILE_SHA256 = "d38b3c0c06ca4e123983cffd2f647d82d14f272e8ee04e13f73f9fd098af8e3d"
BODY_SHA256 = "66d87ea8ea7aa284b7b824fb4525ab72e31d3f1ed3bdd244c7c0a0fc56a847eb"

NAME = "t1_nmc_f2d"
P_MIN = 0.01             # tolerance: p(chi2, ndf) >= 0.01, DIS cut, diagonal
Q2_CUT = 1.0             # GeV^2 -- the declared DIS cut ...
W2_CUT = 4.0             # GeV^2 -- ... and its W^2 half
M_NUCLEON = 0.938272     # GeV, for W^2 only
CONFIGS = ("toy", "mstw", "ct18nlo")
HEADLINE = "toy"         # the shipped default carries the row's verdict
CT18 = "CT18NLO"
REFERENCE = ("NMC F2d per nucleon, NPB 483 (1997) 3, HEPData ins424154 "
             "tables 17-32 via nnpdf-data 4.1.5 (NMC_NC_NOTFIXED_D)")


# ---------------------------------------------------------------- reference

def read_vendored(path=None):
    """The vendored document -- RuntimeError unless the file's bytes AND its
    canonical body hash to the recorded digests."""
    path = path or DATA
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != FILE_SHA256:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, FILE_SHA256))
    doc = json.loads(raw.decode("ascii"))
    canon = json.dumps(doc["body"], sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False)
    body = hashlib.sha256(canon.encode("ascii")).hexdigest()
    if not body == doc["_body_sha256"] == BODY_SHA256:
        raise RuntimeError("%s: body sha256 %s; file records %s, harness %s"
                           % (path, body, doc["_body_sha256"], BODY_SHA256))
    return doc


def load_table(path=None, variant=None):
    """The points as dicts: x, q2, value, err (every component in
    quadrature: the DIAGONAL error), var_uncorr (sum of UNCORR^2) and corr
    ({name: signed absolute value} of the CORR components)."""
    body = read_vendored(path)["body"]
    variant = variant or body["primary_uncertainty_variant"]
    defs = body["uncertainty_variants"][variant]["definitions"]
    for k, d in defs.items():
        if d["type"] not in ("UNCORR", "CORR"):
            raise RuntimeError("%s: uncertainty %s has type %r (only UNCORR "
                               "and CORR are understood)" % (NAME, k, d["type"]))
    out = []
    for i, p in enumerate(body["points"]):
        u = {k: float(p["unc"][variant][k]) for k in defs}
        x, q2 = float(p["x"]), float(p["Q2"])
        out.append(dict(
            i=i, x=x, q2=q2, value=float(p["value"]),
            w2=M_NUCLEON ** 2 + q2 * (1.0 - x) / x,
            err=math.sqrt(sum(v * v for v in u.values())),
            var_uncorr=sum(v * v for k, v in u.items()
                           if defs[k]["type"] == "UNCORR"),
            corr={k: v for k, v in u.items() if defs[k]["type"] == "CORR"},
            unc=u))
    return out


# --------------------------------------------------------------- statistics

def chi2_sf_fallback(chi2, ndf):
    """Q(ndf/2, chi2/2), the regularised upper incomplete gamma, without
    scipy (Numerical Recipes' gammq): the series for P and Q = 1 - P when
    chi2/2 < ndf/2 + 1; otherwise the continued fraction for Q ITSELF
    (modified Lentz), so a small p is not lost to 1 - P cancellation (until
    2026-09-27 the fallback always took 1 - P: p ~ 1e-15 came out as noise
    or 0, and chi2 = 0 raised a math domain error)."""
    chi2, ndf = float(chi2), int(ndf)
    if ndf <= 0 or not chi2 >= 0.0:
        raise ValueError("chi2_sf: need ndf > 0 and chi2 >= 0 (got %r, %r)"
                         % (chi2, ndf))
    if chi2 == 0.0:
        return 1.0
    if math.isinf(chi2):
        return 0.0
    a, xx = ndf / 2.0, chi2 / 2.0
    log_pref = -xx + a * math.log(xx) - math.lgamma(a)
    if xx < a + 1.0:
        term = s = 1.0 / a
        k = 1
        while abs(term) > 1e-17 * abs(s):
            term *= xx / (a + k)
            s += term
            k += 1
        return min(1.0, max(0.0, 1.0 - math.exp(log_pref) * s))
    tiny = 1e-300
    b = xx + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 100000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    q = math.exp(log_pref) * h
    # below the smallest normal double scipy's sf reads 0; so does this
    return 0.0 if q < 2.2250738585072014e-308 else min(1.0, q)


def chi2_sf(chi2, ndf):
    """Upper-tail probability of chi2 with ndf dof (scipy if present, else
    `chi2_sf_fallback`)."""
    chi2, ndf = float(chi2), int(ndf)
    try:
        from scipy.stats import chi2 as _c
    except ImportError:                             # pragma: no cover
        return chi2_sf_fallback(chi2, ndf)
    return float(_c.sf(chi2, ndf))


def chi2_diagonal(points, theory):
    return float(sum(((p["value"] - t) / p["err"]) ** 2
                     for p, t in zip(points, theory)))


def chi2_covariance(points, theory):
    """r^T C^-1 r with C = diag(var_uncorr) + sum_k s_k s_k^T; None when the
    points carry no CORR component (then it IS the diagonal chi2)."""
    names = sorted({k for p in points for k in p["corr"]})
    if not names:
        return None
    import numpy as np
    r = np.array([p["value"] - t for p, t in zip(points, theory)])
    s = np.array([[p["corr"].get(k, 0.0) for p in points] for k in names])
    c = np.diag([p["var_uncorr"] for p in points]) + s.T @ s
    z = np.linalg.solve(np.linalg.cholesky(c), r)
    return float(z @ z)


def stats(points, theory):
    n = len(points)
    cd = chi2_diagonal(points, theory)
    cc = chi2_covariance(points, theory)
    pulls = [(t - p["value"]) / p["err"] for p, t in zip(points, theory)]
    worst = max(range(n), key=lambda j: abs(pulls[j]))
    return dict(
        n=n, chi2=cd, p=chi2_sf(cd, n),
        chi2_cov=cc, p_cov=(None if cc is None else chi2_sf(cc, n)),
        mean_ratio=float(sum(t / p["value"] for p, t in zip(points, theory)) / n),
        worst_pull=float(pulls[worst]), worst_x=points[worst]["x"],
        worst_q2=points[worst]["q2"])


# ---------------------------------------------------------------- the tree

def isoscalar_nucleon(sf, x, q2):
    """RULE 2: the tree's per-nucleon F2 of the deuteron, (F2p + F2n)/2 --
    the ONLY place the theory side is built."""
    return 0.5 * (sf.f2p(x, q2) + sf.f2n(x, q2))


def pythia_version():
    """PYTHIA's version string as the xmldoc next to the pdfdata MstwSF read
    says it (the project pins 8.317; record what actually ran)."""
    import lipolgen._lipolgen as _l
    d = _l.pythia8_pdfdata_dir() if hasattr(_l, "pythia8_pdfdata_dir") else ""
    path = os.path.join(os.path.dirname(d.rstrip("/")), "xmldoc", "Version.xml")
    try:
        with open(path) as f:
            for line in f:
                if 'name="Pythia:versionNumber"' in line:
                    return line.split('default="')[1].split('"')[0]
    except OSError:
        pass
    return "unknown"


def build_backend(name):
    """(sf, q2_floor, label, kernel, blocked_reason).  Blocked ONLY for a
    named environment reason; any other error propagates (broken harness)."""
    import lipolgen as lg
    import lipolgen._lipolgen as _l
    if name == "toy":
        return (lg.ToyF2(), None, "ToyF2 (shipped default)",
                lg.default_inclusive_kernel(lg.deuteron()), None)
    if name == "mstw":
        if not getattr(_l, "HAVE_PYTHIA8", False):
            return (None, None, "MstwSF", None, "environment: lipolgen built "
                    "without the PYTHIA 8 tier (MstwSF)")
        d = _l.pythia8_pdfdata_dir()
        if not (d and os.path.isfile(os.path.join(d, "mstw2008lo.00.dat"))):
            return (None, None, "MstwSF", None, "environment: PYTHIA 8 "
                    "pdfdata/mstw2008lo.00.dat not found (MstwSF)")
        sf = lg.MstwSF()
        label = ("MstwSF (MSTW2008 LO, PYTHIA %s pdfdata, i_fit %d)"
                 % (pythia_version(), sf.i_fit))
        return sf, 1.0, label, _kernel_on(sf), None
    if name == "ct18nlo":
        if not getattr(_l, "HAVE_LHAPDF", False):
            return (None, None, "LhapdfSF(%s)" % CT18, None, "environment: "
                    "lipolgen built without the LHAPDF tier (HAVE_LHAPDF "
                    "false)")
        try:
            sf = lg.LhapdfSF(CT18)
        except RuntimeError as exc:
            if "Info file not found for PDF set" not in str(exc):
                raise
            return (None, None, "LhapdfSF(%s)" % CT18, None,
                    "environment: LHAPDF set %s not installed" % CT18)
        return sf, float(sf.q2_min), "LhapdfSF(%s)" % CT18, _kernel_on(sf), None
    raise ValueError("unknown configuration %r" % name)


def _kernel_on(sf):
    import lipolgen as lg
    o = lg.InclusiveKernel.Options()
    o.f2_source = sf
    return lg.InclusiveKernel(lg.deuteron(), o)


def dis_cut(p, q2_floor=None):
    """THE declared DIS cut (module docstring)."""
    lo = Q2_CUT if q2_floor is None else max(Q2_CUT, q2_floor)
    return p["q2"] >= lo and p["w2"] >= W2_CUT


# --------------------------------------------------------------- the row

def measure(path=None):
    pts = load_table(path)
    configs = {}
    for name in CONFIGS:
        sf, floor, label, kern, why = build_backend(name)
        if why is not None:
            configs[name] = dict(status="blocked", reason=why, label=label)
            continue
        theory = [isoscalar_nucleon(sf, p["x"], p["q2"]) for p in pts]
        kdev = max(abs(kern.tables(p["x"], p["q2"]).f2 / t - 1.0)
                   for p, t in zip(pts, theory))
        sel = [j for j, p in enumerate(pts) if dis_cut(p, floor)]
        dis = stats([pts[j] for j in sel], [theory[j] for j in sel])
        configs[name] = dict(
            label=label, q2_floor=max(Q2_CUT, floor or 0.0),
            all=stats(pts, theory), dis=dis, theory=theory,
            kernel_max_rel_diff=float(kdev),
            status="pass" if dis["p"] >= P_MIN else "fail")
    return pts, configs


def _fmt_p(p):
    return "%.3g" % p


def _line(c):
    d, a = c["dis"], c["all"]
    cov = ("n/a" if d["chi2_cov"] is None else
           "%.1f/%d (p = %s)" % (d["chi2_cov"], d["n"], _fmt_p(d["p_cov"])))
    return ("%s: chi2/ndf = %.1f/%d (p = %s) on the DIS cut [diagonal, the "
            "verdict]; full cov %s; all points %.1f/%d (p = %s); <th/data> = "
            "%.3f; worst pull %+.1f at x = %g, Q2 = %g"
            % (c["label"], d["chi2"], d["n"], _fmt_p(d["p"]), cov, a["chi2"],
               a["n"], _fmt_p(a["p"]), d["mean_ratio"], d["worst_pull"],
               d["worst_x"], d["worst_q2"]))


def run(path=None, verbose=True):
    t0 = time.perf_counter()
    pts, configs = measure(path)
    h = configs[HEADLINE]
    others = "; ".join(
        "%s %s" % (n, "BLOCKED (%s)" % configs[n]["reason"]
                   if configs[n]["status"] == "blocked" else
                   "%.1f/%d (p = %s)" % (configs[n]["dis"]["chi2"],
                                         configs[n]["dis"]["n"],
                                         _fmt_p(configs[n]["dis"]["p"])))
        for n in CONFIGS if n != HEADLINE)
    report = dict(
        name=NAME,
        generator="toy (shipped ToyF2): chi2/ndf = %.1f/%d, p = %s on Q2 >= 1, "
                  "W2 >= 4 [diagonal]; sub-rows: %s"
                  % (h["dis"]["chi2"], h["dis"]["n"], _fmt_p(h["dis"]["p"]),
                     others),
        reference="%s, %d of %d points in the DIS cut"
                  % (REFERENCE, h["dis"]["n"], len(pts)),
        tolerance="p(chi2, ndf) >= %.2f on the DIS cut, stat (+) sys diagonal"
                  % P_MIN,
        status=h["status"], headline=HEADLINE, configs=configs,
        npoints=len(pts), points=pts,
        dis_cut="Q2 >= %g GeV2 and W2 >= %g GeV2 (and >= the backend's grid "
                "floor)" % (Q2_CUT, W2_CUT),
    )
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        print("NMC F2d per nucleon vs (F2p + F2n)/2 -- %d points, DIS cut %s"
              % (len(pts), report["dis_cut"]))
        cols = [n for n in CONFIGS if configs[n]["status"] != "blocked"]
        print("     x       Q2     NMC F2d +- diag   " +
              "  ".join("%10s" % n for n in cols))
        for j, p in enumerate(pts):
            print("%8.4f %7.2f   %.4f +- %.4f %s " % (
                p["x"], p["q2"], p["value"], p["err"],
                " " if dis_cut(p) else "*")
                + "  ".join("%10.4f" % configs[n]["theory"][j] for n in cols))
        print("(* = outside the DIS cut)")
        for n in cols:
            print("kernel check [%s]: max |InclusiveKernel(deuteron).tables().f2"
                  " / isoscalar - 1| = %.3g" % (n, configs[n]["kernel_max_rel_diff"]))
        for n in CONFIGS:
            c = configs[n]
            if c["status"] == "blocked":
                print("SUBROW | %s[%s] | not evaluated (%s) | %s | %s | BLOCKED "
                      "(%s)" % (NAME, n, c["label"], REFERENCE,
                               report["tolerance"], c["reason"]))
            else:
                print("SUBROW | %s[%s] | %s | %s | %s | %s" % (
                    NAME, n, _line(c), REFERENCE, report["tolerance"],
                    c["status"].upper()))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, report["generator"], report["reference"],
            report["tolerance"], report["status"].upper()))
    return report


def main(argv=None):
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    try:
        run()
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
