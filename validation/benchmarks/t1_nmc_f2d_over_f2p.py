#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T1 -- NMC F2d/F2p (survey row F-2) against the tree's per-nucleon
deuteron-to-proton ratio (tier T1, BENCHMARK_PLAN.md sec. 4 "then, by tier").

WHAT IS COMPARED.  NMC's ratio F2d/F2p of PER-NUCLEON structure functions --
M. Arneodo et al., "Accurate measurement of F2d/F2p and Rd - Rp", Nucl.
Phys. B 487 (1997) 3, hep-ex/9611022; HEPData ins426595 tables 2-21 (CC0),
260 points, 0.0015 <= x <= 0.675, 0.16 <= Q2 <= 99 GeV2 -- vendored VERBATIM
from the nnpdf-data 4.1.5 wheel (set NMC_NC_NOTFIXED) by
`vendor_nnpdf_commondata.py` into `data/nnpdf_NMC_NC_NOTFIXED.json` (its
`_provenance` block: wheel URL and sha256, paths inside the wheel, HEPData /
INSPIRE record, licences).  Against it, at each point's OWN (x, Q2):

    R_tree(x, Q2) = [(F2p + F2n)/2] / F2p = (1 + F2n/F2p)/2     (THE VERDICT)

and, RECORDED only, the same ratio through the backend's single-argument
hook `f2n_over_f2p(x)`, which is Q2-FROZEN (06_critic.md sec. 5.1: "the trap
is elsewhere: ToyF2::f2n_over_f2p is Q2-frozen at Q2 = 10"): (1 + hook(x))/2.
For ToyF2 the two are the same function (its F2n IS F2p times the hook); for
MstwSF / LhapdfSF they differ, and the harness prices the difference.
Three backends (the `--unpol-sf` choices):

  toy      ToyF2() -- the SHIPPED default; the numerator is exactly
           `default_inclusive_kernel(deuteron()).tables(x, Q2).f2`.
           THE REPORT ROW'S VERDICT IS THIS CONFIGURATION'S.
  mstw     MstwSF() -- MSTW2008 LO from PYTHIA 8's pdfdata.  Sub-row.
  ct18nlo  LhapdfSF("CT18NLO").  Sub-row; BLOCKED (environment) when LHAPDF
           reports the set not installed -- reported, never dropped.

ZERO free parameters, NO tuning.

WHICH POINTS.  ALL points, and THE DIS CUT, declared here before any chi2
of this harness was computed: Q2 >= 1 GeV2 and W2 = M^2 + Q2 (1 - x)/x >=
4 GeV2 (M = 0.938272 GeV), and Q2 >= the backend's grid floor where it has
one (CT18NLO 1.677; MSTW2008's grid starts at 1 GeV2, so it coincides; ToyF2
has no grid).  The W2 cut removes no point here (min W2 = 4.27 GeV2 above
Q2 = 1); 211 of 260 points pass for toy and mstw.

ERRORS AND WHICH CHI2 IS THE VERDICT.  The file carries two uncertainty
variants, both vendored as given:
  'hepdata' (PRIMARY): stat (ADD, UNCORR) + one 'sys' (MULT, CORR) -- the
      HEPData tables.
  'legacy':  the same stat + five ADD/CORR components sys_corr_1..5, NNPDF's
      port of its pre-HEPData commondata; their quadrature sum matches the
      hepdata 'sys' to <= 6.3e-4 absolute.  The five sources are not named
      in the file (UNVERIFIED).
THE VERDICT is the DIAGONAL chi2 on the 'hepdata' variant, stat (+) sys in
quadrature per point (the systematic treated as uncorrelated).  RECORDED,
never in the verdict: the full-covariance chi2, C = diag(UNCORR^2) +
sum_k s_k s_k^T over CORR components, on 'hepdata' (one component labelled
CORR, i.e. 100 % correlated over all points -- NNPDF's label, not NMC's
statement) AND on 'legacy' (five components; the more physical covariance,
but of unverified origin).  Uncertainties are the ABSOLUTE values given (no
t0 rescaling of the MULT component).

TOLERANCE.  p(chi2, ndf) >= 0.01 on the DIS-cut set, diagonal 'hepdata'
errors, ndf = number of points.  The row's STATUS is the shipped default's
(toy).  The threshold follows row 4 (`t3_nmc_li6_over_d.py`) and was chosen
AFTER the plan scout's previews had been seen (scout, 2026-09-26, all 260
points, diagonal: ToyF2 chi2 = 2252, MSTW 298 Q2-dependent and 310 with the
frozen hook); stated here so that it can be argued with, not tuned.  A FAIL
of the toy is RECORDED with its distance and moves nothing.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations / A-scaling -- NOT APPLICABLE: no b1; no A-scaling
     assumption (no 0.5, no 2/6) anywhere.
  2. isoscalar denominator -- APPLIED AND ASSERTED, in its ratio form: the
     NUMERATOR is the ISOSCALAR nucleon (F2p + F2n)/2 of the free nucleon,
     built in ONE place (`isoscalar_nucleon`) and identical to the kernel's
     per-nucleon deuteron F2 (`run()` records the largest relative
     difference; python/tests/test_bench_t1_nucleon.py asserts both); the
     denominator is the free proton, because that is what NMC divided by.
  3. per nucleus / per nucleon -- a RATIO of PER-NUCLEON structure functions
     on both sides: ~1 at small x (0.98 at x = 0.0015), which a per-deuteron
     numerator would put near 2.  Not like for like in one respect, stated
     rather than corrected: NMC's numerator is the DEUTERON's F2 (nuclear
     effects included; not corrected -- UNVERIFIED offline), the tree's the
     free isoscalar nucleon.

WHAT THIS DOES NOT VALIDATE.  F2n/F2p itself at x -> 1 (the deuteron's
nuclear effects are not unfolded, and the data stop at x = 0.675), Rd - Rp
(the paper's other result is not in the package), and anything polarized.

INTEGRITY.  The sha256 of every byte of the vendored file (`FILE_SHA256`)
and of its canonical body (`BODY_SHA256`, also recorded inside the file)
are re-checked on EVERY read: an altered table is refused with a
RuntimeError, so `__main__` exits 2, never 0 on moved numbers.

Run:  source env.sh && python3 validation/benchmarks/t1_nmc_f2d_over_f2p.py

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
DATA = os.path.join(HERE, "data", "nnpdf_NMC_NC_NOTFIXED.json")
#: sha256 of every byte of the vendored file and of its canonical body, as
#: written by vendor_nnpdf_commondata.py from the pinned wheel (2026-09-26).
FILE_SHA256 = "47e833cfb399eb984455e09ed37d8833edd20445607cefe0a103545cd8e43900"
BODY_SHA256 = "fc903cccb05c5d714851af25c6dc2eacd45526d7dda72a7248978156d73da5cf"

NAME = "t1_nmc_f2d_over_f2p"
P_MIN = 0.01             # tolerance: p(chi2, ndf) >= 0.01, DIS cut, diagonal
Q2_CUT = 1.0             # GeV^2 -- the declared DIS cut ...
W2_CUT = 4.0             # GeV^2 -- ... and its W^2 half
M_NUCLEON = 0.938272     # GeV, for W^2 only
CONFIGS = ("toy", "mstw", "ct18nlo")
HEADLINE = "toy"
PRIMARY, LEGACY = "hepdata", "legacy"
CT18 = "CT18NLO"
REFERENCE = ("NMC F2d/F2p (per nucleon), NPB 487 (1997) 3, HEPData ins426595 "
             "tables 2-21 via nnpdf-data 4.1.5 (NMC_NC_NOTFIXED)")


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
    ({name: signed absolute value} of the CORR components) -- of ONE
    uncertainty variant (default: the file's primary, 'hepdata')."""
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


def stats(points, theory, legacy=None):
    n = len(points)
    cd = chi2_diagonal(points, theory)
    cc = chi2_covariance(points, theory)
    pulls = [(t - p["value"]) / p["err"] for p, t in zip(points, theory)]
    worst = max(range(n), key=lambda j: abs(pulls[j]))
    out = dict(
        n=n, chi2=cd, p=chi2_sf(cd, n),
        chi2_cov=cc, p_cov=(None if cc is None else chi2_sf(cc, n)),
        mean_diff=float(sum(t - p["value"] for p, t in zip(points, theory)) / n),
        worst_pull=float(pulls[worst]), worst_x=points[worst]["x"],
        worst_q2=points[worst]["q2"])
    if legacy is not None:
        cl = chi2_covariance(legacy, theory)
        out.update(chi2_legacy_diag=chi2_diagonal(legacy, theory),
                   chi2_legacy_cov=cl, p_legacy_cov=chi2_sf(cl, n))
    return out


# ---------------------------------------------------------------- the tree

def isoscalar_nucleon(sf, x, q2):
    """RULE 2: the tree's per-nucleon F2 of the deuteron, (F2p + F2n)/2 --
    the ONLY place the numerator is built."""
    return 0.5 * (sf.f2p(x, q2) + sf.f2n(x, q2))


def ratio_q2(sf, x, q2):
    """THE VERDICT's theory: the isoscalar nucleon over the free proton at
    the point's own (x, Q2)."""
    return isoscalar_nucleon(sf, x, q2) / sf.f2p(x, q2)


def ratio_frozen(sf, x):
    """RECORDED only: the same ratio through the Q2-frozen single-argument
    hook `f2n_over_f2p(x)`."""
    return 0.5 * (1.0 + sf.f2n_over_f2p(x))


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


def _kernel_on(sf):
    import lipolgen as lg
    o = lg.InclusiveKernel.Options()
    o.f2_source = sf
    return lg.InclusiveKernel(lg.deuteron(), o)


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


def dis_cut(p, q2_floor=None):
    """THE declared DIS cut (module docstring)."""
    lo = Q2_CUT if q2_floor is None else max(Q2_CUT, q2_floor)
    return p["q2"] >= lo and p["w2"] >= W2_CUT


# --------------------------------------------------------------- the row

def measure(path=None):
    pts = load_table(path, PRIMARY)
    leg = load_table(path, LEGACY)
    configs = {}
    for name in CONFIGS:
        sf, floor, label, kern, why = build_backend(name)
        if why is not None:
            configs[name] = dict(status="blocked", reason=why, label=label)
            continue
        th = [ratio_q2(sf, p["x"], p["q2"]) for p in pts]
        thf = [ratio_frozen(sf, p["x"]) for p in pts]
        kdev = max(abs(kern.tables(p["x"], p["q2"]).f2
                       / isoscalar_nucleon(sf, p["x"], p["q2"]) - 1.0)
                   for p in pts)
        sel = [j for j, p in enumerate(pts) if dis_cut(p, floor)]
        sp, sl = [pts[j] for j in sel], [leg[j] for j in sel]
        dis = stats(sp, [th[j] for j in sel], sl)
        dis_frozen = stats(sp, [thf[j] for j in sel], sl)
        configs[name] = dict(
            label=label, q2_floor=max(Q2_CUT, floor or 0.0),
            all=stats(pts, th, leg), dis=dis,
            all_frozen=stats(pts, thf, leg), dis_frozen=dis_frozen,
            theory=th, theory_frozen=thf,
            max_abs_frozen_minus_q2=float(max(abs(a - b) for a, b in zip(thf, th))),
            kernel_max_rel_diff=float(kdev),
            status="pass" if dis["p"] >= P_MIN else "fail")
    return pts, configs


def _fmt_p(p):
    return "%.3g" % p


def _line(c):
    d, a, f = c["dis"], c["all"], c["dis_frozen"]
    return ("%s: chi2/ndf = %.1f/%d (p = %s) on the DIS cut [diagonal, "
            "hepdata, Q2-dependent: the verdict]; full cov hepdata %.1f (p = "
            "%s), legacy 5-component %.1f (p = %s); frozen f2n_over_f2p hook "
            "%.1f/%d (p = %s); all points %.1f/%d (p = %s); <th - data> = "
            "%+.4f; worst pull %+.1f at x = %g, Q2 = %g"
            % (c["label"], d["chi2"], d["n"], _fmt_p(d["p"]), d["chi2_cov"],
               _fmt_p(d["p_cov"]), d["chi2_legacy_cov"],
               _fmt_p(d["p_legacy_cov"]), f["chi2"], f["n"], _fmt_p(f["p"]),
               a["chi2"], a["n"], _fmt_p(a["p"]), d["mean_diff"],
               d["worst_pull"], d["worst_x"], d["worst_q2"]))


def run(path=None, verbose=True):
    t0 = time.perf_counter()
    pts, configs = measure(path)
    h = configs[HEADLINE]
    others = "; ".join(
        "%s %s" % (n, "BLOCKED (%s)" % configs[n]["reason"]
                   if configs[n]["status"] == "blocked" else
                   "%.1f/%d (p = %s; frozen hook %.1f)"
                   % (configs[n]["dis"]["chi2"], configs[n]["dis"]["n"],
                      _fmt_p(configs[n]["dis"]["p"]),
                      configs[n]["dis_frozen"]["chi2"]))
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
        print("NMC F2d/F2p vs (1 + F2n/F2p)/2 -- %d points, DIS cut %s"
              % (len(pts), report["dis_cut"]))
        cols = [n for n in CONFIGS if configs[n]["status"] != "blocked"]
        print("     x       Q2     NMC d/p +- diag   " +
              "  ".join("%10s %10s" % (n, n + "@hook") for n in cols))
        for j, p in enumerate(pts):
            print("%8.4f %7.2f   %.4f +- %.4f %s " % (
                p["x"], p["q2"], p["value"], p["err"],
                " " if dis_cut(p) else "*")
                + "  ".join("%10.4f %10.4f" % (configs[n]["theory"][j],
                                               configs[n]["theory_frozen"][j])
                            for n in cols))
        print("(* = outside the DIS cut; @hook = through the Q2-frozen "
              "f2n_over_f2p(x))")
        for n in cols:
            print("kernel check [%s]: max |InclusiveKernel(deuteron).tables().f2"
                  " / isoscalar - 1| = %.3g; max |hook - Q2-dependent| = %.4f"
                  % (n, configs[n]["kernel_max_rel_diff"],
                     configs[n]["max_abs_frozen_minus_q2"]))
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
