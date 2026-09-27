#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T1 -- world g1 of the deuteron (survey rows D-2, D-1, D-6; SMC recorded)
against the tree's per-nucleon deuteron g1 (tier T1, BENCHMARK_PLAN.md
sec. 4 "then, by tier").

WHAT IS COMPARED.  Four published g1^d sets, each vendored VERBATIM from the
nnpdf-data 4.1.5 wheel by `vendor_nnpdf_commondata.py` (every file's
`_provenance` block carries the wheel URL and sha256, the paths inside the
wheel, the HEPData / INSPIRE record and the licences -- HEPData CC0, package
GPL-3.0-or-later):

  compass  COMPASS final, PLB 769 (2017) 34 (survey D-2), HEPData ins1501480
           table 1, 15 points, 0.0046 <= x <= 0.567, 1.1 <= Q2 <= 60.8;
           `data/nnpdf_COMPASS15_NC_NOTFIXED_MUD.json`.
  hermes   HERMES, PRD 75 (2007) 012007 (survey D-1), HEPData ins726689
           table 13, 15 points, 0.0264 <= x <= 0.7248, 1.12 <= Q2 <= 12.21;
           `data/nnpdf_HERMES_NC_7GEV_ED.json`.
  e143     E143 (survey D-6: PRD 58 (1998) 112003; the metadata's arXiv
           field says hep-ex/9705012 -- a conflict not resolvable offline),
           HEPData ins467140 table 18, 28 points, 0.031 <= x <= 0.749,
           1.27 <= Q2 <= 9.52; `data/nnpdf_E143_NC_NOTFIXED_ED.json`.
  smc      SMC, HEPData ins471981 table 7, 13 points -- its paper is NOT
           identified offline and the survey has no SMC row, so it is
           RECORDED ONLY and never enters the verdict;
           `data/nnpdf_SMC_NC_NOTFIXED_MUD.json`.

The tree side is the generator's OWN per-nucleon deuteron g1,

    g1^d_tree(x, Q2) = PolSF.g1_nucleus(deuteron(), x, Q2) / 2
                     = (1 - 1.5 P_D) (g1p + g1n)/2,   P_D = P_D_DEUTERON,

with P_D_DEUTERON = 0.045 (beams.hpp) -- exactly what
`InclusiveKernel(deuteron()).tables(x, Q2).g1` returns
(xsec.cpp `g1a` = g1_nucleus / A; beams.hpp: the deuteron's eff_pol_p =
eff_pol_n = DEUTERON_VECTOR_POLARIZATION = vector_dilution_of(P_D_DEUTERON)
= 0.9325).  THIS is the verdict's theory.  RECORDED beside it: the
NO-D-STATE variant (g1p + g1n)/2, i.e. the same nucleons without the
(1 - 1.5 P_D) depolarisation -- the 7 % that 06_critic.md sec. 5.1 (line
524) warns is applied twice or not at all.  Four polarized backends:

  toy       ToyG1() on ToyF2 -- the SHIPPED default (`--pol-sf toy`,
            `--unpol-sf toy`).  THE REPORT ROW'S VERDICT IS THIS ONE'S.
  mstw      ToyG1(MstwSF()) -- the g1 of a `--unpol-sf mstw --pol-sf toy`
            run (the kernel builds its ToyG1 on the run's F2 source): the
            toy A1 times F1 from MSTW2008 LO.  Sub-row.
  ct18nlo   ToyG1(LhapdfSF("CT18NLO")) -- the g1 of `--unpol-sf ct18nlo`.
            Sub-row; BLOCKED (environment) when the set is not installed.
  nnpdfpol  LhapdfG1("NNPDFpol11_100") (`--pol-sf nnpdfpol`).  Sub-row;
            BLOCKED (environment) when the set is not installed.

A blocked configuration is reported with its reason, never dropped, never
faked.  ZERO free parameters, NO tuning.

THE D-STATE CONVENTION OF THE DATA IS UNVERIFIED.  The metadata of none of
the four sets says whether the published g1d is corrected for the D state,
and the papers cannot be read here.  06_critic.md sec. 5.1: "published g1d
is conventionally NOT corrected for the D state; the (1 - 1.5 omega_D)
factor is applied when extracting g1n."  IF that holds, g1d_tree (with the
factor) is the like-for-like comparison, which is why it carries the
verdict; the no-D-state variant prices the other reading.  The tree's
P_D = 0.045 is not necessarily the omega_D each experiment used (not read
offline).

WHICH POINTS.  ALL points, and THE DIS CUT, declared here before any chi2
of this harness was computed: Q2 >= 1 GeV2 and W2 = M^2 + Q2 (1 - x)/x >=
4 GeV2 (M = 0.938272 GeV), and Q2 >= the backend's grid floor where it has
one (`LhapdfSF.q2_min()` of CT18NLO / NNPDFpol11_100; MSTW2008 starts at 1).
The cut removes only SMC's x = 0.002, Q2 = 0.5 point; COMPASS 15, HERMES 15,
E143 28 and SMC 12 points remain (min W2 = 4.07 GeV2, E143).

ERRORS AND WHICH CHI2 IS THE VERDICT.  PRIMARY (the verdict): the DIAGONAL
chi2, every uncertainty component of a point in quadrature (stat (+) sys,
correlations dropped).  RECORDED beside it, never in the verdict: the
FULL-COVARIANCE chi2, C = diag(sum of UNCORR^2) + sum_k s_k s_k^T over the
components NNPDF labels CORR -- for HERMES the 15 'artificial correlated
statistical' vectors sys_0..sys_14 (the statistical covariance; stat = 0 in
every bin, so the per-point statistical error is sqrt(sum_k sys_k^2)) plus
'param' and 'evol'; for E143 the 4.9 % beam-normalisation 'sys_beam' (MULT,
CORR).  COMPASS and SMC carry no CORR component: their full-covariance chi2
IS the diagonal one (reported as n/a).  For HERMES the full covariance is
the more faithful number -- the harness flags whenever the two treatments
disagree on the verdict (`verdict_differs_under_full_cov`) rather than pick
the kinder one.  Uncertainties are the ABSOLUTE values given (no t0
rescaling of MULT components).

TOLERANCE.  p(chi2, ndf) >= 0.01 on the DIS-cut set, diagonal errors,
ndf = number of points (nothing fitted), for EACH of compass, hermes and
e143; the row passes only if all three do.  SMC is recorded, never in the
verdict (its paper is not identified).  The combined three-experiment chi2
is recorded, not asserted.  The row's STATUS is the shipped default's
(toy); the other configurations are sub-rows with their own status.  The
threshold follows row 4 (`t3_nmc_li6_over_d.py`) and was chosen AFTER the
plan scout's previews had been seen (scout, 2026-09-26: ToyG1 chi2 ~ 108/15
on COMPASS, driven by x < 0.01 where the toy gives -0.96 against
-0.13 +- 0.20); stated here so that it can be argued with, not tuned.
ToyG1 is labelled a TOY in sf.hpp ("toy A1(x) shapes ... replace with
JAM/DSSV"): its failure is a RECORDED FAIL with its distance, and moves
nothing.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. normalisations and A-scaling, both sides named (the rule is written for
     b1; g1d has the same trap, 06_critic.md sec. 5.1).  Data: PER NUCLEON
     (NNPDF's y_label g_{1,N}; 06 sec. 5.1; the papers not re-read).  Tree:
     PER NUCLEON, g1_nucleus(deuteron)/A with A = 2 -- a per-DEUTERON g1
     would be twice it.  No A-scaling assumption: no 2/6 (the 6Li inert-alpha
     assumption is NOT applied), no 0.5 beyond the explicit 1/A.  The D-state
     factor is named above and varied.
  2. isoscalar denominator -- APPLIED AND ASSERTED in its g1 form: the tree
     side is the isoscalar nucleon combination (g1p + g1n)/2 times the
     deuteron's vector dilution, built in ONE place (`g1d_tree`) through
     the tree's own `g1_nucleus`; python/tests/test_bench_t1_nucleon.py
     asserts g1d_tree = 0.9325 (g1p + g1n)/2 = kernel tables().g1.  No
     nuclear PDF grid enters.
  3. per nucleus / per nucleon -- structure functions PER NUCLEON on both
     sides.  Not like for like, stated rather than corrected: beyond the
     D-state depolarisation the tree has no deuteron-structure correction
     for g1 (Fermi motion, shadowing), and COMPASS and SMC measured on 6LiD
     with their own dilution models (survey 02 sec. 2).

WHAT THIS DOES NOT VALIDATE.  The tensor sector (g1d is blind to b1), the
effective polarisations of 6Li / 7Li (`eff_pol_p` / `eff_pol_n` of LI6(),
LI7()), g2 (no g2 data is vendored), the proton or neutron g1 separately
(only their D-state-weighted sum enters), and the Bjorken sum.

INTEGRITY.  The sha256 of every byte of each vendored file and of its
canonical body (also recorded inside the file) are re-checked on EVERY
read: an altered table is refused with a RuntimeError, so `__main__` exits
2, never 0 on moved numbers.

Run:  source env.sh && python3 validation/benchmarks/t1_g1d_world.py

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
DATA_DIR = os.path.join(HERE, "data")

NAME = "t1_g1d_world"
P_MIN = 0.01             # tolerance: p(chi2, ndf) >= 0.01, DIS cut, diagonal
Q2_CUT = 1.0             # GeV^2 -- the declared DIS cut ...
W2_CUT = 4.0             # GeV^2 -- ... and its W^2 half
M_NUCLEON = 0.938272     # GeV, for W^2 only
CONFIGS = ("toy", "mstw", "ct18nlo", "nnpdfpol")
HEADLINE = "toy"
CT18 = "CT18NLO"
NNPDFPOL = "NNPDFpol11_100"

#: name -> (file, sha256 of every byte, sha256 of the canonical body,
#: in the verdict?, reference label).  Digests as written by
#: vendor_nnpdf_commondata.py from the pinned wheel (2026-09-26).
EXPERIMENTS = {
    "compass": ("nnpdf_COMPASS15_NC_NOTFIXED_MUD.json",
                "35ef52704ee46f4f7cd49803fc452e3a759a10fa91547f4378307e16b45c7892",
                "d961b0a9f6708d9f2e074c92047517793ffba0bcfd964ab1e2275d016caf05bd",
                True, "COMPASS PLB 769 (2017) 34, HEPData ins1501480 t1"),
    "hermes": ("nnpdf_HERMES_NC_7GEV_ED.json",
               "5b397977b8937866a92bf02693d1e7103b1c871933be0b9df32e3e4f5e0cc0e8",
               "1f6cb0d83fb513cf93dfd756f473d9575b0cd431789f359c2b120008f9d62970",
               True, "HERMES PRD 75 (2007) 012007, HEPData ins726689 t13"),
    "e143": ("nnpdf_E143_NC_NOTFIXED_ED.json",
             "5b5d7faefef49e3b8c3c17cdac3f277921799c2b3d42c721eb7366eccb623915",
             "6509786def09ecb318d45cb750e27763ba19c9a3968dd48811806f9cceb5861d",
             True, "E143 (PRD 58 (1998) 112003 per survey D-6), HEPData "
                   "ins467140 t18"),
    "smc": ("nnpdf_SMC_NC_NOTFIXED_MUD.json",
            "8543a2919d8d9b6be0d31b77bf701679498499e054b876149655e01a0911c4fb",
            "78683115a9c868fc424e2074c918a79a0369a3622b8fa2800c4adbe5391d4673",
            False, "SMC, HEPData ins471981 t7 (paper not identified; "
                   "recorded only)"),
}
VERDICT_EXPERIMENTS = tuple(k for k, v in EXPERIMENTS.items() if v[3])
REFERENCE = ("g1d per nucleon: COMPASS (D-2), HERMES (D-1), E143 (D-6) via "
             "nnpdf-data 4.1.5; SMC recorded only")


# ---------------------------------------------------------------- reference

def data_path(exp):
    return os.path.join(DATA_DIR, EXPERIMENTS[exp][0])


def read_vendored(exp, path=None):
    """The vendored document of one experiment -- RuntimeError unless the
    file's bytes AND its canonical body hash to the recorded digests."""
    _, file_sha, body_sha, _, _ = EXPERIMENTS[exp]
    path = path or data_path(exp)
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != file_sha:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, file_sha))
    doc = json.loads(raw.decode("ascii"))
    canon = json.dumps(doc["body"], sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False)
    body = hashlib.sha256(canon.encode("ascii")).hexdigest()
    if not body == doc["_body_sha256"] == body_sha:
        raise RuntimeError("%s: body sha256 %s; file records %s, harness %s"
                           % (path, body, doc["_body_sha256"], body_sha))
    return doc


def load_table(exp, path=None):
    """The points as dicts: x, q2, value, err (every component in
    quadrature: the DIAGONAL error), var_uncorr (sum of UNCORR^2) and corr
    ({name: signed absolute value} of the CORR components)."""
    body = read_vendored(exp, path)["body"]
    variant = body["primary_uncertainty_variant"]
    defs = body["uncertainty_variants"][variant]["definitions"]
    for k, d in defs.items():
        if d["type"] not in ("UNCORR", "CORR"):
            raise RuntimeError("%s: %s uncertainty %s has type %r (only UNCORR"
                               " and CORR are understood)"
                               % (NAME, exp, k, d["type"]))
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
        worst_pull=float(pulls[worst]), worst_x=points[worst]["x"],
        worst_q2=points[worst]["q2"])


# ---------------------------------------------------------------- the tree

def g1d_tree(pol, deuteron, x, q2):
    """THE VERDICT's theory: the tree's per-nucleon deuteron g1,
    g1_nucleus(deuteron)/A -- (1 - 1.5 P_D) (g1p + g1n)/2."""
    return pol.g1_nucleus(deuteron, x, q2) / float(deuteron.A)


def g1d_no_dstate(pol, x, q2):
    """RECORDED only: the isoscalar nucleon (g1p + g1n)/2 WITHOUT the
    deuteron's D-state depolarisation."""
    return 0.5 * (pol.g1p(x, q2) + pol.g1n(x, q2))


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


def _lhapdf(setname, make):
    """make() or the blocked reason when LHAPDF reports `setname` missing;
    any other error propagates (a broken harness, not a blocked row)."""
    import lipolgen._lipolgen as _l
    if not getattr(_l, "HAVE_LHAPDF", False):
        return None, ("environment: lipolgen built without the LHAPDF tier "
                      "(HAVE_LHAPDF false)")
    try:
        return make(), None
    except RuntimeError as exc:
        if "Info file not found for PDF set" not in str(exc):
            raise
        return None, "environment: LHAPDF set %s not installed" % setname


def _kernel(unpol=None, pol=None):
    import lipolgen as lg
    o = lg.InclusiveKernel.Options()
    if unpol is not None:
        o.f2_source = unpol
    if pol is not None:
        o.g1_model = pol
    return lg.InclusiveKernel(lg.deuteron(), o)


def build_backend(name):
    """(pol, q2_floor, label, kernel, blocked_reason)."""
    import lipolgen as lg
    import lipolgen._lipolgen as _l
    if name == "toy":
        return (lg.ToyG1(), None, "ToyG1 on ToyF2 (shipped default)",
                lg.default_inclusive_kernel(lg.deuteron()), None)
    if name == "mstw":
        if not getattr(_l, "HAVE_PYTHIA8", False):
            return (None, None, "ToyG1(MstwSF)", None, "environment: lipolgen "
                    "built without the PYTHIA 8 tier (MstwSF)")
        d = _l.pythia8_pdfdata_dir()
        if not (d and os.path.isfile(os.path.join(d, "mstw2008lo.00.dat"))):
            return (None, None, "ToyG1(MstwSF)", None, "environment: PYTHIA 8 "
                    "pdfdata/mstw2008lo.00.dat not found (MstwSF)")
        sf = lg.MstwSF()
        return (lg.ToyG1(sf), 1.0, "ToyG1 on MstwSF (MSTW2008 LO, PYTHIA %s "
                "pdfdata; the --unpol-sf mstw g1)" % pythia_version(),
                _kernel(unpol=sf), None)
    if name == "ct18nlo":
        sf, why = _lhapdf(CT18, lambda: lg.LhapdfSF(CT18))
        if why:
            return None, None, "ToyG1(LhapdfSF(%s))" % CT18, None, why
        return (lg.ToyG1(sf), float(sf.q2_min), "ToyG1 on LhapdfSF(%s)" % CT18,
                _kernel(unpol=sf), None)
    if name == "nnpdfpol":
        pol, why = _lhapdf(NNPDFPOL, lambda: lg.LhapdfG1(NNPDFPOL))
        if why:
            return None, None, "LhapdfG1(%s)" % NNPDFPOL, None, why
        floor = float(lg.LhapdfSF(NNPDFPOL).q2_min)    # the grid's own edge
        return (pol, floor, "LhapdfG1(%s)" % NNPDFPOL, _kernel(pol=pol), None)
    raise ValueError("unknown configuration %r" % name)


def dis_cut(p, q2_floor=None):
    """THE declared DIS cut (module docstring)."""
    lo = Q2_CUT if q2_floor is None else max(Q2_CUT, q2_floor)
    return p["q2"] >= lo and p["w2"] >= W2_CUT


# --------------------------------------------------------------- the row

def measure(paths=None):
    import lipolgen as lg
    paths = paths or {}
    deut = lg.deuteron()
    tables = {e: load_table(e, paths.get(e)) for e in EXPERIMENTS}
    configs = {}
    for name in CONFIGS:
        pol, floor, label, kern, why = build_backend(name)
        if why is not None:
            configs[name] = dict(status="blocked", reason=why, label=label)
            continue
        per_exp, kdev = {}, 0.0
        for e, pts in tables.items():
            th = [g1d_tree(pol, deut, p["x"], p["q2"]) for p in pts]
            th0 = [g1d_no_dstate(pol, p["x"], p["q2"]) for p in pts]
            for p, t in zip(pts, th):
                kdev = max(kdev, abs(kern.tables(p["x"], p["q2"]).g1 - t)
                           / max(abs(t), 1e-300))
            sel = [j for j, p in enumerate(pts) if dis_cut(p, floor)]
            dis = stats([pts[j] for j in sel], [th[j] for j in sel])
            dis0 = stats([pts[j] for j in sel], [th0[j] for j in sel])
            ok = dis["p"] >= P_MIN
            ok_cov = ok if dis["p_cov"] is None else dis["p_cov"] >= P_MIN
            per_exp[e] = dict(
                all=stats(pts, th), dis=dis, dis_no_dstate=dis0,
                theory=th, theory_no_dstate=th0, in_verdict=EXPERIMENTS[e][3],
                status="pass" if ok else "fail",
                verdict_differs_under_full_cov=(ok != ok_cov))
        world = [per_exp[e]["dis"] for e in VERDICT_EXPERIMENTS]
        wc, wn = sum(s["chi2"] for s in world), sum(s["n"] for s in world)
        configs[name] = dict(
            label=label, q2_floor=max(Q2_CUT, floor or 0.0),
            experiments=per_exp, kernel_max_rel_diff=float(kdev),
            world_chi2=wc, world_ndf=wn, world_p=chi2_sf(wc, wn),
            status="pass" if all(per_exp[e]["status"] == "pass"
                                 for e in VERDICT_EXPERIMENTS) else "fail")
    dstate = dict(p_d=float(lg.P_D_DEUTERON),
                  factor=float(lg.DEUTERON_VECTOR_POLARIZATION),
                  vector_dilution_of_p_d=float(lg.vector_dilution_of(
                      lg.P_D_DEUTERON)))
    return tables, configs, dstate


def _fmt_p(p):
    return "%.3g" % p


def _exp_line(c, e):
    s = c["experiments"][e]
    d, a, d0 = s["dis"], s["all"], s["dis_no_dstate"]
    cov = ("n/a (no CORR component)" if d["chi2_cov"] is None else
           "%.1f (p = %s)" % (d["chi2_cov"], _fmt_p(d["p_cov"])))
    return ("%s: chi2/ndf = %.1f/%d (p = %s) on the DIS cut [diagonal, "
            "g1_nucleus(d)/2: the verdict%s]; full cov %s%s; no-D-state "
            "(g1p+g1n)/2 %.1f (p = %s); all points %.1f/%d; worst pull %+.1f "
            "at x = %g, Q2 = %g"
            % (c["label"], d["chi2"], d["n"], _fmt_p(d["p"]),
               "" if s["in_verdict"] else " -- RECORDED ONLY", cov,
               " -- VERDICT DIFFERS UNDER FULL COV"
               if s["verdict_differs_under_full_cov"] else "",
               d0["chi2"], _fmt_p(d0["p"]), a["chi2"], a["n"],
               d["worst_pull"], d["worst_x"], d["worst_q2"]))


def run(paths=None, verbose=True):
    t0 = time.perf_counter()
    tables, configs, dstate = measure(paths)
    h = configs[HEADLINE]

    def summary(c):
        return ", ".join("%s %.1f/%d (p = %s)" % (
            e, c["experiments"][e]["dis"]["chi2"],
            c["experiments"][e]["dis"]["n"],
            _fmt_p(c["experiments"][e]["dis"]["p"]))
            for e in VERDICT_EXPERIMENTS)

    others = "; ".join(
        "%s %s" % (n, "BLOCKED (%s)" % configs[n]["reason"]
                   if configs[n]["status"] == "blocked" else
                   "%s %s" % (configs[n]["status"].upper(),
                              summary(configs[n])))
        for n in CONFIGS if n != HEADLINE)
    report = dict(
        name=NAME,
        generator="toy (shipped ToyG1, g1_nucleus(d)/2): %s on Q2 >= 1, "
                  "W2 >= 4 [diagonal]; sub-rows: %s" % (summary(h), others),
        reference="%s; %s points in the DIS cut" % (
            REFERENCE, " + ".join(
                "%d" % h["experiments"][e]["dis"]["n"]
                for e in VERDICT_EXPERIMENTS)),
        tolerance="p(chi2, ndf) >= %.2f on the DIS cut for each of %s, stat "
                  "(+) sys diagonal" % (P_MIN, ", ".join(VERDICT_EXPERIMENTS)),
        status=h["status"], headline=HEADLINE, configs=configs,
        dstate=dstate, tables=tables,
        dis_cut="Q2 >= %g GeV2 and W2 >= %g GeV2 (and >= the backend's grid "
                "floor)" % (Q2_CUT, W2_CUT),
    )
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        print("g1d per nucleon vs g1_nucleus(deuteron)/2 = %.4f (g1p + g1n)/2 "
              "(P_D = %.3f); DIS cut %s" % (dstate["factor"], dstate["p_d"],
                                            report["dis_cut"]))
        cols = [n for n in CONFIGS if configs[n]["status"] != "blocked"]
        for e, pts in tables.items():
            print("-- %s: %s" % (e, EXPERIMENTS[e][4]))
            print("     x       Q2     g1d +- diag        " +
                  "  ".join("%9s %9s" % (n, n + "-noD") for n in cols))
            for j, p in enumerate(pts):
                print("%8.4f %7.2f  %+.4f +- %.4f %s " % (
                    p["x"], p["q2"], p["value"], p["err"],
                    " " if dis_cut(p) else "*")
                    + "  ".join("%+9.4f %+9.4f" % (
                        configs[n]["experiments"][e]["theory"][j],
                        configs[n]["experiments"][e]["theory_no_dstate"][j])
                        for n in cols))
        print("(* = outside the DIS cut; -noD = (g1p + g1n)/2 without the "
              "D-state factor)")
        for n in cols:
            c = configs[n]
            print("kernel check [%s]: max |InclusiveKernel(deuteron).tables().g1"
                  " / g1d_tree - 1| = %.3g; combined %s chi2 = %.1f/%d (p = %s, "
                  "recorded)" % (n, c["kernel_max_rel_diff"],
                                 "+".join(VERDICT_EXPERIMENTS), c["world_chi2"],
                                 c["world_ndf"], _fmt_p(c["world_p"])))
        for n in CONFIGS:
            c = configs[n]
            for e in EXPERIMENTS:
                if c["status"] == "blocked":
                    print("SUBROW | %s[%s/%s] | not evaluated (%s) | %s | %s | "
                          "BLOCKED (%s)" % (NAME, n, e, c["label"],
                                            EXPERIMENTS[e][4],
                                            report["tolerance"], c["reason"]))
                    continue
                s = c["experiments"][e]
                print("SUBROW | %s[%s/%s] | %s | %s | %s | %s" % (
                    NAME, n, e, _exp_line(c, e), EXPERIMENTS[e][4],
                    report["tolerance"],
                    s["status"].upper() if s["in_verdict"]
                    else "RECORDED (%s)" % s["status"].upper()))
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
