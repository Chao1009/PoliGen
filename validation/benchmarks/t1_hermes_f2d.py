#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T1 -- HERMES F2 of the deuteron (and of the proton) against the tree's
per-nucleon F2 (survey row F-7 of docs/benchmarking/02_data_nucleon_deuteron.md
sec. 3; sec. 10 rank 2; run 2026-10-26 phase C1).

WHAT IS COMPARED.  HERMES's F2^d PER NUCLEON and F2^p -- A. Airapetian et
al., JHEP 05 (2011) 126, arXiv:1103.5704; HEPData ins894309 Tables 1 (p) and
2 (d), CC0 -- 81 (x, Q2) bins each, 0.0069 <= <x> <= 0.6641, 0.291 <= <Q2>
<= 12.778 GeV2, vendored from the record's yaml with every row checked
against the paper's Tables 5 and 6 (data/hepdata_ins894309_hermes_f2.json,
whose _provenance carries the record, the per-table sha256 and a record
defect found while vendoring: Table 2's cross-section column repeats the
proton's; only F2 is used).  Against it, at each point's OWN (<x>, <Q2>),
the tree's per-nucleon F2 of the deuteron, (F2p + F2n)/2, and F2p, for the
three unpolarized backends (the `--unpol-sf` choices), exactly as
t1_nmc_f2d.py builds them (its structure re-used, nothing imported from it):

  toy      ToyF2() -- the SHIPPED default (what the generator's
           `default_inclusive_kernel(deuteron()).tables(x, Q2).f2` returns).
           THE REPORT ROW'S VERDICT IS THIS CONFIGURATION'S, ON THE DEUTERON.
  mstw     MstwSF() -- MSTW2008 LO from PYTHIA 8's pdfdata (`--unpol-sf
           mstw`).  Sub-row; BLOCKED (environment) only if the PYTHIA 8 tier
           or its pdfdata grid is absent.
  ct18nlo  LhapdfSF("CT18NLO") (`--unpol-sf ct18nlo`), F2 at leading order
           as the tree evaluates it.  Sub-row; BLOCKED (environment) when
           LHAPDF reports the set not installed -- never dropped, never
           faked.

Each configuration carries two targets: `deuteron` (its status is the
configuration's) and `proton` (recorded, its own pass/fail).  ZERO free
parameters, NO tuning.

WHICH POINTS.  THE DIS CUT, as t1_nmc_f2d.py's and declared before any chi2
of this harness was computed: Q2 >= 1 GeV2 and W2 = M^2 + Q2 (1 - x)/x >= 4
GeV2 (M = 0.938272 GeV), and Q2 >= the backend's own grid floor where it has
one (CT18NLO's `LhapdfSF.q2_min`, 1.677; MSTW2008's grid starts at 1).  The
chi2 is computed on the cut set; the points below Q2 = 1 (HERMES's low-Q2
reach, 0.006 < x < 0.04) are listed and not compared.

ERRORS, AND WHICH CHI2 IS THE VERDICT -- DIFFERENT FROM t1_nmc_f2d.py, AND
WHY.  The record gives, per bin and IN PERCENT of the value, stat and four
systematics (PID, MOD = model dependence outside the acceptance, MIS =
misalignment, RAD = Bethe-Heitler efficiencies), and ONE overall
normalisation uncertainty per target, 7.5 % (d) / 7.6 % (p), which the paper
states is "dominated by the uncertainty of the year-dependent luminosity
constant" (sec. 5.1): common to every bin.  t1_nmc_f2d.py's verdict is the
diagonal chi2 with every component in quadrature because NMC's file carries
no statement of correlation that NMC makes; here the reference states one.
So THE VERDICT is
    chi2 = r^T C^-1 r,  C = diag(stat^2 + PID^2 + MOD^2 + MIS^2 + RAD^2)
                            + n n^T,   n_i = 0.075 F2_i (0.076 for p),
the per-bin components uncorrelated (the paper states no correlation for
them) and the normalisation ONE fully correlated component (absolute values
from the data central values, no t0 rescaling, as t1_nmc_f2d.py).
RECORDED beside it, never in the verdict: (a) the diagonal chi2 with the
per-bin components only (no normalisation) -- the strictest reading; (b) the
diagonal chi2 with the normalisation added in quadrature to every point --
t1_nmc_f2d.py's literal "every component in quadrature", the most lenient;
(c) the normalisation pull b = the shift, in units of the 7.5 %, that the
covariance chi2 implicitly applies (chi2 = min_b [sum (r_i - b n_i)^2/D_i +
b^2], b-hat = sum(r n/D)/(1 + sum(n^2/D))).

TOLERANCE, DECLARED BEFORE ANY COMPARISON WAS COMPUTED: p(chi2, ndf) >= 0.01
on the DIS-cut set, ndf = the number of points (nothing is fitted), with the
covariance above.  The threshold is the T1 rows' and row 4's; it was chosen
for those rows after their previews (t1_nmc_f2d.py says so) and is applied
here unchanged.  WHAT WAS SEEN BEFORE THIS HEADER WAS WRITTEN: the HERMES
table (every row, while it was vendored and checked against the paper);
t1_nmc_f2d's committed verdicts on NMC's F2d (ToyF2 8307.9/155, MSTW
416.5/155, CT18NLO 3401.8/146 in phase 0) -- a strong prior that the toy
fails here too; and one bin of the b1 inversion below, estimated by hand
(bin 1: F2d_implied / F2d_HERMES about 1.05).  No chi2 of this harness had
been computed.

THE b1 INVERSION (RECORDED sub-row `b1_inversion`; 02 sec. 10 rank 2 asks for
exactly this datum).  WHAT A PER-NUCLEON F2d MEANS FOR HERMES's b1: HERMES's
b1 paper (PRL 95 (2005) 242001, Eq. (5)) extracted b1^d = -(3/2) A_zz^d F1^d
with F1^d = (1 + gamma^2) F2^d/(2x(1 + R)), F2^d = F2^p (1 + F2^n/F2^p)/2
(ALLM97 and NMC) and R = R1990 -- a per-NUCLEON F2^d, so its b1^d is per
nucleon (the arbiter of constants.hpp's B1_PER_DEUTERON_TO_PER_NUCLEON
block).  This record measures that F2^d itself, per nucleon, in the same
apparatus.  Inverting the six rows of Table II (data/hermes2005_b1d_table2.json,
sha256-checked here) gives the F1 HERMES divided by, F1_implied = -(2/3)
b1/A_zz (exact up to the printed digits), and F2d_implied = 2x(1 + R1990)
F1_implied/(1 + gamma^2), gamma^2 = 4 M^2 x^2/Q2, with R1990 the form of
Whitlow et al. as HERMES used it (data/whitlow1990_r1990.json, the three
parameters t2_hermes_b1_table2.py reads, sha256-checked here).  Each is
divided by THIS record's F2d at the nearest HERMES bin in (ln x, ln Q2).
RULE, declared before computing it: "per nucleon" is CONFIRMED if all six
ratios lie in [0.8, 1.25]; "per deuteron" if all six lie in [1.6, 2.5];
otherwise UNDETERMINED.  The window absorbs the nearest-bin offset (F2d
changes by a few % between neighbouring HERMES bins) and the 2005 F2d model
(ALLM97 x NMC) against the 2011 measurement.  What it can and cannot settle:
it decides whether HERMES's PUBLISHED b1^d is per nucleon -- the factor of 2
between HERMES and a per-deuteron b1 -- against a published HERMES table; it
says NOTHING about Miller's own normalisation (B1_MILLER_TABLE_TO_PER_NUCLEON
= 0.5 rests on reading Miller's equations, an author decision recorded in
docs/OPEN_ITEMS_SOLUTIONS.md sec. 10), and nothing about b1 itself.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations / A-scaling -- the F2 comparison: NOT APPLICABLE (A = 2
     both sides, no 0.5 or 2/6 anywhere).  The b1_inversion sub-row EXISTS
     to test rule 1's premise for HERMES: whether its b1^d is per nucleon.
  2. isoscalar denominator -- APPLIED AND ASSERTED.  The tree side of the
     deuteron is the ISOSCALAR nucleon (F2p + F2n)/2 built in ONE place
     (`isoscalar_nucleon`), and run() records the largest relative
     difference to what the generator's kernel returns for the deuteron
     (`InclusiveKernel(deuteron).tables().f2`); the pytest asserts it.  The
     proton is compared to F2p.  Reference side: the deuteron is isoscalar;
     HERMES applies no non-isoscalarity correction.
  3. per nucleus / per nucleon -- both sides PER NUCLEON: HERMES normalises to
     "the integrated luminosity L per nucleon" (its Eq. (4.3), the factor
     A/Z), and its F2d is 0.036-0.354 (about 0.33 at x ~ 0.1), one nucleon's
     size.  NOT like for like in one respect, stated rather than corrected:
     HERMES's F2d is the deuteron's (Fermi motion, binding, shadowing
     included), the tree's the free isoscalar nucleon (no deuteron structure
     correction exists for F2 in the tree).  HERMES derived F2 from its Born
     cross section with R = R1998 (Tables 5, 6 captions); the tree's F2 does
     not depend on R.

WHAT THIS DOES NOT VALIDATE.  R (t1_r_world.py); F2 below Q2 = 1 GeV2;
anything nuclear beyond A = 2; the polarized and tensor sectors; b1 (the
inversion recovers the F1 HERMES used, nothing about b1's size or shape);
HERMES's sigma_d/sigma_p (Table 3, not compared).

INTEGRITY.  The sha256 of every byte of the vendored F2 file and of its
canonical body, and of every byte of the two older files the b1 inversion
reads, are constants below, re-checked on EVERY read: an altered file raises
RuntimeError, so `__main__` exits 2.

Run:  source env.sh && python3 validation/benchmarks/t1_hermes_f2d.py [--plot]

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
DATA = os.path.join(HERE, "data", "hepdata_ins894309_hermes_f2.json")
#: sha256 of every byte of the vendored file and of its canonical body
#: (written 2026-10-26 from the HEPData yaml; the file's _provenance names
#: each download and its sha256)
FILE_SHA256 = "d32f7ea09a18ee59bf2515071ea149c29dad9a3def3490e64bf46384a7b26875"
BODY_SHA256 = "3b2c456a7628c5b0e7404bab91633d8508864f5daffb1b4187a33ecfa4039af5"
#: the two files the b1 inversion reads -- vendored 2026-09-23/26 for
#: t2_hermes_b1_table2.py, whose DATA_SHA256 / R1990_SHA256 are these values
B1_DATA = os.path.join(HERE, "data", "hermes2005_b1d_table2.json")
B1_SHA256 = "20b9104edf52414822e88305f7144566568ee2979bd531f77dd79f3a92006fc5"
R1990_DATA = os.path.join(HERE, "data", "whitlow1990_r1990.json")
R1990_SHA256 = "ce330eeef6305315a78d66422ff81fe8307bf71fdca4b8ef6e781e71698591d9"

NAME = "t1_hermes_f2d"
P_MIN = 0.01             # tolerance: p(chi2, ndf) >= 0.01, DIS cut, covariance
Q2_CUT = 1.0             # GeV^2 -- the declared DIS cut ...
W2_CUT = 4.0             # GeV^2 -- ... and its W^2 half
M_NUCLEON = 0.938272     # GeV, for W^2 and gamma^2 only
CONFIGS = ("toy", "mstw", "ct18nlo")
TARGETS = ("deuteron", "proton")
HEADLINE = "toy"         # the shipped default carries the row's verdict ...
HEADLINE_TARGET = "deuteron"   # ... on the deuteron
CT18 = "CT18NLO"
#: the b1-inversion rule (header): per nucleon if every ratio is in
#: PER_NUCLEON, per deuteron if every ratio is in PER_DEUTERON
PER_NUCLEON = (0.8, 1.25)
PER_DEUTERON = (1.6, 2.5)
REFERENCE = ("HERMES F2d (and F2p) per nucleon, JHEP 05 (2011) 126, HEPData "
             "ins894309 Tables 2 (d) and 1 (p)")


# ---------------------------------------------------------------- reference

def _check_bytes(path, sha):
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != sha:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, sha))
    return raw


def read_vendored(path=None):
    """The vendored document -- RuntimeError unless the file's bytes AND its
    canonical body hash to the recorded digests."""
    path = path or DATA
    doc = json.loads(_check_bytes(path, FILE_SHA256).decode("ascii"))
    canon = json.dumps(doc["body"], sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False)
    body = hashlib.sha256(canon.encode("ascii")).hexdigest()
    if not body == doc["_body_sha256"] == BODY_SHA256:
        raise RuntimeError("%s: body sha256 %s; file records %s, harness %s"
                           % (path, body, doc["_body_sha256"], BODY_SHA256))
    return doc


def load_table(target, path=None):
    """The points of one target as dicts: x, q2, value, the per-bin
    components (absolute), `err_pp` (their quadrature sum: the point-to-point
    error), `norm` (the overall normalisation, absolute, 100 % correlated)."""
    body = read_vendored(path)["body"]
    t = body[target]
    norm_pct = float(t["norm_pct"])
    out = []
    for i, p in enumerate(t["points"]):
        x, q2, v = float(p["x"]), float(p["q2"]), float(p["value"])
        comps = {k: float(p[k + "_pct"]) / 100.0 * v
                 for k in ("stat", "pid", "mod", "mis", "rad")}
        out.append(dict(
            i=i, x=x, q2=q2, value=v, x_bin=p["x_bin"], q2_bin=p["q2_bin"],
            w2=M_NUCLEON ** 2 + q2 * (1.0 - x) / x, comps=comps,
            err_pp=math.sqrt(sum(c * c for c in comps.values())),
            norm=norm_pct / 100.0 * v, norm_pct=norm_pct))
    return out


# --------------------------------------------------------------- statistics

def chi2_sf_fallback(chi2, ndf):
    """Q(ndf/2, chi2/2) without scipy -- t1_nmc_f2d.py's `chi2_sf_fallback`
    (Numerical Recipes' gammq), copied unchanged."""
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


def chi2_norm(points, theory):
    """(chi2, b-hat) for C = diag(err_pp^2) + n n^T (Sherman-Morrison), with
    b-hat the normalisation shift in units of the normalisation error."""
    r = [p["value"] - t for p, t in zip(points, theory)]
    d = [p["err_pp"] ** 2 for p in points]
    n = [p["norm"] for p in points]
    a = sum(ri * ri / di for ri, di in zip(r, d))
    b = sum(ri * ni / di for ri, ni, di in zip(r, n, d))
    c = sum(ni * ni / di for ni, di in zip(n, d))
    return a - b * b / (1.0 + c), b / (1.0 + c)


def stats(points, theory):
    nn = len(points)
    c, bhat = chi2_norm(points, theory)
    pulls_pp = [(t - p["value"]) / p["err_pp"] for p, t in zip(points, theory)]
    c_pp = float(sum(u * u for u in pulls_pp))
    c_quad = float(sum(((t - p["value"]) ** 2) / (p["err_pp"] ** 2 + p["norm"] ** 2)
                       for p, t in zip(points, theory)))
    worst = max(range(nn), key=lambda j: abs(pulls_pp[j]))
    return dict(
        n=nn, chi2=float(c), p=chi2_sf(c, nn), norm_pull=float(bhat),
        chi2_pp=c_pp, p_pp=chi2_sf(c_pp, nn),
        chi2_quad=c_quad, p_quad=chi2_sf(c_quad, nn),
        mean_ratio=float(sum(t / p["value"] for p, t in zip(points, theory)) / nn),
        worst_pull_pp=float(pulls_pp[worst]), worst_x=points[worst]["x"],
        worst_q2=points[worst]["q2"])


# ---------------------------------------------------------------- the tree

def isoscalar_nucleon(sf, x, q2):
    """RULE 2: the tree's per-nucleon F2 of the deuteron, (F2p + F2n)/2 --
    the ONLY place the deuteron's theory side is built."""
    return 0.5 * (sf.f2p(x, q2) + sf.f2n(x, q2))


def theory_of(sf, target, x, q2):
    return isoscalar_nucleon(sf, x, q2) if target == "deuteron" else sf.f2p(x, q2)


def pythia_version():
    """PYTHIA's version string as the xmldoc next to MstwSF's pdfdata says
    it (the project pins 8.317; record what actually ran)."""
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
    """(sf, q2_floor, label, deuteron kernel, blocked_reason) -- t1_nmc_f2d's
    logic: blocked ONLY for a named environment reason; any other error
    propagates (broken harness)."""
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


# -------------------------------------------------------- the b1 inversion

def r1990_factory(path=None):
    """R1990 (Whitlow et al. 1990, Eqs. (6)-(7) as data/whitlow1990_r1990.json
    transcribes them): the R inside HERMES's b1."""
    p = json.loads(_check_bytes(path or R1990_DATA, R1990_SHA256).decode("utf-8"))

    def r1990(x, q2):
        th = 1.0 + 12.0 * (q2 / (q2 + 1.0)) * (0.125 ** 2 / (0.125 ** 2 + x * x))
        return (p["b1"] / math.log(q2 / p["Lambda2_GeV2"]) * th + p["b2"] / q2
                + p["b3"] / (q2 * q2 + 0.09))
    return r1990


def b1_inversion(deuteron_points, b1_path=None, r1990_path=None):
    d = json.loads(_check_bytes(b1_path or B1_DATA, B1_SHA256).decode("utf-8"))
    r1990 = r1990_factory(r1990_path)
    bins = []
    for row in d["rows"]:
        x, q2, azz, _, _, b1, _, _ = row
        f1 = -(2.0 / 3.0) * b1 / azz            # the scale 1e-2 cancels
        g2 = 4.0 * M_NUCLEON ** 2 * x * x / q2
        r = r1990(x, q2)
        f2_implied = 2.0 * x * (1.0 + r) * f1 / (1.0 + g2)
        near = min(deuteron_points, key=lambda p: math.hypot(
            math.log(p["x"] / x), math.log(p["q2"] / q2)))
        bins.append(dict(x=x, q2=q2, azz=azz, b1=b1, f1_implied=f1,
                         r1990=r, gamma2=g2, f2d_implied=f2_implied,
                         hermes_x=near["x"], hermes_q2=near["q2"],
                         hermes_f2d=near["value"],
                         ratio=f2_implied / near["value"]))
    ratios = [b["ratio"] for b in bins]
    geo = math.exp(sum(math.log(v) for v in ratios) / len(ratios))
    return dict(bins=bins, reading=classify(ratios), geometric_mean=geo,
                min_ratio=min(ratios), max_ratio=max(ratios))


def classify(ratios):
    """THE declared b1-inversion rule (module docstring)."""
    if all(PER_NUCLEON[0] <= v <= PER_NUCLEON[1] for v in ratios):
        return "per nucleon"
    if all(PER_DEUTERON[0] <= v <= PER_DEUTERON[1] for v in ratios):
        return "per deuteron"
    return "undetermined"


# --------------------------------------------------------------- the row

def measure(path=None, b1_path=None, r1990_path=None):
    pts = {t: load_table(t, path) for t in TARGETS}
    configs = {}
    for name in CONFIGS:
        sf, floor, label, kern, why = build_backend(name)
        if why is not None:
            configs[name] = dict(name=name, status="blocked", reason=why,
                                 label=label)
            continue
        exps = {}
        for t in TARGETS:
            theory = [theory_of(sf, t, p["x"], p["q2"]) for p in pts[t]]
            sel = [j for j, p in enumerate(pts[t]) if dis_cut(p, floor)]
            s = stats([pts[t][j] for j in sel], [theory[j] for j in sel])
            exps[t] = dict(s, theory=theory, selected=sel,
                           status="pass" if s["p"] >= P_MIN else "fail")
        kdev = max(abs(kern.tables(p["x"], p["q2"]).f2 / th - 1.0)
                   for p, th in zip(pts["deuteron"], exps["deuteron"]["theory"]))
        configs[name] = dict(
            name=name, label=label, q2_floor=max(Q2_CUT, floor or 0.0),
            kernel_max_rel_diff=float(kdev), experiments=exps,
            status=exps[HEADLINE_TARGET]["status"])
    inv = b1_inversion(pts["deuteron"], b1_path, r1990_path)
    return pts, configs, inv


def _fmt_p(p):
    return "%.3g" % p


def _line(c, t):
    e = c["experiments"][t]
    return ("%s %s: chi2/ndf = %.1f/%d (p = %s) [per-bin diagonal + 1 correlated "
            "normalisation, the verdict]; normalisation pull %+.2f sigma_norm; "
            "recorded: per-bin only %.1f/%d (p = %s), normalisation in quadrature "
            "%.1f/%d (p = %s); <th/data> = %.3f; worst per-bin pull %+.1f at x = %g, "
            "Q2 = %g" % (c["label"], t, e["chi2"], e["n"], _fmt_p(e["p"]),
                         e["norm_pull"], e["chi2_pp"], e["n"], _fmt_p(e["p_pp"]),
                         e["chi2_quad"], e["n"], _fmt_p(e["p_quad"]), e["mean_ratio"],
                         e["worst_pull_pp"], e["worst_x"], e["worst_q2"]))


def run(path=None, verbose=True, b1_path=None, r1990_path=None):
    t0 = time.perf_counter()
    pts, configs, inv = measure(path, b1_path, r1990_path)
    h = configs[HEADLINE]
    hd = h["experiments"][HEADLINE_TARGET]
    others = "; ".join(
        "%s %s" % (n, "BLOCKED (%s)" % configs[n]["reason"]
                   if configs[n]["status"] == "blocked" else
                   "%.1f/%d (p = %s)" % (configs[n]["experiments"]["deuteron"]["chi2"],
                                         configs[n]["experiments"]["deuteron"]["n"],
                                         _fmt_p(configs[n]["experiments"]["deuteron"]["p"])))
        for n in CONFIGS if n != HEADLINE)
    subrows = []
    for n in CONFIGS:
        c = configs[n]
        if c["status"] == "blocked":
            subrows.append(dict(name=n, status="blocked", reason=c["reason"],
                                label=c["label"]))
            continue
        subrows.append(dict(
            name=n, status=c["status"], label=c["label"], q2_floor=c["q2_floor"],
            kernel_max_rel_diff=c["kernel_max_rel_diff"],
            # `experiments` is the nested container _report.py reads: here the
            # two TARGETS, each with its own status
            experiments={t: dict((k, v) for k, v in c["experiments"][t].items()
                                 if k not in ("theory", "selected"))
                         for t in TARGETS}))
    subrows.append(dict(
        name="b1_inversion", status="recorded", reading=inv["reading"],
        geometric_mean=inv["geometric_mean"], min_ratio=inv["min_ratio"],
        max_ratio=inv["max_ratio"],
        ratios=[b["ratio"] for b in inv["bins"]],
        rule="per nucleon iff every ratio in [%g, %g]; per deuteron iff in [%g, %g]"
             % (PER_NUCLEON + PER_DEUTERON)))
    report = dict(
        name=NAME,
        generator=("toy (shipped ToyF2) on F2d: chi2/ndf = %.1f/%d, p = %s on Q2 >= 1, "
                   "W2 >= 4 [per-bin diagonal + correlated 7.5 %% normalisation]; "
                   "sub-rows: %s; HERMES b1 inversion: F2d_implied/F2d_HERMES = %s "
                   "-> %s" % (hd["chi2"], hd["n"], _fmt_p(hd["p"]), others,
                              " / ".join("%.3f" % b["ratio"] for b in inv["bins"]),
                              inv["reading"].upper())),
        reference="%s, %d of %d deuteron points in the DIS cut"
                  % (REFERENCE, hd["n"], len(pts["deuteron"])),
        tolerance=("p(chi2, ndf) >= %.2f on the DIS cut, per-bin errors diagonal + the "
                   "overall normalisation fully correlated" % P_MIN),
        status=h["status"], headline="%s/%s" % (HEADLINE, HEADLINE_TARGET),
        subrows=subrows, configs_detail=configs, b1_inversion=inv,
        npoints={t: len(pts[t]) for t in TARGETS}, points=pts,
        dis_cut="Q2 >= %g GeV2 and W2 >= %g GeV2 (and >= the backend's grid floor)"
                % (Q2_CUT, W2_CUT),
        config=dict(data_file_sha256=FILE_SHA256, b1_file_sha256=B1_SHA256,
                    r1990_file_sha256=R1990_SHA256),
    )
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        _print(report)
    return report


def _print(rep):
    configs, pts, inv = rep["configs_detail"], rep["points"], rep["b1_inversion"]
    cols = [n for n in CONFIGS if configs[n]["status"] != "blocked"]
    for t in TARGETS:
        print("\nHERMES F2 %s per nucleon vs the tree -- %d points, DIS cut %s"
              % (t, len(pts[t]), rep["dis_cut"]))
        print("    <x>   <Q2>     W2   F2 +- per-bin (norm %.1f %%)   %s"
              % (pts[t][0]["norm_pct"], "  ".join("%10s" % n for n in cols)))
        for j, p in enumerate(pts[t]):
            print("%7.4f %6.3f %6.2f   %.4f +- %.4f %s  %s" % (
                p["x"], p["q2"], p["w2"], p["value"], p["err_pp"],
                " " if dis_cut(p) else "*",
                "  ".join("%10.4f" % configs[n]["experiments"][t]["theory"][j]
                          for n in cols)))
        print("(* = outside Q2 >= 1, W2 >= 4; CT18NLO also needs Q2 >= its grid floor)")
    for n in cols:
        print("kernel check [%s]: max |InclusiveKernel(deuteron).tables().f2 / "
              "isoscalar - 1| = %.3g" % (n, configs[n]["kernel_max_rel_diff"]))
    print("\nHERMES b1 inversion (PRL 95 (2005) 242001 Table II; R1990 as HERMES used):")
    for b in inv["bins"]:
        print("  x = %.3f Q2 = %.2f: A_zz %+.2f, b1 %+.2f (1e-2) -> F1_implied %.4f, "
              "R1990 %.3f, F2d_implied %.4f; HERMES 2011 F2d %.4f at (%.4f, %.3f): "
              "ratio %.3f" % (b["x"], b["q2"], b["azz"], b["b1"], b["f1_implied"],
                              b["r1990"], b["f2d_implied"], b["hermes_f2d"],
                              b["hermes_x"], b["hermes_q2"], b["ratio"]))
    print("  -> %s (geometric mean %.3f; range %.3f ... %.3f; rule: per nucleon iff "
          "every ratio in [%g, %g], per deuteron iff in [%g, %g])"
          % (inv["reading"].upper(), inv["geometric_mean"], inv["min_ratio"],
             inv["max_ratio"], PER_NUCLEON[0], PER_NUCLEON[1], PER_DEUTERON[0],
             PER_DEUTERON[1]))
    print()
    for n in CONFIGS:
        c = configs[n]
        if c["status"] == "blocked":
            print("SUBROW | %s[%s] | not evaluated (%s) | %s | %s | BLOCKED (%s)"
                  % (NAME, n, c["label"], REFERENCE, rep["tolerance"], c["reason"]))
            continue
        for t in TARGETS:
            e = c["experiments"][t]
            print("SUBROW | %s[%s/%s] | %s | %s | %s | %s%s" % (
                NAME, n, t, _line(c, t), REFERENCE, rep["tolerance"],
                e["status"].upper(), "" if t == HEADLINE_TARGET else " (recorded)"))
    print("SUBROW | %s[b1_inversion] | F2d_implied/F2d_HERMES = %s | HERMES PRL 95 "
          "(2005) 242001 Table II and this record's F2d | per nucleon iff all in "
          "[%g, %g] | RECORDED: %s" % (
              NAME, " / ".join("%.3f" % b["ratio"] for b in inv["bins"]),
              PER_NUCLEON[0], PER_NUCLEON[1], inv["reading"].upper()))
    print("REPORT | %s | %s | %s | %s | %s" % (
        NAME, rep["generator"], rep["reference"], rep["tolerance"],
        rep["status"].upper()))


# ------------------------------------------------------------------ figure

def plot(report, out_dir=None):
    """F2d per nucleon vs x (the DIS-cut points, ink; each backend's value at
    the point's own (x, Q2) as a tick in its slot colour), over a theory/data
    ratio panel with the 7.5 % overall normalisation as a grey band."""
    sys.path.insert(0, HERE)
    from _plots import band, new_figure, save, status_title, style_series
    fig, (ax, axr) = new_figure(two_panel=True, height=4.0)
    configs, pts = report["configs_detail"], report["points"]["deuteron"]
    floor_free = [j for j, p in enumerate(pts) if dis_cut(p)]
    xs = [pts[j]["x"] for j in floor_free]
    ax.errorbar(xs, [pts[j]["value"] for j in floor_free],
                yerr=[pts[j]["err_pp"] for j in floor_free],
                **style_series("reference", label="HERMES F2d per nucleon (per-bin errors)"))
    slots = {"toy": ("default", "ToyF2 (shipped)"), "mstw": ("mstw", "MSTW2008 LO"),
             "ct18nlo": ("ct18nlo", "CT18NLO (at LO)")}
    for n in CONFIGS:
        c = configs[n]
        if c["status"] == "blocked":
            continue
        e = c["experiments"]["deuteron"]
        sel = e["selected"]
        slot, lab = slots[n]
        kw = style_series(slot, label="%s: %.0f/%d" % (lab, e["chi2"], e["n"]))
        kw.update(marker="_", linestyle="none", markersize=9, markeredgewidth=2)
        ax.plot([pts[j]["x"] for j in sel], [e["theory"][j] for j in sel], **kw)
        rk = style_series(slot, marker="o")
        rk.update(markersize=3, label="_nolegend_")
        axr.plot([pts[j]["x"] for j in sel],
                 [e["theory"][j] / pts[j]["value"] for j in sel], **rk)
    xl = sorted(xs)
    band(axr, [xl[0] * 0.9, xl[-1] * 1.1], [1 - pts[0]["norm_pct"] / 100.0] * 2,
         [1 + pts[0]["norm_pct"] / 100.0] * 2, label="_nolegend_")
    axr.axhline(1.0, color="#c3c2b7", linewidth=1.0, zorder=0)
    axr.set_ylim(0.5, 1.5)
    ax.set_xscale("log")
    ax.set_ylabel("F2 per nucleon (deuteron)")
    ax.legend(loc="lower center", ncol=2, fontsize=7)
    axr.set_xlabel("<x> (each point at its own <Q2>, 1 <= Q2 <= 12.8 GeV2)")
    axr.set_ylabel("theory / data")
    h = configs[HEADLINE]["experiments"]["deuteron"]
    status_title(fig, NAME, report["status"],
                 "ToyF2 on F2d: %.1f/%d (p = %.2g); tolerance p >= 0.01, Q2 >= 1, W2 >= 4, "
                 "7.5 %% norm. correlated (band)" % (h["chi2"], h["n"], h["p"]))
    return save(fig, NAME, out_dir)


def main(argv=None):
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    argv = sys.argv[1:] if argv is None else argv
    try:
        rep = run()
        if "--plot" in argv:
            print("figure: %s" % plot(rep))
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
