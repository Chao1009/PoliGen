#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T1 -- R = sigma_L/sigma_T against the measured R of the world (survey rows
R-1, R-2, R-3, R-5 of docs/benchmarking/02_data_nucleon_deuteron.md sec. 4),
with NMC's R at the letter's sweet spots as the PRIMARY sub-row (run
2026-10-26 phase C1; docs/benchmarking/08_channel_priority.md ingredient A).

WHAT IS COMPARED.  The tree's two R backends, evaluated at each data point's
OWN (x, Q2):

  r_sigma_lt  `lg.r_sigma_lt(x, q2)` = 0.18/(1 + Q2/50) -- the SHIPPED
              default R of the whole program (sf.hpp: "a placeholder of the
              right magnitude ... not a fit"; a null `RFunc` hook resolves to
              it).  THE ROW'S VERDICT IS THIS BACKEND'S, ON THE PRIMARY
              (nmc) SUB-ROW.
  r1998       `lg.r1998(x, q2)` -- the E143 world fit R1998 (Abe et al., PLB
              452 (1999) 194, Table II; the average of its three forms), with
              `clip = True` AS SHIPPED (every generation path calls it so):
              outside its support 0.005 <= x <= 0.86, 0.5 <= Q2 <= 130 GeV2
              (sf.hpp R1998_X_MIN ... R1998_Q2_MAX; this harness copies the
              four numbers and checks them against the clip at run time) it
              is FROZEN at the boundary.  Recorded beside the verdict, with
              its own pass/fail.
  r1998_unclipped  `lg.r1998(x, q2, clip=False)` -- the bare formula
              continued outside the support.  RECORDED only (no generation
              path calls it).

against four data sets, vendored under data/ from HEPData (CC0) with their
provenance, each re-checked by sha256 on every read.  The SUB-ROWS, in this
order and with these LABELS:

  nmc    PRIMARY -- NMC, Nucl. Phys. B 483 (1997) 3, HEPData ins424154 Table
         33 (= the paper's Table 7): ten x bins 0.0045 <= x <= 0.11 at their
         mean Q2 1.31-20.6 GeV2.  The only measured R that covers the
         letter's sweet-spot x values (0.011 / 0.028 / 0.071 / 0.141) -- in
         x, NOT in Q2: at x = 0.0125 NMC's R sits at <Q2> = 2.2 GeV2, the
         cells at 1.14-14.3.  NMC's R is "the average R for the proton and
         the deuteron", one Q2-averaged value per x bin (paper sec. 4).
         x = 0.0045 is below R1998's support (marked extrapolated).
  h1     EXTRAPOLATION -- H1, Phys. Lett. B 665 (2008) 139, HEPData
         ins786161: F_L of the proton at 0.00024 <= x <= 0.0036, 12 <= Q2 <=
         90 GeV2.  EVERY H1 point is below R1998's support (x < 0.005): r1998
         there is its value frozen at x = 0.005 -- this sub-row measures what
         that freeze costs.  F_L IS COMPARED DIRECTLY, not converted to R: the
         record carries F_L only -- no F2 and no covariance -- so R =
         F_L/(F2 - F_L) cannot be formed from it.  The tree's F_L is its own
         `UnpolSF.flp` = F2p R/(1 + R) (sf.hpp `fl_from_f2`) with R from the
         backend under test through the `set_r_func` hook (a null hook =
         r_sigma_lt) and F2p from MstwSF (MSTW2008 LO, a fit that includes
         HERA's F2 data; the opt-in `--unpol-sf mstw`) -- F2 is the nuisance
         here and MSTW is the tree's backend closest to the proton F2 at these
         (x, Q2).  RECORDED beside it: the same with the SHIPPED ToyF2's F2p
         (the generator's own proton F_L) and with CT18NLO's (evaluated at
         leading order).  The verdict set is the record's Table 9 (= the
         paper's Table 2: eight Q2-averaged points, stat / uncorrelated /
         correlated); the 25 per-bin points of Tables 1-8 (total and stat
         only) are recorded.  For the FIGURE only, H1 is shown in R as
         F_L/(F2p_MSTW - F_L) (the verdict's own F2, so the figure and the
         chi2 agree).  If MstwSF is unavailable (no PYTHIA 8 tier) the h1
         verdict cells are BLOCKED (environment) and the ToyF2 variant is
         still recorded.
  r1990  CLOSURE -- Whitlow et al., Phys. Lett. B 250 (1990) 193, HEPData
         ins296980 Tables 3-13: the global SLAC extraction (89 bins, 0.1 <=
         x <= 0.86; R_p and R_d averaged), independent of E140's deuterium
         cross sections (paper p. 6).  These data are R1998's OWN fit input
         (Abe et al. 1999: "a new fit using the present data [3,11-14]", [3]
         = Whitlow et al.): for r1998 a pass here tests the TRANSCRIPTION of
         the fit, not the physics.
  e140   CLOSURE -- Dasu et al. (SLAC E140), Phys. Rev. D 49 (1994) 5641,
         HEPData ins360765 Tables 23-25: R on deuterium at x = 0.2, 0.35, 0.5
         (11 points).  Also in R1998's input (through [3], whose record
         carries E140's R): a CLOSURE for r1998 as well.

AND A FACT FOUND WHILE READING THE PAPERS, STATED BEFORE ANY NUMBER: NMC's R
is ALSO an R1998 fit input ([12] of Abe et al. = "NMC, M. Arneodo et al.,
Nucl. Phys. B483, 3 (1997)", the very table of the nmc sub-row).  So for
r1998 every sub-row but h1 is a closure, and the only R1998-independent data
here (h1) lie outside its support.  For r_sigma_lt -- which is not a fit at
all -- every sub-row is a genuine test.  The verdict is r_sigma_lt's on nmc.

TOLERANCE, DECLARED BEFORE ANY COMPARISON WAS COMPUTED.  For every (sub-row,
backend) cell:  p(chi2, ndf) >= 0.01, ndf = the number of points (nothing is
fitted), chi2 DIAGONAL with every published uncertainty component of a point
in quadrature (nmc: stat (+) sys; h1: stat (+) uncorrelated sys (+)
correlated sys; r1990, e140: stat (+) sys), on THE DIS CUT Q2 >= 1 GeV2 and
W2 = M^2 + Q2 (1 - x)/x >= 4 GeV2 (M = 0.938272 GeV) -- the T1 rows' rule
(t1_nmc_f2d.py) and row 4's threshold, applied unchanged.  Points outside
the cut are listed and not compared.  THE ROW'S STATUS = the cell
(nmc, r_sigma_lt).  RECORDED beside the diagonal chi2, never in a verdict:
  * nmc: the chi2 with the systematic 100 % correlated over the ten bins
    (C = diag(stat^2) + s s^T) -- the extreme the paper's "largely
    correlated" (secs. 4.2, 6) allows; NOT a covariance the record gives;
  * h1: the chi2 with the correlated component 100 % correlated over the
    eight Q2 points, and the diagonal chi2 of the 25 per-bin points;
  * e140: the chi2 with the paper's "additional normalization error of 0.03"
    on R (its Table VII caption; NOT in the record's numbers) added as one
    fully correlated component;
  * r1990: none (the record gives no correlation; the paper's correlated
    radiative +-0.025 and normalisations are not separable from 'sys').
WHAT WAS SEEN BEFORE THIS HEADER WAS WRITTEN, said so that the tolerance can
be argued with and not suspected: every data table was read while it was
vendored (so a reader of the tables, as the author of this header, could see
that NMC's R = 0.537 at x = 0.0045 sits far above r_sigma_lt's 0.18/(1 +
1.38/50) = 0.175); r_sigma_lt and r1998 were each evaluated at ONE point while
probing the API (r1998(0.01, 1.14) = 0.394, r_sigma_lt = 0.176), and the
kernel's cos 2phi amplitude ratio at one sweet-spot cell (0.849 at x =
0.0282, Q2 = 1.14) while checking the bias-table code path.  NO chi2, pull or
p-value against any data set had been computed.  docs/USAGE.md's r1998 vs
r_sigma_lt prices at x = 0.05 / 0.10 / 0.30 (Q2 = 2.5) were known.

THE BIAS TABLE (08_channel_priority.md ingredient A asks for it; RECORDED,
no tolerance): (1 + r1998)/(1 + r_sigma_lt) at the four sweet-spot cells of
the letter -- PolarizedLithiumSim/plans/07 sec. 7.4 WP1 addendum (the toy ->
grid drift table, 10 x 6Li 99.5 GeV/u): spot 2 (x 0.0112, Q2 1.14), spot 1
(x 0.0282, Q2 1.14), spot 3 (x 0.0708, Q2 3.13), spot 4 (x 0.141, Q2 14.3);
sec. 7.1's YR addendum quotes the same cells rounded, (0.011, 1.14),
(0.028, 1.14), (0.071, 3.13), (0.141, 14.3) -- printed as their own SUBROW
lines; and over a grid spanning the generator window (x = 0.001 ... 0.7,
Q2 = 1 ... 100 GeV2), the cells below R1998's x support marked
'extrapolated' and those with W2 < 4 GeV2 marked.  At the four cells the
GENERATOR'S OWN FACTOR is read from the kernel as well: the cos 2phi
amplitude a_2 of `default_inclusive_kernel(li6())` (the shipped R) against
`default_inclusive_kernel(li6(), Li6Convolution, r_source=R1998)` (the only
path that accepts an R source; the Delta slot is the same for every b1
model, which run() asserts bit for bit), at y = Q2/(s x) of the letter's mid
configuration (`default_configs('6Li')[1]`, e 10 x 6Li 99.5 GeV/u).  In the
kernel a_2 = -(1/2)(1-y)/y^2 c_eff sin^2(theta_S) Delta/D_phi (Jaffe-Manohar since D4, 2026-10-27; no 1/2 before: a common factor, the ratio is bit for bit the same) with D_phi = F1 +
(1-y)/(x y^2) F2 and the shipped Delta = 1e-2 F1 x^0.3 (1-x)^4
(`toy_delta_gluon` in pipeline.cpp), so a_2 is proportional to F1/D_phi =
1/[1 + 2 (1-y)(1+R)/y^2]: F2 cancels, and the ratio is
[1 + 2(1-y)(1+R_lt)/y^2]/[1 + 2(1-y)(1+R_98)/y^2] = (1+R_lt)/(1+R_98) up to
O(y^2) -- run() prints both and their difference.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations / A-scaling -- NOT APPLICABLE: no b1, no A-scaling
     assumption anywhere (the bias table's kernel ratio is a ratio of two
     6Li kernels that differ in R only; its normalisation cancels).
  2. isoscalar denominator -- NOT APPLICABLE to R itself: R is a ratio of
     two cross sections of the same target, and the tree's R takes (x, Q2)
     only (no target, no isospin).  The data: NMC's R is p and d combined
     (R_d - R_p = 0 within 0.02, the paper's sec. 4), Whitlow's p and d
     averaged (R_d - R_p = -0.001 +- 0.009 +- 0.009, its Eq. (5)), E140's on
     deuterium, H1's on the proton -- all read as the nucleon's R.  For h1
     the tree side is the PROTON's F_L (F2p), like for like.
  3. per nucleus / per nucleon -- R is convention-free (a ratio); H1's F_L
     is per proton and the tree's flp is the proton's.  The kernel ratio of
     the bias table is per nucleon on both sides (InclusiveKernel is per
     nucleon).

WHAT THIS DOES NOT VALIDATE.  R on a nucleus (E140's R_A - R_D average,
0.001 +- 0.018 +- 0.016, is quoted from the paper; the tree's R has no A
dependence, 6Li included, and no Li data exist); the TENSOR L/T separation
(Cosyn Eqs. (16)-(17) -- never measured); R's Q2 dependence within an NMC x
bin (NMC averages over Q2); F2 (rows F-1, F-2; for h1 the F2 is a nuisance
taken from MSTW, so a pass or fail there carries MSTW2008 LO's F2 error);
the sweet-spot cells' Q2 (NMC covers their x, not their Q2); the Delta
model itself (the bias table prices R at fixed Delta/F1).

INTEGRITY.  Each vendored file's sha256 (every byte) and its canonical-body
sha256 (also recorded inside the file) are constants below, re-checked on
EVERY read: an altered file raises RuntimeError, so `__main__` exits 2.

Run:  source env.sh && python3 validation/benchmarks/t1_r_world.py [--plot]

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

NAME = "t1_r_world"
P_MIN = 0.01             # tolerance: p(chi2, ndf) >= 0.01, DIS cut, diagonal
Q2_CUT = 1.0             # GeV^2 -- the declared DIS cut ...
W2_CUT = 4.0             # GeV^2 -- ... and its W^2 half
M_NUCLEON = 0.938272     # GeV, for W^2 only (t1_nmc_f2d.py's value)

#: sf.hpp's R1998 support, copied; `check_r1998_support()` verifies each
#: edge against the shipped clip at run time (RuntimeError if they differ)
R1998_X_MIN, R1998_X_MAX = 0.005, 0.86
R1998_Q2_MIN, R1998_Q2_MAX = 0.5, 130.0

BACKENDS = ("r_sigma_lt", "r1998", "r1998_unclipped")
SUBROWS = ("nmc", "h1", "r1990", "e140")
LABELS = {"nmc": "PRIMARY", "h1": "EXTRAPOLATION", "r1990": "CLOSURE",
          "e140": "CLOSURE"}
HEADLINE = ("nmc", "r_sigma_lt")
#: the F2p used to build the tree's F_L for the h1 sub-row: the verdict's
#: first, then the recorded variants
H1_F2 = ("mstw", "toy", "ct18nlo")
CT18 = "CT18NLO"

#: (file name, sha256 of every byte, sha256 of the canonical body) --
#: written 2026-10-26 from the HEPData yaml downloads (each file's
#: _provenance names them and their sha256)
FILES = {
    "nmc": ("hepdata_ins424154_nmc_r.json",
            "9a6120c15e817127ef9c9458e26609bd076793f96c8bbc163c20c88f940be952",
            "4fac8c13fd6d65336b28fd5968714f56df126809809fe4287a4b914b2f6e32a0"),
    "h1": ("hepdata_ins786161_h1_fl.json",
           "3ea1b6f997ff4f23193b8618520474462a9d7cbe694cddf0d2e701ec9ea76ca6",
           "2ff5cc7ebe4bc5a2e7d0bdf432cfe0889aa509d914ee6a22af74854c6e425c13"),
    "r1990": ("hepdata_ins296980_whitlow_r.json",
              "e24ed0929d20da172e5a0e14a283b7d57767a458f08613c9165c25c8aa1435a7",
              "86c0cffc28d017a2a0845378ab4ca59a4b33a0266116de768eed04fe801cdd29"),
    "e140": ("hepdata_ins360765_e140_r.json",
             "e5326261b6d604a5ac167293b5577141d2afda6bf7c7b66068e357285c7410f9",
             "96498c131899688acbcca9218669d4849952522c48c5256a6e4191200ddc8b58"),
}

REFERENCES = {
    "nmc": "NMC R, NPB 483 (1997) 3, HEPData ins424154 Table 33",
    "h1": "H1 F_L, PLB 665 (2008) 139, HEPData ins786161 Table 9 (+ Tables 1-8 recorded)",
    "r1990": "Whitlow et al. R (global SLAC), PLB 250 (1990) 193, HEPData ins296980 Tables 3-13",
    "e140": "SLAC E140 R (deuterium), PRD 49 (1994) 5641, HEPData ins360765 Tables 23-25",
}

#: the letter's four sweet-spot cells, ordered by x, with plans/07's labels
#: (sec. 7.4 WP1 addendum, the toy -> grid drift table)
SWEET_SPOTS = (
    dict(label="spot 2", x=0.0112, q2=1.14),
    dict(label="spot 1", x=0.0282, q2=1.14),
    dict(label="spot 3", x=0.0708, q2=3.13),
    dict(label="spot 4", x=0.141, q2=14.3),
)
MID_CONFIG = 1           # lg.default_configs('6Li')[1]: e(10) x 6Li(99.5/u)
BIAS_GRID_X = (0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02, 0.03, 0.05,
               0.07, 0.1, 0.2, 0.3, 0.5, 0.7)
BIAS_GRID_Q2 = (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0)


# ---------------------------------------------------------------- reference

def read_vendored(key, path=None):
    """The vendored document of sub-row `key` -- RuntimeError unless the
    file's bytes AND its canonical body hash to the recorded digests."""
    fname, file_sha, body_sha = FILES[key]
    path = path or os.path.join(DATA_DIR, fname)
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


def w2_of(x, q2):
    return M_NUCLEON ** 2 + q2 * (1.0 - x) / x


def _point(key, i, x, q2, value, comps, **extra):
    """One data point: the diagonal error is every component in quadrature."""
    x, q2 = float(x), float(q2)
    p = dict(set=key, i=i, x=x, q2=q2, value=float(value),
             comps={k: float(v) for k, v in comps.items()},
             err=math.sqrt(sum(float(v) ** 2 for v in comps.values())),
             w2=w2_of(x, q2))
    p["dis"] = p["q2"] >= Q2_CUT and p["w2"] >= W2_CUT
    p["in_r1998_support"] = (R1998_X_MIN <= x <= R1998_X_MAX
                             and R1998_Q2_MIN <= q2 <= R1998_Q2_MAX)
    p.update(extra)
    return p


def load_points(key, path=None):
    """The points of sub-row `key` (h1: the verdict set, Table 9)."""
    body = read_vendored(key, path)["body"]
    if key == "nmc":
        return [_point(key, i, p["x"], p["q2"], p["value"],
                       dict(stat=p["stat"], sys=p["sys"]),
                       y_range=(p["y_low"], p["y_high"]))
                for i, p in enumerate(body["points"])]
    if key == "h1":
        return [_point(key, i, p["x"], p["q2"], p["value"],
                       dict(stat=p["stat"], sys_uncorr=p["sys_uncorr"],
                            sys_corr=p["sys_corr"]))
                for i, p in enumerate(body["averaged"])]
    if key in ("r1990", "e140"):
        return [_point(key, i, p["x"], p["q2"], p["value"],
                       dict(stat=p["stat"], sys=p["sys"]), table=p["table"])
                for i, p in enumerate(body["points"])]
    raise ValueError("unknown sub-row %r" % key)


def load_h1_per_bin(path=None):
    body = read_vendored("h1", path)["body"]
    return [_point("h1_per_bin", i, p["x"], p["q2"], p["value"],
                   dict(total=p["total"]), stat=p["stat"], table=p["table"])
            for i, p in enumerate(body["per_bin"])]


# --------------------------------------------------------------- statistics

def chi2_sf_fallback(chi2, ndf):
    """Q(ndf/2, chi2/2) without scipy -- t1_nmc_f2d.py's `chi2_sf_fallback`
    (Numerical Recipes' gammq: the series for P when chi2/2 < ndf/2 + 1,
    else the continued fraction for Q itself), copied unchanged."""
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


def chi2_rank1(resid, diag_var, vec):
    """r^T C^-1 r for C = diag(diag_var) + v v^T (Sherman-Morrison; one fully
    correlated component on top of a diagonal)."""
    a = sum(r * r / d for r, d in zip(resid, diag_var))
    b = sum(r * v / d for r, v, d in zip(resid, vec, diag_var))
    c = sum(v * v / d for v, d in zip(vec, diag_var))
    return a - b * b / (1.0 + c)


def stats(points, theory):
    """chi2 (diagonal), p, pulls (theory - data)/err, the worst pull."""
    n = len(points)
    pulls = [(t - p["value"]) / p["err"] for p, t in zip(points, theory)]
    chi2 = float(sum(u * u for u in pulls))
    worst = max(range(n), key=lambda j: abs(pulls[j]))
    return dict(n=n, chi2=chi2, p=chi2_sf(chi2, n), pulls=pulls,
                mean_diff=float(sum(t - p["value"] for p, t in zip(points, theory)) / n),
                worst_pull=float(pulls[worst]), worst_x=points[worst]["x"],
                worst_q2=points[worst]["q2"])


# ---------------------------------------------------------------- the tree

def r_backend(name):
    """The R function of backend `name` (x, q2) -> R."""
    import lipolgen as lg
    if name == "r_sigma_lt":
        return lambda x, q2: lg.r_sigma_lt(x, q2)
    if name == "r1998":
        return lambda x, q2: lg.r1998(x, q2)
    if name == "r1998_unclipped":
        return lambda x, q2: lg.r1998(x, q2, clip=False)
    raise ValueError("unknown R backend %r" % name)


def check_r1998_support():
    """The four support edges copied above are the ones the shipped clip
    uses: inside an edge the clipped value moves, beyond it it is frozen."""
    import lipolgen as lg
    eps = 1e-9
    probes = [  # (frozen beyond the edge, at the edge)
        (lambda: lg.r1998(R1998_X_MIN * 0.5, 5.0), lambda: lg.r1998(R1998_X_MIN, 5.0),
         lambda: lg.r1998(R1998_X_MIN * (1 + 1e-3), 5.0)),
        (lambda: lg.r1998(0.95, 5.0), lambda: lg.r1998(R1998_X_MAX, 5.0),
         lambda: lg.r1998(R1998_X_MAX * (1 - 1e-3), 5.0)),
        (lambda: lg.r1998(0.1, 0.3), lambda: lg.r1998(0.1, R1998_Q2_MIN),
         lambda: lg.r1998(0.1, R1998_Q2_MIN * (1 + 1e-3))),
        (lambda: lg.r1998(0.1, 500.0), lambda: lg.r1998(0.1, R1998_Q2_MAX),
         lambda: lg.r1998(0.1, R1998_Q2_MAX * (1 - 1e-3))),
    ]
    for beyond, edge, inside in probes:
        if beyond() != edge() or abs(inside() - edge()) <= eps * abs(edge()):
            raise RuntimeError("R1998 support constants of %s no longer match "
                               "the shipped clip (sf.hpp)" % NAME)
    return dict(x=(R1998_X_MIN, R1998_X_MAX), q2=(R1998_Q2_MIN, R1998_Q2_MAX))


def pythia_version():
    """PYTHIA's version as the xmldoc next to MstwSF's pdfdata says it."""
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


def f2_backend(name):
    """(factory, label, blocked_reason) for the F2p of the h1 sub-row; the
    factory returns a FRESH backend object (its R hook is set per R backend).
    Blocked ONLY for a named environment reason."""
    import lipolgen as lg
    import lipolgen._lipolgen as _l
    if name == "toy":
        return (lambda: lg.ToyF2()), "ToyF2 (shipped default)", None
    if name == "mstw":
        if not getattr(_l, "HAVE_PYTHIA8", False):
            return None, "MstwSF", ("environment: lipolgen built without the "
                                    "PYTHIA 8 tier (MstwSF)")
        d = _l.pythia8_pdfdata_dir()
        if not (d and os.path.isfile(os.path.join(d, "mstw2008lo.00.dat"))):
            return None, "MstwSF", ("environment: PYTHIA 8 pdfdata/"
                                    "mstw2008lo.00.dat not found (MstwSF)")
        return ((lambda: lg.MstwSF()),
                "MstwSF (MSTW2008 LO, PYTHIA %s pdfdata)" % pythia_version(), None)
    if name == "ct18nlo":
        if not getattr(_l, "HAVE_LHAPDF", False):
            return None, "LhapdfSF(%s)" % CT18, ("environment: lipolgen built "
                                                  "without the LHAPDF tier")
        try:
            lg.LhapdfSF(CT18)
        except RuntimeError as exc:
            if "Info file not found for PDF set" not in str(exc):
                raise
            return None, "LhapdfSF(%s)" % CT18, ("environment: LHAPDF set %s "
                                                  "not installed" % CT18)
        return (lambda: lg.LhapdfSF(CT18)), "LhapdfSF(%s) at LO" % CT18, None
    raise ValueError("unknown F2 backend %r" % name)


def tree_fl(factory, rname, points):
    """The tree's proton F_L at each point: a fresh F2 backend whose R hook
    is the backend under test (none for r_sigma_lt: the null hook IS
    r_sigma_lt, as shipped), read through `flp`; the F2p used is returned
    too, and the identity F_L = F2p R/(1+R) is checked."""
    sf = factory()
    if rname != "r_sigma_lt":
        sf.set_r_func(r_backend(rname))
    rf = r_backend(rname)
    fl, f2 = [], []
    for p in points:
        v = sf.flp(p["x"], p["q2"])
        f = sf.f2p(p["x"], p["q2"])
        r = rf(p["x"], p["q2"])
        if abs(v - f * r / (1.0 + r)) > 1e-12 * max(1.0, abs(v)):
            raise RuntimeError("flp != F2p R/(1+R) for %s at x = %g" % (rname, p["x"]))
        fl.append(v)
        f2.append(f)
    return fl, f2


# ------------------------------------------------------------ the sub-rows

def recorded_correlated(key, pts, theory, body=None):
    """The recorded chi2 with one fully correlated component (see header)."""
    resid = [t - p["value"] for p, t in zip(pts, theory)]
    if key == "nmc":
        diag = [p["comps"]["stat"] ** 2 for p in pts]
        vec = [p["comps"]["sys"] for p in pts]
        what = "sys 100 % correlated over the bins (the extreme of 'largely correlated')"
    elif key == "h1":
        diag = [p["comps"]["stat"] ** 2 + p["comps"]["sys_uncorr"] ** 2 for p in pts]
        vec = [p["comps"]["sys_corr"] for p in pts]
        what = "the correlated component 100 % correlated over the Q2 points"
    elif key == "e140":
        norm = float(body["paper_normalisation_error_r"])
        diag = [p["err"] ** 2 for p in pts]
        vec = [norm for _ in pts]
        what = "the paper's additional normalisation error 0.03 on R, fully correlated"
    else:
        return None
    c = chi2_rank1(resid, diag, vec)
    return dict(chi2=float(c), n=len(pts), p=chi2_sf(c, len(pts)), what=what)


def measure_r_subrow(key, path=None):
    pts = load_points(key, path)
    body = read_vendored(key, path)["body"]
    sel = [p for p in pts if p["dis"]]
    out = dict(name=key, label=LABELS[key], reference=REFERENCES[key],
               n_points=len(pts), n_dis=len(sel),
               n_extrapolated=sum(1 for p in sel if not p["in_r1998_support"]),
               points=pts, backends={})
    for b in BACKENDS:
        rf = r_backend(b)
        theory_all = [rf(p["x"], p["q2"]) for p in pts]
        theory = [rf(p["x"], p["q2"]) for p in sel]
        s = stats(sel, theory)
        cell = dict(s, theory_all=theory_all,
                    status="pass" if s["p"] >= P_MIN else "fail",
                    correlated=recorded_correlated(key, sel, theory, body))
        out["backends"][b] = cell
    if key == "e140":
        out["paper_quotes"] = read_vendored(key, path)["_provenance"]["paper_quotes"]
    if key == "r1990":
        out["rd_minus_rp"] = body["rd_minus_rp"]
    return out


def measure_h1(path=None):
    pts = load_points("h1", path)
    per_bin = load_h1_per_bin(path)
    sel = [p for p in pts if p["dis"]]
    out = dict(name="h1", label=LABELS["h1"], reference=REFERENCES["h1"],
               n_points=len(pts), n_dis=len(sel),
               n_extrapolated=sum(1 for p in sel if not p["in_r1998_support"]),
               points=pts, f2_backends={}, backends={})
    for f2name in H1_F2:
        factory, label, why = f2_backend(f2name)
        if why is not None:
            out["f2_backends"][f2name] = dict(label=label, status="blocked",
                                              reason=why)
            continue
        cells = {}
        for b in BACKENDS:
            fl, f2 = tree_fl(factory, b, sel)
            s = stats(sel, fl)
            fl_bin, _ = tree_fl(factory, b, per_bin)
            cells[b] = dict(s, theory=fl, f2p=f2,
                            status="pass" if s["p"] >= P_MIN else "fail",
                            correlated=recorded_correlated("h1", sel, fl),
                            per_bin=dict((k, v) for k, v in stats(per_bin, fl_bin).items()
                                         if k != "pulls"))
        out["f2_backends"][f2name] = dict(label=label, status="ran", cells=cells)
    verdict_f2 = out["f2_backends"][H1_F2[0]]
    if verdict_f2["status"] == "blocked":
        for b in BACKENDS:
            out["backends"][b] = dict(status="blocked", reason=verdict_f2["reason"])
    else:
        out["backends"] = verdict_f2["cells"]
        # the figure's R view of H1, with the verdict's own F2p
        f2 = verdict_f2["cells"]["r_sigma_lt"]["f2p"]
        out["r_display"] = [
            dict(x=p["x"], q2=p["q2"], f2p=f, r=p["value"] / (f - p["value"]),
                 err=p["err"] * f / (f - p["value"]) ** 2)
            for p, f in zip(sel, f2)]
    out["f2_verdict"] = H1_F2[0]
    return out


# --------------------------------------------------------- the bias table

def bias_cell(x, q2):
    import lipolgen as lg
    r_lt = lg.r_sigma_lt(x, q2)
    r98 = lg.r1998(x, q2)
    r98u = lg.r1998(x, q2, clip=False)
    return dict(x=x, q2=q2, r_sigma_lt=r_lt, r1998=r98, r1998_unclipped=r98u,
                r1998_spread=lg.r1998_spread(x, q2),
                r1998_fit_error=lg.r1998_fit_error(x, q2),
                ratio=(1.0 + r98) / (1.0 + r_lt),
                ratio_unclipped=(1.0 + r98u) / (1.0 + r_lt),
                extrapolated=not (R1998_X_MIN <= x <= R1998_X_MAX
                                  and R1998_Q2_MIN <= q2 <= R1998_Q2_MAX),
                w2=w2_of(x, q2), w2_below_4=w2_of(x, q2) < W2_CUT)


def kernel_factor(cells):
    """a_2(R1998)/a_2(shipped) of the 6Li inclusive kernel at each cell, at
    y = Q2/(s x) of the mid configuration, with the closed form beside it."""
    import lipolgen as lg
    s = lg.default_configs("6Li")[MID_CONFIG].s_per_nucleon()
    label = lg.default_configs("6Li")[MID_CONFIG].label()
    k_ship = lg.default_inclusive_kernel(lg.li6())
    k_unset = lg.default_inclusive_kernel(lg.li6(), model=lg.B1Model.Li6Convolution,
                                          r_source=lg.RSource.Unset)
    k_r98 = lg.default_inclusive_kernel(lg.li6(), model=lg.B1Model.Li6Convolution,
                                        r_source=lg.RSource.R1998)
    st = lg.EventSpinState()
    st.j, st.m, st.theta_s, st.phi_s = 1.0, 1.0, math.pi / 2.0, 0.0
    out = []
    for c in cells:
        x, q2 = c["x"], c["q2"]
        y = q2 / (s * x)
        a_ship = k_ship.amplitudes_at(x, q2, s, st).a2
        a_unset = k_unset.amplitudes_at(x, q2, s, st).a2
        a_r98 = k_r98.amplitudes_at(x, q2, s, st).a2
        if a_unset != a_ship:
            raise RuntimeError("the Delta slot differs between b1 models at x = %g" % x)
        rl, r9 = c["r_sigma_lt"], c["r1998"]
        closed = ((1.0 + 2.0 * (1.0 - y) * (1.0 + rl) / (y * y))
                  / (1.0 + 2.0 * (1.0 - y) * (1.0 + r9) / (y * y)))
        out.append(dict(y=y, a2_shipped=a_ship, a2_r1998=a_r98,
                        kernel_ratio=a_r98 / a_ship, closed_form=closed,
                        inverse_bias=1.0 / c["ratio"],
                        kernel_minus_inverse_bias=a_r98 / a_ship - 1.0 / c["ratio"]))
    return dict(s=s, config=label), out


def nearest_nmc(cell, nmc_points):
    p = min(nmc_points, key=lambda q: abs(math.log(q["x"] / cell["x"])))
    return dict(x=p["x"], q2=p["q2"], r=p["value"], err=p["err"])


def bias_table(nmc_points):
    cells = [dict(bias_cell(c["x"], c["q2"]), label=c["label"]) for c in SWEET_SPOTS]
    meta, kf = kernel_factor(cells)
    for c, k in zip(cells, kf):
        c.update(k)
        c["nmc_nearest"] = nearest_nmc(c, nmc_points)
    grid = [bias_cell(x, q2) for q2 in BIAS_GRID_Q2 for x in BIAS_GRID_X]
    inside = [g["ratio"] for g in grid if not g["extrapolated"] and not g["w2_below_4"]]
    extrap = [g["ratio"] for g in grid if g["extrapolated"] and not g["w2_below_4"]]
    return dict(cells=cells, kernel=meta, grid=grid,
                grid_x=list(BIAS_GRID_X), grid_q2=list(BIAS_GRID_Q2),
                summary=dict(in_support_min=min(inside), in_support_max=max(inside),
                             extrapolated_min=min(extrap), extrapolated_max=max(extrap)))


# --------------------------------------------------------------- the row

def _cell_text(c):
    if c.get("status") == "blocked":
        return "BLOCKED (%s)" % c["reason"]
    return "%.2f/%d (p = %.3g)" % (c["chi2"], c["n"], c["p"])


def _subrow_dict(m):
    """The REPORT.json view of a sub-row: scalars per backend."""
    d = dict(name=m["name"], label=m["label"], n_points=m["n_points"],
             n_dis=m["n_dis"], n_extrapolated_r1998=m["n_extrapolated"])
    st = None
    for b in BACKENDS:
        c = m["backends"].get(b, {})
        if c.get("status") == "blocked":
            d[b] = dict(status="blocked", reason=c["reason"])
            st = st or "blocked"
            continue
        d[b] = dict(chi2=c["chi2"], ndf=c["n"], p=c["p"], status=c["status"],
                    worst_pull=c["worst_pull"], worst_x=c["worst_x"],
                    worst_q2=c["worst_q2"], mean_diff=c["mean_diff"])
        if c.get("correlated"):
            d[b]["chi2_correlated"] = c["correlated"]["chi2"]
            d[b]["p_correlated"] = c["correlated"]["p"]
    d["status"] = (m["backends"]["r_sigma_lt"].get("status")
                   if m["backends"].get("r_sigma_lt") else "blocked")
    if m["name"] == "h1":
        d["f2_verdict"] = m["f2_verdict"]
        for f2name, fb in m["f2_backends"].items():
            if fb["status"] == "blocked":
                d["f2_" + f2name] = dict(status="blocked", reason=fb["reason"])
                continue
            d["f2_" + f2name] = {b: "%.2f/%d (p = %.3g)" % (c["chi2"], c["n"], c["p"])
                                 for b, c in fb["cells"].items()}
    return d


def run(paths=None, verbose=True):
    paths = paths or {}
    t0 = time.perf_counter()
    support = check_r1998_support()
    subs = {}
    for key in SUBROWS:
        subs[key] = (measure_h1(paths.get(key)) if key == "h1"
                     else measure_r_subrow(key, paths.get(key)))
    bias = bias_table(subs["nmc"]["points"])
    head = subs[HEADLINE[0]]["backends"][HEADLINE[1]]
    status = head["status"]
    import lipolgen as lg
    subrows = [_subrow_dict(subs[k]) for k in SUBROWS]
    for c in bias["cells"]:
        subrows.append(dict(
            name="bias[%s x=%g Q2=%g]" % (c["label"], c["x"], c["q2"]),
            status="recorded", x=c["x"], q2=c["q2"], y_mid_config=c["y"],
            r_sigma_lt=c["r_sigma_lt"], r1998=c["r1998"],
            r1998_spread=c["r1998_spread"], r1998_fit_error=c["r1998_fit_error"],
            ratio_1pR=c["ratio"], kernel_a2_ratio=c["kernel_ratio"],
            closed_form_a2_ratio=c["closed_form"], in_r1998_support=not c["extrapolated"],
            nmc_nearest_r=c["nmc_nearest"]["r"], nmc_nearest_x=c["nmc_nearest"]["x"],
            nmc_nearest_q2=c["nmc_nearest"]["q2"], nmc_nearest_err=c["nmc_nearest"]["err"]))
    others = "; ".join(
        "%s[%s] %s" % (k, b, _cell_text(subs[k]["backends"][b]))
        for k in SUBROWS for b in ("r_sigma_lt", "r1998")
        if k != HEADLINE[0])
    report = dict(
        name=NAME,
        generator=("r_sigma_lt (shipped R) on NMC R: chi2/ndf = %.2f/%d, p = %.3g "
                   "[diagonal, DIS cut]; r1998 (clip) on NMC %s; %s; bias "
                   "(1+r1998)/(1+r_sigma_lt) at the sweet spots %s"
                   % (head["chi2"], head["n"], head["p"],
                      _cell_text(subs["nmc"]["backends"]["r1998"]), others,
                      " / ".join("%.4f" % c["ratio"] for c in bias["cells"]))),
        reference="; ".join("%s (%s, %d of %d points in the DIS cut)"
                            % (REFERENCES[k], LABELS[k], subs[k]["n_dis"],
                               subs[k]["n_points"]) for k in SUBROWS),
        tolerance=("p(chi2, ndf) >= %.2f per data set on Q2 >= 1, W2 >= 4 GeV2, "
                   "stat (+) sys diagonal; status = r_sigma_lt on the PRIMARY nmc set"
                   % P_MIN),
        status=status, headline="%s[%s]" % HEADLINE,
        subrows=subrows, detail=subs, bias=bias, r1998_support=support,
        dis_cut="Q2 >= %g GeV2 and W2 >= %g GeV2" % (Q2_CUT, W2_CUT),
        config=dict(lipolgen=getattr(lg, "__version__", None),
                    python=sys.version.split()[0],
                    data_file_sha256={k: v[1] for k, v in FILES.items()},
                    kernel_config=bias["kernel"]["config"],
                    kernel_s=bias["kernel"]["s"]),
    )
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        _print(report)
    return report


def _print(rep):
    subs, bias = rep["detail"], rep["bias"]
    print("R = sigma_L/sigma_T: r_sigma_lt (shipped), r1998 (clip, shipped), "
          "r1998 unclipped (recorded) -- DIS cut %s; R1998 support x in [%g, %g], "
          "Q2 in [%g, %g]" % (rep["dis_cut"], R1998_X_MIN, R1998_X_MAX,
                               R1998_Q2_MIN, R1998_Q2_MAX))
    for key in SUBROWS:
        m = subs[key]
        print("\n[%s] %s -- %s; %d points, %d in the DIS cut, %d of them outside "
              "R1998's support" % (key, m["label"], m["reference"], m["n_points"],
                                   m["n_dis"], m["n_extrapolated"]))
        if key == "h1":
            hb = m["backends"]
            print("     x        Q2     F_L +- diag     F2p(%s)  tree F_L: %s"
                  % (m["f2_verdict"], "  ".join("%10s" % b for b in BACKENDS)))
            sel = [p for p in m["points"] if p["dis"]]
            for j, p in enumerate(sel):
                f2 = hb["r_sigma_lt"]["f2p"][j] if "f2p" in hb.get("r_sigma_lt", {}) else float("nan")
                print("%9.5f %6.1f   %6.3f +- %.3f   %7.4f            %s  extrapolated"
                      % (p["x"], p["q2"], p["value"], p["err"], f2,
                         "  ".join("%10.4f" % hb[b]["theory"][j] if "theory" in hb.get(b, {})
                                   else "%10s" % "blocked" for b in BACKENDS)))
            for f2name, fb in m["f2_backends"].items():
                if fb["status"] == "blocked":
                    print("  F2p from %s: BLOCKED (%s)" % (f2name, fb["reason"]))
                    continue
                print("  F2p from %-32s %s" % (fb["label"] + (" [VERDICT]" if f2name == m["f2_verdict"]
                                                               else " [recorded]") + ":",
                                              "; ".join("%s %s, corr-sys %.2f, per-bin %.2f/%d"
                                                        % (b, _cell_text(c), c["correlated"]["chi2"],
                                                           c["per_bin"]["chi2"], c["per_bin"]["n"])
                                                        for b, c in fb["cells"].items())))
        else:
            print("     x        Q2      W2     R +- diag       %s" %
                  "  ".join("%14s" % b for b in BACKENDS))
            for j, p in enumerate(m["points"]):
                print("%8.4f %7.2f %7.2f   %.3f +- %.3f %s%s  %s" % (
                    p["x"], p["q2"], p["w2"], p["value"], p["err"],
                    " " if p["dis"] else "*", " " if p["in_r1998_support"] else "e",
                    "  ".join("%14.4f" % m["backends"][b]["theory_all"][j] for b in BACKENDS)))
            print("  (* = outside the DIS cut, not compared; e = outside R1998's support)")
        for b in BACKENDS:
            c = m["backends"][b]
            if c.get("status") == "blocked":
                print("  %-16s BLOCKED (%s)" % (b, c["reason"]))
                continue
            cor = c.get("correlated")
            print("  %-16s chi2/ndf = %7.2f/%-3d p = %-10.3g %s  worst pull %+.2f at x = %g, "
                  "Q2 = %g; <theory - data> = %+.4f%s"
                  % (b, c["chi2"], c["n"], c["p"], c["status"].upper(), c["worst_pull"],
                     c["worst_x"], c["worst_q2"], c["mean_diff"],
                     "" if not cor else "; recorded (%s): %.2f/%d, p = %.3g"
                     % (cor["what"], cor["chi2"], cor["n"], cor["p"])))
        if key == "e140":
            q = m["paper_quotes"]["ra_minus_rd_average"]
            print("  recorded: E140's average R_A - R_D = %.3f +- %.3f (stat) +- %.3f (sys) "
                  "(paper sec. IV.B); the tree's R has no A dependence (0 by construction)"
                  % (q["value"], q["stat"], q["sys"]))
        if key == "r1990":
            q = m["rd_minus_rp"]
            print("  recorded: Whitlow's R_d - R_p = %.3f +- %.3f +- %.3f (record Table 1); "
                  "the tree's R has no target dependence" % (q["value"], q["stat"], q["sys"]))
    print("\nBIAS TABLE -- (1 + r1998)/(1 + r_sigma_lt); the kernel's cos 2phi amplitude "
          "ratio a2(R1998)/a2(shipped) at %s (s = %.1f GeV2)" % (bias["kernel"]["config"],
                                                                 bias["kernel"]["s"]))
    for c in bias["cells"]:
        n = c["nmc_nearest"]
        print("  %s x = %-7g Q2 = %-5g y = %.5f: R_lt %.4f, R98 %.4f (spread %.4f, fit err "
              "%.4f) -> (1+R98)/(1+R_lt) = %.4f; kernel a2 ratio %.6f, closed form %.6f, "
              "1/bias %.6f; nearest NMC R %.3f +- %.3f at x = %g, <Q2> = %g"
              % (c["label"], c["x"], c["q2"], c["y"], c["r_sigma_lt"], c["r1998"],
                 c["r1998_spread"], c["r1998_fit_error"], c["ratio"], c["kernel_ratio"],
                 c["closed_form"], c["inverse_bias"], n["r"], n["err"], n["x"], n["q2"]))
    print("  grid (rows Q2, columns x; e = x below R1998's support, frozen; w = W2 < 4):")
    print("  Q2 \\ x " + " ".join("%7g" % x for x in BIAS_GRID_X))
    for q2 in BIAS_GRID_Q2:
        row = [g for g in bias["grid"] if g["q2"] == q2]
        print("  %6g " % q2 + " ".join(
            "%6.3f%s" % (g["ratio"], "e" if g["extrapolated"] else ("w" if g["w2_below_4"] else " "))
            for g in row))
    s = bias["summary"]
    print("  in support (W2 >= 4): %.4f ... %.4f; extrapolated x < %g: %.4f ... %.4f"
          % (s["in_support_min"], s["in_support_max"], R1998_X_MIN,
             s["extrapolated_min"], s["extrapolated_max"]))
    print()
    for sr in rep["subrows"]:
        if sr["name"].startswith("bias["):
            print("SUBROW | %s[%s] | (1+r1998)/(1+r_sigma_lt) = %.4f (R98 %.4f, R_lt %.4f); "
                  "kernel cos 2phi a2(R1998)/a2(shipped) = %.6f at y = %.5f | plans/07 "
                  "sec. 7.4 WP1 cell; nearest NMC R %.3f +- %.3f (x %g, <Q2> %g) | none "
                  "(a price, not a test) | RECORDED"
                  % (NAME, sr["name"], sr["ratio_1pR"], sr["r1998"], sr["r_sigma_lt"],
                     sr["kernel_a2_ratio"], sr["y_mid_config"], sr["nmc_nearest_r"],
                     sr["nmc_nearest_err"], sr["nmc_nearest_x"], sr["nmc_nearest_q2"]))
            continue
        cells = "; ".join("%s %s" % (b, _cell_text(subs[sr["name"]]["backends"][b]))
                          for b in BACKENDS)
        print("SUBROW | %s[%s] | %s: %s | %s | %s | %s" % (
            NAME, sr["name"], sr["label"], cells, REFERENCES[sr["name"]],
            rep["tolerance"], str(sr["status"]).upper()))
    print("REPORT | %s | %s | %s | %s | %s" % (
        NAME, rep["generator"], rep["reference"], rep["tolerance"],
        rep["status"].upper()))


# ------------------------------------------------------------------ figure

def plot(report, out_dir=None):
    """R vs x: the data (ink, one marker per set; H1 shown in R with the
    verdict's own F2p), r_sigma_lt and r1998 at each point's own Q2 (lines
    through the NMC and H1 points in x order, ticks at the SLAC points), a
    pull panel (the verdict chi2's own pulls)."""
    sys.path.insert(0, HERE)
    from _plots import new_figure, save, status_title, style_series, MUTED, INK
    fig, (ax, axp) = new_figure(two_panel=True, height=4.2)
    subs = report["detail"]
    markers = {"nmc": "o", "h1": "D", "r1990": "s", "e140": "^"}
    series = (("r_sigma_lt", "default", "r_sigma_lt (shipped)"),
              ("r1998", "r1998", "r1998 (clip, as shipped)"))
    for key in SUBROWS:
        m = subs[key]
        if key == "h1":
            if "r_display" not in m:
                continue
            xs = [d["x"] for d in m["r_display"]]
            ys = [d["r"] for d in m["r_display"]]
            es = [d["err"] for d in m["r_display"]]
            lab = "H1 F_L -> R (F2p MSTW; %s)" % LABELS[key]
        else:
            sel = [p for p in m["points"] if p["dis"]]
            xs = [p["x"] for p in sel]
            ys = [p["value"] for p in sel]
            es = [p["err"] for p in sel]
            lab = {"nmc": "NMC R (PRIMARY)", "r1990": "Whitlow R1990 data (CLOSURE)",
                   "e140": "E140 R, D (CLOSURE)"}[key]
        kw = style_series("reference", label=lab, marker=markers[key])
        if key == "r1990":
            kw.update(color=MUTED, markersize=3, zorder=4)
        ax.errorbar(xs, ys, yerr=es, **kw)
    for b, slot, lab in series:
        first = True
        for key in SUBROWS:
            m = subs[key]
            c = m["backends"].get(b, {})
            if c.get("status") == "blocked":
                continue
            if key == "h1":
                sel = [p for p in m["points"] if p["dis"]]
                th = [r_backend(b)(p["x"], p["q2"]) for p in sel]
            else:
                sel = [p for p in m["points"] if p["dis"]]
                th = [t for p, t in zip(m["points"], c["theory_all"]) if p["dis"]]
            order = sorted(range(len(sel)), key=lambda j: sel[j]["x"])
            xs = [sel[j]["x"] for j in order]
            ys = [th[j] for j in order]
            kw = style_series(slot, label=lab)
            if not first:
                kw["label"] = "_nolegend_"
            if key in ("nmc", "h1"):
                ax.plot(xs, ys, **kw)
            else:
                kw.update(marker="_", linestyle="none", markersize=9)
                ax.plot(xs, ys, **kw)
            first = False
            pk = style_series(slot, marker=markers[key])
            pk.update(markersize=4, label="_nolegend_")
            axp.plot([p["x"] for p in sel], c["pulls"], **pk)
    for c in report["bias"]["cells"]:
        ax.axvline(c["x"], color=MUTED, linewidth=0.8, linestyle=":", zorder=0)
    ax.axvline(R1998_X_MIN, color=INK, linewidth=0.8, linestyle="--", zorder=0)
    ax.text(R1998_X_MIN * 0.93, 0.03, "r1998 support edge", rotation=90, ha="right",
            va="bottom", fontsize=7, color=MUTED, transform=ax.get_xaxis_transform())
    ax.text(SWEET_SPOTS[0]["x"] * 0.95, 0.03, "sweet spots", rotation=90, ha="right",
            va="bottom", fontsize=7, color=MUTED, transform=ax.get_xaxis_transform())
    ax.set_xscale("log")
    ax.set_ylabel("R = sigma_L / sigma_T")
    ax.set_ylim(-0.45, 1.85)
    ax.legend(loc="upper center", ncol=3, fontsize=7)
    axp.set_xlabel("x (each point at its own Q2)")
    axp.set_ylabel("pull (th - data)/sigma")
    h = report["detail"]["nmc"]["backends"]["r_sigma_lt"]
    status_title(fig, NAME, report["status"],
                 "r_sigma_lt on NMC R: %.1f/%d (p = %.2g); tolerance p >= 0.01 per set, "
                 "Q2 >= 1, W2 >= 4 GeV2" % (h["chi2"], h["n"], h["p"]))
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
