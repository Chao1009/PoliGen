#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T5 MC closures -- the polarized sum rules of 05_chain.md sec. 1.2 (CH-25)
and the FSI unitarity sum of 05 sec. 1.1 (tier T5).  Three sub-rows:

  bc_ww                Burkhardt-Cottingham on the tree's `g2_ww`
                       (asserted; PASS)
  bjorken[toyg1]       Bjorken sum on the shipped `ToyG1`
                       (asserted; FAIL, recorded)
  bjorken[nnpdfpol11]  the same on `LhapdfG1("NNPDFpol11_100")`
                       (BLOCKED: environment, the set is not installed)
  fsi_unitarity        Strikman-Weiss int dGamma S[FSI + FSI^2] = 0
                       (BLOCKED: not definable in the tree as built)

WHAT IS COMPARED.

(a) bc_ww -- AN IMPLEMENTATION CLOSURE, NOT A PHYSICS TEST.  The tree builds
g2 by Wandzura-Wilczek, g2^WW(x) = -g1(x) + int_x^1 g1(u)/u du
(`lipolgen.g2_ww`, src/core/sf.cpp; the substitution u = x^(1-t) and a
96-point trapezoid in t).  For ANY g1, exchanging the order of integration
gives the exact truncated identity

    int_{x_min}^1 g2^WW(x) dx  =  T(x_min)  =  -x_min int_{x_min}^1 g1(u)/u du ,

whose x_min -> 0 limit is the Burkhardt-Cottingham sum rule int_0^1 g2 dx = 0
-- BC HOLDS IDENTICALLY BY CONSTRUCTION ON g2^WW (05 sec. 1.2), provided
int_0^1 g1 dx exists.  The harness integrates the TREE's g2_ww over
[x_min, 1] (scipy `quad` in ln x, piecewise between successive x_min) and
compares it with T(x_min), computed independently of g2_ww (quad of
g1(u)/u).  It therefore catches a
wrong WW integral -- a lost Jacobian -ln x, a wrong range, a u^-2 kernel --
and nothing about physics: the twist-3 part of g2 that BC is sensitive to is
zero by construction here.  Targets: `ToyG1()` proton (`g2p`), neutron
(`g2n`) and deuteron (`lipolgen.g2_ww` on `ToyG1().g1_nucleus(deuteron())`,
the free function through its Python binding), at Q2 = 2, 5, 10 GeV^2 and
x_min = 1e-2 ... 1e-6 (45 cells).

    TOLERANCE (fixed before any residual was computed).  Per cell,
    tol = B(x_min) + e_I + e_T, where e_I, e_T are quad's own error
    estimates of the two integrals, and B is the RIGOROUS Peano-kernel bound
    of the tree's inner trapezoid: on each panel the trapezoid error is
    -int f''(t) (t - t_i)(t_i+1 - t)/2 dt with a kernel <= h^2/8, so with
    f(t) = L g1(x^(1-t)), L = -ln x, h = 1/(npts - 1):
        |g2_tree(x) - g2^WW(x)| <= (h^2/8) int_0^1 |f''| dt
                                 = (h^2/8) L^2 int_{ln x}^0 |d^2 g1/ds^2| ds ,
    integrated over x in [x_min, 1].  d^2 g1/ds^2 (s = ln u) is taken by
    second differences on a uniform s grid of >= 4000 points per unit of s,
    so the bound's own discretisation error is O(1e-7) relative.  The
    Euler-Maclaurin leading term (h^2/12) L^2 [Dg1(1) - Dg1(x)] (D = d/ds)
    is RECORDED as the predicted residual; |predicted|/bound ~ 2/3 for a
    power-law g1, so the bound is tight, not padded.  The bound follows
    from `npts` alone: it is not tuned to 96 (pytest re-runs the closure at
    npts = 8).

    THE LIMIT, recorded and not asserted.  T(x_min) -> 0 for the toy PROTON
    (g1p ~ x^-0.5).  For the toy NEUTRON and hence the DEUTERON it does NOT:
    ToyG1's a1n(0) = -0.07 times F1n ~ x^(-1-lambda(Q2)) (ToyF2's sea,
    lambda = 0.045 ln(max(Q2, 1.1)/0.04) = 0.176 / 0.217 / 0.248 at
    Q2 = 2 / 5 / 10) makes g1n ~ x^(-1-lambda) NOT INTEGRABLE, so
    int_0^1 g2n^WW dx does not exist and T(x_min) grows like x_min^-lambda.
    The harness reports the effective exponent from the last decade.  This
    is a property of the toy's small-x shape, not of the WW code.

(b) bjorken -- A MEASUREMENT OF THE SHIPPED MODEL.  Gamma_1^{p-n}(Q2) =
int_0^1 (g1p - g1n) dx against the leading-twist Bjorken sum rule
    Gamma_ref = (g_A/6) C_Bj(a),  a = alpha_s(Q2)/pi (MS-bar),
    C_Bj = 1 - a - (55/12 - nf/3) a^2 - c3(nf) a^3
(Gorishny-Larin 1986; Larin-Vermaseren, PLB 259 (1991) 345; c3 = 20.215 at
nf = 3, 13.850 at nf = 4), nf = 4 above m_c.  alpha_s(Q) is run here, by the
3-loop MS-bar beta function (RK4 in ln mu^2), from alpha_s(M_Z), nf = 5 -> 4
at mu = m_b, continuous at the threshold (the O(alpha_s^3) matching step is
neglected); target-mass and higher-twist (1/Q2) terms are NOT included.  The
reference's own uncertainty (g_A, alpha_s(M_Z), and the size of the last
included term) is reported; it is 2.6 / 1.3 / 0.8 % at Q2 = 2 / 5 / 10,
the truncation term dominating.  In 05's words the sum rule is a
property of the PDF SET, not of the generator: on `ToyG1` this sub-row
measures what the shipped default g1 does; on `LhapdfG1` it would check the
isovector flavour plumbing.

    TOLERANCE (declared before computing; see the honesty note below).
    For each Q2, over the x_min scan (toy: 1e-2 ... 1e-6; LHAPDF: 1e-2 ...
    1e-5, not below a typical grid edge) PASS iff (i) the scan CONVERGES:
    the last decade adds <= 1 % of Gamma_ref, and (ii) |Gamma(x_min,last) -
    Gamma_ref| <= 10 % of Gamma_ref -- about the precision to which world
    data test the sum rule and generous for a toy, well above the reference's
    own 0.8-2.6 %.  The sub-row passes iff every Q2 passes.  A divergent
    integral fails any finite tolerance.

    THE REFERENCE CONSTANTS are embedded below (`REFERENCE`), not vendored
    under data/, because this harness's task named three files only; they
    carry a sha256 of their canonical JSON that `load_reference()` re-checks
    on every read (an altered constant -> RuntimeError -> exit 2).
    Provenance: PDG 2026 (F. Takahashi et al. (Particle Data Group), Int. J.
    Mod. Phys. A 41, 2630011 (2026)), read 2026-09-26 from the `pdg.sqlite`
    inside the PyPI wheel pdg-2026.0-py3-none-any.whl (PDG's own Python API,
    https://files.pythonhosted.org via https://pypi.org/project/pdg/, wheel
    sha256 681e11f8c9a5accb1cb41ccb87bf7efd398adec8ea39ddb894630b17b45c8e8a;
    data licence CC BY 4.0 per its pdginfo table, package Modified BSD;
    pdg.lbl.gov itself answers 403 through this machine's proxy):
      S017AV  lambda = g_A/g_V = -1.2753 +- 0.0013 (S = 2.7; sign dropped
              here, Gamma_1^{p-n} > 0 is g_A/6 with g_A = |lambda|);
      S044M   M_Z = 91.1879 +- 0.0020 GeV;
      Q005M   m_b(m_b) = 4.186 +- 0.006 GeV;  Q004M  m_c(m_c) = 1.2729 +-
              0.0045 GeV (MS-bar; summary-table values).
    alpha_s(M_Z) = 0.1180 +- 0.0009 is NOT in that database (it lives in the
    QCD review): it is the PDG 2024 world average as the author knows it,
    NOT re-read, and is the least verified input of the row.

    LhapdfG1 NOTE, for when the set is installed (never run here): LhapdfG1
    is LO g1 (no coefficient functions), so its first moment is a3/6 WITHOUT
    C_Bj; if NNPDFpol1.1 carries a3 ~ g_A it will sit ~ (1/C_Bj - 1), i.e.
    ~ +14 % at Q2 = 5, above Gamma_ref -- a convention of the LO g1, not a
    flavour bug.  Not verified here.

(c) fsi_unitarity -- BLOCKED: "not definable in the tree as built".
Strikman-Weiss's sum rule (04_theory.md T-29) says the FSI terms integrate
to zero over the spectator phase space.  No option of `GlauberFsiWeight`
(include/lipolgen/fsi.hpp) realises such a sum: (1) the default keeps the
absorptive loss on purpose -- the linear term carries sigma_tot(X-cluster),
the quadratic sigma_el, ~73 % of rescatterings destroy the alpha tag, and
int w dGamma / int dGamma is a SURVIVAL probability (`survival()`, ~0.52),
not 1; (2) `elastic_gain` rescales the quadratic term alone, and a value that
zeroes int (w - 1) dGamma for one channel at one sigma_XN would be a TUNE,
not a realisation; (3) `weight_normalised` integrates to 1 by DIVISION
(tests/test_fsi.cpp), a normalisation, not a unitarity sum; (4) the weight is
tabulated only on k <= k_max = 1.2 GeV and clipped at w_max = 50, so the full
spectator integral the sum rule needs is not available.  The tree side is
printed for the record (survival, sigma_tot and sigma_el of the default
6Li alpha channel) and compared with nothing.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: no b1 anywhere in this row, and no
     A-scaling assumption (no 0.5, no 2/6) is used.
  2. isoscalar denominator -- NOT APPLICABLE: no nuclear-to-nucleon ratio;
     the Bjorken sum is the isovector p - n difference of FREE-nucleon g1.
  3. per nucleon / per nucleus -- bc_ww: g1p, g1n and g2p, g2n are PER
     NUCLEON; the deuteron's `g1_nucleus` is PER NUCLEUS, Z P_p g1p + N P_n
     g1n with P_p = P_n = 0.9325 = 1 - 1.5 P_D (NOT divided by A = 2; the
     inclusive kernel's g1a divides by A).  The closure is linear and
     homogeneous in g1, so the normalisation cancels between its two sides;
     it is stated, and pytest asserts the composition.  bjorken: per nucleon,
     free p and n.  fsi_unitarity: a per-event weight ratio, dimensionless.

WHAT THIS DOES NOT VALIDATE: the physics of g2 (twist-3 is zero by
construction; measured g2 -- E155x, RSS, SANE -- would test the WW
approximation, a structure-function-model question, 05 sec. 1.2); the tensor
sector (neither sum rule touches b1..b4); any nucleus beyond A = 2; the
LhapdfG1 flavour plumbing (blocked here); FSI unitarity (not definable).

HONESTY NOTE.  Before these tolerances were fixed the author had seen: the
scout's preview "BC on ToyG1 proton at Q2 = 5: int g2^WW = -6.3e-4 against
int g1 = 0.128" (a truncated integral, x_min not stated); ToyG1's g1/g2 at
nine x values at Q2 = 5 (g1n(1e-5) = -7344); hence that the toy's neutron
is not integrable and that the Bjorken sub-row MUST fail on the toy at any
finite tolerance; and the FSI survival 0.520.  No residual of the bc_ww
closure and no Bjorken integral was computed before the tolerances above
were written down.  The reference constants were first typed from memory
(PDG 2024: 1.2754, 91.1876, 4.18, 1.27) and then replaced by the PDG 2026
database values above; that moved Gamma_ref by 2e-5 and no verdict.

WHAT WAS MEASURED (2026-09-26; 1.6 s wall).  bc_ww PASS: 45/45 cells
within tolerance, max |resid|/tol = 0.667 -- every residual is the tree's
96-point trapezoid error, matching its Euler-Maclaurin prediction to
resid/EM = 0.9994 ... 1.0002 (e.g. proton, Q2 = 5, x_min = 1e-6:
int g2_tree = -2.035413e-4 against T = -2.111946e-4, resid +7.65e-6,
tol 1.15e-5).  Limit: |T| ~ x_min^+0.45..0.53 -> 0 for the proton; for the
neutron |T| grows as x_min^-0.175 / -0.217 / -0.248 at Q2 = 2 / 5 / 10 (the
toy's lambda to 1e-3), the deuteron likewise -- no BC limit.  bjorken[toyg1]
FAIL (recorded): int_{x_min}^1 (g1p - g1n) = 0.139 / 0.203 / 0.289 / 0.423 /
0.643 at x_min = 1e-2 ... 1e-6, Q2 = 5, against Gamma_ref = 0.18538 +- 0.00233
(alpha_s = 0.2850, nf = 4, C_Bj = 0.8722) -- the last decade alone adds
0.220 (119 % of the reference); Gamma^p converges to 0.128 (the scout's
number), Gamma^n does not.  bjorken[nnpdfpol11] BLOCKED (environment);
fsi_unitarity BLOCKED (not definable); survival recorded 0.520239.  The row
as a whole is therefore FAIL, recorded: the status rule is "fail if any
evaluated sub-row fails, pass if all evaluated pass, blocked if none was
evaluated", and the failing sub-row is the shipped toy's Bjorken integral,
not the WW implementation.

Run:  source env.sh && python3 validation/benchmarks/t5_sum_rules.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT rows -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (an altered
reference constant, a missing dependency, any exception).  A harness is a
measurement, not a CI gate: the pytests are the gate.
"""

import hashlib
import json
import math
import sys
import time
import traceback

NAME = "t5_sum_rules"

Q2_POINTS = (2.0, 5.0, 10.0)                       # GeV^2
X_MIN_SCAN = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6)        # toy scan (bc_ww, bjorken)
X_MIN_SCAN_LHAPDF = (1e-2, 1e-3, 1e-4, 1e-5)       # not below a grid edge
G2_NPTS = 96              # g2_ww's default = InclusiveKernel::Options::g2_npts
BOUND_POINTS_PER_UNIT_S = 4000                     # s = ln u grid of the bound
QUAD_EPSABS = 1e-15
QUAD_EPSREL = 1e-12

BJ_TOL_REL = 0.10         # |Gamma - Gamma_ref| <= 10 % of Gamma_ref
BJ_CONV_REL = 0.01        # last decade adds <= 1 % of Gamma_ref
LHAPDF_G1_SET = "NNPDFpol11_100"

FSI_BLOCKED_REASON = "not definable in the tree as built"

#: The Bjorken reference inputs (see the module docstring for provenance).
REFERENCE = {
    "g_A": 1.2753, "g_A_err": 0.0013,
    "alpha_s_mz": 0.1180, "alpha_s_mz_err": 0.0009,
    "m_z_gev": 91.1879, "m_b_gev": 4.186, "m_c_gev": 1.2729,
}
#: sha256 of json.dumps(REFERENCE, sort_keys=True), recorded 2026-09-26.
REFERENCE_SHA256 = "a91f1488481fdba9c3f2a3219e878a217e6e9e9479ad746b34159d1b007567bd"

ZETA3 = 1.2020569031595942
ZETA5 = 1.0369277551433699


def _canonical(ref):
    return json.dumps(ref, sort_keys=True).encode("utf-8")


def load_reference(ref=None, sha256=None):
    """The Bjorken reference constants -- RuntimeError unless the sha256 of
    their canonical JSON is the recorded one, so a moved number is refused."""
    ref = REFERENCE if ref is None else ref
    sha256 = REFERENCE_SHA256 if sha256 is None else sha256
    got = hashlib.sha256(_canonical(ref)).hexdigest()
    if got != sha256:
        raise RuntimeError("%s: reference constants sha256 %s != recorded %s"
                           % (NAME, got, sha256))
    return dict(ref)


# ------------------------------------------------------------ quadrature

def _quad(f, a, b):
    from scipy.integrate import quad
    v, e = quad(f, a, b, epsabs=QUAD_EPSABS, epsrel=QUAD_EPSREL, limit=400)
    return float(v), float(e)


def _decades(x_mins):
    """[(hi, lo), ...] from 1 down through the sorted scan."""
    edges = [1.0] + sorted(x_mins, reverse=True)
    return list(zip(edges[:-1], edges[1:]))


def _s_grid(x_lo, per_unit=BOUND_POINTS_PER_UNIT_S, x_mins=X_MIN_SCAN):
    """Uniform-per-piece grid in s = ln u on [ln x_lo, 0] whose nodes include
    every ln x_min, so cumulative integrals are exact at the scan points."""
    import numpy as np
    knots = sorted({0.0, math.log(x_lo)} | {math.log(x) for x in x_mins
                                            if x >= x_lo})
    pieces = []
    for a, b in zip(knots[:-1], knots[1:]):
        n = int(math.ceil((b - a) * per_unit)) + 1
        pieces.append(np.linspace(a, b, n)[:-1])
    pieces.append(np.array([0.0]))
    return np.concatenate(pieces)


def _cumtrapz_from_top(y, s):
    """C[j] = int_{s_j}^{s_last} y ds by the trapezoid rule."""
    import numpy as np
    seg = 0.5 * (y[1:] + y[:-1]) * np.diff(s)
    c = np.zeros_like(y)
    c[:-1] = np.cumsum(seg[::-1])[::-1]
    return c


def trapezoid_bound(g1, x_mins=X_MIN_SCAN, npts=G2_NPTS,
                    per_unit=BOUND_POINTS_PER_UNIT_S):
    """For each x_min: (B, EM) where B = int_{x_min}^1 (h^2/8) L^2
    int_{ln x}^0 |D^2 g1| ds dx is the rigorous bound on the tree's inner
    trapezoid integrated over x, and EM the Euler-Maclaurin leading
    prediction int_{x_min}^1 (h^2/12) L^2 [Dg1(1) - Dg1(x)] dx of the
    residual  int g2_tree - T."""
    import numpy as np
    s = _s_grid(min(x_mins), per_unit, x_mins)
    g = np.array([g1(math.exp(v)) for v in s])
    ds = np.diff(s)
    # second derivative in s by second differences on the (piecewise)
    # uniform grid; three-point formula valid for unequal spacing
    d2 = np.empty_like(g)
    h0, h1 = ds[:-1], ds[1:]
    d2[1:-1] = 2.0 * (h0 * g[2:] - (h0 + h1) * g[1:-1] + h1 * g[:-2]) \
        / (h0 * h1 * (h0 + h1))
    d2[0], d2[-1] = d2[1], d2[-2]
    d1 = np.gradient(g, s, edge_order=2)
    h = 1.0 / (npts - 1)
    c_abs = _cumtrapz_from_top(np.abs(d2), s)           # int_s^0 |D^2 g1|
    xs = np.exp(s)
    b_int = (h * h / 8.0) * s * s * c_abs * xs          # per ds of the x integral
    em_int = (h * h / 12.0) * s * s * (d1[-1] - d1) * xs
    b_cum = _cumtrapz_from_top(b_int, s)
    em_cum = _cumtrapz_from_top(em_int, s)
    out = {}
    for xm in x_mins:
        j = int(np.argmin(np.abs(s - math.log(xm))))
        out[xm] = (float(b_cum[j]), float(em_cum[j]))
    return out


def bc_closure(g1, g2, x_mins=X_MIN_SCAN, npts=G2_NPTS):
    """The bc_ww closure for one target at one Q2.  `g1`, `g2` are callables
    of x alone; `g2` is the TREE's g2 (built with `npts` trapezoid points).
    Returns one dict per x_min, largest first."""
    bound = trapezoid_bound(g1, x_mins, npts)
    rows = []
    i_tot = i_err = gg = g_err = 0.0
    for hi, lo in _decades(x_mins):
        a, b = math.log(lo), math.log(hi)
        v, e = _quad(lambda t: g2(math.exp(t)) * math.exp(t), a, b)
        i_tot += v
        i_err += e
        v, e = _quad(lambda t: g1(math.exp(t)), a, b)   # int g1(u)/u du
        gg += v
        g_err += e
        trunc = -lo * gg
        resid = i_tot - trunc
        bnd, em = bound[lo]
        tol = bnd + i_err + lo * g_err
        rows.append(dict(
            x_min=lo, integral=i_tot, truncation=trunc, residual=resid,
            bound=bnd, quad_err=i_err + lo * g_err, tol=tol,
            em_prediction=em, ratio=abs(resid) / tol,
            ok=abs(resid) <= tol))
    return rows


def _eff_exponent(rows):
    """-d ln|T| / d ln x_min over the last decade: > 0 means T grows as
    x_min -> 0 (no BC limit), < 0 that it vanishes."""
    t1, t0 = rows[-1]["truncation"], rows[-2]["truncation"]
    x1, x0 = rows[-1]["x_min"], rows[-2]["x_min"]
    if t1 == 0.0 or t0 == 0.0 or (t1 > 0) != (t0 > 0):
        return float("nan")
    return -math.log(abs(t1) / abs(t0)) / math.log(x1 / x0)


def measure_bc(npts=G2_NPTS, q2s=Q2_POINTS, x_mins=X_MIN_SCAN):
    import lipolgen as lg
    toy = lg.ToyG1()
    d = lg.deuteron()
    out = {}
    for q2 in q2s:
        targets = {
            "p": (lambda x, q=q2: toy.g1p(x, q),
                  lambda x, q=q2: toy.g2p(x, q, npts)),
            "n": (lambda x, q=q2: toy.g1n(x, q),
                  lambda x, q=q2: toy.g2n(x, q, npts)),
            "d": (lambda x, q=q2: toy.g1_nucleus(d, x, q),
                  lambda x, q=q2: lg.g2_ww(
                      lambda xx, qq: toy.g1_nucleus(d, xx, qq), x, q, npts)),
        }
        for tgt, (g1, g2) in targets.items():
            rows = bc_closure(g1, g2, x_mins, npts)
            out[(tgt, q2)] = dict(rows=rows, exponent=_eff_exponent(rows),
                                  toy_lambda=0.045 * math.log(max(q2, 1.1) / 0.04))
    return out


# ------------------------------------------------------------------ Bjorken

def _beta(nf):
    pi = math.pi
    b0 = (33.0 - 2.0 * nf) / (12.0 * pi)
    b1 = (153.0 - 19.0 * nf) / (24.0 * pi ** 2)
    b2 = (2857.0 - 5033.0 * nf / 9.0 + 325.0 * nf * nf / 27.0) / (128.0 * pi ** 3)
    return b0, b1, b2


def _run_alpha(alpha, mu_from, mu_to, nf, steps=4000):
    """3-loop MS-bar: mu^2 d alpha/d mu^2 = -(b0 a^2 + b1 a^3 + b2 a^4), RK4
    in ln mu^2."""
    b0, b1, b2 = _beta(nf)
    f = lambda a: -(b0 * a * a + b1 * a ** 3 + b2 * a ** 4)   # noqa: E731
    t0, t1 = 2.0 * math.log(mu_from), 2.0 * math.log(mu_to)
    h = (t1 - t0) / steps
    a = alpha
    for _ in range(steps):
        k1 = f(a)
        k2 = f(a + 0.5 * h * k1)
        k3 = f(a + 0.5 * h * k2)
        k4 = f(a + h * k3)
        a += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return a


def alpha_s(q, ref, alpha_mz=None):
    """alpha_s(Q) [Q in GeV], nf = 5 above m_b, 4 down to m_c, 3 below;
    continuous at the thresholds."""
    a = ref["alpha_s_mz"] if alpha_mz is None else alpha_mz
    mz, mb, mc = ref["m_z_gev"], ref["m_b_gev"], ref["m_c_gev"]
    if q >= mb:
        return _run_alpha(a, mz, q, 5)
    a = _run_alpha(a, mz, mb, 5)
    if q >= mc:
        return _run_alpha(a, mb, q, 4)
    a = _run_alpha(a, mb, mc, 4)
    return _run_alpha(a, mc, q, 3)


def nf_at(q, ref):
    return 5 if q >= ref["m_b_gev"] else (4 if q >= ref["m_c_gev"] else 3)


def c_bjorken_coeffs(nf):
    c2 = 55.0 / 12.0 - nf / 3.0
    c3 = (13841.0 / 216.0 + 44.0 / 9.0 * ZETA3 - 55.0 / 2.0 * ZETA5
          + nf * (-10339.0 / 1296.0 - 61.0 / 54.0 * ZETA3 + 5.0 / 3.0 * ZETA5)
          + nf * nf * 115.0 / 648.0)
    return 1.0, c2, c3


def bjorken_reference(q2, ref=None):
    """Gamma_1^{p-n}(Q2) at leading twist with its uncertainty."""
    ref = load_reference() if ref is None else ref
    q = math.sqrt(q2)
    nf = nf_at(q, ref)
    c1, c2, c3 = c_bjorken_coeffs(nf)

    def gam(ga, amz):
        a = alpha_s(q, ref, amz) / math.pi
        return ga / 6.0 * (1.0 - c1 * a - c2 * a * a - c3 * a ** 3), a

    val, a = gam(ref["g_A"], ref["alpha_s_mz"])
    d_ga = ref["g_A_err"] / 6.0 * (val / (ref["g_A"] / 6.0))
    up, _ = gam(ref["g_A"], ref["alpha_s_mz"] + ref["alpha_s_mz_err"])
    dn, _ = gam(ref["g_A"], ref["alpha_s_mz"] - ref["alpha_s_mz_err"])
    d_as = 0.5 * abs(up - dn)
    d_tr = ref["g_A"] / 6.0 * abs(c3 * a ** 3)
    return dict(q2=q2, nf=nf, alpha_s=a * math.pi, c_bj=val / (ref["g_A"] / 6.0),
                gamma=val, err=math.sqrt(d_ga ** 2 + d_as ** 2 + d_tr ** 2),
                err_g_a=d_ga, err_alpha_s=d_as, err_truncation=d_tr)


def bjorken_scan(g1p, g1n, x_mins):
    """Truncated first moments over the scan: Gamma^{p-n}, Gamma^p, Gamma^n."""
    rows = []
    gp = gn = ep = en = 0.0
    for hi, lo in _decades(x_mins):
        a, b = math.log(lo), math.log(hi)
        v, e = _quad(lambda t: g1p(math.exp(t)) * math.exp(t), a, b)
        gp += v
        ep += e
        v, e = _quad(lambda t: g1n(math.exp(t)) * math.exp(t), a, b)
        gn += v
        en += e
        rows.append(dict(x_min=lo, gamma_pn=gp - gn, gamma_p=gp, gamma_n=gn,
                         quad_err=ep + en))
    return rows


def bjorken_verdict(rows, ref_row):
    """The declared rule: converged over the last decade AND within 10 %."""
    last, prev = rows[-1], rows[-2]
    g_ref = ref_row["gamma"]
    step = abs(last["gamma_pn"] - prev["gamma_pn"])
    dist = last["gamma_pn"] - g_ref
    converged = step <= BJ_CONV_REL * g_ref
    within = abs(dist) <= BJ_TOL_REL * g_ref
    return dict(last_decade_step=step, distance=dist,
                distance_rel=dist / g_ref, converged=converged, within=within,
                ok=converged and within)


def measure_bjorken(polsf, x_mins, q2s=Q2_POINTS, ref=None):
    ref = load_reference() if ref is None else ref
    per_q2 = {}
    for q2 in q2s:
        rows = bjorken_scan(lambda x, q=q2: polsf.g1p(x, q),
                            lambda x, q=q2: polsf.g1n(x, q), x_mins)
        r = bjorken_reference(q2, ref)
        per_q2[q2] = dict(rows=rows, reference=r, verdict=bjorken_verdict(rows, r))
    status = "pass" if all(v["verdict"]["ok"] for v in per_q2.values()) else "fail"
    return dict(per_q2=per_q2, status=status)


def lhapdf_g1():
    """(LhapdfG1 object, None) or (None, blocked reason).  Any error other
    than a missing tier / a missing set propagates (the harness is broken)."""
    import lipolgen as lg
    if not getattr(lg, "HAVE_LHAPDF", False) or not hasattr(lg, "LhapdfG1"):
        return None, ("environment: lipolgen built without the LHAPDF tier "
                      "(HAVE_LHAPDF false); needs LHAPDF set %s" % LHAPDF_G1_SET)
    try:
        return lg.LhapdfG1(LHAPDF_G1_SET, 0), None
    except RuntimeError as exc:
        if "Info file not found for PDF set" not in str(exc):
            raise
        return None, "environment: LHAPDF set %s not installed" % LHAPDF_G1_SET


# ---------------------------------------------------------------------- FSI

def fsi_record():
    """The tree side of the would-be unitarity sub-row, recorded only."""
    import lipolgen as lg
    f = lg.GlauberFsiWeight(lg.li6_alpha_channel())
    o = f.options
    sxn = o.sigma_xn_mb
    return dict(
        channel="li6_alpha_channel() defaults (Hulthen beta 0.3)",
        sigma_xn_mb=float(sxn), survival=float(f.survival()),
        int_w_minus_1=float(f.survival()) - 1.0,
        sigma_tot_cluster_mb=float(f.sigma_cluster_mb(sxn)),
        sigma_el_cluster_mb=float(f.sigma_cluster_el_mb(sxn)),
        elastic_gain=float(o.elastic_gain), k_max_gev=float(o.k_max),
        w_max=float(o.w_max), clipped_grid_fraction=float(f.clipped_grid_fraction()),
    )


# ---------------------------------------------------------------------- run

def run(verbose=True):
    t_start = time.perf_counter()
    ref = load_reference()

    # (a) Burkhardt-Cottingham on g2_ww
    bc = measure_bc()
    cells = [(k, r) for k, v in bc.items() for r in v["rows"]]
    bc_ok = all(r["ok"] for _, r in cells)
    worst_key, worst = max(cells, key=lambda kr: kr[1]["ratio"])
    limits = {}
    for (tgt, q2), v in bc.items():
        limits.setdefault(tgt, {})[q2] = v["exponent"]
    bc_row = dict(
        status="pass" if bc_ok else "fail", n_cells=len(cells),
        n_within=sum(1 for _, r in cells if r["ok"]),
        max_ratio=worst["ratio"], worst_cell=dict(target=worst_key[0],
                                                  q2=worst_key[1],
                                                  x_min=worst["x_min"]),
        limit_exponent=limits, cells=bc)

    # (b) Bjorken on ToyG1, and LhapdfG1 if present
    import lipolgen as lg
    bj_toy = measure_bjorken(lg.ToyG1(), X_MIN_SCAN, ref=ref)
    lh, why = lhapdf_g1()
    if lh is None:
        bj_lh = dict(status="blocked", reason=why)
    else:                                            # never run in this tree's env
        bj_lh = measure_bjorken(lh, X_MIN_SCAN_LHAPDF, ref=ref)

    # (c) FSI unitarity
    fsi = dict(status="blocked",
               reason=FSI_BLOCKED_REASON + " (no GlauberFsiWeight option "
                      "realises int dGamma S[FSI + FSI^2] = 0: the default is "
                      "absorptive by design, elastic_gain would be a tune, "
                      "weight_normalised is a division, the table stops at "
                      "k_max)",
               tree_side=fsi_record())

    subrows = {"bc_ww": bc_row, "bjorken[toyg1]": bj_toy,
               "bjorken[nnpdfpol11]": bj_lh, "fsi_unitarity": fsi}
    evaluated = [v["status"] for v in subrows.values() if v["status"] != "blocked"]
    if not evaluated:
        status = "blocked"
    else:
        status = "pass" if all(s == "pass" for s in evaluated) else "fail"

    q5 = bj_toy["per_q2"][5.0]
    g_last = q5["rows"][-1]
    generator = ("bc_ww closure max|resid|/tol = %.3f over %d cells (%s); "
                 "Bjorken ToyG1 Q2=5: int_{1e-6}^1 (g1p-g1n) = %.4f, last decade "
                 "+%.4f (divergent: g1n ~ x^-%.3f)"
                 % (worst["ratio"], len(cells), bc_row["status"].upper(),
                    g_last["gamma_pn"], q5["verdict"]["last_decade_step"],
                    1.0 + bc[("n", 5.0)]["toy_lambda"]))
    reference = ("BC: -x_min int g1/u du (exact truncation of int g2^WW); "
                 "Bjorken: g_A/6 C_Bj(alpha_s) = %.4f +- %.4f at Q2=5"
                 % (q5["reference"]["gamma"], q5["reference"]["err"]))
    tolerance = ("BC: Peano trapezoid bound (h^2/8)L^2 int|D^2 g1| + quad errors; "
                 "Bjorken: converged (last decade <= 1 %) and within 10 %")
    report = dict(name=NAME, generator=generator, reference=reference,
                  tolerance=tolerance, status=status, subrows=subrows,
                  reference_constants=ref, runtime_s=None)

    if verbose:
        print("(a) BC closure  int_{x_min}^1 g2_ww dx  vs  T = -x_min int_{x_min}^1 g1/u du"
              "   (ToyG1, npts = %d)" % G2_NPTS)
        print("  tgt  Q2    x_min    int g2_tree      T(x_min)        resid"
              "        tol      resid/EM   |r|/tol")
        for (tgt, q2), v in bc.items():
            for r in v["rows"]:
                em = r["em_prediction"]
                print("  %s  %4.1f  %7.0e  %+.6e  %+.6e  %+.3e  %.3e  %7.4f  %6.3f"
                      % (tgt, q2, r["x_min"], r["integral"], r["truncation"],
                         r["residual"], r["tol"],
                         (r["residual"] / em) if em else float("nan"),
                         r["ratio"]))
        print("  x_min -> 0: effective exponent of |T| (> 0: T grows, no BC limit):")
        for tgt in ("p", "n", "d"):
            print("    %s: %s" % (tgt, ", ".join(
                "Q2=%g: %+.4f" % (q2, float(e)) for q2, e in sorted(limits[tgt].items()))))
        print("  toy lambda(Q2) = %s  (g1n ~ x^(-1-lambda))" % ", ".join(
            "%.4f" % bc[("n", q2)]["toy_lambda"] for q2 in Q2_POINTS))
        print("(b) Bjorken  int_{x_min}^1 (g1p - g1n) dx, ToyG1:")
        for q2, v in bj_toy["per_q2"].items():
            r = v["reference"]
            print("  Q2 = %4.1f: Gamma_ref = %.5f +- %.5f  (alpha_s = %.4f, nf = %d, "
                  "C_Bj = %.4f)" % (q2, r["gamma"], r["err"], r["alpha_s"], r["nf"],
                                    r["c_bj"]))
            for row in v["rows"]:
                print("     x_min %7.0e  Gamma^{p-n} = %.5f  Gamma^p = %.5f  "
                      "Gamma^n = %+.5f" % (row["x_min"], row["gamma_pn"],
                                          row["gamma_p"], row["gamma_n"]))
            vd = v["verdict"]
            print("     last decade +%.5f (%.1f %% of ref), distance %+.5f (%+.1f %%)"
                  " -> %s" % (vd["last_decade_step"],
                              100 * vd["last_decade_step"] / r["gamma"],
                              vd["distance"], 100 * vd["distance_rel"],
                              "PASS" if vd["ok"] else "FAIL"))
        fr = fsi["tree_side"]
        print("(c) FSI, tree side recorded (NOT a unitarity sum): %s, sigma_XN = %g mb: "
              "survival = int w dGamma/int dGamma = %.6f; sigma_tot(X-alpha) = %.2f mb, "
              "sigma_el = %.2f mb; k_max = %g GeV, w_max = %g"
              % (fr["channel"], fr["sigma_xn_mb"], fr["survival"],
                 fr["sigma_tot_cluster_mb"], fr["sigma_el_cluster_mb"],
                 fr["k_max_gev"], fr["w_max"]))
        print("REPORT | %s[bc_ww] | max |resid|/tol = %.3f, %d/%d cells within "
              "(worst %s Q2=%g x_min=%g) | exact truncation -x_min int g1/u du "
              "| Peano trapezoid bound + quad errors | %s"
              % (NAME, worst["ratio"], bc_row["n_within"], len(cells),
                 worst_key[0], worst_key[1], worst["x_min"],
                 bc_row["status"].upper()))
        print("REPORT | %s[bjorken[toyg1]] | %s | g_A/6 C_Bj(alpha_s), leading "
              "twist | converged (last decade <= 1 %%) and |dist| <= 10 %% | %s"
              % (NAME, "; ".join(
                  "Q2=%g: %.4f (ref %.4f, %+.0f %%, last decade +%.4f)"
                  % (q2, v["rows"][-1]["gamma_pn"], v["reference"]["gamma"],
                     100 * v["verdict"]["distance_rel"],
                     v["verdict"]["last_decade_step"])
                  for q2, v in bj_toy["per_q2"].items()),
                 bj_toy["status"].upper()))
        if bj_lh["status"] == "blocked":
            print("REPORT | %s[bjorken[nnpdfpol11]] | BLOCKED: %s"
                  % (NAME, bj_lh["reason"]))
        else:
            print("REPORT | %s[bjorken[nnpdfpol11]] | %s | g_A/6 C_Bj(alpha_s) | "
                  "converged and within 10 %% | %s" % (NAME, "; ".join(
                      "Q2=%g: %.4f (%+.0f %%)" % (q2, v["rows"][-1]["gamma_pn"],
                                                  100 * v["verdict"]["distance_rel"])
                      for q2, v in bj_lh["per_q2"].items()),
                      bj_lh["status"].upper()))
        print("REPORT | %s[fsi_unitarity] | BLOCKED: %s" % (NAME, fsi["reason"]))
    report["runtime_s"] = time.perf_counter() - t_start
    if verbose:
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, generator, reference, tolerance, status.upper()))
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
