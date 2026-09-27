#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Scout item 5 / 02_data_nucleon_deuteron.md sec. 5.2 -- the deuteron's two
static observables, eta_d and Q_d, against every deuteron wave function the
tree ships (tier T2).

WHAT IS COMPARED.  Two measured numbers, vendored in
`data/deuteron_static_measured.json` as a SECOND-HAND transcription (from
02_data_nucleon_deuteron.md sec. 5.2 and include/lipolgen/b1_nuclear.hpp's
"100 % band" block; neither paper was re-read):

  eta_d = A_D/A_S = 0.0256(4)       Rodning & Knutson, PRC 41 (1990) 898
  Q_d   = +0.2859(3) fm^2           Bishop & Cheung, PRA 20 (1979) 381

against the three deuteron wave functions the tree ships, each at its
SHIPPED settings, nothing moved:

  cdbonn   `lg.cdbonn_wave()`, Machleidt's Appendix-D closed form (C_j, D_j,
           m_j of Table XX, C_11 and D_9..D_11 solved for).  Reached by
           `DeuteronConvolutionB1` with `wave = kCdBonn` (not a default).
           eta = the tree's own `eta()` = D_1/C_1; Q_d from the closed-form
           u(r), w(r) of Eqs. (D19)/(D20), integrated here.
  av18     `data/vmc/deuteron/fdeut.av18` (`lg.VMC_DEUTERON_WAVE`), the file
           the tree reads for `DeuteronConvolutionB1`'s default and for
           `deuteron_channel(source=VmcAV18)`.  The tree reads its k-space
           block; eta needs the asymptotic r-space ratio, which a k table on
           a 0.1 fm^-1 grid resolves only to a few 1e-3 at r ~ 10-20 fm and
           not at all beyond (kappa = 0.23 fm^-1 is ~2 grid steps; at 40 fm
           the spline's transform has the wrong sign), so eta and Q_d are
           computed here from the SAME file's r-space block (u, w to
           r = 100 fm), and the k-space path the generator uses is carried
           as a cross-check on both (Q_d via the tree's
           `alpha_d_quadrupole_fm2`, eta via the tree's spline transformed
           to r = 20 fm: -0.24 %).
  hulthen  `deuteron_channel()` -- THE DEFAULT deuteron control of the
           tagged channel: the analytic pair psi_0 = 1/(k^2+kappa^2) -
           1/(k^2+beta^2), psi_2 = k^2/((k^2+kappa^2)(k^2+beta^2)^2) at
           beta = BETA_DEFAULT = 0.30 GeV, P_D = P_D_DEUTERON = 0.045, kappa
           from DEUTERON_P_TAG's separation energy, normalised exactly as
           `TaggedModel` normalises them (the constants N_L are READ from
           `TaggedModel(...).radial_table(L)`, not recomputed).

HOW eta AND Q_d ARE COMPUTED (impulse approximation, no meson-exchange
currents, no relativistic corrections -- the textbook one-body operator):

  Q_d = [ (sqrt2/10) INT u w r^2 dr - (1/20) INT w^2 r^2 dr ] / INT (u^2+w^2) dr
  eta = A_D/A_S,  u -> A_S e^{-kappa r},  w -> A_D e^{-kappa r} h(kappa r),
                  h(x) = 1 + 3/x + 3/x^2.

  u, w are in the `fdeut.av18` sign convention (the bare-j_L Fourier pair,
  w = -psi_2^a for CD-Bonn; cluster.hpp `CdBonnWave`), in which Q_d > 0.
  For the Hulthen pair, partial fractions give the r-space forms in closed
  form: u = N_0 sqrt(pi/2) (e^{-kappa r} - e^{-beta r}) and
  w = N_2 [A (T_kappa - T_beta) + C T2_beta], with A = -kappa^2/(beta^2 -
  kappa^2)^2, C = beta^2/(beta^2 - kappa^2), T_m(r) = -sqrt(pi/2) [e^{-mr}
  h(mr) - 3/(mr)^2] the bare-j_2 transform of 1/(k^2+m^2) and T2_m =
  -dT_m/d(m^2); hence eta = -N_2 A/N_0 exactly (the pole residues at
  k^2 = -kappa^2).  The forms take ONE beta for both waves; `hulthen()`
  refuses (RuntimeError) a channel whose D wave carries a beta different
  from its S wave's (both 0.30 GeV today).
  `python/tests/test_bench_static.py` checks w(r) against a direct
  numerical transform and eta against the large-r ratio.
  Quadrature: composite Simpson on a uniform r grid; its error is MEASURED
  by halving the grid and carried in the row.

TOLERANCE, DECLARED BEFORE THE VERDICT.
  * eta: ASSERTED.  |eta_model - 0.0256| <= 2 sigma (sigma = 0.0004) for
    EVERY shipped wave function; the row passes only if all three do.  The
    2 sigma was fixed by the task that commissioned this harness before the
    AV18 and Hulthen values were looked at (AV18's is printed in
    fdeut.av18's own header, 0.025045, and was read only afterwards;
    Hulthen's had never been computed).  CD-Bonn's 0.07 sigma WAS known (02
    sec. 5.2 prints it, as does the scout's preview), so for CD-Bonn the
    tolerance was chosen after its number was seen -- stated here, as
    t3_nmc_li6_over_d.py states its own.  Nothing was re-chosen after the
    Hulthen number (+17 sigma) came out.
  * Q_d: RECORDED, NOT ASSERTED.  The measured Q_d contains two-body
    currents and relativistic corrections an impulse-approximation Q_d
    omits; realistic wave functions are expected to fall short by about
    5 % (the scout's preview, before this harness: CD-Bonn 0.2705, -5.4 %).
    The row carries each model's distance (absolute and relative); the
    status reflects ONLY the asserted eta part.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: no b1 on either side, and no
     A-scaling assumption (no 0.5, no 2/6): A = 2 on both sides.
  2. isoscalar denominator -- NOT APPLICABLE: no structure function enters.
  3. per nucleus / per nucleon -- both sides are PER DEUTERON: eta and Q_d
     are properties of the bound n-p state (Q_d in e fm^2 with e = 1; the
     1/4 from the proton sitting at r/2 is inside the operator's sqrt2/10
     and 1/20).  eta is a ratio and independent of the wave function's
     normalisation; Q_d is quadratic in it, so every model is divided by
     its own INT (u^2 + w^2) dr (the Hulthen tables integrate to 1.000144
     over all r, because TaggedModel normalises on its own finite k grid;
     that is divided out, not tuned away).

WHAT THIS DOES NOT VALIDATE.  P_D (not an observable; CD-Bonn 4.86 %, AV18
5.76 %, the scenario 4.5 % differ by more than any error bar and nothing
here constrains them); w(k) at the k >~ 300 MeV/c that dominates the tagged
channel's high-k tail (eta is an asymptotic, r -> infinity quantity); the
light-cone convolution or b1; the 6Li alpha-d D wave or its eta (row C-1 of
00_in_tree_checks.md, a consistency band on the quadrupole dial); the
deuteron elastic form factors of the RC tail (their F_q(0) is NORMALISED on
Q_d = 0.2859, an input, so comparing it would be circular).  The S-wave-only
Hulthen of `spectator.hpp`'s `MomentumSampler` (the P_D = 0 limit of the
same pair) is not a separate row: eta = Q_d = 0 there by construction.  The
survey's sentence (02 sec. 5.2) that eta validates "the S-D mixture of the
deuteron wave function the tagged channel ... rides on" holds for CD-Bonn
and AV18 only: the tagged channel's DEFAULT deuteron is the Hulthen pair.

Run:  source env.sh && python3 validation/benchmarks/t2_deuteron_static.py

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

try:
    import numpy as np
except ImportError:  # a missing dependency is a broken harness: exit 2
    if __name__ == "__main__":
        traceback.print_exc()
        sys.exit(2)
    raise

# numpy < 2.0 spells it np.trapz; >= 2.0 spells it np.trapezoid (and 2.4
# removed np.trapz), as in validation/vmc_reconcile.py.  pyproject.toml
# declares numpy >= 1.20.
_trapz = getattr(np, "trapezoid", None) or np.trapz

HERE =os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "deuteron_static_measured.json")
#: sha256 of every byte of the vendored file, and of its canonical body
#: (json.dumps(body, sort_keys=True, separators=(",", ":"))), which the file
#: also records.  Both are re-checked on EVERY read.
FILE_SHA256 = "6433c31e3c837bd2a88c71417eca5fa4eba51f17495486dd2c21bc99d8f692d5"
BODY_SHA256 = "18f805e529544681f4b793771f0b2011ba540f1333c30a77f8c0979ad6ede849"

NAME = "t2_deuteron_static"
ETA_MAX_PULL = 2.0          # tolerance: |eta - eta_d| <= 2 sigma, every model
MODELS = ("cdbonn", "av18", "hulthen")
R_MAX_FM = 80.0             # closed-form quadrature: e^{-2 kappa 80} ~ 1e-16
N_INTERVALS = 100000        # Simpson intervals on (0, R_MAX_FM]
ETA_WINDOW_FM = (20.0, 40.0)  # AV18: the r range the asymptotic ratio is read
HULTHEN_BETA_BAND = (0.20, 0.40)  # cluster.hpp's scanned band -- RECORDED only
SQ = math.sqrt(math.pi / 2.0)


# ---------------------------------------------------------------- reference

def read_vendored(path=DATA):
    """The vendored document -- RuntimeError unless the file's bytes AND its
    canonical body hash to the recorded digests."""
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


def load_reference(path=DATA):
    """{quantity: point} for the two vendored points, in file order."""
    pts = read_vendored(path)["body"]["points"]
    names = [p["quantity"] for p in pts]
    if names != ["eta_d", "Q_d"]:
        raise RuntimeError("%s: expected points eta_d, Q_d; found %s"
                           % (path, names))
    return {p["quantity"]: dict(p, value=float(p["value"]),
                                sigma=float(p["sigma"])) for p in pts}


# ---------------------------------------------------------------- numerics

def h_asym(x):
    """The D-wave asymptotic factor 1 + 3/x + 3/x^2."""
    return 1.0 + 3.0 / x + 3.0 / (x * x)


def simpson(y, h):
    """Composite Simpson on a uniform grid with an even number of intervals."""
    n = len(y) - 1
    if n < 2 or n % 2:
        raise ValueError("simpson needs an even number of intervals, got %d" % n)
    return float(h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum()
                            + 2.0 * y[2:-1:2].sum()))


def quadrupole_ia(r, u, w, h):
    """Impulse-approximation Q_d [fm^2 when r is in fm], P_D and the norm of
    a (u, w) pair on a uniform grid, normalised by INT (u^2 + w^2) dr."""
    norm = simpson(u * u + w * w, h)
    q = ((math.sqrt(2.0) / 10.0) * simpson(u * w * r * r, h)
         - simpson(w * w * r * r, h) / 20.0) / norm
    return q, simpson(w * w, h) / norm, norm


def quadrupole_with_error(r, u, w, h):
    """`quadrupole_ia` on the grid and on every second point of it; the
    difference is the measured quadrature error (Simpson on h vs 2h)."""
    q, pd, norm = quadrupole_ia(r, u, w, h)
    if (len(r) - 1) % 4 == 0:
        q2, _, _ = quadrupole_ia(r[::2], u[::2], w[::2], 2.0 * h)
        err = abs(q - q2)
    else:                                            # pragma: no cover
        err = float("nan")
    return q, pd, norm, err


def r_grid(r_max=R_MAX_FM, n=N_INTERVALS):
    """[h, r_max + h] in n equal steps, r = 0 dropped: u ~ r and w ~ r^3
    there, so every integrand vanishes at least like r^2 and [0, h]
    contributes O(h^3) ~ 1e-10 of the norm at h = 8e-4 fm."""
    h = r_max / n
    r = h * np.arange(1, n + 2)            # n + 1 points: h .. r_max + h
    return r, h


# ---------------------------------------------------------------- the models

def cdbonn():
    import lipolgen as lg  # noqa: WPS433 (env.sh puts it on the path)
    wf = lg.cdbonn_wave()
    m = np.asarray(wf.m, dtype=float)
    c = np.asarray(wf.c, dtype=float)
    d = np.asarray(wf.d, dtype=float)

    def uw(r):
        u = np.zeros_like(r)
        w = np.zeros_like(r)
        for mj, cj, dj in zip(m, c, d):    # (D19), (D20): r in fm
            e = np.exp(-mj * r)
            u += cj * e
            w += dj * e * h_asym(mj * r)
        return u, w

    r, h = r_grid()
    u, w = uw(r)
    q, pd_r, norm, qerr = quadrupole_with_error(r, u, w, h)
    r_far = np.array([60.0])
    u_far, w_far = uw(r_far)
    eta_far = float(w_far[0] / (u_far[0] * h_asym(m[0] * r_far[0])))
    # the tree's own quadrupole routine -- the number tests/test_b1_nuclear.cpp
    # item 5 prints (0.270178) -- at its default numerics, and converged
    t = lg.cdbonn_fdeut_table()
    q_tree = lg.alpha_d_quadrupole_fm2(
        lg.ClusterPartialWave.from_uw(t["k_gev"], t["u"], 0),
        lg.ClusterPartialWave.from_uw(t["k_gev"], t["w"], 2))
    tf = lg.cdbonn_fdeut_table(20.0, 0.02)
    q_tree_fine = lg.alpha_d_quadrupole_fm2(
        lg.ClusterPartialWave.from_uw(tf["k_gev"], tf["u"], 0),
        lg.ClusterPartialWave.from_uw(tf["k_gev"], tf["w"], 2),
        r_max_fm=60.0, n_r=6000)
    return dict(
        label="CD-Bonn closed form (lg.cdbonn_wave(), Machleidt PRC 63 "
              "(2001) 024001 App. D)",
        eta=float(wf.eta()), eta_source="tree: CdBonnWave.eta() = D_1/C_1",
        eta_numerical_error=0.0,
        eta_check_r60=eta_far,
        q_d=q, q_d_quad_error=qerr, p_d=pd_r, norm=norm,
        p_d_tree=float(wf.norm_d()), a_s=float(wf.a_s()),
        kappa_fm=float(m[0]),
        cross_checks=dict(
            q_tree_default=float(q_tree),
            q_tree_default_note="lg.alpha_d_quadrupole_fm2 on "
                                "cdbonn_fdeut_table(), r_max 20 fm, dk 0.1 "
                                "fm^-1 (tests/test_b1_nuclear.cpp item 5)",
            q_tree_converged=float(q_tree_fine),
            q_tree_converged_note="the same routine on dk = 0.02 fm^-1, "
                                  "r_max = 60 fm, n_r = 6000"),
        config=dict(r_max_fm=R_MAX_FM, n_intervals=N_INTERVALS),
    )


def read_fdeut_rspace(path):
    """Header numbers and the r-space block (r, u, w) of an ANL `fdeut.*`."""
    with open(path) as f:
        lines = f.read().splitlines()
    header = None
    start = None
    for i, line in enumerate(lines):
        tok = line.split()
        if tok[:2] == ["ebind", "dstate"]:
            header = dict(zip(tok, (float(v) for v in lines[i + 1].split())))
        if tok[:5] == ["r", "u", "du/dr", "w", "dw/dr"]:
            start = i + 1
            break
    if header is None or start is None:
        raise RuntimeError("%s: header or r-space block not found" % path)
    rows = []
    for line in lines[start:]:
        tok = line.split()
        if not tok:
            if rows:
                break
            continue
        if len(tok) != 5:
            break
        rows.append([float(v) for v in tok])
    a = np.asarray(rows)
    return header, a[:, 0], a[:, 1], a[:, 3]


def av18():
    import lipolgen as lg  # noqa: WPS433
    path = lg.data_path(lg.VMC_DEUTERON_WAVE)
    header, r, u, w = read_fdeut_rspace(path)
    n = len(r)
    h = r[1] - r[0]
    if n != 10000 or abs(r[0] - h) > 1e-9 or np.max(np.abs(np.diff(r) - h)) > 1e-9:
        raise RuntimeError("%s: r-space block is not the uniform 0.01 .. 100 fm "
                           "grid this harness reads (%d rows)" % (path, n))
    # r = 0 (u = w = 0) prepended: 10000 intervals, Simpson on the file grid
    r0 = np.concatenate([[0.0], r])
    u0 = np.concatenate([[0.0], u])
    w0 = np.concatenate([[0.0], w])
    q, pd_r, norm, qerr = quadrupole_with_error(r0, u0, w0, h)
    lo, hi = ETA_WINDOW_FM
    i_lo, i_hi = int(np.searchsorted(r, lo)), int(np.searchsorted(r, hi))
    kappa = -math.log(u[i_hi] / u[i_lo]) / (r[i_hi] - r[i_lo])
    sel = slice(i_lo, i_hi + 1)
    eta_r = w[sel] / (u[sel] * h_asym(kappa * r[sel]))
    eta = float(np.mean(eta_r))
    a_s = float(np.mean(u[sel] * np.exp(kappa * r[sel]))) / math.sqrt(norm)
    ch = lg.deuteron_channel(source=lg.ClusterWaveSource.VmcAV18)
    v0, v2 = ch.waves[0].vmc, ch.waves[1].vmc
    p0 = lg.ClusterPartialWave.from_uw(v0.k, v0.psi, 0)
    p2 = lg.ClusterPartialWave.from_uw(v2.k, v2.psi, 2)     # phi_2 = -w(k)
    q_tree = lg.alpha_d_quadrupole_fm2(p0, p2)
    # eta from the k block the generator reads, through the tree's own spline
    # and the bare-j_L transform at r = ETA_WINDOW_FM[0]: a cross-check only
    hc = lg.HBARC_GEV_FM
    kf = np.linspace(0.0, float(v0.k[-1]) / hc, 20001)          # fm^-1
    s0 = np.array([p0(float(x * hc)) for x in kf])
    s2 = np.array([p2(float(x * hc)) for x in kf])
    rk = ETA_WINDOW_FM[0]
    x = kf[1:] * rk
    j0 = np.concatenate([[0.0], np.sin(x) / x])
    j2 = np.concatenate([[0.0], (3.0 / x ** 3 - 1.0 / x) * np.sin(x)
                         - 3.0 * np.cos(x) / x ** 2])
    dk = kf[1] - kf[0]
    u_k = float(_trapz(kf * kf * j0 * s0, dx=dk))
    w_k = -float(_trapz(kf * kf * j2 * s2, dx=dk))
    eta_kblock = w_k / (u_k * h_asym(kappa * rk))
    return dict(
        label="AV18 table (%s, r-space block; the tree reads its k block)"
              % lg.VMC_DEUTERON_WAVE,
        eta=eta, eta_source="this harness: w/(u h(kappa r)) averaged over "
                            "r in [%g, %g] fm of the file's r block" % ETA_WINDOW_FM,
        eta_numerical_error=float(np.max(eta_r) - np.min(eta_r)),
        kappa_fm=kappa, a_s=a_s,
        q_d=q, q_d_quad_error=qerr, p_d=pd_r, norm=norm,
        p_d_tree=float(lg.deuteron_av18_p_d()),
        header=header,
        cross_checks=dict(
            header_eta=header["eta"], header_qm=header["qm"],
            header_dstate=header["dstate"], header_as=header["as"],
            q_tree_kblock=float(q_tree),
            q_tree_kblock_note="lg.alpha_d_quadrupole_fm2 on the k block "
                               "the generator reads (deuteron_channel("
                               "source=VmcAV18)), default numerics "
                               "(tests/test_b1_nuclear.cpp item 5: 0.269362)",
            eta_kblock_r=eta_kblock,
            eta_kblock_note="w/(u h(kappa r)) at r = %g fm from the k block "
                            "through ClusterPartialWave's spline and a "
                            "20001-point bare-j_L transform" % rk),
        config=dict(file=path, eta_window_fm=ETA_WINDOW_FM),
    )


def hulthen_forms(kappa, beta):
    """The Hulthen pair's r-space closed forms (r in GeV^-1), per unit N_L,
    and the D-wave residue coefficient A.  See the module docstring."""
    b2, k2 = beta * beta, kappa * kappa
    a_res = -k2 / (b2 - k2) ** 2
    c_dbl = b2 / (b2 - k2)

    def t1(m, r):
        x = m * r
        return -SQ * (np.exp(-x) * h_asym(x) - 3.0 / (x * x))

    def t2(m, r):           # -d t1 / d(m^2) = sqrt(pi/2) (1/2m) df/dm
        e = np.exp(-m * r)
        dfdm = (-r * e * h_asym(m * r)
                + e * (-3.0 / (m * m * r) - 6.0 / (m ** 3 * r * r))
                + 6.0 / (m ** 3 * r * r))
        return SQ * dfdm / (2.0 * m)

    def u(r):
        return SQ * (np.exp(-kappa * r) - np.exp(-beta * r))

    def w(r):
        return a_res * (t1(kappa, r) - t1(beta, r)) + c_dbl * t2(beta, r)

    return u, w, a_res


def hulthen(beta=None):
    import lipolgen as lg  # noqa: WPS433
    ch = lg.deuteron_channel() if beta is None else lg.deuteron_channel(beta=beta)
    kappa = float(ch.base.kappa())
    s_wave, d_wave = ch.waves[0], ch.waves[1]
    if (s_wave.l, d_wave.l) != (0, 2) or s_wave.vmc is not None:
        raise RuntimeError("deuteron_channel() is no longer the analytic "
                           "Hulthen S + D pair")
    # the closed forms (hulthen_forms) and the D-wave residue eta = -N_2 A/N_0
    # take ONE beta for both waves: refuse a pair whose D wave carries its own
    if float(d_wave.beta) != float(s_wave.beta):
        raise RuntimeError("deuteron_channel()'s D wave has beta = %r GeV, its S "
                           "wave %r: the Hulthen closed forms here assume one "
                           "beta for both" % (float(d_wave.beta),
                                              float(s_wave.beta)))
    beta_used = float(s_wave.beta)
    tm = lg.TaggedModel(ch)
    k = np.asarray(tm.k, dtype=float)
    norms = []
    for wave in (s_wave, d_wave):
        tab = np.asarray(tm.radial_table(wave.l), dtype=float)
        ana = np.array([wave.radial(float(x), kappa) for x in k])
        ratio = tab / ana
        spread = float((ratio.max() - ratio.min()) / abs(ratio.mean()))
        if spread > 1e-12:
            raise RuntimeError("TaggedModel's L = %d table is not a constant "
                               "times Wave.radial (spread %.3g)" % (wave.l, spread))
        norms.append(float(ratio.mean()))
    n0, n2 = norms
    u_f, w_f, a_res = hulthen_forms(kappa, beta_used)
    eta = -n2 * a_res / n0
    hc = lg.HBARC_GEV_FM
    r_fm, h = r_grid()
    rg = r_fm / hc                                     # GeV^-1
    # u, w [GeV^1/2] -> [fm^-1/2], so INT (u^2 + w^2) dr is the same number
    # in either unit system (Q_d, a ratio, does not care)
    to_fm = 1.0 / math.sqrt(hc)
    q, pd_r, norm, qerr = quadrupole_with_error(r_fm, to_fm * n0 * u_f(rg),
                                                to_fm * n2 * w_f(rg), h)
    r_far = np.array([80.0 / hc])
    eta_far = float(n2 * w_f(r_far)[0] / (n0 * u_f(r_far)[0]
                                          * h_asym(kappa * r_far[0])))
    # the continuum normalisation (INT_0^inf), for comparison only
    p_d_nom = float(d_wave.prob)
    n0_an = math.sqrt(1.0 - p_d_nom) / math.sqrt(
        (math.pi / 4.0) * (beta_used - kappa) ** 2
        / (kappa * beta_used * (kappa + beta_used)))
    ratio_an = n0_an / n0
    out = dict(
        label="Hulthen pair (deuteron_channel() default: beta = %.2f GeV, "
              "P_D = %.3f, kappa = %.6f GeV), TaggedModel normalisation"
              % (beta_used, p_d_nom, kappa),
        eta=eta, eta_source="this harness: -N_2 A/N_0, the k^2 = -kappa^2 "
                            "pole residues, N_L read from TaggedModel",
        eta_numerical_error=abs(eta - eta_far),
        eta_check_r80=eta_far,
        kappa_fm=kappa / hc, beta_gev=beta_used, p_d_nominal=p_d_nom,
        n0=n0, n2=n2, n0_continuum_over_tree=ratio_an,
        q_d=q, q_d_quad_error=qerr, p_d=pd_r, norm=norm,
        config=dict(beta_gev=beta_used, p_d=p_d_nom, k_grid_gev=(
            float(k[0]), float(k[-1]), len(k)), r_max_fm=R_MAX_FM,
            n_intervals=N_INTERVALS),
    )
    if beta is None:
        p0 = lg.ClusterPartialWave.from_uw(k, tm.radial_table(0), 0)
        p2 = lg.ClusterPartialWave.from_uw(k, tm.radial_table(2), 2)
        out["cross_checks"] = dict(
            q_tree_tables=float(lg.alpha_d_quadrupole_fm2(p0, p2)),
            q_tree_tables_note="lg.alpha_d_quadrupole_fm2 on the TaggedModel "
                               "tables (k <= 1.2 GeV), default numerics")
    return out


# ---------------------------------------------------------------- verdict

def verdict(models, ref):
    """Per model: the eta pull and pass/fail at ETA_MAX_PULL; the Q_d
    distance, RECORDED.  The row passes iff every model's eta passes."""
    e, q = ref["eta_d"], ref["Q_d"]
    for m in models.values():
        m["eta_pull"] = (m["eta"] - e["value"]) / e["sigma"]
        m["eta_status"] = "pass" if abs(m["eta_pull"]) <= ETA_MAX_PULL else "fail"
        m["q_d_distance"] = m["q_d"] - q["value"]
        m["q_d_rel_distance"] = m["q_d_distance"] / q["value"]
    return "pass" if all(m["eta_status"] == "pass" for m in models.values()) \
        else "fail"


def measure(path=DATA):
    ref = load_reference(path)
    models = dict(cdbonn=cdbonn(), av18=av18(), hulthen=hulthen())
    band = {b: hulthen(beta=b) for b in HULTHEN_BETA_BAND}
    status = verdict(models, ref)
    for b in band.values():                 # recorded, never in the verdict
        b["eta_pull"] = (b["eta"] - ref["eta_d"]["value"]) / ref["eta_d"]["sigma"]
    return ref, models, band, status


def run(path=DATA, verbose=True):
    t0 = time.perf_counter()
    ref, models, band, status = measure(path)
    e, q = ref["eta_d"], ref["Q_d"]
    failing = [n for n in MODELS if models[n]["eta_status"] == "fail"]
    generator = ("eta: %s; Q_d (impulse approx., recorded): %s fm^2" % (
        ", ".join("%s %.6f (%+.2f sigma)" % (n, models[n]["eta"],
                                             models[n]["eta_pull"])
                  for n in MODELS),
        " / ".join("%.6f (%+.2f %%)" % (models[n]["q_d"],
                                        100 * models[n]["q_d_rel_distance"])
                   for n in MODELS)))
    reference = ("eta_d = %s (Rodning-Knutson PRC 41 (1990) 898), Q_d = %s "
                 "(Bishop-Cheung PRA 20 (1979) 381); second-hand transcription"
                 % (e["as_printed"], q["as_printed"]))
    tolerance = ("abs(eta - eta_d) <= %g sigma for each of %s; Q_d recorded, not "
                 "asserted" % (ETA_MAX_PULL, ", ".join(MODELS)))
    import lipolgen as lg  # noqa: WPS433
    report = dict(name=NAME, generator=generator, reference=reference,
                  tolerance=tolerance, status=status, models=models,
                  hulthen_beta_band=band, reference_points=ref,
                  failing=failing,
                  config=dict(lipolgen=getattr(lg, "__version__", None),
                              numpy=np.__version__,
                              python=sys.version.split()[0],
                              data_file_sha256=FILE_SHA256))
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        print("  model     eta          pull     Q_d(IA) fm^2  vs %.4f     P_D      "
              "norm" % q["value"])
        for n in MODELS:
            m = models[n]
            print("  %-8s  %.7f  %+7.2f   %.6f      %+.5f (%+.2f %%)  %.5f  %.6f"
                  % (n, m["eta"], m["eta_pull"], m["q_d"], m["q_d_distance"],
                     100 * m["q_d_rel_distance"], m["p_d"], m["norm"]))
        c = models["cdbonn"]
        print("cdbonn: eta from the closed-form u, w at r = 60 fm %.12f (tree "
              "eta() %.12f); Q_d quadrature error %.1e; tree "
              "alpha_d_quadrupole_fm2 %.6f (default numerics), %.7f (converged)"
              % (c["eta_check_r60"], c["eta"], c["q_d_quad_error"],
                 c["cross_checks"]["q_tree_default"],
                 c["cross_checks"]["q_tree_converged"]))
        a = models["av18"]
        x = a["cross_checks"]
        print("av18:   r block: eta %.7f (spread over the window %.1e), kappa "
              "%.7f fm^-1, A_S %.6f, P_D %.6f, Q_d %.6f (quadrature error %.1e); "
              "file header: eta %.6f, as %.6f, dstate %.6f, qm %.6f; the k "
              "block the generator reads: Q_d %.6f (%+.2f %% of the r block), "
              "P_D %.6f, eta at r = %g fm %.6f (%+.2f %%)"
              % (a["eta"], a["eta_numerical_error"], a["kappa_fm"],
                 a["a_s"], a["p_d"], a["q_d"], a["q_d_quad_error"],
                 x["header_eta"], x["header_as"], x["header_dstate"],
                 x["header_qm"], x["q_tree_kblock"],
                 100 * (x["q_tree_kblock"] / a["q_d"] - 1),
                 a["p_d_tree"], ETA_WINDOW_FM[0], x["eta_kblock_r"],
                 100 * (x["eta_kblock_r"] / a["eta"] - 1)))
        hm = models["hulthen"]
        print("hulthen: N_0 = %.9f, N_2 = %.9f GeV-units (TaggedModel); eta "
              "from the large-r ratio at 80 fm %.10f; continuum-normalised "
              "N_0 differs by %+.2e; Q_d quadrature error %.1e; tree "
              "alpha_d_quadrupole_fm2 on the tables %.6f"
              % (hm["n0"], hm["n2"], hm["eta_check_r80"],
                 hm["n0_continuum_over_tree"] - 1, hm["q_d_quad_error"],
                 hm["cross_checks"]["q_tree_tables"]))
        print("hulthen beta band (RECORDED, not asserted, nothing moved): %s"
              % ", ".join("beta = %.2f: eta %.5f (%+.1f sigma), Q_d %.4f"
                          % (b, v["eta"], v["eta_pull"], v["q_d"])
                          for b, v in sorted(band.items())))
        for n in MODELS:
            m = models[n]
            print("SUBROW | %s[%s] | eta = %.7f (%+.2f sigma); Q_d = %.6f fm^2 "
                  "(%+.2f %%, recorded) | %s | abs(pull) <= %g | %s"
                  % (NAME, n, m["eta"], m["eta_pull"], m["q_d"],
                     100 * m["q_d_rel_distance"], m["label"], ETA_MAX_PULL,
                     m["eta_status"].upper()))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, generator, reference, tolerance,
            status.upper() + ("" if status == "pass" else
                              " (recorded: %s outside; nothing moved)"
                              % ", ".join("%s %+.1f sigma"
                                          % (n, models[n]["eta_pull"])
                                          for n in failing))))
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
