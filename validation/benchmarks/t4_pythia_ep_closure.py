#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T4 (BENCHMARK_PLAN.md sec. 2 T4 "PYTHIA 8 ep"; sec. 3 "T4 beyond row 7 ...
not run"; 01_generators.md sec. 4 D-1 "broaden, don't replace: add HFS
energy/pT flow and a second kinematic point") -- LiPolGen's PYTHIA bridge
(the T2 tier) against STOCK PYTHIA 8 e p DIS, at two (x, Q^2) windows, on
four hadronic-final-state observables.

WHAT IS COMPARED.  Generator side: `lipolgen.Pipeline` on a PROTON target
(`make_config(isotope="p", channel="inclusive", config=1)`: e 10 GeV x p 100
GeV/c, head-on, proton along +z), unpolarised (`helicity_flip_plan` at
J = 1/2, P_z = P_e = 0), its sampler grid restricted to the window, with the
T2 tier bound exactly as `lipolgen-run --hadronize` binds it
(`set_pythia_hadronizer(cfg, PythiaBridge(beams, options))`, default
`PythiaBridgeOptions` but the seed).  Reference side: stock PYTHIA 8
through its own Python module, `WeakBosonExchange:ff2ff(t:gmZ)` with the
settings of the tree's own D-1 gate (tests/test_pythia.cpp, = PolarizedLithium
Sim's gen_dis_hfs.py: frameType 2, dipoleRecoil, pTmaxMatch 2, PDF:lepton off,
QEDshowerByL off, the two silent-cut fixes pTHatMinDiverge = mHatMin = 0.5)
PLUS `SpaceShower:QEDshowerByQ = off`, the one shower setting the bridge
applies that D-1's stock run does not (docs/PYTHIA_BRIDGE.md sec. 7: QED ISR
off the quark recoils against the lepton), and `PhaseSpace:Q2Min` at the
window's lower Q^2 edge.  Same beams (Beams:eA is the energy of a 100 GeV/c
proton at PYTHIA's own m_p), same library (both sides map ONE libpythia8, read
from /proc/self/maps and recorded), same PDF: PYTHIA's default proton set on
both (PDF:pSet recorded; the bridge draws its struck flavour from
`getPDFPtr(2212)` of that same default and sets no proton-PDF key -- its
only PDF:* lines are PDF:lepton = off, as the stock's, and the Pomeron
instance's PDF:PomSet / PDF:PomRescale; asserted from `applied_settings()`).  This is the configuration in which the
two SHOULD agree: a free proton at rest, no Fermi motion, no off-shellness,
no nucleus, no QED on either lepton leg.  What still differs, and is what
this row measures: the hard process (stock: PYTHIA's own LO gamma*/Z
matrix element, its own scale choice and phase space; bridge: gamma*-only
flavour draw e_q^2 x f_q at scale Q^2 through LHAup, SCALUP = Q), the surrogate
beam + Lorentz frame map, and the massless LHA lepton.

NO VENDORED DATA.  The reference is a generator run live in the same
process, not a published table, so there is nothing under data/ to vendor or
sha256-check (README.md convention 2 does not apply); what makes the row
reproducible is the configuration it records -- the stock settings verbatim,
the bridge's `applied_settings()`, PDF:pSet, the seeds, the PYTHIA version
and library path -- and both samples are bit-reproducible at fixed seeds.

WINDOWS (fixed before any comparison was run).  W1 = x in [0.005, 0.10],
Q^2 in [4, 30] GeV^2: D-1's own window.  W2 = x in [0.10, 0.50], Q^2 in
[30, 300] GeV^2: the second kinematic point D-1 asks for (valence x, an order
of magnitude higher Q^2).  Both lie inside the generator's Scenario (y in
[0.004, 0.985], W^2 >= 8, E' >= 0.3, eta_e in [-3.8, 3.8]); the SAME
Scenario cuts are applied to the stock events, read from the config.

OBSERVABLES, per event, over the hadronic final state = every final-state
particle except the scattered electron (stock: the bottom copy of the hard
process's outgoing lepton; bridge: `Role.Hadron`, `Status.Final`, via
`lipolgen.hfs_arrays`), in the head-on lab frame, no pT or eta cut:
  n_ch        the number of charged particles;
  sum_empz    sum(E - p_z) over all of them (charged and neutral);
  sum_pt_ch   the scalar sum of charged-particle p_T;
  n_ch_eta[a,b)  the charged count in lab pseudorapidity bins with edges
              -inf, -2, 0, 2, 4, +inf (the charged eta flow).
(x, Q^2) MATCHING.  Each generator fills the window with its own (x, Q^2)
density (PYTHIA's LO cross section vs LiPolGen's kernel on ToyF2), and every
observable above moves with (x, Q^2) -- sum_empz IS 2 E_e y.  So the bridge's
mean is taken in 4 x 4 log-uniform (x, Q^2) sub-bins of the window and
re-weighted to the STOCK sample's sub-bin populations; stock sub-bins with no
bridge event are dropped from both sides (counted).  The un-matched ratio is
recorded next to it and never enters the verdict.

TOLERANCE (declared before any number was produced): D-1's +-20 %.  Every
observable in every window: |bridge_matched / stock - 1| <= 0.20.  An eta
bin is COMPARED only if the stock sample has >= 100 charged particles in it
in that window (declared with the bins, before the run); otherwise it is
reported and not compared.  The row passes iff every compared sub-row
passes; otherwise FAIL with the worst distance.  The statistical error of
each ratio is printed and does not enter the verdict (N = 20 000 events per
window per generator puts it far below 0.20).

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: unpolarised e p on both sides
     (P_z = P_e = 0, a spin-1/2 target has no rank-2 sector); no A-scaling
     assumption anywhere, A = 1.
  2. isoscalar denominator -- NOT APPLICABLE: a FREE PROTON on both sides,
     PYTHIA's own proton PDF on both; no nuclear grid, no average nucleon to
     divide by.  Asserted instead: the target has A = 1 and every bridge event
     struck a proton (pdg 2212).
  3. per nucleus / per nucleon -- A = 1, so per nucleus == per nucleon; and
     nothing compared is a cross section or a structure function: all four
     observables are per-EVENT averages of the hadronic final state at matched
     (x, Q^2).  Neither generator's rate enters.

WHAT THIS DOES NOT VALIDATE.  PYTHIA itself: both sides run the same
shower and string code, so this checks the BRIDGE (surrogate, frame map,
flavour draw, LHA record) against PYTHIA's own DIS, not PYTHIA against HERA
(that is Rivet's, T6).  Nothing nuclear: no Fermi motion, no off-shell
nucleon, no neutron instance, no tagged or coherent (Pomeron) record, no
nuclear effect on the final state (none is modelled -- PYTHIA_BRIDGE.md
sec. 9).  No cross section, no spin, no QED radiation (off on both sides).
PYTHIA version: this row ran PYTHIA 8.312 (conda-forge); the project pins
8.317 -- the version actually used is recorded in `config`.

MEASURED 2026-09-26 (PYTHIA 8.312, after the tolerance, windows, bins,
threshold and matching above were fixed; nothing was changed after the first
numbers were seen): PASS, max |bridge/stock - 1| = 0.023 over 15 compared
sub-rows (W2 n_ch_eta[-2,0), 0.977 +- 0.068 statistical).  <n_ch> matched
ratio 1.0034 +- 0.0037 (W1: stock 8.027, bridge 8.054) and 0.9976 +- 0.0038
(W2: 7.261 / 7.243); sum_empz 1.0065 / 1.0024, sum_pt_ch 1.0058 / 1.0006;
every other compared ratio within 0.016 of 1.  W2's n_ch_eta[-inf,-2) holds
3 stock particles and is not compared.  Un-matched, sum_empz reads 1.067
(W1) / 1.029 (W2): the two (x, Q^2) densities differ, which is what the
matching is for.  Wall time 9.7-10.0 s for the whole row (20 000 + 20 000
events per window; stock 2.5 + 3.6 s, bridge 1.8 + 1.7 s).  After the first
full run ONE harness detail was corrected: the stock side's (x, Q^2) now use
PYTHIA's own massive beam electron (event entry 2), not a massless
(10, 0, 0, -10); it moved the stock means in the fourth digit (n_ch W1
8.0261 -> 8.0267) and changed no tolerance, window, bin or threshold.

EXPENSIVE = True: a generator-vs-generator sample is an expensive row by
kind (BENCHMARK_PLAN.md sec. 3), so its full-size pytest runs only with
LIPOLGEN_BENCH=1; a small-N smoke run is always on.  (Measured, the full row
takes ~10 s, under README.md's 30 s line: the flag is set by the row's KIND,
and whether to keep it opt-in at this cost is the README owner's call.)

Run:  source env.sh && python3 validation/benchmarks/t4_pythia_ep_closure.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (any
exception, lipolgen not importable included).  BLOCKED, with an
`environment:` reason, is only the PYTHIA 8 tier off (HAVE_PYTHIA8) or the
pythia8 module missing.  A harness is a measurement, not a CI gate: the
pytests are the gate.
"""

import contextlib
import io
import math
import sys
import time
import traceback

NAME = "t4_pythia_ep_closure"
EXPENSIVE = True            # full N costs a generator-vs-generator sample
TOL = 0.20                  # D-1's +-20 %, declared before any run
N_EVENTS = 20000            # per window, per generator
SEED = 20260926
BEAM_CONFIG = 1             # default_configs("p")[1]: e 10 GeV x p 100 GeV/c
PYTHIA_PINNED = "8.317"     # docs/PYTHIA_BRIDGE.md, env.sh
N_SUB = 4                   # (x, Q^2) matching sub-bins per axis, log-uniform
ETA_EDGES = (-math.inf, -2.0, 0.0, 2.0, 4.0, math.inf)
MIN_STOCK_COUNTS = 100      # an eta bin is compared only above this
MAX_TRY_FACTOR = 100        # stock: at most this many next() per wanted event

WINDOWS = (
    dict(label="W1", x=(0.005, 0.10), q2=(4.0, 30.0),
         note="D-1's own window (tests/test_pythia.cpp)"),
    dict(label="W2", x=(0.10, 0.50), q2=(30.0, 300.0),
         note="the second kinematic point: valence x, higher Q2"),
)

# tests/test_pythia.cpp D-1, verbatim, + the bridge's QEDshowerByQ = off.
STOCK_SETTINGS_D1 = (
    "Beams:frameType = 2",
    "Beams:idA = 2212",
    "Beams:idB = 11",
    "WeakBosonExchange:ff2ff(t:gmZ) = on",
    "PhaseSpace:pTHatMinDiverge = 0.5",
    "PhaseSpace:mHatMin = 0.5",
    "SpaceShower:dipoleRecoil = on",
    "SpaceShower:pTmaxMatch = 2",
    "PDF:lepton = off",
    "TimeShower:QEDshowerByL = off",
    "Random:setSeed = on",
    "Next:numberCount = 0",
    "Next:numberShowInfo = 0",
    "Next:numberShowProcess = 0",
    "Next:numberShowEvent = 0",
    "Print:quiet = on",
)
STOCK_SETTINGS_ADDED = ("SpaceShower:QEDshowerByQ = off",)


def eta_labels():
    def f(v):
        return "-inf" if v == -math.inf else "+inf" if v == math.inf else "%g" % v
    return ["n_ch_eta[%s,%s)" % (f(a), f(b)) for a, b in zip(ETA_EDGES[:-1], ETA_EDGES[1:])]


OBSERVABLES = ["n_ch", "sum_empz", "sum_pt_ch"] + eta_labels()


# ------------------------------------------------------------ environment

def environment():
    """None when both sides can run here, else the blocking reason.

    lipolgen itself is NOT an environment reason: a missing or broken
    lipolgen (numpy absent, an undefined symbol) is the harness broken, so
    its ImportError propagates and `__main__` exits 2, like every other
    harness (validation/benchmarks/README.md convention 6; until
    2026-09-27 it was a BLOCKED row with exit 0).  BLOCKED is kept for the
    PYTHIA 8 tier being off and for the pythia8 module."""
    import lipolgen as lg                                  # noqa: WPS433
    if not getattr(lg, "HAVE_PYTHIA8", False):
        return "environment: lipolgen built without PYTHIA 8 (HAVE_PYTHIA8 false)"
    try:
        import pythia8                                     # noqa: F401,WPS433
    except Exception:
        return "environment: pythia8 Python module not importable (stock side)"
    return None


def pythia_libraries():
    """Every libpythia8 mapped into this process (Linux), so the row says
    whether the stock module and the bridge ran ONE library."""
    try:
        with open("/proc/self/maps") as f:
            return sorted({ln.split()[-1] for ln in f if "libpythia8" in ln})
    except OSError:                                        # pragma: no cover
        return None


# ------------------------------------------------------------ kinematics

def _dot(a, b):
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def dis_kinematics(k, p, kp):
    """(x, Q2, y, W2) from beam e k, target p, scattered e kp (E, px, py, pz)."""
    q = tuple(a - b for a, b in zip(k, kp))
    q2 = -_dot(q, q)
    pq = _dot(p, q)
    x = q2 / (2.0 * pq)
    y = pq / _dot(p, k)
    w2 = _dot(p, p) + 2.0 * pq - q2
    return x, q2, y, w2


def eta_of(px, py, pz):
    pt = math.hypot(px, py)
    if pt == 0.0:
        return math.inf if pz > 0 else -math.inf
    return math.asinh(pz / pt)


def eta_bin(eta):
    for i in range(len(ETA_EDGES) - 1):
        if ETA_EDGES[i] <= eta < ETA_EDGES[i + 1] or (
                i == len(ETA_EDGES) - 2 and eta == math.inf):
            return i
    return None                                            # pragma: no cover


def hfs_row(parts):
    """The per-event observables from [(E, px, py, pz, charge), ...]."""
    n_ch, empz, ptch = 0, 0.0, 0.0
    etas = [0] * (len(ETA_EDGES) - 1)
    for e, px, py, pz, q in parts:
        empz += e - pz
        if q != 0.0:
            n_ch += 1
            ptch += math.hypot(px, py)
            etas[eta_bin(eta_of(px, py, pz))] += 1
    return [float(n_ch), empz, ptch] + [float(c) for c in etas]


# ---------------------------------------------------------------- config

def make_config(window, n_events, seed):
    """The generator-side PipelineConfig: a proton, the grid = the window."""
    import lipolgen as lg                                  # noqa: WPS433
    cfg = lg.make_config(isotope="p", channel="inclusive", config=BEAM_CONFIG,
                         events=n_events, seed=seed)
    g = cfg.grid
    g.x_min, g.x_max = window["x"]
    g.q2_min, g.q2_max = window["q2"]
    cfg.grid = g
    return cfg


def make_plan():
    import lipolgen as lg                                  # noqa: WPS433
    return lg.make_plan("helicity-flip", j=0.5, pz=0.0, pe=0.0)


def scenario_cuts(cfg):
    sc = cfg.scenario
    return dict(y_min=sc.y_min, y_max=sc.y_max, w2_min=sc.w2_min,
                e_prime_min=sc.e_prime_min, eta_min=sc.eta_min,
                eta_max=sc.eta_max, x_max=sc.x_max)


def in_window(window, x, q2):
    return (window["x"][0] <= x <= window["x"][1]
            and window["q2"][0] <= q2 <= window["q2"][1])


def stock_settings(window, seed, e_a, e_b):
    return (list(STOCK_SETTINGS_D1) + list(STOCK_SETTINGS_ADDED)
            + ["Random:seed = %d" % (1 + seed % 900000000),
               "PhaseSpace:Q2Min = %.12g" % window["q2"][0],
               "Beams:eA = %.12g" % e_a, "Beams:eB = %.12g" % e_b])


# ------------------------------------------------------------ the samples

def stock_sample(window, n_events, seed, cuts, beams):
    """Stock PYTHIA e p, selected after the fact into the window."""
    import pythia8                                         # noqa: WPS433
    with contextlib.redirect_stdout(io.StringIO()):
        py = pythia8.Pythia("", False)
    p_u, e_e = beams["p_u"], beams["e_e"]
    m_p = py.particleData.m0(2212)
    e_a = math.sqrt(p_u * p_u + m_p * m_p)
    settings = stock_settings(window, seed, e_a, e_e)
    for s in settings:
        if not py.readString(s):
            raise RuntimeError("stock PYTHIA refused setting %r" % s)
    if not py.init():
        raise RuntimeError("stock PYTHIA init() failed")
    info = dict(settings=settings, pdf_pset=py.settings.word("PDF:pSet"),
                version=float(py.settings.parm("Pythia:versionNumber")),
                m_p=m_p, e_a=e_a)
    rows, xs, q2s = [], [], []
    n_try = n_fail = n_fallback = 0
    while len(rows) < n_events and n_try < MAX_TRY_FACTOR * n_events:
        n_try += 1
        if not py.next():
            n_fail += 1
            continue
        ev = py.event
        # the hard process's outgoing lepton, followed to its bottom copy
        i_e = -1
        for j in range(ev.size()):
            if ev[j].id() == 11 and abs(ev[j].status()) == 23:
                i_e = ev[j].iBotCopyId()
                break
        if i_e < 0 or not ev[i_e].isFinal():
            n_fallback += 1
            best = -1.0
            for j in range(ev.size()):
                if ev[j].isFinal() and ev[j].id() == 11 and ev[j].e() > best:
                    best, i_e = ev[j].e(), j
            if i_e < 0:
                continue
        pe = ev[i_e]
        kp = (pe.e(), pe.px(), pe.py(), pe.pz())
        pb, eb = ev[1], ev[2]            # PYTHIA's own beams A (p) and B (e)
        p = (pb.e(), pb.px(), pb.py(), pb.pz())
        k = (eb.e(), eb.px(), eb.py(), eb.pz())
        x, q2, y, w2 = dis_kinematics(k, p, kp)
        if not in_window(window, x, q2):
            continue
        eta_e = eta_of(kp[1], kp[2], kp[3])
        if not (cuts["y_min"] <= y <= cuts["y_max"] and w2 >= cuts["w2_min"]
                and kp[0] >= cuts["e_prime_min"]
                and cuts["eta_min"] <= eta_e <= cuts["eta_max"]
                and x <= cuts["x_max"]):
            continue
        parts = []
        for j in range(ev.size()):
            if j == i_e:
                continue
            a = ev[j]
            if a.isFinal():
                parts.append((a.e(), a.px(), a.py(), a.pz(), a.charge()))
        rows.append(hfs_row(parts))
        xs.append(x)
        q2s.append(q2)
    info.update(n_try=n_try, n_fail=n_fail, n_fallback=n_fallback,
                n_accepted=len(rows), efficiency=len(rows) / max(1, n_try))
    return dict(x=xs, q2=q2s, rows=rows, info=info)


def bridge_sample(window, n_events, seed):
    """LiPolGen's T2 tier, bound the way `lipolgen-run --hadronize` binds it."""
    import lipolgen as lg                                  # noqa: WPS433
    cfg = make_config(window, n_events, seed)
    beams = lg.default_configs(cfg.isotope)[cfg.beam_config]
    opts = lg.PythiaBridgeOptions()
    opts.seed = seed
    opts.f2_source = cfg.unpol_sf_obj           # what cli.py does, verbatim
    bridge = lg.PythiaBridge(beams, opts)
    lg.set_pythia_hadronizer(cfg, bridge)
    pipe = lg.Pipeline(cfg, make_plan())
    cols = pipe.generate(0, events=True, nthreads=1)   # bridge not re-entrant
    evs = cols["events"]
    h = lg.hfs_arrays(evs)
    offs, p4, ch = h["offsets"], h["p4"], h["charge"]
    rows, xs, q2s = [], [], []
    n_veto = 0
    for i in range(len(evs)):
        lo, hi = int(offs[i]), int(offs[i + 1])
        if hi == lo:
            n_veto += 1                  # PYTHIA vetoed: no T2 record
            continue
        parts = [(float(p4[j, 0]), float(p4[j, 1]), float(p4[j, 2]),
                  float(p4[j, 3]), float(ch[j])) for j in range(lo, hi)]
        rows.append(hfs_row(parts))
        xs.append(float(h["x"][i]))
        q2s.append(float(h["q2"][i]))
    st = bridge.stats
    struck = [int(v) for v in cols["struck_pdg"]]
    applied = list(bridge.applied_settings())
    info = dict(
        isotope=cfg.isotope, A=int(beams.ion.A), e_e=beams.electron_energy,
        p_u=beams.ion_momentum_per_nucleon, n_generated=len(evs),
        n_vetoed=n_veto, n_ok=int(st.n_ok), n_failed=int(st.n_failed),
        n_retries=int(st.n_retries), n_proton=int(st.n_proton),
        n_neutron=int(st.n_neutron),
        max_rescale_dev=float(st.max_rescale_dev),
        struck_all_proton=all(s == 2212 for s in struck),
        applied_settings=applied,
        pdf_setting_applied=[s for s in applied if s.startswith("PDF:")
                             and not s.startswith("PDF:lepton")
                             and not s.startswith("PDF:Pom")],
        unpol_sf=str(cfg.unpol_sf), scenario=scenario_cuts(cfg),
        x_range=(min(xs), max(xs)) if xs else None,
        q2_range=(min(q2s), max(q2s)) if q2s else None)
    return dict(x=xs, q2=q2s, rows=rows, info=info, cfg=cfg, beams=beams)


# -------------------------------------------------------------- matching

def sub_bin(window, x, q2):
    lx0, lx1 = (math.log(v) for v in window["x"])
    lq0, lq1 = (math.log(v) for v in window["q2"])
    ix = min(N_SUB - 1, max(0, int(N_SUB * (math.log(x) - lx0) / (lx1 - lx0))))
    iq = min(N_SUB - 1, max(0, int(N_SUB * (math.log(q2) - lq0) / (lq1 - lq0))))
    return ix * N_SUB + iq


def _mean_var(vals):
    n = len(vals)
    m = sum(vals) / n
    v = sum((a - m) ** 2 for a in vals) / (n - 1) if n > 1 else 0.0
    return m, v


def compare(window, stock, bridge):
    """Per observable: stock mean, bridge mean re-weighted to the stock's
    (x, Q2) sub-bin populations, the raw bridge mean, ratios and errors."""
    nb = N_SUB * N_SUB
    sb = [[] for _ in range(nb)]
    bb = [[] for _ in range(nb)]
    for x, q2, r in zip(stock["x"], stock["q2"], stock["rows"]):
        sb[sub_bin(window, x, q2)].append(r)
    for x, q2, r in zip(bridge["x"], bridge["q2"], bridge["rows"]):
        bb[sub_bin(window, x, q2)].append(r)
    keep = [b for b in range(nb) if sb[b] and bb[b]]
    dropped = sum(len(sb[b]) for b in range(nb) if sb[b] and not bb[b])
    n_s = sum(len(sb[b]) for b in keep)
    out = []
    for k, name in enumerate(OBSERVABLES):
        s_vals = [r[k] for b in keep for r in sb[b]]
        s_mean, s_var = _mean_var(s_vals)
        s_err = math.sqrt(s_var / len(s_vals))
        b_mean, b_var = 0.0, 0.0
        for b in keep:
            w = len(sb[b]) / n_s
            m, v = _mean_var([r[k] for r in bb[b]])
            b_mean += w * m
            b_var += w * w * v / len(bb[b])
        b_err = math.sqrt(b_var)
        raw_mean, _ = _mean_var([r[k] for r in bridge["rows"]])
        raw_stock, _ = _mean_var([r[k] for r in stock["rows"]])
        ratio = b_mean / s_mean if s_mean else float("nan")
        rerr = (abs(ratio) * math.hypot(s_err / s_mean, b_err / b_mean)
                if s_mean and b_mean else float("nan"))
        is_eta = name.startswith("n_ch_eta")
        n_stock_particles = sum(s_vals) if is_eta else None
        compared = (not is_eta) or n_stock_particles >= MIN_STOCK_COUNTS
        ok = compared and math.isfinite(ratio) and abs(ratio - 1.0) <= TOL
        out.append(dict(
            window=window["label"], observable=name, stock=s_mean,
            stock_err=s_err, bridge=b_mean, bridge_err=b_err,
            ratio=ratio, ratio_err=rerr, distance=abs(ratio - 1.0),
            bridge_raw=raw_mean, stock_raw=raw_stock,
            ratio_raw=raw_mean / raw_stock if raw_stock else float("nan"),
            n_stock_particles=n_stock_particles, compared=compared,
            status=("pass" if ok else "fail") if compared else "not compared"))
    return out, dict(n_sub_kept=len(keep), n_sub_total=nb,
                     n_stock_dropped=dropped, n_stock_matched=n_s,
                     n_bridge=len(bridge["rows"]))


def mean_kinematics(sample):
    n = len(sample["x"])
    return (sum(sample["x"]) / n, sum(sample["q2"]) / n) if n else (None, None)


# ------------------------------------------------------------------- run

def measure(n_events=N_EVENTS, seed=SEED, windows=WINDOWS):
    blocked = environment()
    if blocked:
        return dict(blocked=blocked)
    res = dict(windows=[], subrows=[])
    for w_i, window in enumerate(windows):
        t0 = time.perf_counter()
        br = bridge_sample(window, n_events, seed + 100 + w_i)
        t1 = time.perf_counter()
        beams = dict(p_u=br["info"]["p_u"], e_e=br["info"]["e_e"])
        stk = stock_sample(window, n_events, seed + w_i,
                           br["info"]["scenario"], beams)
        t2 = time.perf_counter()
        rows, match = compare(window, stk, br)
        res["subrows"] += rows
        res["windows"].append(dict(
            label=window["label"], x=window["x"], q2=window["q2"],
            note=window["note"], match=match,
            stock_mean_x_q2=mean_kinematics(stk),
            bridge_mean_x_q2=mean_kinematics(br),
            stock=stk["info"], bridge=br["info"],
            seconds_bridge=t1 - t0, seconds_stock=t2 - t1))
    res["pythia_libraries"] = pythia_libraries()
    return res


def run(verbose=True, n_events=N_EVENTS, seed=SEED, windows=WINDOWS):
    t0 = time.perf_counter()
    m = measure(n_events=n_events, seed=seed, windows=windows)
    tol_txt = ("|bridge/stock - 1| <= %.2f per observable per window, at "
               "matched (x, Q2) (D-1)" % TOL)
    if m.get("blocked"):
        report = dict(name=NAME, generator="not run -- %s" % m["blocked"],
                      reference="stock PYTHIA 8 e p ff2ff(t:gmZ)",
                      tolerance=tol_txt, status="blocked", reason=m["blocked"],
                      subrows=[], windows=[], config=dict(expensive=EXPENSIVE))
        if verbose:
            print("REPORT | %s | %s | %s | %s | BLOCKED" % (
                NAME, report["generator"], report["reference"], tol_txt))
        report["runtime_s"] = time.perf_counter() - t0
        return report
    subs = m["subrows"]
    comp = [s for s in subs if s["compared"]]
    fails = [s for s in comp if s["status"] != "pass"]
    status = "fail" if fails else "pass"
    worst = max(comp, key=lambda s: s["distance"])
    w0 = m["windows"][0]
    version = w0["stock"]["version"]
    by = {(s["window"], s["observable"]): s for s in subs}
    nch = "; ".join("<n_ch> %s %.3f" % (w["label"], by[(w["label"], "n_ch")]["ratio"])
                    for w in m["windows"])
    report = dict(
        name=NAME,
        generator=("max |bridge/stock - 1| = %.3f (%s %s) over %d compared "
                   "sub-rows, %d fail; %s"
                   % (worst["distance"], worst["window"], worst["observable"],
                      len(comp), len(fails), nch)),
        reference=("stock PYTHIA %.3f e p WeakBosonExchange:ff2ff(t:gmZ), "
                   "D-1 settings, e 10 x p 100, PDF:pSet %s, %d + %d events"
                   % (version, w0["stock"]["pdf_pset"],
                      m["windows"][0]["match"]["n_stock_matched"],
                      m["windows"][-1]["match"]["n_stock_matched"])),
        tolerance=tol_txt,
        status=status,
        worst=dict(window=worst["window"], observable=worst["observable"],
                   ratio=worst["ratio"], distance=worst["distance"]),
        subrows=subs,
        windows=m["windows"],
        config=dict(
            expensive=EXPENSIVE, n_events=n_events, seed=seed,
            beam_config=BEAM_CONFIG, pythia_version_used=version,
            pythia_version_pinned=PYTHIA_PINNED,
            pythia_libraries=m["pythia_libraries"],
            pdf_pset=w0["stock"]["pdf_pset"], n_sub=N_SUB,
            eta_edges=list(ETA_EDGES), min_stock_counts=MIN_STOCK_COUNTS,
            stock_settings_d1=list(STOCK_SETTINGS_D1),
            stock_settings_added=list(STOCK_SETTINGS_ADDED),
            bridge_applied_settings=w0["bridge"]["applied_settings"]),
    )
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        print("PYTHIA used %.3f (pinned %s); libraries mapped: %s; PDF:pSet %s"
              % (version, PYTHIA_PINNED, m["pythia_libraries"], w0["stock"]["pdf_pset"]))
        for w in m["windows"]:
            s, b, mt = w["stock"], w["bridge"], w["match"]
            print("%s x in [%g, %g], Q2 in [%g, %g] -- %s" % (
                w["label"], w["x"][0], w["x"][1], w["q2"][0], w["q2"][1], w["note"]))
            print("   stock : %d events of %d tries (eff %.3f), %d next() "
                  "failures, %d e' fallbacks, <x> = %.4f, <Q2> = %.2f, %.1f s"
                  % (s["n_accepted"], s["n_try"], s["efficiency"], s["n_fail"],
                     s["n_fallback"], w["stock_mean_x_q2"][0],
                     w["stock_mean_x_q2"][1], w["seconds_stock"]))
            print("   bridge: %d events, %d vetoed, %d retries, all struck p: %s, "
                  "max|lambda-1| %.2g, <x> = %.4f, <Q2> = %.2f, %.1f s"
                  % (b["n_generated"], b["n_vetoed"], b["n_retries"],
                     b["struck_all_proton"], b["max_rescale_dev"],
                     w["bridge_mean_x_q2"][0], w["bridge_mean_x_q2"][1],
                     w["seconds_bridge"]))
            print("   matching: %d/%d sub-bins, %d stock events dropped"
                  % (mt["n_sub_kept"], mt["n_sub_total"], mt["n_stock_dropped"]))
            print("   %-20s %10s %10s %8s %7s %8s  %s" % (
                "observable", "stock", "bridge", "ratio", "+-", "raw", "status"))
            for r in subs:
                if r["window"] != w["label"]:
                    continue
                print("   %-20s %10.4f %10.4f %8.4f %7.4f %8.4f  %s%s" % (
                    r["observable"], r["stock"], r["bridge"], r["ratio"],
                    r["ratio_err"], r["ratio_raw"], r["status"].upper(),
                    "" if r["n_stock_particles"] is None
                    else "  (%d stock particles)" % r["n_stock_particles"]))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, report["generator"], report["reference"], tol_txt,
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
