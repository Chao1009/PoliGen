#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T4 row 7 (BENCHMARK_PLAN.md sec. 4): DJANGOH's Rad/noRad four-bin table
against LiPolGen's radiative tail, per Q^2 bin, on a PROTON.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), FOR THIS ROW
----------------------------------------------------------------------
1. b1.  Neither side carries a b1.  DJANGOH is unpolarised e-p; LiPolGen runs
   an UNPOLARISED proton (spin 1/2, `helicity_flip_plan(0.5, 0, 0)`: P_z = 0,
   P_e = 0), where every rank-2 slot is identically zero.  No A-scaling
   assumption enters, because A = 1 on both sides.
2. Nuclear PDF grids.  None.  DJANGOH uses the FREE-proton CTEQ6.1M set
   (LHAPDF5 10150, LO); LiPolGen's proton Born is `InclusiveKernel::
   dsigma_unpol` of the proton ion with the SHIPPED default unpolarised
   structure functions (`UnpolSfSource.Toy`).  There is no average-nucleon
   grid to divide by 1/2(F2p + F2n); the harness asserts A == 1 instead.
3. Cross sections vs structure functions.  DJANGOH's table is sigma in
   microbarn per e-p collision.  LiPolGen's `dsigma_unpol` is per nucleon.
   With A = 1 per nucleus == per nucleon, and the compared quantity is in any
   case a RATIO (sigma_RC / sigma_Born over one Q^2 bin), which carries no
   normalisation convention at all.  The one absolute number printed (the
   Born cross section per bin) is per e-p collision on both sides.

WHAT IS COMPARED, AND WHY THE ROW IS **BLOCKED** (measured 2026-09-23)
---------------------------------------------------------------------
Reference: sigma(Rad=1)/sigma(Rad=0) of DJANGOH 4.6.21 (HERACLES), e- 18 GeV x
p 275 GeV, per Q^2 bin [1,10], [10,100], [100,1000], [1000,10000] GeV^2,
x_l in [1e-5, 1], y_l in [1e-4, 1], W_h > 3 GeV.  Vendored in
`data/djangoh_rad_norad_ep18x275.json` with provenance.

Generator: the bin-integrated sigma_RC/sigma_Born = int Born (1 + r) / int
Born, with r = `RcModel.tail_ratio_at(x, Q2, q_n = 0)` = sigma_tail/sigma_Born,
on LiPolGen's proton at the SAME beams (config 2 = 18 x 275 GeV), the same x
floor, W^2 >= 9 GeV^2, and y in [1e-4, 0.985] (0.985 = the shipped
`Scenario.y_max`, which is the top of the tail table's y support -- `RcModel`
refuses to extrapolate past it; DJANGOH's y runs to 1).  Two tail models:
`t-peak` (the shipped default) and `polrad-full` (opt-in).

WHAT THE TWO RC MODELS SHARE AND DO NOT -- the tolerance comes from here.
  DJANGOH Rad=1 carries: the non-radiative O(alpha) piece (virtual + soft,
    INC2 with LPAR(2) = 1), the inelastic hard-photon radiation from the
    lepton (initial and final state) and the Compton term (INC31..INC33),
    with the HADRONIC vertex; and NOT the elastic radiative tail: all eight
    ePIC logs read IEL2 = IEL31 = IEL32 = IEL33 = 0, which the DJANGOH manual
    (MZ-TH/05-15 sec. 3.1.1) defines as elastic ep -> ep and the
    "quasielastic tail ep -> ep gamma" from ISR, FSR and Compton.
  LiPolGen's `rc_tail` carries: ONLY the radiative tails under the DIS bin --
    POLRAD Eq. (38) (t-peak) or Eq. (18) (all three peaks) of the coherent
    elastic vertex, plus the Eq. (44) nucleon-sum quasi-elastic tail.  On a
    proton the coherent spin-1 vertex does not exist (a zero form factor is
    supplied, the only way to build `RcModel` on an ion that
    `HoSpin1FF::for_ion` refuses) and the Eq. (44) slot with Z = 1, N = 0 and
    no Pauli factor (`qe_kf_gev = 0`, a free proton) IS the proton's elastic
    radiative tail.  It carries NO lepton-vertex / inelastic correction as a
    shift (rc.hpp: that is priced by the tensor BAND, never applied).
  => The two share NO O(alpha) term.  DJANGOH's ratio is the inelastic
    correction with the elastic tail switched off; LiPolGen's is the elastic
    tail and nothing else.  They are ADDITIVE, disjoint pieces of one O(alpha)
    correction, so no value of either can contradict the other: the
    "must-not-contradict" bound holds vacuously and no tolerance is
    definable.  A comparison that cannot fail is not a benchmark, so the row
    reports BLOCKED, prints both sides and their distance, and tunes nothing.
  It unblocks with a DJANGOH run at IEL31..IEL33 != 0 (the elastic tail ON)
  at these cuts -- which needs the DJANGOH source (not public, 01_generators.md
  sec. 5.6) or a new ePIC production.  A LiPolGen-side inelastic RC shift
  would not unblock it: none exists and adding one is a physics change.

Run:  source env.sh && python3 validation/benchmarks/t4_djangoh_rad_noRad.py
(17.3-19.6 s measured 2026-09-23; the polrad-full tail table is ~16 s of it).

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import json
import math
import pathlib
import sys
import traceback
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data" / "djangoh_rad_norad_ep18x275.json"

NAME = "t4_djangoh_rad_noRad"
#: LiPolGen beam configuration index for e 18 GeV x p 275 GeV/u.
CONFIG = 2
#: The DJANGOH run's leptonic cuts (runcards KINEM-CUTS), and LiPolGen's own
#: y ceiling -- the RC tail tables refuse y above `RC_TAIL_Y_CEILING`.
X_MIN = 1.0e-5
Y_MIN = 1.0e-4
W2_MIN = 9.0
#: Gauss-Legendre nodes per panel, and panels, in ln Q^2 and in each of the
#: two y pieces (ln y below Y_SPLIT, -ln(1 - y) above it).  Converged: see
#: `convergence()`.
N_GL = 12
PANELS_Q2 = 4
PANELS_Y = 6
Y_SPLIT = 0.5
TAIL_MODELS = ("t-peak", "polrad-full")


def load_reference(path=DATA):
    """The four vendored rows -> list of dicts, with the ratio computed."""
    doc = json.loads(pathlib.Path(path).read_text())
    rows = []
    for lo, hi, rad1, rad0 in doc["rows"]:
        rows.append({"q2_lo": lo, "q2_hi": hi, "sigma_rad1_ub": rad1,
                     "sigma_rad0_ub": rad0, "ratio": rad1 / rad0})
    return doc, rows


def _gl(n):
    t, w = np.polynomial.legendre.leggauss(n)
    return t, w


def _panel_nodes(a, b, panels, n):
    """Composite Gauss-Legendre nodes/weights on [a, b]."""
    t, w = _gl(n)
    edges = np.linspace(a, b, panels + 1)
    nodes, weights = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        h = 0.5 * (hi - lo)
        nodes.append(lo + h * (t + 1.0))
        weights.append(h * w)
    return np.concatenate(nodes), np.concatenate(weights)


def build_pipeline(tail_model, qe_kf_gev=0.0):
    """LiPolGen's proton, RC on, at 18 x 275 with DJANGOH's x/W cuts.

    The ONLY non-default RC settings: a zero coherent form factor (the proton
    has no spin-1 coherent vertex; `HoSpin1FF::for_ion` throws for it) and,
    unless `qe_kf_gev` is passed, no Pauli factor (a free proton).  The
    `tail_model` is the knob under test.  Nothing here is a physics default
    of the library; this is a harness configuration."""
    import lipolgen as lg
    from lipolgen import _lipolgen as _l

    class _NoCoherentVertex(_l.Spin1ElasticFF):
        def fc(self, t):
            return 0.0

        def fm(self, t):
            return 0.0

        def fq(self, t):
            return 0.0

        def provenance(self):
            return ("t4_djangoh_rad_noRad: zero coherent spin-1 vertex (a "
                    "proton has none); the proton elastic tail is the "
                    "Eq. (44) Z = 1, N = 0 slot")

    cfg = lg.make_config(isotope="p", channel="inclusive", config=CONFIG,
                         events=10, seed=1, rc="tensor-band",
                         rc_tail_model=tail_model)
    g = cfg.grid
    g.x_min = X_MIN
    g.q2_min = 1.0
    g.q2_max = 1.0e4
    cfg.grid = g
    sc = cfg.scenario
    sc.w2_min = W2_MIN
    sc.y_min = Y_MIN
    sc.q2_min = 1.0
    # DJANGOH has no scattered-electron acceptance: open LiPolGen's.
    sc.eta_min, sc.eta_max, sc.e_prime_min = -50.0, 50.0, 0.0
    cfg.scenario = sc
    ff = _NoCoherentVertex()
    cfg.rc_options.ff = ff
    cfg.rc_options.qe_kf_gev = qe_kf_gev
    p = lg.Pipeline(cfg, lg.helicity_flip_plan(0.5, 0.0, 0.0))
    assert p.dis_sampler.kernel.ion.A == 1, "rule 2/3: A must be 1"
    # `ff` is returned so the caller holds the Python trampoline for as long
    # as the model lives.
    return p, cfg, ff


def _y_nodes(y_lo, y_hi):
    """Nodes/weights in y on [y_lo, y_hi]: ln y below Y_SPLIT (the Born's
    low-x peak) and -ln(1 - y) above it (the tail's 1/(1 - y) growth)."""
    ys, ws = [], []
    a, b = y_lo, min(y_hi, Y_SPLIT)
    if b > a:
        v, w = _panel_nodes(math.log(a), math.log(b), PANELS_Y, N_GL)
        ys.append(np.exp(v))
        ws.append(w * np.exp(v))
    a, b = max(y_lo, Y_SPLIT), y_hi
    if b > a:
        t, w = _panel_nodes(-math.log(1.0 - a), -math.log(1.0 - b),
                            PANELS_Y, N_GL)
        ys.append(1.0 - np.exp(-t))
        ws.append(w * np.exp(-t))
    if not ys:
        return np.zeros(0), np.zeros(0)
    return np.concatenate(ys), np.concatenate(ws)


def integrate_bin(p, q2_lo, q2_hi, y_max):
    """int Born and int Born * r over one Q^2 bin, in (ln Q^2, y).

    Born = d2sigma/(dx dQ2) [pb/GeV^2], so dx dQ2 = (Q^2/(s y^2)) dy Q^2 dlnQ^2.
    The cuts are DJANGOH's (x >= X_MIN, y >= Y_MIN, W^2 >= W2_MIN, x <= 1)
    with LiPolGen's y ceiling `y_max`."""
    import lipolgen as lg
    ds = p.dis_sampler
    k, m, s = ds.kernel, p.rc_model, ds.s
    m_p2 = lg.PROTON_MASS ** 2
    uq, wq = _panel_nodes(math.log(q2_lo), math.log(q2_hi), PANELS_Q2, N_GL)
    born = born_r = uncovered = 0.0
    for u, wu in zip(uq, wq):
        q2 = math.exp(u)
        y_lo = max(Y_MIN, q2 / s, (W2_MIN - m_p2 + q2) / s)
        y_hi = min(y_max, q2 / (s * X_MIN))
        for y, wy in zip(*_y_nodes(y_lo, y_hi)):
            x = q2 / (s * y)
            b = (k.dsigma_unpol(x, q2, s) * q2 * q2 / (s * y * y)) * wu * wy
            born += b
            try:
                born_r += b * m.tail_ratio_at(x, q2, 0.0)
            except Exception:                 # outside the tail table
                uncovered += b
    return born, born_r, uncovered


def convergence(p, rows, y_max, factor=2):
    """The bin ratios at the shipped node counts and at `factor` x the panels
    in both variables -- the quadrature error bound quoted in the record."""
    global PANELS_Q2, PANELS_Y
    base = (PANELS_Q2, PANELS_Y)
    out = []
    try:
        for pq, py in (base, (base[0] * factor, base[1] * factor)):
            PANELS_Q2, PANELS_Y = pq, py
            vals = []
            for r in rows:
                b, br, u = integrate_bin(p, r["q2_lo"], r["q2_hi"], y_max)
                vals.append(1.0 + br / (b - u))
            out.append(vals)
    finally:
        PANELS_Q2, PANELS_Y = base
    return out


def lipolgen_side(tail_models=TAIL_MODELS, rows=None, y_max=None):
    import lipolgen as lg
    if rows is None:
        rows = load_reference()[1]
    out = {}
    for tm in tail_models:
        t0 = time.time()
        p, cfg, _ff = build_pipeline(tm)
        build_s = time.time() - t0
        # The shipped Scenario.y_max, read and never retyped: it is the top of
        # the tail table's y support (`RcModel` refuses to extrapolate).
        ym = cfg.scenario.y_max if y_max is None else y_max
        per_bin = []
        for r in rows:
            born, born_r, unc = integrate_bin(p, r["q2_lo"], r["q2_hi"], ym)
            per_bin.append({
                "ratio": 1.0 + born_r / (born - unc),
                # dsigma_unpol is d2sigma/(dx dQ2) in pb/GeV^2 (rc.hpp:
                # born_pb_at = x s dsigma_unpol [pb]); pb -> microbarn.
                "born_ub": born * 1e-6,
                "uncovered_born_fraction": unc / born,
            })
        out[tm] = {"bins": per_bin, "build_s": build_s, "y_max": ym,
                   "clipped_node_fraction": p.rc_model.clipped_cell_fraction,
                   "ff_provenance": p.rc_model.ff_provenance,
                   "unpol_sf": str(cfg.unpol_sf)}
    return out


def log_slope(q2_mid, excess):
    """d ln(excess) / d ln Q^2 by least squares -- the Q^2 TREND, as a sign
    and a size, of (ratio - 1)."""
    xs = np.log(np.asarray(q2_mid))
    ys = np.log(np.asarray(excess))
    return float(np.polyfit(xs, ys, 1)[0])


def run(tail_models=TAIL_MODELS, verbose=True):
    doc, rows = load_reference()
    gen = lipolgen_side(tail_models, rows)
    q2_mid = [math.sqrt(r["q2_lo"] * r["q2_hi"]) for r in rows]
    ref_excess = [r["ratio"] - 1.0 for r in rows]
    detail = []
    for i, r in enumerate(rows):
        d = {"bin": (r["q2_lo"], r["q2_hi"]), "djangoh_ratio": r["ratio"],
             "djangoh_sigma_rad0_ub": r["sigma_rad0_ub"]}
        for tm in tail_models:
            b = gen[tm]["bins"][i]
            d[tm] = b["ratio"]
            d[tm + ":excess_over_djangoh_excess"] = (
                (b["ratio"] - 1.0) / (r["ratio"] - 1.0))
            d[tm + ":born_ub"] = b["born_ub"]
            d[tm + ":uncovered"] = b["uncovered_born_fraction"]
        detail.append(d)
    trend = {"djangoh": log_slope(q2_mid, ref_excess)}
    for tm in tail_models:
        trend[tm] = log_slope(q2_mid, [gen[tm]["bins"][i]["ratio"] - 1.0
                                       for i in range(len(rows))])
    # Internal sanity -- the only thing that CAN fail here: a tail is >= 0.
    sane = all(gen[tm]["bins"][i]["ratio"] >= 1.0 and
               math.isfinite(gen[tm]["bins"][i]["ratio"])
               for tm in tail_models for i in range(len(rows)))
    status = "blocked" if sane else "fail"
    reason = (
        "no shared O(alpha) term: DJANGOH Rad=1 ran with IEL2=IEL31=IEL32="
        "IEL33=0 (elastic radiative tail OFF; all eight ePIC logs), LiPolGen's "
        "rc_tail is the elastic(+QE) tail ONLY and applies no inelastic/vertex "
        "shift -- disjoint additive pieces, so must-not-contradict holds "
        "vacuously; unblocks with a DJANGOH run at IEL31..33 != 0"
        if sane else "a LiPolGen tail ratio below 1 or non-finite")
    row = {
        "name": NAME,
        "generator_value": {tm: [round(gen[tm]["bins"][i]["ratio"], 8)
                                 for i in range(len(rows))]
                            for tm in tail_models},
        "reference_value": [round(r["ratio"], 8) for r in rows],
        "tolerance": "none definable (the two RC models share no term)",
        "status": status,
        "reason": reason,
        "trend_dlnexcess_dlnQ2": trend,
        "detail": detail,
        "config": {"beams": "e 18 x p 275 (config 2)", "x_min": X_MIN,
                   "y_min": Y_MIN, "y_max": gen[tail_models[0]]["y_max"],
                   "W2_min": W2_MIN,
                   "rc_ff": "zero coherent vertex", "qe_kf_gev": 0.0,
                   "clipped_node_fraction": {tm: gen[tm]["clipped_node_fraction"]
                                             for tm in tail_models},
                   "tail_build_s": {tm: round(gen[tm]["build_s"], 1)
                                    for tm in tail_models},
                   "unpol_sf": {tm: gen[tm]["unpol_sf"] for tm in tail_models}},
    }
    if verbose:
        print_row(row, tail_models)
    return row


def print_row(row, tail_models):
    w = sys.stdout.write
    w("%-22s %-12s %s\n" % ("Q2 bin [GeV2]", "DJANGOH", "  ".join(
        "%-26s" % ("LiPolGen " + tm) for tm in tail_models)))
    for d in row["detail"]:
        w("%-22s %-12.6f %s\n" % (
            "[%g, %g]" % d["bin"], d["djangoh_ratio"], "  ".join(
                "%-26s" % ("%.6f (x%.4f of DJ excess)" % (
                    d[tm], d[tm + ":excess_over_djangoh_excess"]))
                for tm in tail_models)))
    for d in row["detail"]:
        w("  Born sigma [ub], bin [%g, %g]: DJANGOH Rad=0 %.6g (CTEQ6.1M, y<=1)"
          " | LiPolGen %s (shipped Toy SF, y<=%g) -- informational\n" % (
              d["bin"][0], d["bin"][1], d["djangoh_sigma_rad0_ub"],
              " / ".join("%.6g" % d[tm + ":born_ub"] for tm in tail_models),
              row["config"]["y_max"]))
    w("trend d ln(ratio-1)/d ln Q2: %s\n" % ", ".join(
        "%s %+.4f" % (k, v) for k, v in row["trend_dlnexcess_dlnQ2"].items()))
    w("ROW | %s | generator %s | reference %s | tolerance: %s | %s -- %s\n" % (
        row["name"], row["generator_value"], row["reference_value"],
        row["tolerance"], row["status"].upper(), row["reason"]))


if __name__ == "__main__":
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    try:
        run()
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        sys.exit(2)
    sys.exit(0)
