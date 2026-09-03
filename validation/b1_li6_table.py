#!/usr/bin/env python3
"""Regenerate design D's number tables and its A = 2 validation-gate report.

Design: `docs/open_items/run_2026-09-02/design_D_b1_li6.md` (sections 5 and 7).
Outputs (with `--write`, and it writes NOTHING outside `docs/`):
    docs/open_items/run_2026-09-02/phase_D_gate.md      -- the gate, section 5
    docs/open_items/run_2026-09-02/phase_D_numbers.md   -- the numbers, section 7

Without `--write` everything goes to stdout, which is how you check a change
before it lands in a document.

WHY THIS SCRIPT EXISTS RATHER THAN A DOCTEST.  Checklist item 4 of design
section 5.4 -- the nucleon PDF, which is the DOMINANT term in the gate's
residual -- has to be run with a real global fit, and the core must never link
LHAPDF (`sf.hpp`).  `LhapdfSF` is reachable only from a build that has it and
from outside the kernel, so the one measurement that closes the gate's budget
lives here.  Everything else this script prints is also pinned in
`tests/test_b1_nuclear.cpp`; the script is the reproducibility handle, the
doctest is the regression gate.

READ THE GATE BEFORE QUOTING ANY 6Li NUMBER.  G3b (magnitude) FAILS at the
default settings: a factor 2.27 low at CDKS Eq. (21)'s delta-function, 3.68
at Eq. (17)'s kappa = 1 form.  Under design section 5.4's Escalation clause no
6Li number from `Li6ConvolutionB1` may be published while that stands; the
tables are recorded so the phase is reproducible and its regressions
catchable.

Run:
    source env.sh
    python3 validation/b1_li6_table.py            # to stdout
    python3 validation/b1_li6_table.py --write    # into docs/
"""

import argparse
import pathlib
import sys

import lipolgen as lg

HERE = pathlib.Path(__file__).resolve().parent
DOCS = HERE.parent / "docs" / "open_items" / "run_2026-09-02"

Q2 = 2.5
X_TABLE = (0.05, 0.10, 0.20, 0.30, 0.50)
X_GATE = (0.10, 0.20, 0.30, 0.50, 0.70, 0.80)
#: x points per gate scan.  Small on purpose: see `gate_report`'s cost note.
#: design section 5.4's own table is the 300-point one and lives in the doctest.
N_X = 100


def x_grid(lo=0.01, hi=1.0, n=300):
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


# --------------------------------------------------------------- the numbers

def four_terms(model=None):
    """Design section 7: the four terms of x*b1 per nucleon, 6Li, Q2 = 2.5."""
    b = lg.Li6ConvolutionB1() if model is None else model
    miller = lg.Li6B1(lg.MillerB1())
    cdks = lg.Li6B1(lg.CdksB1())
    rows = []
    for x in X_TABLE:
        t1 = x * b.b1_embedded_s(x, Q2)
        t3 = x * b.b1_cg_dwave(x, Q2)
        t2d = x * b.b1_alpha_d_dwave_d(x, Q2)
        t2a = x * b.b1_alpha_d_dwave_alpha(x, Q2)
        rows.append((x, t1, t3, t2d, t2a, x * b.b1(x, Q2, 0.0),
                     x * miller.b1(x, Q2, 0.0), x * cdks.b1(x, Q2, 0.0)))
    return rows


def bands(b):
    """The band rows and the term-(2) knob rows -- never quote one alone."""
    out = {"band": [], "w_alpha_d": []}
    for scale in (0.0, 1.0, 2.0):
        out["band"].append((scale, 0.30 * b.banded(scale).b1(0.30, Q2, 0.0)))
    for w in (0.0, 1.0, 2.0):
        o = lg.Li6ConvolutionOptions()
        o.w_alpha_d_dwave = w
        m = lg.Li6ConvolutionB1(o)
        out["w_alpha_d"].append(
            (w, [x * m.b1(x, Q2, 0.0) for x in (0.10, 0.30, 0.50)]))
    return out


def densities(b):
    d, a = b.densities(), b.densities_alpha()
    return {
        "norm_d": d.norm(), "norm_a": a.norm(),
        "mean_y_d": d.mean_y(), "mean_y_a": a.mean_y(),
        "p_d": d.p_d(), "p_d_momentum": d.p_d_momentum(),
        "vmc_p_d": lg.VMC_P_D_LI6, "n_alpha_d": lg.VMC_N_ALPHA_D_LI6,
        "y_max_d": d.kinematics.y_max(), "y_max_a": a.kinematics.y_max(),
    }


# ------------------------------------------------------------------ the gate

def gate_report(unpol=None, finite_q=True, label="", n_x=N_X):
    """One row of the gate table: landmarks + the G3b peak ratio.

    THE GRID IS THIS SCRIPT'S OWN and it is NOT the doctest's: the x points
    are the union of `x_grid(0.01, 1.0, n_x)` and `x_grid(0.10, 0.80, n_x)`
    (199 points at the default --n-x 100, denser inside the G3b window),
    while `tests/test_b1_nuclear.cpp` scans `linspace(0.01, 1.59, 300)`.  The
    two therefore disagree in the 4th digit of a landmark, which is a grid
    difference and nothing else -- label any table you copy out of here.

    COST NOTE.  `finite_q_delta = true` -- the DEFAULT of
    `DeuteronConvolutionB1` since 2026-09-03 -- rebuilds the light-cone
    density set at every x (kappa = kappa(x, Q2) is inside the
    delta-function), which is ~0.24 s per point against ~0.7 ms for the
    kappa = 1 path -- a factor 300.
    So x*b1 is evaluated ONCE per grid point here and memoized; `b1_landmarks`
    and the [0.10, 0.80] peak scan of G3b then share the same evaluations, and
    `--n-x` is small by default.  The doctest (T1) is the one that runs the
    300-point grid design section 5.4 quotes.
    """
    o = lg.DeuteronConvolutionB1.Options()
    if unpol is not None:
        o.unpol = unpol
    o.finite_q_delta = finite_q
    gate = lg.DeuteronConvolutionB1(o)

    cache = {}

    def xb1(x):
        if x not in cache:
            cache[x] = x * gate.b1(x, Q2, 0.0)
        return cache[x]

    grid = sorted(set(x_grid(0.01, 1.0, n_x)) | set(x_grid(0.10, 0.80, n_x)))
    lm = lg.b1_landmarks(xb1, grid)
    table = lg.b1_landmarks_of_table()
    peak = max(abs(xb1(x)) for x in x_grid(0.10, 0.80, n_x))
    return {
        "label": label or ("ToyF2" if unpol is None else "custom"),
        "finite_q": finite_q,
        "zeros": list(lm.zeros), "zero_slope": list(lm.zero_slope),
        "min": (lm.x_min, lm.xb1_min), "max": (lm.x_max, lm.xb1_max),
        "integral": lm.integral_b1,
        "peak_window": peak,
        "g3b_ratio": peak / abs(table.xb1_max),
        "table": {"zeros": list(table.zeros),
                  "min": (table.x_min, table.xb1_min),
                  "max": (table.x_max, table.xb1_max),
                  "integral": table.integral_b1},
    }


def checklist_item4(n_x=N_X):
    """The DOMINANT gate item: ToyF2 against a real global fit.

    Falls back to the two ToyF2 rows when this build has no LHAPDF -- the
    measurement is then simply unavailable, which is a fact about the build
    and not a failure.
    """
    rows = [gate_report(None, True,
                        "ToyF2, kappa = sqrt(1+gamma^2)  [Eq. 21, the DEFAULT]",
                        n_x),
            gate_report(None, False, "ToyF2, kappa = 1  [Eq. 17]", n_x)]
    if not lg.HAVE_LHAPDF:
        return rows
    lg.lhapdf_quiet()
    ct18 = lg.LhapdfSF("CT18NLO")
    rows += [gate_report(ct18, True, "CT18NLO, kappa  [Eq. 21]", n_x),
             gate_report(ct18, False, "CT18NLO, kappa = 1", n_x)]
    return rows


def sign_gate():
    w0, w2 = lg.li6_alpha_d_partial_waves()
    hbar = 0.19733
    return {
        "q_alpha_d_fm2": lg.alpha_d_quadrupole_fm2(w0, w2),
        "phi0phi2": [(k, w0(k * hbar) * w2(k * hbar))
                     for k in (0.2, 1.0, 3.0)],
    }


# ----------------------------------------------------------------- reporting

def emit(fh, n_x=N_X):
    b = lg.Li6ConvolutionB1()
    w = fh.write
    w("### section 7 -- the four terms, x*b1 per nucleon, 6Li, Q2 = %g\n\n" % Q2)
    w("| x | (1) | (3) | (2d) | (2a) | total | Li6B1(MillerB1) | Li6B1(CdksB1) |\n")
    w("|---|---|---|---|---|---|---|---|\n")
    for r in four_terms(b):
        w("| %.2f | %+.4e | %+.4e | %+.4e | %+.4e | **%+.4e** | %+.4e | %+.4e |\n" % r)
    w("\nRATIOS: (3)/(1) and [(2d)+(2a)]/(1), and (2a)/(2d) against the\n"
      "analytic 2 (M_d/M_alpha)^2 = 0.5064:\n\n")
    for x, t1, t3, t2d, t2a, tot, _m, _c in four_terms(b):
        w("  x = %.2f   (3)/(1) = %+.4f %%   [(2d)+(2a)]/(1) = %+.4f   "
          "(2a)/(2d) = %.4f\n"
          % (x, 100.0 * t3 / t1, (t2d + t2a) / t1, t2a / t2d))

    bd = bands(b)
    w("\n### the MANDATORY 100 %% band (x*b1 at x = 0.30) -- never quote one row\n\n")
    for scale, v in bd["band"]:
        w("  --b1-band-scale %g   %+.4e\n" % (scale, v))
    w("\n### the term-(2) knob, w_alpha_d_dwave (x*b1 at x = 0.10/0.30/0.50)\n\n")
    for wv, vals in bd["w_alpha_d"]:
        w("  w = %g   " % wv + "   ".join("%+.4e" % v for v in vals) + "\n")

    d = densities(b)
    w("\n### densities\n\n")
    for k in sorted(d):
        w("  %-14s %.7f\n" % (k, d[k]))
    w("\n  int b1 dx [0.01, 1.2] = %+.5e\n" % b.close_kumano_integral())

    w("\n### the A = 2 gate -- layer 2, the alpha-d SIGN gate\n\n")
    sg = sign_gate()
    w("  Q(alpha-d relative motion) = %+.6f fm^2   MUST BE < 0\n"
      % sg["q_alpha_d_fm2"])
    for k, v in sg["phi0phi2"]:
        w("  phi0*phi2 at k = %.1f fm^-1: %+.4e\n" % (k, v))

    w("\n### the A = 2 gate -- layer 3, CDKS Fig. 4 (checklist item 4 included)\n\n")
    w("Grid: THIS SCRIPT's own, the union of x_grid(0.01, 1.0, %d) and\n"
      "x_grid(0.10, 0.80, %d) = %d points -- NOT the doctest's\n"
      "linspace(0.01, 1.59, 300).  Landmarks differ in the 4th digit between\n"
      "the two for that reason alone.\n\n"
      % (n_x, n_x, len(sorted(set(x_grid(0.01, 1.0, n_x))
                              | set(x_grid(0.10, 0.80, n_x)))))) 
    reports = checklist_item4(n_x)
    if not lg.HAVE_LHAPDF:
        w("  (this build has no LHAPDF: checklist item 4, the DOMINANT term,\n"
          "   could not be run -- see phase_D_gate.md for the measured values)\n\n")
    w("| configuration | zeros | minimum | maximum | G3b ratio |\n")
    w("|---|---|---|---|---|\n")
    for r in reports:
        w("| %s | %s | %+.4e @ %.4f | %+.4e @ %.4f | **%.4f** |\n"
          % (r["label"], ", ".join("%.5f" % z for z in r["zeros"]),
             r["min"][1], r["min"][0], r["max"][1], r["max"][0],
             r["g3b_ratio"]))
    t = reports[0]["table"]
    w("| *digitized CDKS (b1_landmarks_of_table)* | %s | %+.4e @ %.4f | "
      "%+.4e @ %.4f | *1* |\n"
      % (", ".join("%.5f" % z for z in t["zeros"]), t["min"][1], t["min"][0],
         t["max"][1], t["max"][0]))
    w("\nG3b needs the ratio inside a factor 2, i.e. in [0.5, 2].  At the "
      "DEFAULT it is %.4f\n-- a factor %.2f LOW -- so THE GATE FAILS and no "
      "6Li number may be published\n(design section 5.4, Escalation).  The "
      "ratio is quoted BOTH ways throughout the phase D\ndocuments: 'ratio "
      "R' and 'a factor 1/R low' are the same statement.\n"
      % (reports[0]["g3b_ratio"], 1.0 / reports[0]["g3b_ratio"]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--write", action="store_true",
                    help="also refresh the docs/ reports (writes nothing "
                         "outside docs/)")
    ap.add_argument("--n-x", type=int, default=N_X,
                    help="x points per gate scan (default %d).  The "
                         "finite_q_delta rows rebuild both density sets at "
                         "every point, so this is ~0.14 s per point there"
                         % N_X)
    args = ap.parse_args(argv)
    emit(sys.stdout, args.n_x)
    if args.write:
        # The two documents are hand-written prose around these numbers, so
        # this appends a machine-regenerated appendix rather than overwriting
        # an author's argument with a table.
        out = DOCS / "phase_D_regenerated.md"
        with open(out, "w") as fh:
            fh.write("# Phase D -- regenerated by validation/b1_li6_table.py\n\n"
                     "Machine output.  The prose lives in `phase_D_gate.md` "
                     "and `phase_D_numbers.md`;\nthis file exists so a "
                     "reviewer can diff the numbers against them.\n\n")
            emit(fh, args.n_x)
        print("\nwrote %s" % out, file=sys.stderr)


if __name__ == "__main__":
    main()
