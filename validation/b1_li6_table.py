#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regenerate design D's number tables and its A = 2 validation-gate report.

Design: `docs/open_items/run_2026-09-02/design_D_b1_li6.md` (sections 5 and 7).
Outputs (with `--write`, and it writes NOTHING outside `docs/`):
    docs/open_items/run_2026-09-02/phase_D_regenerated.md
It APPENDS nothing to `phase_D_gate.md` / `phase_D_numbers.md` and overwrites
neither: those two are hand-written prose around these numbers, and the
regenerated file exists so a reviewer can diff the numbers against them.

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

READ THE GATE BEFORE QUOTING ANY 6Li NUMBER.  G3b (magnitude) is read off
CDKS's OWN nucleon PDF, MSTW2008 LO, at CDKS Eq. (21)'s delta-function:
design section 5.4's Escalation clause reads the peak ratio AFTER the
checklist, and checklist item 4 IS the nucleon PDF.  Measured 2026-09-03
(docs/open_items/run_2026-09-03/phase_A_numbers.md): 0.843, a factor 1.19
low, INSIDE the factor of 2 -- G3b passes there.  It does NOT pass on the
SHIPPED DEFAULT unpolarised backend, the library ToyF2 (0.440, a factor 2.27
low, OUTSIDE the [0.5, 2] window), and it clears by only 4 % at Eq. (17)'s
kappa = 1 (0.520).  So the lift is about a CONFIGURATION: quote numbers made
with --b1-unpol mstw (PipelineConfig.b1_unpol = B1UnpolSource.Mstw); the
default toy backend's numbers are not covered by it.  With the real CD-Bonn
wave function as well (checklist item 5) the ratio is 1.000338 -- read that as
"the residual is below the error of digitizing a published figure", NEVER as
three-digit agreement with CDKS, and note it is specific to CD-Bonn AND MSTW
together and is not the default.

THE 6Li PUBLICATION BAN WAS LIFTED ON 2026-09-03 on the strength of that,
design section 5.4's Escalation clause not being triggered
(docs/OPEN_ITEMS_SOLUTIONS.md section 10).  What did NOT lift is the mandatory
+-100 % band: it comes from Q(6Li) = -0.0818(17) fm^2 against
Q_d = +0.2859(3) fm^2, not from the gate, and the gate is A = 2 and tests
nothing about the alpha-d step.  Quote the band and the configuration, never
one row.

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
    are the union of `x_grid(0.001, 1.0, n_x)` and `x_grid(0.10, 0.80, n_x)`
    (200 points at the default --n-x 100, denser inside the G3b window),
    while `tests/test_b1_nuclear.cpp` scans `linspace(0.001, 1.59, 300)`.  The
    two therefore disagree in the 4th digit of a landmark, which is a grid
    difference and nothing else -- label any table you copy out of here.

    THE SCAN FLOOR IS 0.001, not 0.01, since 2026-09-03: a sign change below
    the floor is invisible to `b1_landmarks` no matter what G3a's counting
    window says, and with CT18NLO the low-x zero is at 0.0098.  Design
    section 5.4 now states the floor as part of G3a.  `integral` is the one
    quantity that must NOT follow it -- see below.

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

    grid = sorted(set(x_grid(0.001, 1.0, n_x)) | set(x_grid(0.10, 0.80, n_x)))
    lm = lg.b1_landmarks(xb1, grid)
    # G3c on the x >= 0.01 SUB-GRID (design section 5.4's G3c, amended
    # 2026-09-03).  `integral_b1` is a trapezoid of b1 = (x*b1)/x, so running
    # it from the 0.001 scan floor moves it ~+20 % on 1/x weighting alone,
    # while the reference `close_kumano_integral(true)` is over the digitized
    # table's own [0.0100, 1.590].  The evaluations are memoized above, so
    # this second landmark pass costs nothing.
    lm_g3c = lg.b1_landmarks(xb1, [x for x in grid if x >= 0.01])
    table = lg.b1_landmarks_of_table()
    peak = max(abs(xb1(x)) for x in x_grid(0.10, 0.80, n_x))
    return {
        "label": label or ("ToyF2" if unpol is None else "custom"),
        "finite_q": finite_q,
        "zeros": list(lm.zeros), "zero_slope": list(lm.zero_slope),
        "min": (lm.x_min, lm.xb1_min), "max": (lm.x_max, lm.xb1_max),
        "integral": lm_g3c.integral_b1,          # G3c, x >= 0.01
        "integral_full_scan": lm.integral_b1,    # NOT G3c -- see above
        "peak_window": peak,
        "g3b_ratio": peak / abs(table.xb1_max),
        "table": {"zeros": list(table.zeros),
                  "min": (table.x_min, table.xb1_min),
                  "max": (table.x_max, table.xb1_max),
                  "integral": table.integral_b1},
    }


def mstw_sf():
    """CDKS's own nucleon PDF -- MSTW2008 LO -- or None if this build lacks it.

    CDKS build the Fig. 4 curve this gate compares against on MSTW2008 LO, so
    this is the LIKE-FOR-LIKE version of checklist item 4; the CT18NLO row is
    a modern cross-check, not the same measurement.  `MstwSF` wraps
    `Pythia8::MSTWpdf`, so it rides the PYTHIA tier and NOT the LHAPDF one and
    has to be guarded independently -- a build with PYTHIA and no LHAPDF still
    gets these rows.  The pdfdata grids are a separate install artefact from
    the PYTHIA library itself, so construction is attempted and a missing grid
    is reported as an absent row rather than as an exception.
    """
    if not getattr(lg, "HAVE_PYTHIA8", False) or not hasattr(lg, "MstwSF"):
        return None
    try:
        return lg.MstwSF()
    except RuntimeError:
        return None


def checklist_item4(n_x=N_X):
    """The DOMINANT gate item: ToyF2 against the real global fits.

    Three PDF families, each at both delta-functions: the library ToyF2, the
    modern CT18NLO, and MSTW2008 LO -- the PDF CDKS themselves used, which is
    why its row is the one G3b is finally read off.  A family whose structure
    function this build cannot supply is skipped; that is a fact about the
    build and not a failure.
    """
    rows = [gate_report(None, True,
                        "ToyF2, kappa = sqrt(1+gamma^2)  [Eq. 21, the DEFAULT]",
                        n_x),
            gate_report(None, False, "ToyF2, kappa = 1  [Eq. 17]", n_x)]
    if lg.HAVE_LHAPDF:
        lg.lhapdf_quiet()
        ct18 = lg.LhapdfSF("CT18NLO")
        rows += [gate_report(ct18, True, "CT18NLO, kappa  [Eq. 21]", n_x),
                 gate_report(ct18, False, "CT18NLO, kappa = 1", n_x)]
    mstw = mstw_sf()
    if mstw is not None:
        rows += [gate_report(mstw, True,
                             "MSTW2008 LO, kappa  [Eq. 21]  <- CDKS's own PDF",
                             n_x),
                 gate_report(mstw, False, "MSTW2008 LO, kappa = 1", n_x)]
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

def gate_footer(verdict_ratio, toy_ratio):
    """The report's closing paragraph -- CONDITIONAL on the verdict row PASSING.

    It used to announce the ban lift UNCONDITIONALLY, including in the branch
    that had just printed "MSTW2008 LO is not available in this build, so
    G3b's VERDICT row was not measured here".  A report that states a
    conclusion its own body does not contain is exactly the failure the gate
    documents warn about, and this one is printed INTO every regenerated
    report, where it outlives the run that made it.

    Keying it on "was the row measured?" was still not enough: a MEASURED but
    FAILING MSTW row would print "OUTSIDE [0.5, 2]: G3b FAILS" in the body and
    then announce the lift in the footer of the same report.  It is keyed on
    the ratio itself, so the footer cannot contradict the table above it.

    A pure function of (the verdict row's ratio or None, the ToyF2 ratio), so
    all three branches are testable without running the 200-point scan.
    """
    if verdict_ratio is not None and not (0.5 <= verdict_ratio <= 2.0):
        head = ("\nTHE VERDICT ROW WAS MEASURED AND IS OUTSIDE [0.5, 2]: G3b "
                "FAILS at %.6f.\nThis report does NOT support the 2026-09-03 "
                "ban lift -- it contradicts it.  The\nlift on record in "
                "docs/OPEN_ITEMS_SOLUTIONS.md section 10 rests on a measured "
                "0.843243\nat CDKS Eq. (21); if this run disagrees, the "
                "configurations differ and the\ndifference has to be found "
                "before anything here is quoted.\n" % verdict_ratio)
    elif verdict_ratio is None:
        head = ("\nTHIS REPORT DID NOT MEASURE THE VERDICT ROW, so it does "
                "not establish the\n2026-09-03 6Li publication ban lift and "
                "must not be read as doing so.  The lift\nis on record in "
                "docs/OPEN_ITEMS_SOLUTIONS.md section 10 and rests on the "
                "MSTW2008 LO\nrow at CDKS Eq. (21), which is MISSING above.  "
                "Every 6Li row in this report was\nmade with the ToyF2 "
                "default, whose own G3b is %.4f -- OUTSIDE [0.5, 2] -- so "
                "none\nof them is covered by the lift.  Rerun in a build "
                "with MSTW2008 LO on disk\nbefore quoting anything here.\n"
                % toy_ratio)
    else:
        head = ("\nTHE 6Li PUBLICATION BAN WAS LIFTED ON 2026-09-03 (design "
                "section 5.4's Escalation\nclause is not triggered after the "
                "checklist; docs/OPEN_ITEMS_SOLUTIONS.md section 10),\nfor "
                "the MSTW2008 LO row above and NOT for the ToyF2 default row "
                "(%.4f), which is\noutside the window.  Quote the "
                "configuration with the number.\n" % toy_ratio)
    return head + (
        "\nWhat did NOT lift, on either branch, is the MANDATORY +-100 % "
        "BAND -- it comes from\nQ(6Li) against Q_d, not from this gate, and "
        "this gate is A = 2 and tests nothing\nabout the alpha-d step.  "
        "Quote the band rows above and the configuration that made\nthem, "
        "never one row.\n")


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
    w("\n### the MANDATORY 100 % band (x*b1 at x = 0.30) -- never quote one row\n\n")
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
    w("Grid: THIS SCRIPT's own, the union of x_grid(0.001, 1.0, %d) and\n"
      "x_grid(0.10, 0.80, %d) = %d points -- NOT the doctest's\n"
      "linspace(0.001, 1.59, 300).  Landmarks differ in the 4th digit between\n"
      "the two for that reason alone.  The int-b1-dx column is over\n"
      "x in [0.01, 1.0] -- the x >= 0.01 sub-grid of THIS grid, which stops at\n"
      "1.0.  G3c proper is over the digitized range [0.0100, 1.590] and lives\n"
      "in the doctest; the two are NOT the same number.  Neither includes the\n"
      "0.001 scan floor: the trapezoid is of b1 = (x*b1)/x and the extension\n"
      "would move it ~20 %% on 1/x weighting alone.\n\n"
      % (n_x, n_x, len(sorted(set(x_grid(0.001, 1.0, n_x))
                              | set(x_grid(0.10, 0.80, n_x))))))
    reports = checklist_item4(n_x)
    if not lg.HAVE_LHAPDF:
        w("  (this build has no LHAPDF: the CT18NLO cross-check of checklist\n"
          "   item 4 could not be run -- see phase_D_gate.md for its values)\n\n")
    if not any(r["label"].startswith("MSTW") for r in reports):
        w("  (this build has no MSTW2008 LO grid: checklist item 4's\n"
          "   LIKE-FOR-LIKE row -- CDKS's own PDF -- could not be run, and it\n"
          "   is the row G3b is read off; see phase_A_numbers.md)\n\n")
    w("| configuration | zeros (G3a) | minimum | maximum | int b1 dx "
      "(x in [0.01, 1.0]) | G3b ratio |\n")
    w("|---|---|---|---|---|---|\n")
    for r in reports:
        w("| %s | %s | %+.4e @ %.4f | %+.4e @ %.4f | %+.4e | **%.4f** |\n"
          % (r["label"], ", ".join("%.5f" % z for z in r["zeros"]),
             r["min"][1], r["min"][0], r["max"][1], r["max"][0],
             r["integral"], r["g3b_ratio"]))
    t = reports[0]["table"]
    w("| *digitized CDKS (b1_landmarks_of_table)* | %s | %+.4e @ %.4f | "
      "%+.4e @ %.4f | %+.4e | *1* |\n"
      % (", ".join("%.5f" % z for z in t["zeros"]), t["min"][1], t["min"][0],
         t["max"][1], t["max"][0], t["integral"]))
    # G3b's verdict is READ OFF the MSTW2008 LO row at Eq. (21): design
    # section 5.4's Escalation clause reads the peak ratio AFTER the
    # checklist, checklist item 4 is the nucleon PDF, and CDKS's own PDF is
    # MSTW2008 LO.  The ToyF2 row is the library toy and is NOT the verdict.
    verdict = next((r for r in reports
                    if r["label"].startswith("MSTW2008 LO, kappa  [Eq. 21]")),
                   None)
    toy = reports[0]["g3b_ratio"]
    w("\nG3b needs the ratio inside a factor 2, i.e. in [0.5, 2].  The row it "
      "is READ OFF\nis CDKS's own PDF at CDKS's own Eq. (21) delta-function: "
      "design section 5.4's\nEscalation clause reads the peak ratio AFTER the "
      "checklist, and checklist item 4\nIS the nucleon PDF.\n\n")
    if verdict is None:
        w("  MSTW2008 LO is not available in this build, so G3b's VERDICT row "
          "was not\n  measured here.  The ToyF2 default row is %.4f (a factor "
          "%.2f low) and is\n  the library toy, NOT the verdict.\n"
          % (toy, 1.0 / toy))
    else:
        g = verdict["g3b_ratio"]
        w("  MSTW2008 LO, kappa [Eq. 21]   ratio = %.4f   factor %.2f low   "
          "-> %s\n" % (g, 1.0 / g, "INSIDE [0.5, 2]: G3b PASSES"
                        if 0.5 <= g <= 2.0 else "OUTSIDE [0.5, 2]: G3b FAILS"))
        w("  ToyF2, kappa [Eq. 21]         ratio = %.4f   factor %.2f low   "
          "-- the library\n  toy, kept so the PDF dependence is visible; NOT "
          "the verdict.\n" % (toy, 1.0 / toy))
    w("\nThe ratio is quoted BOTH ways throughout the phase D documents: "
      "'ratio R' and\n'a factor 1/R low' are the same statement.\n")
    w(gate_footer(verdict["g3b_ratio"] if verdict is not None else None, toy))


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
