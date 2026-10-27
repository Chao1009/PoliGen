#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""docs/benchmarking/04_theory.md T-20 -- Wiringa & Schiavilla's 6Li elastic
form factors as a THEORY anchor for the tree's C0 zero (tier T3).

New 2026-10-26 (run 2026-10-26 phase C4, PLAN.md "C4"; under
docs/benchmarking/08_channel_priority.md sec. 4 ingredient F -- the 6Li
elastic form factor's one consumer is the elastic radiative tail of the
inclusive 6Li channels, ranks 1 and 4).  The DATA level of the same form
factor stays where it is: row 5 (t3_li6_charge_ff_fb.py, the UVa FB density
and Li et al. 1971's minimum) and row 6 (t3_li_magnetization_rfy.py, RFY's
F_M points).  This row reads a CALCULATION.

THE DOCUMENT.  R. B. Wiringa and R. Schiavilla, "Microscopic calculation of
6Li elastic and transition form factors", Phys. Rev. Lett. 81 (1998) 4317,
arXiv:nucl-th/9807037v2 -- AV18 + Urbana IX variational Monte Carlo, charge
and current operators in impulse approximation (IA) and with two-body terms
(IA+MEC).  Fig. 1 (PDF page 10) is a VECTOR drawing: its lines, symbols,
error bars and data points are drawing operators, read back from the PDF
(sha256 fc51cb93...91d20e8, the author's corpus copy, NOT in this
repository) into data/ws1998_li6_formfactors.csv -- 503 record rows (every
line vertex, clip point, symbol and bar of both panels) and 68 SLOT rows in
data/ff/li6_elastic.csv's columns (q, F_C0, F_C2, F_M1; IA+MEC; promotion =
one copy, a registry decision).  The file's provenance block gives the frame
boxes, the axis calibration and its check, the device resolution, the
reading of every column and the commands; `--redigitize PDF` re-reads the
PDF (dev time, PyMuPDF) and compares.  The sha256 of every byte of the file
is DATA_SHA256 below, re-checked on every read; q and F2 of every row are
re-derived from its PDF coordinates and every slot row from the record rows.

WHAT THE FIGURE LOCATES (a fact of the reference, stated before the design).
The left panel plots F_L^2, F_C0^2 and F_C2^2 on a log axis that stops at
1e-7.  The C0 lines are drawn from (q, F^2) = (0, 1) DOWN TO THAT FLOOR and
not beyond: the IA+MEC line's last computed vertex is at q ~ 2.92 fm^-1, the
IA line's at ~ 3.20 fm^-1, and no C0 line re-enters the frame up to 4 fm^-1
(the C2 lines, drawn to 4 fm^-1, dip under the floor at 2.54 .. 2.68 fm^-1
and come back: their zero).  So the figure does NOT draw the C0 zero; it
shows F_C0 != 0 below the last vertex.

THE TREE SIDE.  The SHIPPED HoSpin1FF.for_ion(li6()) (C0 shape `ho`,
fold_nucleon = True, nothing moved) and its other C0 edge, C0Shape.VmcFt
(`vmc-ft`, the j0 transform of the VMC point-proton density, continued
exponentially above 3 fm^-1 with NO zero).  Read as WS98 plot them, per
nucleus: F_C0 = fc(t)/Z, F_C2 = (2 sqrt 2/3) eta fq(t)/Z (rc.hpp's multipole
conversion, eta = t/(4 M_A^2)), F_T^2 = (4/3) eta (1 + eta) fm(t)^2/Z^2,
t = (q hbar c)^2, M_A = lg.li6().mass(), Z = 3.

WHAT IS COMPARED, AND THE TOLERANCE -- fixed in this header on 2026-10-27
(UTC) before this harness computed any comparison (PLAN.md ground rule 4).
The block between the two marker lines is hashed (DESIGN_BLOCK_SHA256
below); `design_block_sha256()` re-checks it, and the numerical rules
against its sentences, on every run; the harness exits 2 if either moved.
==== DESIGN BLOCK BEGIN ====
THE VERDICT, one per C0 edge; the row's status is the SHIPPED edge's (`ho`),
`vmc-ft` is a sub-row by the same rule.  Reference: the IA+MEC C0 line, WS98's
full calculation (the one the paper says reproduces the F_L data); its
computed vertices are the vendored rows (curve C0, calc IA+MEC, kind vertex).
q_last is the largest q among them.  Each has F_C0^2 >= 1e-7 (drawn above the
frame floor), so F_C0 != 0 on [0, q_last]; the line meets the floor after
q_last and no C0 line is drawn above it up to q = 4 fm^-1.  The figure bounds
the first C0 zero from BELOW only: q0_WS > q_last.
The planned rule (PLAN.md phase C4 and the task: the tree's C0 zero within
the VMC curve's zero +- (digitization error + 0.1 fm^-1)) is applied to the
SET of zeros the figure allows, q0_WS in (q_last, infinity): the distance of
q0_tree from that set must not exceed delta + 0.1 fm^-1, i.e.
    PASS iff q0_tree >= q_last - (delta + 0.1 fm^-1),
with delta = one device unit of the drawing in q = 0.0047619 fm^-1 (the
vendored dq_fm).  q0_tree is the first sign change of the edge's fc(t),
t = (q hbar c)^2, on (0, 5] fm^-1 (steps of 0.01 fm^-1, bisected to 1e-12 --
row 5's first_zero).  An edge with NO sign change below 5 fm^-1 has no zero
where the VMC C0 is drawn non-zero, so it meets the rule; it is recorded as
"no zero below 5 fm^-1" and never given a distance.  THE UPPER SIDE IS NOT
TESTED: the figure does not locate the zero above q_last, so a PASS means
"not contradicted by this calculation's figure", no more.
RECORDED, NEVER IN THE VERDICT -- each rule fixed here, before computing:
 (r1) where the zero might be, as information: the IA+MEC drawn branch as the
      signed F_C0 = +sqrt(F2), continued beyond q_last by four fits in F (not
      log F) to its last computed vertices: (e1) the straight line through the
      last two; (e2) the least-squares line through the last five; (e3) the
      least-squares quadratic through the last five; (e4) the least-squares
      quadratic through the last ten.  Each gives its smallest root above
      q_last, or "no root below 4 fm^-1".  Recorded with their spread; and the
      IA line's q_last, with the verdict's rule applied to it, for both edges.
 (r2) the C0 ceiling: above q_last the IA+MEC C0 is under the frame floor,
      F_C0^2 < 1e-7 up to 4 fm^-1, IF the line was computed that far (an
      inference: the C2 lines run to 4 fm^-1 on the same 0.04 fm^-1 grid, the
      C0 line is not drawn there).  Each edge's maximum (fc/Z)^2 on
      (q_last, 4] fm^-1 is recorded against 1e-7.
 (r3) the slope radius: weighted least squares of ln F2 = a1 q^2 + a2 q^4 (no
      constant: F_C0(0) = 1 is the frame corner) over the computed C0
      vertices with 0 < q <= 1.0 fm^-1, sigma = one device unit in ln F^2
      (0.010440) per vertex; r = sqrt(-3 a1).  Model spread: q_max in
      {0.6, 1.0, 1.4} fm^-1 times the orders (a1, a2) and (a1, a2, a3 q^6).
      The SAME fit of each edge's (fc/Z)^2 at the same q is recorded beside it
      (like for like: the truncation bias is common), and each edge's slope
      radius by Richardson extrapolation (as t3_li6_rms_devries.py).  For
      IA+MEC and IA.
 (r4) the C0 shape with ONE free scale: chi2 = sum over the computed C0
      vertices with 0 < q <= q_last of (ln F2_i - ln s - ln m_i)^2/sigma^2,
      m_i = (fc(t_i)/Z)^2, sigma = 0.010440; ln s in closed form (the mean log
      ratio); chi2 and ndf = N - 1, s, and F2_i/(s m_i) at the vertices
      nearest q = 1.0, 2.0, 2.5 fm^-1 and at q_last.  Both edges, IA+MEC and
      IA.  No tolerance: the digitization error is not a theory error.
 (r5) F_M1 against fm: the 20 IA+MEC F_T^2 symbols, sigma_i = the upper
      half-bar and the digitization error (1.05 % of F2) in quadrature,
      against s (4/3) eta_i (1 + eta_i) fm(t_i)^2/Z^2 (the shipped F_mag, q_z =
      1.30 fm^-1, b = 1.85 fm; one F_m on both C0 edges), s free in closed
      form: chi2, ndf = 19, p; and s = 1 (each side normalised to its own
      moment, the tree's measured 0.822 and WS98's 0.832 mu_N): chi2, ndf =
      20.  The F_T zero: the symbol of least F_T^2 with 0.5 < q < 2.0 fm^-1
      and its two neighbours, against the tree's q_z = 1.30 fm^-1.
 (r6) F_C2 against fq: the tree's F_C2 = (2 sqrt 2/3) eta fq(t)/Z at every
      computed IA+MEC C2 vertex; the ratio |F_C2^WS|/|F_C2^tree| at the
      vertices with q <= 0.5 fm^-1 (mean, min, max) -- the low-q end, set by
      WS98's VMC quadrupole and by the tree's MEASURED Q = -0.0818 fm^2; WS98's
      C2 zero (the last vertex before and the first after the floor gap)
      against the tree's (its fq shares the C0 zero); the second-lobe maximum.
 (r7) closures of the READING against numbers the paper prints (they flag the
      normalisation read into F2 and the slot columns; they move no verdict):
      mu: the two lowest IA+MEC F_T^2 symbols as mu_i^2 = Z^2 F2_i/((4/3)
        eta_i (1 + eta_i) (M_A/m_p)^2), extrapolated linearly in q^2 to q = 0:
        HOLDS iff |mu/0.832 - 1| <= 0.02 (the IA symbols against 0.829, the
        same rule).
      Q: the three lowest computed IA+MEC C2 vertices as c_i = 3 sqrt(2) Z
        sqrt(F2_i)/q_i^2 (the Rosenbluth reading F_C2 -> q^2 Q/(3 sqrt 2 Z)),
        extrapolated linearly in q^2 to q = 0: HOLDS iff |c_0/0.23 - 1| <=
        0.05 (the printed -0.23(9) fm^2's central value).  One that does not
        hold says the figure's C2 is not the Rosenbluth F_C2/Z of the printed
        Q, and the slot's F_C2 column inherits whatever factor it is.
      lines vs symbols: at each IA+MEC FL2 symbol with 0.25 <= q <= 2.75
        fm^-1 bracketed by computed vertices of both lines, C0^2 + C2^2 (each
        interpolated linearly in ln q between its bracketing vertices) against
        the symbol's F2: HOLDS iff within 2 sigma_i, sigma_i = the upper
        half-bar and 2.1 % of F2 in quadrature (three readings of one device
        unit each, rounded up).
      RFY: the seven right-panel Stanford points against RFY Table I
        (data/rfy1966_li_magnetic_ff.csv, row 6's file, its sha256 re-checked
        here): the q's HOLD iff |q_dig - sqrt(q2_printed)| <= 0.0095 fm^-1
        (two device units) at every point; the ratio F2_dig/[F_M^2 (4/3) eta
        (1 + eta) (M_A/m_p)^2 0.822^2/Z^2] per point, its mean and max/min
        recorded -- a CONSTANT ratio is a normalisation difference between
        RFY's F_M^2 and WS98's plotted F_T^2, not a digitization error; the
        ratio is constant iff max/min <= 1.042 (four device units).
WHAT WAS SEEN BEFORE THIS BLOCK WAS WRITTEN (by the worker, 2026-10-26/27):
the page rendered at 150 dpi; every drawing operator of page 10, the vertex
coordinates of the C0 and C2 lines among them, and by hand from them: the
IA+MEC C0 line's last computed vertex at q = 2.919 fm^-1 (F^2 = 1.34e-7), its
floor crossing at 2.933 and the IA line's at 3.214, the C2 floor gap
2.543 .. 2.676 fm^-1, the signed F of the last six IA+MEC C0 vertices, a
straight line through the last two (root ~ 3.02 fm^-1) and a quadratic by
finite differences that turns up before reaching zero, F_C0^2(2.5) ~ 2.6e-5;
on the tree side, by hand, the ho edge's (fc/Z)^2 ~ 4.7e-10 at q = 2.93 and
~ 4.9e-7 at 2.5 fm^-1 (about 50x under WS98's IA+MEC there), besides the
committed q0 = 3.0998 fm^-1 and r_ch = 2.5710 fm and rc.hpp's statement that
the vmc-ft edge has no C0 zero and is x6.15 the ho edge at q = 2.5 fm^-1; on
the reference side, by hand, the C2 line's low-q coefficient read as
|Q| ~ 0.48-0.49 fm^2 through the Rosenbluth convention (so the (r7) Q closure
is expected NOT to hold; a factor 3/sqrt(2) in the convention would reconcile
it), the q = 0.05 F_T^2 symbol read as mu ~ 0.84, and the first Stanford F_T^2
point against RFY Table I row 1 read as a ratio ~ 0.71; the digitized symbols
with their bars; and the tool's C2 extraction agreeing with the vendored rows
to 1.1e-4.  So the verdict is foreseeable: q0_tree = 3.0998 against
q_last - 0.1048 ~ 2.814 fm^-1.  The tolerance is the plan's and the task's,
fixed before any of this; what the worker chose after seeing the figure is
how the rule reads a zero the figure does not locate (the one-sided form
above) and the recorded rules (r1)-(r7).  No rule here was changed after a
number it governs was computed by this harness.
==== DESIGN BLOCK END ====

POST-HOC ADDENDUM (2026-10-27, written after the first run of 00:25:46Z and
OUTSIDE the hashed block, which is unchanged).  (a) The (r7) mu closure does
NOT hold as declared: 0.8952 (IA+MEC, +7.6 %) and 0.8790 (IA, +6.0 %).  The
cause is the digitization, not the reading: the theory F_T^2 symbols stand
on a grid of q = 0.05 n (n = 1 .. 5) and 0.25 n fm^-1, and the device put the
odd multiples of 0.05 half a unit low (q = 0.0476 for 0.05: 4.8 % of q,
9.8 % in eta at the lowest symbol), which the declared extrapolation through
q = 0.0476 and 0.1000 turns into ~ 5 % in mu.  A DIAGNOSTIC, added after the
run and never a verdict: the same closure at each symbol's NOMINAL q (its
digitized q rounded to the nearest 0.05 fm^-1; all 40 theory F_T^2 symbols
lie within one device unit of it -- checked on every run) gives mu = 0.8414
(IA+MEC, +1.1 %) and 0.8268 (IA, -0.3 %), both inside the declared 2 %.  The
declared closure stays FAIL as recorded; the diagnostic says the B/Z^2
reading of F_T^2 holds to about 1 %.  (b) For the same reason the (r5)
chi2 carries the quantization of the lowest symbols (the model goes as q^2
there); at the nominal q it is 77.70/19 instead of 88.91/19 -- the shape
mismatch (pulls -2.4 .. -4.4 across the first peak, q = 0.5 .. 1.25 fm^-1,
and +2.0 .. +3.1 across 2.25 .. 3.0 fm^-1) is not the quantization.
(c) The block's sigma "(0.010440)" is the device unit in ln F^2, 0.0104392,
rounded to four significant figures; the code uses the block's number
(0.008 % larger) and its self-check against the device unit allows that
rounding.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations / A-scaling -- NOT APPLICABLE: no b1, no A-scaling.
  2. isoscalar denominator -- NOT APPLICABLE: no structure function, no
     nuclear PDF grid.  Its analogue is the nucleon folding: the tree's fc is
     the point shape times the ISOSCALAR G_E^p + G_E^n (N = Z); WS98's charge
     operator carries the nucleon form factors in its one-body term (the
     parametrization is not printed) -- both sides are CHARGE form factors.
  3. per nucleus / per nucleon -- both sides PER NUCLEUS: WS98's F_L^2 is
     normalised to 1 at q = 0 (the tree's (fc^2 + F_C2^2)/Z^2, Z = 3) and its
     F_T^2 is B(Q^2)/Z^2; the tree's fc(0) = Z, fm(0) = (M_A/m_p) mu,
     fq(0) = M_A^2 Q (rc.hpp).  The zero and the shapes do not depend on the
     normalisation; the (r7) closures check it.

WHAT THIS DOES NOT VALIDATE.  Data: WS98 is a calculation (a VMC wave
function 10 % underbound, its point radius held near 2.43 fm by
construction) -- rows 5 and 6 are the data level, and a PASS here says only
that the tree's zero is not contradicted by one ab-initio calculation's
figure.  The zero's position above q_last (the figure does not draw it).
The elastic radiative tail itself (row 7) and the form factor where the
t-peak integrand reaches (q ~ 88-278 fm^-1; WS98 stops at 4 fm^-1).  The C2
normalisation (record (r7)).  7Li (HoSpin1FF.for_ion throws).  Nothing here
moves a default: shipping the WS98 table (the slot rows) is a registry
decision, and the HO refit is open item Q1.

Run:  source env.sh && python3 validation/benchmarks/t3_li6_ff_wiringa_schiavilla.py [--plot]
      (dev time: ... --redigitize PATH/TO/nucl-th_9807037.pdf)

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, an edited design block, a missing dependency, any
exception).  A harness is a measurement, not a CI gate: the pytests are the
gate.  `--redigitize` exits 0 when the re-read body is identical, 1 when not.
"""

import csv
import hashlib
import io
import math
import os
import re
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "ws1998_li6_formfactors.csv")
#: sha256 of EVERY byte of the vendored file, as written 2026-10-27T00:19:51Z;
#: re-checked on every read
DATA_SHA256 = "c51d31f9e6e2328dcae7818c4bf1f2a3f759780d6ea5da0d7c0715797dee219b"
#: row 6's RFY Table I/II file (t3_li_magnetization_rfy.py's FF_SHA256), read
#: ONLY for the (r7) RFY closure of the digitized Stanford F_T^2 points
RFY_DATA = os.path.join(HERE, "data", "rfy1966_li_magnetic_ff.csv")
RFY_SHA256 = "a3ada6b099c50e8fd6c90453e6fd964f8436b97675d1bdb17c0a75bc4e593aa8"
#: sha256 of this header's DESIGN BLOCK (its two marker lines included), as
#: written at 2026-10-27T00:21:38Z, before this harness computed any
#: comparison.  `design_block_sha256()` re-checks it, and the rule constants
#: below against the block's own sentences, on every run.
DESIGN_BLOCK_SHA256 = (
    "a9f1b3d79dc8a771ef2dbbc0c4397ec0aff44290de90048e165d23b647556a0b")
DESIGN_BEGIN = "==== DESIGN BLOCK BEGIN ===="
DESIGN_END = "==== DESIGN BLOCK END ===="

NAME = "t3_li6_ff_wiringa_schiavilla"
#: a library evaluation on a vendored table, well under 30 s, no event sample
EXPENSIVE = False

# ------------------------------------------------------------ the digitization
#: the PDF the numbers were read from (the author's corpus copy of
#: arXiv:nucl-th/9807037v2); NOT in this repository -- `--redigitize PDF`
#: re-reads it (dev time, PyMuPDF) and compares with the vendored body
WS_PDF_SHA256 = "fc51cb930e81b16238bb02561e81d4f9160a7d1708e813f8380ba669991d20e8"
WS_PAGE = 10                       # 1-based: "Fig.1 Wiringa & Schiavilla"
#: the two plot boxes, PDF points (x0, y0, x1, y1), y down -- each is ONE
#: 'qu' (quad) path in the page's drawing operators
FRAMES = {"left": (104.40, 214.56, 306.00, 585.12),
          "right": (306.00, 214.56, 507.60, 585.12)}
#: the axis values at the box edges: (q left, q right) fm^-1 and
#: (log10 bottom, log10 top) of the plotted F^2 (left F_L^2, right F_T^2)
AXES = {"left": ((0.0, 4.0), (-7.0, 0.0)),
        "right": ((0.0, 4.0), (-10.0, -3.0))}
#: the drawing's device resolution: every coordinate on the page is a
#: multiple of 0.24 pt (1/300 in); ONE unit is the declared digitization
#: error of every coordinate (the quantization, and the +-1-unit irregularity
#: of the tick marks: the q = 2 tick of the left panel sits one unit low)
UNIT_PT = 0.24
DQ_FM = UNIT_PT / ((FRAMES["left"][2] - FRAMES["left"][0])
                   / (AXES["left"][0][1] - AXES["left"][0][0]))
DLOG10F2 = UNIT_PT / ((FRAMES["left"][3] - FRAMES["left"][1])
                      / (AXES["left"][1][1] - AXES["left"][1][0]))

RED, BLUE, BLACK, WHITE = (1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)
SOLID, DASHED = "[] 0", "[ 4.8 4.8 ] 0"
#: line style -> calculation (the caption: dashed = IA, solid = IA+MEC; the
#: legend: blue = IA, red = IA+MEC)
LINE_CALC = {(RED, SOLID): "IA+MEC", (BLUE, DASHED): "IA"}

COLUMNS = ["q_fm", "F_C0", "F_C2", "F_M1", "source", "curve", "calc", "kind",
           "drawing", "branch", "x_pt", "y_pt", "y_bar_lo_pt", "y_bar_hi_pt",
           "F2", "F2_bar_lo", "F2_bar_hi", "dq_fm", "dlog10F2"]
#: every (curve, calc, kind) the body must hold, with its row count
EXPECTED_COUNTS = {
    ("C0", "IA+MEC", "vertex"): 74, ("C0", "IA+MEC", "clip"): 1,
    ("C0", "IA", "vertex"): 81, ("C0", "IA", "clip"): 1,
    ("C2", "IA+MEC", "vertex"): 95, ("C2", "IA+MEC", "clip"): 3,
    ("C2", "IA", "vertex"): 95, ("C2", "IA", "clip"): 3,
    ("FL2", "IA+MEC", "symbol"): 20, ("FL2", "IA", "symbol"): 20,
    ("FL2", "data:Stanford", "symbol"): 39,
    ("FT2", "IA+MEC", "symbol"): 20, ("FT2", "IA", "symbol"): 20,
    ("FT2", "data:Amsterdam", "symbol"): 9, ("FT2", "data:Saskatoon", "symbol"): 15,
    ("FT2", "data:Stanford", "symbol"): 7,
    ("slot", "IA+MEC", "slot"): 68,
}


def to_q(x_pt, panel):
    x0, _, x1, _ = FRAMES[panel]
    (q0, q1), _ = AXES[panel]
    return q0 + (x_pt - x0) / (x1 - x0) * (q1 - q0)


def to_log10f2(y_pt, panel):
    _, y0, _, y1 = FRAMES[panel]
    _, (lo, hi) = AXES[panel]
    return lo + (y1 - y_pt) / (y1 - y0) * (hi - lo)


def to_f2(y_pt, panel):
    return 10.0 ** to_log10f2(y_pt, panel)


def _r2(v):
    """A drawing coordinate as the device wrote it: a multiple of 0.24 pt,
    exact to two decimals (PyMuPDF hands back its float32 neighbour)."""
    return round(float(v), 2)


def _col(c):
    return None if c is None else tuple(round(float(v), 3) for v in c)


def _shape(d):
    return "".join(it[0] for it in d["items"])


def _panels(ds):
    """Drawing indices of each panel's CONTENT.  PAW/HIGZ draws a panel's
    data, then its frame (one 'qu' quad), then its legend: the left panel is
    everything before the first quad; the left legend follows that quad and
    lies inside the LEFT frame, so the right panel is what lies between the
    two quads with its centre right of the left frame."""
    quads = [i for i, d in enumerate(ds) if _shape(d) == "qu"]
    if len(quads) != 2:
        raise RuntimeError("expected two frame quads on page %d, found %d"
                           % (WS_PAGE, len(quads)))
    for name, qi in zip(("left", "right"), quads):
        q = ds[qi]["items"][0][1]
        box = (_r2(min(p.x for p in q)), _r2(min(p.y for p in q)),
               _r2(max(p.x for p in q)), _r2(max(p.y for p in q)))
        if box != FRAMES[name]:
            raise RuntimeError("frame %s is %r, the declared one %r"
                               % (name, box, FRAMES[name]))
    xr = FRAMES["right"][0] - 3.0
    right = [i for i in range(quads[0] + 1, quads[1])
             if 0.5 * (ds[i]["rect"].x0 + ds[i]["rect"].x1) > xr]
    return {"left": list(range(0, quads[0])), "right": right}


def _polyline(d):
    """The vertices of a stroked polyline in DRAWING order, split where a new
    subpath starts (a clip at the frame floor ends one and starts the next)."""
    branches, cur, last = [], [], None
    for it in d["items"]:
        if it[0] != "l":
            raise RuntimeError("a curve path holds a %r item" % (it[0],))
        p1 = (_r2(it[1].x), _r2(it[1].y))
        p2 = (_r2(it[2].x), _r2(it[2].y))
        if last is None or p1 != last:
            if cur:
                branches.append(cur)
            cur = [p1]
        cur.append(p2)
        last = p2
    branches.append(cur)
    return branches


def _centre(d):
    r = d["rect"]
    return _r2(0.5 * (r.x0 + r.x1)), _r2(0.5 * (r.y0 + r.y1))


def _series(ds, idx):
    """The symbol series of one panel.  PAW draws each as [bars][caps][n
    symbols]: a path of n vertical segments (one per symbol, in the symbols'
    order; zero length = an error below the device resolution), a path of
    the bar caps, then the n symbols.  Returns [(calc, [(drawing, x, y,
    y_bar_lo, y_bar_hi), ...]), ...]."""
    out, k = [], 0
    while k < len(idx):
        d = ds[idx[k]]
        col, fill, shape = _col(d.get("color")), _col(d.get("fill")), _shape(d)
        calc, step, stroke = None, 1, None
        if d.get("type") == "fs" and col in (RED, BLUE) and fill == col:
            calc, stroke = {(RED, "re"): "IA+MEC", (BLUE, "cccc"): "IA"}.get((col, shape)), col
        elif d.get("type") == "fs" and col == BLACK and fill == WHITE:
            calc, stroke = {"re": "data:Stanford", "cccc": "data:Saskatoon"}.get(shape), BLACK
        elif d.get("type") == "f" and col is None and fill == WHITE and shape == "ll":
            calc, step, stroke = "data:Amsterdam", 2, BLACK   # white fill + black outline
        if calc is None:
            k += 1
            continue
        pts, j = [], k
        while j + step - 1 < len(idx):
            if step == 1:
                dj = ds[idx[j]]
                if not (dj.get("type") == d.get("type") and _col(dj.get("color")) == col
                        and _col(dj.get("fill")) == fill and _shape(dj) == shape):
                    break
                pts.append((idx[j],) + _centre(dj))
            else:
                f_, s_ = ds[idx[j]], ds[idx[j + 1]]
                if not (f_.get("type") == "f" and _col(f_.get("fill")) == WHITE
                        and _shape(f_) == "ll" and s_.get("type") == "s"
                        and _col(s_.get("color")) == BLACK and _shape(s_) == "lll"):
                    break
                pts.append((idx[j + 1],) + _centre(s_))
            j += step
        bar = ds[idx[k] - 2] if idx[k] >= 2 else None
        if (bar is None or bar.get("type") != "s" or _col(bar.get("color")) != stroke
                or len(bar["items"]) != len(pts)
                or not all(it[0] == "l" and abs(it[1].x - it[2].x) < 0.01
                           for it in bar["items"])):
            raise RuntimeError("series %s at drawing %d: drawing %d is not its "
                               "%d-bar path" % (calc, idx[k], idx[k] - 2, len(pts)))
        rows = []
        for (di, x, y), it in zip(pts, bar["items"]):
            if abs(_r2(it[1].x) - x) > 0.5:
                raise RuntimeError("a bar at x = %.2f is not on its symbol at %.2f"
                                   % (it[1].x, x))
            ya, yb = _r2(it[1].y), _r2(it[2].y)
            rows.append((di, x, y, max(ya, yb), min(ya, yb)))
        out.append((calc, rows))
        k = j
    return out


def digitize(pdf):
    """Every object of WS98 Fig. 1 as a record row (DEV TIME: PyMuPDF and the
    PDF; nothing at run time needs either)."""
    import pymupdf  # noqa: WPS433 (dev-time only)
    with open(pdf, "rb") as f:
        got = hashlib.sha256(f.read()).hexdigest()
    if got != WS_PDF_SHA256:
        raise RuntimeError("%s: sha256 %s is not the arXiv nucl-th/9807037v2 PDF "
                           "this table was read from (%s)" % (pdf, got, WS_PDF_SHA256))
    ds = pymupdf.open(pdf)[WS_PAGE - 1].get_drawings()
    pan = _panels(ds)
    rows = []
    for panel in ("left", "right"):
        x0, y0, _, y1 = FRAMES[panel]
        for i in pan[panel]:                              # the theory lines
            d = ds[i]
            col = _col(d.get("color"))
            if not (d.get("type") == "s" and col in (RED, BLUE) and len(d["items"]) >= 60):
                continue
            calc = LINE_CALC.get((col, d.get("dashes")))
            if calc is None:
                raise RuntimeError("drawing %d: an unexpected line style" % i)
            br = _polyline(d)
            # C0 starts at (q, F^2) = (0, 1), the frame's top-left corner;
            # C2 starts on the frame floor (it rises as q^4)
            c0 = abs(br[0][0][0] - x0) < 0.01 and abs(br[0][0][1] - y0) < 0.01
            for b, branch in enumerate(br):
                for (x, y) in branch:
                    rows.append(dict(panel=panel, curve="C0" if c0 else "C2",
                                     calc=calc, drawing=i, branch=b, x=x, y=y,
                                     kind="clip" if abs(y - y1) < 0.01 else "vertex",
                                     ylo=None, yhi=None))
        quantity = "FL2" if panel == "left" else "FT2"
        for calc, pts in _series(ds, pan[panel]):          # symbols + bars
            for (di, x, y, ylo, yhi) in pts:
                rows.append(dict(panel=panel, curve=quantity, calc=calc,
                                 drawing=di, branch=0, x=x, y=y, kind="symbol",
                                 ylo=ylo, yhi=yhi))
    return rows


# ------------------------------------------------------------ the slot table

def _signed_m1(sym, m_a_fm, z):
    """The IA+MEC F_T^2 symbols as the tree's Sachs F_m (TabulatedSpin1FF's
    m1 column): Z^2 F_T^2 = B = (4/3) eta (1 + eta) F_m^2, eta = q^2/(4 M_A^2).
    Sign + up to the symbol at the first-zero minimum, - after it; the
    minimum symbol itself takes the sign of the ln-q interpolation of its
    signed neighbours at its q (the reader's own interpolation)."""
    pts = sorted(sym, key=lambda r: r["q"])
    mag = []
    for r in pts:
        eta = r["q"] ** 2 / (4.0 * m_a_fm ** 2)
        mag.append(z * math.sqrt(r["F2"] / (4.0 / 3.0 * eta * (1.0 + eta))))
    inner = [k for k, r in enumerate(pts) if 0.5 < r["q"] < 2.0]
    kz = min(inner, key=lambda k: pts[k]["F2"])
    a, b = mag[kz - 1], -mag[kz + 1]
    t = math.log(pts[kz]["q"] / pts[kz - 1]["q"]) / math.log(pts[kz + 1]["q"] / pts[kz - 1]["q"])
    sz = 1.0 if a + (b - a) * t >= 0.0 else -1.0
    sign = [1.0 if k < kz else (sz if k == kz else -1.0) for k in range(len(pts))]
    return [(r["q"], s * m) for r, s, m in zip(pts, sign, mag)], kz


def _interp_lnq(qs, fs, q):
    """Linear in ln q between the bracketing points -- TabulatedSpin1FF's
    `log_interp_strict`, refusing to extrapolate."""
    if not qs[0] <= q <= qs[-1]:
        raise RuntimeError("q = %g outside [%g, %g]" % (q, qs[0], qs[-1]))
    for i in range(len(qs) - 1):
        if qs[i] <= q <= qs[i + 1]:
            t = math.log(q / qs[i]) / math.log(qs[i + 1] / qs[i])
            return fs[i] + (fs[i + 1] - fs[i]) * t
    return fs[-1]


def slot_rows(record, m_a_fm, z):
    """The SLOT TABLE (data/ff/li6_elastic.csv's columns q, F_C0, F_C2, F_M1),
    IA+MEC, derived from the record rows by the rule stated in the vendored
    file's provenance block: one row per computed C0 vertex (q > 0) that has
    a computed C2 vertex at the same x; F_C0 = +Z |F_C0| (the drawn branch is
    before the first zero); F_C2 = -Z |F_C2| on the C2 line's first branch
    (the sign of WS98's printed Q = -0.23 fm^2) and +Z |F_C2| on its second
    (a simple zero in the floor gap); F_M1 interpolated linearly in ln q of
    the signed Sachs F_m between the IA+MEC F_T^2 symbols (_signed_m1)."""
    c0 = [r for r in record if r["curve"] == "C0" and r["calc"] == "IA+MEC"
          and r["kind"] == "vertex" and r["q"] > 0.0]
    c2 = {r["x"]: r for r in record if r["curve"] == "C2" and r["calc"] == "IA+MEC"
          and r["kind"] == "vertex"}
    m1, _ = _signed_m1([r for r in record if r["curve"] == "FT2"
                        and r["calc"] == "IA+MEC"], m_a_fm, z)
    mq, mf = [p[0] for p in m1], [p[1] for p in m1]
    out = []
    for r in sorted(c0, key=lambda r: r["q"]):
        s = c2.get(r["x"])
        if s is None:
            continue
        out.append(dict(q=r["q"], x=r["x"],
                        F_C0=z * math.sqrt(r["F2"]),
                        F_C2=(-1.0 if s["branch"] == 0 else 1.0) * z * math.sqrt(s["F2"]),
                        F_M1=_interp_lnq(mq, mf, r["q"])))
    return out


def _tree_constants():
    import lipolgen as lg  # noqa: WPS433 (env.sh puts it on the path)
    ion = lg.li6()
    return ion.mass() / lg.HBARC_GEV_FM, float(ion.Z)


def _fmt(v, f):
    return "" if v is None else f % v


def body_text(record, m_a_fm, z):
    """The vendored file's body (the column line and every row), exactly as
    written: record rows (slot columns empty, so TabulatedSpin1FF's reader
    skips them) and then the slot rows (all four numeric)."""
    for r in record:
        r["q"] = to_q(r["x"], r["panel"])
        r["F2"] = to_f2(r["y"], r["panel"])
    buf = io.StringIO()
    buf.write(",".join(COLUMNS) + "\n")
    for r in record:
        p = r["panel"]
        buf.write(",".join([
            "%.6f" % r["q"], "", "", "", "WS98 Fig.1 " + p, r["curve"], r["calc"],
            r["kind"], "%d" % r["drawing"], "%d" % r["branch"],
            "%.2f" % r["x"], "%.2f" % r["y"], _fmt(r["ylo"], "%.2f"),
            _fmt(r["yhi"], "%.2f"), "%.6g" % r["F2"],
            _fmt(None if r["ylo"] is None else to_f2(r["ylo"], p), "%.6g"),
            _fmt(None if r["yhi"] is None else to_f2(r["yhi"], p), "%.6g"),
            "%.7f" % DQ_FM, "%.7f" % DLOG10F2]) + "\n")
    for s in slot_rows(record, m_a_fm, z):
        buf.write(",".join([
            "%.6f" % s["q"], "%.6g" % s["F_C0"], "%.6g" % s["F_C2"],
            "%.6g" % s["F_M1"], "WS98 Fig.1 slot", "slot", "IA+MEC", "slot", "",
            "", "%.2f" % s["x"], "", "", "", "", "", "", "%.7f" % DQ_FM,
            "%.7f" % DLOG10F2]) + "\n")
    return buf.getvalue()


def redigitize(pdf):
    """DEV TIME: re-read the PDF and return the body this harness would
    vendor (compare with the file's body; `--redigitize PDF` does)."""
    m_a_fm, z = _tree_constants()
    return body_text(digitize(pdf), m_a_fm, z)


# ------------------------------------------------------------ the design block

#: the verdict (design block): PASS iff q0_tree >= q_last - (DELTA_FM + ZERO_MARGIN_FM)
ZERO_MARGIN_FM = 0.1
DELTA_FM = 0.0047619
#: row 5's first_zero: steps of Q_SCAN_STEP on (0, Q_SCAN_MAX], bisected to BISECT_TOL
Q_SCAN_STEP = 0.01
Q_SCAN_MAX = 5.0
BISECT_TOL = 1e-12
#: (r3) the slope-radius fit window and its model variants
SLOPE_QMAX = 1.0
SLOPE_QMAX_VARIANTS = (0.6, 1.0, 1.4)
SLOPE_ORDERS = (2, 3)
#: one device unit in ln F^2 (the vendored dlog10F2 times ln 10)
SIGMA_LNF2 = 0.010440
#: (r5) the digitization error of an F_T^2 symbol, relative
M1_DIG_REL = 0.0105
#: (r7) closure rules
MU_TOL = 0.02
MU_PRINTED = {"IA+MEC": 0.832, "IA": 0.829}      # WS98 p. 5, mu_N
Q_TOL = 0.05
Q_PRINTED_FM2 = -0.23                              # WS98 p. 4, -0.23(9) fm^2
Q_PRINTED_ERR_FM2 = 0.09
LINESUM_NSIGMA = 2.0
LINESUM_REL = 0.021
LINESUM_QRANGE = (0.25, 2.75)
RFY_Q_TOL_FM = 0.0095
RFY_RATIO_SPREAD = 1.042
RFY_MU0 = 0.822                                    # RFY Table I caption, nm
FLOOR_F2 = 1e-7

_RULE_RES = (
    (r"PASS iff q0_tree >= q_last - \(delta \+ ([0-9.]+) fm\^-1\)", float, "ZERO_MARGIN_FM"),
    (r"delta = one device unit of the drawing in q = ([0-9.]+) fm\^-1", float, "DELTA_FM"),
    (r"\(steps of ([0-9.]+) fm\^-1, bisected to ([0-9e.+-]+) --", None, None),
    (r"on \(0, ([0-9.]+)\] fm\^-1 \(steps", float, "Q_SCAN_MAX"),
    (r"vertices with 0 < q <= ([0-9.]+) fm\^-1, sigma = one device unit", float, "SLOPE_QMAX"),
    (r"q_max in\s+\{([0-9., ]+)\} fm\^-1", None, None),
    (r"in ln F\^2\s+\(([0-9.]+)\) per vertex", float, "SIGMA_LNF2"),
    (r"the digitization error \(([0-9.]+) % of F2\)", None, None),
    (r"HOLDS iff \|mu/0\.832 - 1\| <= ([0-9.]+)", float, "MU_TOL"),
    (r"HOLDS iff \|c_0/0\.23 - 1\| <=\s+([0-9.]+)", float, "Q_TOL"),
    (r"HOLDS iff within ([0-9.]+) sigma_i", float, "LINESUM_NSIGMA"),
    (r"half-bar and ([0-9.]+) % of F2 in quadrature \(three", None, None),
    (r"with ([0-9.]+) <= q <= ([0-9.]+)\s+fm\^-1 bracketed", None, None),
    (r"<= ([0-9.]+) fm\^-1\s+\(two device units\)", float, "RFY_Q_TOL_FM"),
    (r"constant iff max/min <= ([0-9.]+)", float, "RFY_RATIO_SPREAD"),
)


def design_block(source=None):
    """The DESIGN BLOCK of this file's header, both marker lines included."""
    if source is None:
        with open(os.path.abspath(__file__), encoding="utf-8") as f:
            source = f.read()
    i0 = source.index(DESIGN_BEGIN)
    i1 = source.index(DESIGN_END) + len(DESIGN_END)
    return source[i0:i1]


def rules_from_design_block(block):
    """Every numerical rule as the block's sentences state it -- RuntimeError
    if a sentence is missing."""
    got = {}
    for rx, conv, key in _RULE_RES:
        m = re.search(rx, block)
        if m is None:
            raise RuntimeError("the design block no longer states %r" % rx)
        got[rx] = m.groups()
    out = dict(
        ZERO_MARGIN_FM=float(got[_RULE_RES[0][0]][0]),
        DELTA_FM=float(got[_RULE_RES[1][0]][0]),
        Q_SCAN_STEP=float(got[_RULE_RES[2][0]][0]),
        BISECT_TOL=float(got[_RULE_RES[2][0]][1]),
        Q_SCAN_MAX=float(got[_RULE_RES[3][0]][0]),
        SLOPE_QMAX=float(got[_RULE_RES[4][0]][0]),
        SLOPE_QMAX_VARIANTS=tuple(float(v) for v in got[_RULE_RES[5][0]][0].split(",")),
        SIGMA_LNF2=float(got[_RULE_RES[6][0]][0]),
        M1_DIG_REL=float(got[_RULE_RES[7][0]][0]) / 100.0,
        MU_TOL=float(got[_RULE_RES[8][0]][0]),
        Q_TOL=float(got[_RULE_RES[9][0]][0]),
        LINESUM_NSIGMA=float(got[_RULE_RES[10][0]][0]),
        LINESUM_REL=float(got[_RULE_RES[11][0]][0]) / 100.0,
        LINESUM_QRANGE=tuple(float(v) for v in got[_RULE_RES[12][0]]),
        RFY_Q_TOL_FM=float(got[_RULE_RES[13][0]][0]),
        RFY_RATIO_SPREAD=float(got[_RULE_RES[14][0]][0]),
    )
    return out


def design_block_sha256(source=None):
    """sha256 of the design block -- RuntimeError unless it is the recorded
    one AND the block's rules are the constants this code uses."""
    block = design_block(source)
    got = hashlib.sha256(block.encode("utf-8")).hexdigest()
    if got != DESIGN_BLOCK_SHA256:
        raise RuntimeError("the header's design block sha256 %s != recorded %s: "
                           "the design or the tolerance was edited"
                           % (got, DESIGN_BLOCK_SHA256))
    stated = rules_from_design_block(block)
    used = dict(ZERO_MARGIN_FM=ZERO_MARGIN_FM, DELTA_FM=DELTA_FM,
                Q_SCAN_STEP=Q_SCAN_STEP, BISECT_TOL=BISECT_TOL,
                Q_SCAN_MAX=Q_SCAN_MAX, SLOPE_QMAX=SLOPE_QMAX,
                SLOPE_QMAX_VARIANTS=SLOPE_QMAX_VARIANTS, SIGMA_LNF2=SIGMA_LNF2,
                M1_DIG_REL=M1_DIG_REL, MU_TOL=MU_TOL, Q_TOL=Q_TOL,
                LINESUM_NSIGMA=LINESUM_NSIGMA, LINESUM_REL=LINESUM_REL,
                LINESUM_QRANGE=LINESUM_QRANGE, RFY_Q_TOL_FM=RFY_Q_TOL_FM,
                RFY_RATIO_SPREAD=RFY_RATIO_SPREAD)
    for k, v in used.items():
        w = stated[k]
        same = (all(abs(a - b) <= 1e-12 * max(1.0, abs(b)) for a, b in zip(v, w))
                and len(v) == len(w)) if isinstance(v, tuple) else \
            abs(v - w) <= 1e-12 * max(1.0, abs(w))
        if not same:
            raise RuntimeError("rule %s = %r differs from the hashed design "
                               "block's %r" % (k, v, w))
    # the block's delta and sigma are the drawing's one device unit, to the
    # precision the block prints them (0.0047619 = 0.00476190..., 0.010440 =
    # 0.01043918... rounded to four significant figures; the code uses the
    # block's numbers, as stated)
    if abs(DELTA_FM - DQ_FM) > 5e-8 or abs(SIGMA_LNF2 - DLOG10F2 * math.log(10.0)) > 5e-6:
        raise RuntimeError("the design block's delta / sigma are not the "
                           "drawing's one device unit")
    return got


# ------------------------------------------------------------ the reference

def read_vendored(path, sha256):
    """The bytes of a vendored file -- RuntimeError unless the sha256 of
    every byte is the recorded one."""
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != sha256:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, sha256))
    return raw


def _num(v):
    return None if v == "" else float(v)


def _close(a, b, rel, absol=0.0):
    return abs(a - b) <= max(absol, rel * max(abs(a), abs(b)))


def load_reference(path=DATA, m_a_fm=None, z=None):
    """The vendored rows, checked: sha256 of every byte, the columns, the row
    count of every (curve, calc, kind), q and F2 of every row re-derived from
    its PDF coordinates, and every slot row re-derived from the record rows
    (that needs M_A and Z: the tree's, unless given).  Returns
    dict(record=[...], slot=[...], body=text)."""
    raw = read_vendored(path, DATA_SHA256).decode("utf-8")
    lines = [l for l in raw.splitlines() if l.strip() and not l.startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames != COLUMNS:
        raise RuntimeError("%s: columns %r, expected %r" % (path, reader.fieldnames, COLUMNS))
    record, slot, counts = [], [], {}
    for r in reader:
        key = (r["curve"], r["calc"], r["kind"])
        counts[key] = counts.get(key, 0) + 1
        if r["kind"] == "slot":
            slot.append(dict(q=float(r["q_fm"]), x=float(r["x_pt"]),
                             F_C0=float(r["F_C0"]), F_C2=float(r["F_C2"]),
                             F_M1=float(r["F_M1"])))
            continue
        if r["F_C0"] or r["F_C2"] or r["F_M1"]:
            raise RuntimeError("%s: a record row carries a slot value" % path)
        panel = r["source"].rsplit(" ", 1)[-1]
        if panel not in FRAMES:
            raise RuntimeError("%s: source %r" % (path, r["source"]))
        row = dict(panel=panel, curve=r["curve"], calc=r["calc"], kind=r["kind"],
                   drawing=int(r["drawing"]), branch=int(r["branch"]),
                   x=float(r["x_pt"]), y=float(r["y_pt"]), ylo=_num(r["y_bar_lo_pt"]),
                   yhi=_num(r["y_bar_hi_pt"]), q=float(r["q_fm"]), F2=float(r["F2"]),
                   F2_lo=_num(r["F2_bar_lo"]), F2_hi=_num(r["F2_bar_hi"]))
        # every number of the row from its coordinates, through the calibration
        if abs(to_q(row["x"], panel) - row["q"]) > 6e-7:
            raise RuntimeError("%s: q %.6f is not x_pt %.2f calibrated" % (path, row["q"], row["x"]))
        if not _close(to_f2(row["y"], panel), row["F2"], 6e-6):
            raise RuntimeError("%s: F2 %g is not y_pt %.2f calibrated" % (path, row["F2"], row["y"]))
        for yk, fk in (("ylo", "F2_lo"), ("yhi", "F2_hi")):
            if (row[yk] is None) != (row[fk] is None) or (
                    row[yk] is not None and not _close(to_f2(row[yk], panel), row[fk], 6e-6)):
                raise RuntimeError("%s: a bar end does not match its coordinate" % path)
        if abs(float(r["dq_fm"]) - DQ_FM) > 1e-7 or abs(float(r["dlog10F2"]) - DLOG10F2) > 1e-7:
            raise RuntimeError("%s: a row's digitization error is not one device unit" % path)
        record.append(row)
    if counts != EXPECTED_COUNTS:
        raise RuntimeError("%s: row counts %r, expected %r"
                           % (path, sorted(counts.items()), sorted(EXPECTED_COUNTS.items())))
    if m_a_fm is None:
        m_a_fm, z = _tree_constants()
    exact = [dict(r, q=to_q(r["x"], r["panel"]), F2=to_f2(r["y"], r["panel"])) for r in record]
    derived = slot_rows(exact, m_a_fm, z)
    if len(derived) != len(slot):
        raise RuntimeError("%s: %d slot rows, the rule gives %d" % (path, len(slot), len(derived)))
    for a, b in zip(slot, derived):
        if not (abs(a["q"] - b["q"]) <= 6e-7 and abs(a["x"] - b["x"]) < 0.005
                and all(_close(a[k], b[k], 6e-6) for k in ("F_C0", "F_C2", "F_M1"))):
            raise RuntimeError("%s: the slot row at q = %.6f is not the rule's %r"
                               % (path, a["q"], b))
    return dict(record=record, slot=slot, body="\n".join(lines) + "\n")


def body_of_file(path=DATA):
    raw = read_vendored(path, DATA_SHA256).decode("utf-8")
    return "".join(l + "\n" for l in raw.splitlines() if not l.startswith("#"))


def load_rfy_table1(path=RFY_DATA):
    """RFY Table I (6Li) from row 6's file, its sha256 re-checked."""
    raw = read_vendored(path, RFY_SHA256).decode("utf-8")
    lines = [l for l in raw.splitlines() if l.strip() and not l.startswith("#")]
    rows = [r for r in csv.DictReader(lines) if r["nucleus"] == "6Li" and r["table"] == "I"]
    if len(rows) != 7:
        raise RuntimeError("%s: expected 7 Table I rows, found %d" % (path, len(rows)))
    return [dict(row=int(r["row"]), q2_fm2=float(r["q2_fm2"]), fm2=float(r["fm2"]),
                 fm2_err=float(r["fm2_err"])) for r in rows]


def pick(record, curve, calc, kind="vertex"):
    return sorted((r for r in record if r["curve"] == curve and r["calc"] == calc
                   and r["kind"] == kind), key=lambda r: (r["q"], r["branch"]))


# ------------------------------------------------------------ the tree side

def tree_side():
    """The two C0 edges of the shipped 6Li form factor, as WS98 plot them."""
    import lipolgen as lg  # noqa: WPS433
    ion = lg.li6()
    ho = lg.HoSpin1FF.for_ion(ion)                      # SHIPPED, nothing moved
    o = ho.options
    if o.c0_shape != lg.C0Shape.Ho or not o.fold_nucleon:
        raise RuntimeError("HoSpin1FF.for_ion(li6()) is no longer the shipped "
                           "(ho, fold_nucleon) configuration this row names")
    ov = lg.HoSpin1FFOptions()
    ov.c0_shape = lg.C0Shape.VmcFt
    vmc = lg.HoSpin1FF.for_ion(ion, ov)
    hbarc, m_a, z = lg.HBARC_GEV_FM, ion.mass(), float(ion.Z)
    return dict(edges={"ho": ho, "vmc-ft": vmc}, hbarc=hbarc, m_a_gev=m_a,
                m_a_fm=m_a / hbarc, z=z, m_p_gev=lg.PROTON_MASS,
                options=dict(a_fm=o.a_fm, alpha=o.alpha, fm_qz_fm=o.fm_qz_fm,
                             fm_b_fm=o.fm_b_fm, mu_n=o.mu_n, q_fm2=o.q_fm2),
                provenance={k: v.provenance() for k, v in
                            (("ho", ho), ("vmc-ft", vmc))})


class Edge:
    """One C0 edge read as WS98 plot it, per nucleus, at q in fm^-1."""

    def __init__(self, ff, tree):
        self.ff, self.hbarc, self.m_a_fm, self.z = ff, tree["hbarc"], tree["m_a_fm"], tree["z"]

    def t(self, q):
        return (q * self.hbarc) ** 2

    def eta(self, q):
        return q * q / (4.0 * self.m_a_fm ** 2)

    def c0(self, q):
        return self.ff.fc(self.t(q)) / self.z

    def c2(self, q):
        return 2.0 * math.sqrt(2.0) / 3.0 * self.eta(q) * self.ff.fq(self.t(q)) / self.z

    def ft2(self, q):
        e = self.eta(q)
        return 4.0 / 3.0 * e * (1.0 + e) * self.ff.fm(self.t(q)) ** 2 / self.z ** 2


def bisect(f, lo, hi, tol=BISECT_TOL):
    flo = f(lo)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def first_zero(f, q_max=Q_SCAN_MAX, step=Q_SCAN_STEP):
    """The first sign change of f on (0, q_max], bisected; None if none
    (t3_li6_charge_ff_fb.py's first_zero)."""
    q, fq = step, f(step)
    while q < q_max - 1e-12:
        q2 = q + step
        f2 = f(q2)
        if fq * f2 <= 0:
            return bisect(f, q, q2)
        q, fq = q2, f2
    return None


def richardson_r(f, s0=0.05, levels=5):
    """The slope radius sqrt(-6 dF/ds / F(0)) at s = q^2 -> 0 (fm), Richardson
    on g(s) = 6 (1 - F/F(0))/s (t3_li6_rms_devries.py's scheme)."""
    f0 = f(0.0)
    tab = [[6.0 * (1.0 - f(s0 / 2 ** i) / f0) / (s0 / 2 ** i)] for i in range(levels)]
    for i in range(1, levels):
        for j in range(1, i + 1):
            tab[i].append(tab[i][j - 1] + (tab[i][j - 1] - tab[i - 1][j - 1]) / (2 ** j - 1))
    return math.sqrt(tab[-1][-1])


# ------------------------------------------------------------ statistics

def _gammq(a, x):
    """Regularized upper incomplete gamma Q(a, x) (Numerical Recipes)."""
    if x <= 0.0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(10000):
            ap += 1.0
            d *= x / ap
            s += d
            if abs(d) < abs(s) * 1e-15:
                break
        return 1.0 - s * math.exp(-x + a * math.log(x) - gln)
    b, c, d = x + 1.0 - a, 1.0 / 1e-300, 1.0 / (x + 1.0 - a)
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def chi2_sf(chi2, ndf):
    try:
        from scipy.stats import chi2 as _c  # noqa: WPS433
        return float(_c.sf(chi2, ndf))
    except Exception:                                   # pragma: no cover
        return _gammq(0.5 * ndf, 0.5 * chi2)


def _lstsq(cols, y, w=None):
    """Least squares y ~ sum_k a_k cols[k] (no intercept unless a column of
    ones is given); returns (a, covariance / sigma^2)."""
    import numpy as np  # noqa: WPS433
    a_mat = np.array(cols, dtype=float).T
    yv = np.array(y, dtype=float)
    if w is not None:
        sw = np.sqrt(np.array(w, dtype=float))
        a_mat, yv = a_mat * sw[:, None], yv * sw
    cov = np.linalg.inv(a_mat.T @ a_mat)
    return cov @ (a_mat.T @ yv), cov


# ------------------------------------------------------------ the comparisons

def zero_verdict(q0, q_last):
    """The design block's rule.  q0 None = no sign change below Q_SCAN_MAX."""
    bound = q_last - (DELTA_FM + ZERO_MARGIN_FM)
    if q0 is None:
        return dict(status="pass", q0_tree=None, q_last=q_last, bound=bound,
                    distance_above_bound=None,
                    reason="no C0 zero below %g fm^-1: none where the VMC C0 is "
                           "drawn non-zero" % Q_SCAN_MAX)
    return dict(status="pass" if q0 >= bound else "fail", q0_tree=q0, q_last=q_last,
                bound=bound, distance_above_bound=q0 - bound,
                distance_above_q_last=q0 - q_last, reason=None)


def extrapolations(c0v):
    """(r1): roots above q_last of four fits to the signed F of the last
    computed vertices."""
    import numpy as np  # noqa: WPS433
    q = np.array([r["q"] for r in c0v])
    f = np.sqrt(np.array([r["F2"] for r in c0v]))
    q_last = q[-1]
    out = {}
    for key, n, deg in (("e1_line_last2", 2, 1), ("e2_line_last5", 5, 1),
                        ("e3_quad_last5", 5, 2), ("e4_quad_last10", 10, 2)):
        coef = np.polyfit(q[-n:], f[-n:], deg)
        roots = [float(r.real) for r in np.roots(coef)
                 if abs(r.imag) < 1e-12 and q_last < r.real < 4.0]
        out[key] = dict(root=min(roots) if roots else None, n=n, degree=deg,
                        coefficients=[float(c) for c in coef])
    found = [v["root"] for v in out.values() if v["root"] is not None]
    return dict(fits=out, n_with_root=len(found),
                min_root=min(found) if found else None,
                max_root=max(found) if found else None)


def slope_fit(qs, ln_f2, qmax, order):
    pts = [(q, y) for q, y in zip(qs, ln_f2) if 0.0 < q <= qmax + 1e-9]
    cols = [[q ** (2 * k) for q, _ in pts] for k in range(1, order + 1)]
    a, cov = _lstsq(cols, [y for _, y in pts])
    r2 = -3.0 * a[0]
    r = math.sqrt(r2) if r2 > 0 else float("nan")
    sig_a1 = math.sqrt(cov[0][0]) * SIGMA_LNF2
    return dict(r=r, r_err_dig=1.5 * sig_a1 / r if r == r else None, n=len(pts),
                qmax=qmax, order=order)


def slope_radius(c0v, edges):
    """(r3)."""
    qs = [r["q"] for r in c0v]
    ws = [math.log(r["F2"]) for r in c0v]
    central = slope_fit(qs, ws, SLOPE_QMAX, 2)
    variants = [slope_fit(qs, ws, qm, o) for qm in SLOPE_QMAX_VARIANTS for o in SLOPE_ORDERS]
    rs = [v["r"] for v in variants]
    out = dict(ws=dict(r=central["r"], r_err_dig=central["r_err_dig"], n=central["n"],
                       r_min_variants=min(rs), r_max_variants=max(rs),
                       spread=max(rs) - min(rs), variants=variants))
    for name, e in edges.items():
        same = slope_fit(qs, [math.log(e.c0(q) ** 2) if q > 0 else 0.0 for q in qs],
                         SLOPE_QMAX, 2)
        out[name] = dict(r_same_fit=same["r"],
                         r_richardson=richardson_r(lambda s, e=e: e.c0(math.sqrt(s))))
    return out


def shape_chi2(c0v, edge):
    """(r4): ln F2 against ln(s m), one free scale."""
    pts = [r for r in c0v if r["q"] > 0.0]
    d = [math.log(r["F2"]) - math.log(edge.c0(r["q"]) ** 2) for r in pts]
    c = sum(d) / len(d)
    chi2 = sum((x - c) ** 2 for x in d) / SIGMA_LNF2 ** 2

    def ratio_near(q):
        k = min(range(len(pts)), key=lambda i: abs(pts[i]["q"] - q))
        return dict(q=pts[k]["q"], ratio=math.exp(d[k] - c))

    return dict(chi2=chi2, ndf=len(d) - 1, s=math.exp(c),
                ratio_q1=ratio_near(1.0), ratio_q2=ratio_near(2.0),
                ratio_q25=ratio_near(2.5), ratio_qlast=dict(q=pts[-1]["q"],
                                                            ratio=math.exp(d[-1] - c)),
                max_abs_log_ratio=max(abs(x - c) for x in d))


def m1_compare(sym, edge):
    """(r5)."""
    d, m, sg = [], [], []
    for r in sym:
        up = 0.0 if r["F2_hi"] is None else r["F2_hi"] - r["F2"]
        sg.append(math.hypot(up, M1_DIG_REL * r["F2"]))
        d.append(r["F2"])
        m.append(edge.ft2(r["q"]))
    sw = sum(mi * mi / s ** 2 for mi, s in zip(m, sg))
    s_hat = sum(di * mi / s ** 2 for di, mi, s in zip(d, m, sg)) / sw
    chi2 = sum(((di - s_hat * mi) / s) ** 2 for di, mi, s in zip(d, m, sg))
    chi2_1 = sum(((di - mi) / s) ** 2 for di, mi, s in zip(d, m, sg))
    ndf = len(d) - 1
    inner = [k for k, r in enumerate(sym) if 0.5 < r["q"] < 2.0]
    kz = min(inner, key=lambda k: sym[k]["F2"])
    return dict(chi2=chi2, ndf=ndf, p=chi2_sf(chi2, ndf), s_hat=s_hat,
                s_err=1.0 / math.sqrt(sw), chi2_s1=chi2_1, ndf_s1=len(d),
                p_s1=chi2_sf(chi2_1, len(d)),
                ws_zero_bracket=(sym[kz - 1]["q"], sym[kz]["q"], sym[kz + 1]["q"]),
                ws_zero_symbol_q=sym[kz]["q"], ws_zero_symbol_f2=sym[kz]["F2"],
                pulls=[(di - s_hat * mi) / s for di, mi, s in zip(d, m, sg)])


def c2_compare(c2v, edge, tree_q0):
    """(r6)."""
    low = [r for r in c2v if r["q"] <= 0.5]
    ratios = [math.sqrt(r["F2"]) / abs(edge.c2(r["q"])) for r in low]
    b0 = [r for r in c2v if r["branch"] == 0]
    b1 = [r for r in c2v if r["branch"] == 1]
    lobe = max(b1, key=lambda r: r["F2"])
    q0_fq = first_zero(lambda q: edge.ff.fq(edge.t(q)))
    return dict(ratio_low_q_mean=sum(ratios) / len(ratios), ratio_low_q_min=min(ratios),
                ratio_low_q_max=max(ratios), n_low_q=len(low),
                ws_zero_bracket=(b0[-1]["q"], b1[0]["q"]), tree_c2_zero=q0_fq,
                tree_c2_zero_equals_c0_zero=(q0_fq is not None and tree_q0 is not None
                                             and abs(q0_fq - tree_q0) < 1e-9),
                ws_second_lobe_q=lobe["q"], ws_second_lobe_f2=lobe["F2"],
                tree_f2_at_ws_lobe=edge.c2(lobe["q"]) ** 2)


def closures(rec, tree):
    """(r7)."""
    m_a, z, ratio_m = tree["m_a_fm"], tree["z"], tree["m_a_gev"] / tree["m_p_gev"]
    out = {}
    # POST-HOC DIAGNOSTIC (header addendum (a)): the theory F_T^2 symbols at
    # their nominal q (the 0.05 fm^-1 grid the device quantized)
    theory = pick(rec, "FT2", "IA+MEC", "symbol") + pick(rec, "FT2", "IA", "symbol")
    off = max(abs(r["q"] - nominal_q(r["q"])) for r in theory)
    if off > DQ_FM + 1e-9:
        raise RuntimeError("a theory F_T^2 symbol is %.4f fm^-1 off the 0.05 grid" % off)
    for calc in ("IA+MEC", "IA"):
        s2 = [dict(r, q=nominal_q(r["q"])) for r in pick(rec, "FT2", calc, "symbol")[:2]]
        mu2 = []
        for r in s2:
            e = r["q"] ** 2 / (4.0 * m_a ** 2)
            mu2.append(z * z * r["F2"] / (4.0 / 3.0 * e * (1.0 + e) * ratio_m ** 2))
        q1, q2 = s2[0]["q"] ** 2, s2[1]["q"] ** 2
        mu0 = math.sqrt(mu2[0] - (mu2[1] - mu2[0]) * q1 / (q2 - q1))
        out["mu_nominal_q[%s]" % calc] = dict(
            mu=mu0, printed=MU_PRINTED[calc], rel=mu0 / MU_PRINTED[calc] - 1.0,
            holds=abs(mu0 / MU_PRINTED[calc] - 1.0) <= MU_TOL, post_hoc=True,
            q_used=[s2[0]["q"], s2[1]["q"]], max_offset_from_grid=off)
    for calc in ("IA+MEC", "IA"):
        s2 = pick(rec, "FT2", calc, "symbol")[:2]
        mu2 = []
        for r in s2:
            e = r["q"] ** 2 / (4.0 * m_a ** 2)
            mu2.append(z * z * r["F2"] / (4.0 / 3.0 * e * (1.0 + e) * ratio_m ** 2))
        q1, q2 = s2[0]["q"] ** 2, s2[1]["q"] ** 2
        mu0 = math.sqrt(mu2[0] - (mu2[1] - mu2[0]) * q1 / (q2 - q1))
        out["mu[%s]" % calc] = dict(mu=mu0, printed=MU_PRINTED[calc],
                                    rel=mu0 / MU_PRINTED[calc] - 1.0,
                                    holds=abs(mu0 / MU_PRINTED[calc] - 1.0) <= MU_TOL,
                                    q_used=[s2[0]["q"], s2[1]["q"]])
    c2v = pick(rec, "C2", "IA+MEC")[:3]
    c = [3.0 * math.sqrt(2.0) * z * math.sqrt(r["F2"]) / r["q"] ** 2 for r in c2v]
    a, _ = _lstsq([[1.0] * 3, [r["q"] ** 2 for r in c2v]], c)
    c0q = float(a[0])
    out["Q[IA+MEC]"] = dict(abs_q_fm2=c0q, printed_fm2=Q_PRINTED_FM2,
                            rel=c0q / abs(Q_PRINTED_FM2) - 1.0,
                            holds=abs(c0q / abs(Q_PRINTED_FM2) - 1.0) <= Q_TOL,
                            q_used=[r["q"] for r in c2v],
                            pulls_of_printed_mc_error=(c0q - abs(Q_PRINTED_FM2)) / Q_PRINTED_ERR_FM2,
                            ratio_sqrt2_over_3=c0q * math.sqrt(2.0) / 3.0)
    # lines vs symbols
    c0l, c2l = pick(rec, "C0", "IA+MEC"), pick(rec, "C2", "IA+MEC")

    def at(line, q):
        for br in sorted({r["branch"] for r in line}):
            pts = [r for r in line if r["branch"] == br]
            for u, v in zip(pts, pts[1:]):
                if u["q"] <= q <= v["q"] and u["q"] > 0:
                    t = math.log(q / u["q"]) / math.log(v["q"] / u["q"])
                    return u["F2"] + (v["F2"] - u["F2"]) * t
        return None

    rows = []
    for r in pick(rec, "FL2", "IA+MEC", "symbol"):
        if not LINESUM_QRANGE[0] - 0.005 <= r["q"] <= LINESUM_QRANGE[1] + 0.005:
            continue
        a0, a2 = at(c0l, r["q"]), at(c2l, r["q"])
        if a0 is None or a2 is None:
            continue
        up = 0.0 if r["F2_hi"] is None else r["F2_hi"] - r["F2"]
        sig = math.hypot(up, LINESUM_REL * r["F2"])
        rows.append(dict(q=r["q"], symbol=r["F2"], lines=a0 + a2,
                         pull=(a0 + a2 - r["F2"]) / sig,
                         holds=abs(a0 + a2 - r["F2"]) <= LINESUM_NSIGMA * sig))
    out["lines_vs_symbols[IA+MEC]"] = dict(holds=all(x["holds"] for x in rows),
                                           n=len(rows), points=rows,
                                           max_abs_pull=max(abs(x["pull"]) for x in rows))
    # RFY
    rfy = load_rfy_table1()
    dig = pick(rec, "FT2", "data:Stanford", "symbol")
    pts = []
    for d_, t_ in zip(dig, rfy):
        q_r = math.sqrt(t_["q2_fm2"])
        e = t_["q2_fm2"] / (4.0 * m_a ** 2)
        conv = t_["fm2"] * 4.0 / 3.0 * e * (1.0 + e) * ratio_m ** 2 * RFY_MU0 ** 2 / z ** 2
        pts.append(dict(row=t_["row"], q_dig=d_["q"], q_rfy=q_r, dq=d_["q"] - q_r,
                        ratio=d_["F2"] / conv,
                        rel_err_dig=(d_["F2_hi"] - d_["F2"]) / d_["F2"],
                        rel_err_rfy=t_["fm2_err"] / t_["fm2"]))
    rs = [p["ratio"] for p in pts]
    out["rfy_table1[Stanford F_T^2]"] = dict(
        q_holds=all(abs(p["dq"]) <= RFY_Q_TOL_FM for p in pts),
        max_abs_dq=max(abs(p["dq"]) for p in pts),
        ratio_mean=sum(rs) / len(rs), ratio_min=min(rs), ratio_max=max(rs),
        constant=max(rs) / min(rs) <= RFY_RATIO_SPREAD, points=pts)
    return out


def _word(ok):
    return "pass" if ok else "fail"


def nominal_q(q):
    """POST HOC (header addendum (a)): a theory symbol's q on the 0.05 fm^-1
    grid it was computed on (the device put odd multiples of 0.05 half a
    unit low)."""
    return round(q / 0.05) * 0.05


def measure(path=DATA):
    design = design_block_sha256()
    tree = tree_side()
    ref = load_reference(path, tree["m_a_fm"], tree["z"])
    rec = ref["record"]
    edges = {k: Edge(v, tree) for k, v in tree["edges"].items()}
    c0 = {calc: pick(rec, "C0", calc) for calc in ("IA+MEC", "IA")}
    q_last = {calc: v[-1]["q"] for calc, v in c0.items()}
    q0 = {k: first_zero(e.c0) for k, e in edges.items()}
    verdict = {k: zero_verdict(q0[k], q_last["IA+MEC"]) for k in edges}
    verdict_ia = {k: zero_verdict(q0[k], q_last["IA"]) for k in edges}
    ceiling = {}
    for k, e in edges.items():
        grid = [q_last["IA+MEC"] + (4.0 - q_last["IA+MEC"]) * i / 1000.0 for i in range(1, 1001)]
        v = [e.c0(q) ** 2 for q in grid]
        i = max(range(len(v)), key=lambda j: v[j])
        ceiling[k] = dict(max_f2=v[i], at_q=grid[i], above_floor=v[i] > FLOOR_F2)
    shape = {(k, calc): shape_chi2(c0[calc], e) for k, e in edges.items() for calc in c0}
    return dict(
        design_block_sha256=design, tree=tree, ref=ref, edges=edges, q_last=q_last,
        q0=q0, verdict=verdict, verdict_ia=verdict_ia,
        status=verdict["ho"]["status"],
        extrapolation={calc: extrapolations(v) for calc, v in c0.items()},
        ceiling=ceiling,
        slope={calc: slope_radius(v, edges) for calc, v in c0.items()},
        shape=shape,
        m1=m1_compare(pick(rec, "FT2", "IA+MEC", "symbol"), edges["ho"]),
        m1_nominal_q=m1_compare([dict(r, q=nominal_q(r["q"])) for r in
                                 pick(rec, "FT2", "IA+MEC", "symbol")], edges["ho"]),
        c2=c2_compare(pick(rec, "C2", "IA+MEC"), edges["ho"], q0["ho"]),
        closures=closures(rec, tree),
        slot=dict(n=len(ref["slot"]), q_min=ref["slot"][0]["q"], q_max=ref["slot"][-1]["q"]))


def run(path=DATA, verbose=True):
    t0 = time.perf_counter()
    m = measure(path)
    v, q0, ql = m["verdict"], m["q0"], m["q_last"]
    ex, cl, sl, m1, c2 = (m["extrapolation"]["IA+MEC"], m["closures"],
                          m["slope"]["IA+MEC"], m["m1"], m["c2"])

    def qs(x):
        return "no zero below %g fm^-1" % Q_SCAN_MAX if x is None else "%.4f fm^-1" % x

    generator = ("HoSpin1FF C0 zero q0 = %s (shipped ho); vmc-ft: %s"
                 % (qs(q0["ho"]), qs(q0["vmc-ft"])))
    reference = ("Wiringa-Schiavilla PRL 81 (1998) 4317 Fig. 1, digitized from the "
                 "vector PDF: the IA+MEC C0 line is drawn non-zero to q_last = %.4f "
                 "fm^-1 (first zero > q_last; not located above it); IA to %.4f"
                 % (ql["IA+MEC"], ql["IA"]))
    tolerance = ("q0 >= q_last - (%.7f + %g) = %.4f fm^-1 (the plan's |q0 - q0_WS| <= "
                 "delta + %g on the zeros the figure allows; the upper side is not "
                 "tested; declared in the hashed design block)"
                 % (DELTA_FM, ZERO_MARGIN_FM, v["ho"]["bound"], ZERO_MARGIN_FM))
    subrows = []
    for k in ("ho", "vmc-ft"):
        subrows.append(dict(name="c0_zero[%s]" % k, status=v[k]["status"],
                            in_verdict=(k == "ho"), q0_tree=v[k]["q0_tree"],
                            q_last=v[k]["q_last"], bound=v[k]["bound"],
                            distance_above_bound=v[k]["distance_above_bound"],
                            reason=v[k]["reason"]))
    for k in ("ho", "vmc-ft"):
        w = m["verdict_ia"][k]
        subrows.append(dict(name="c0_zero_vs_IA[%s]" % k, status="recorded",
                            in_verdict=False, q0_tree=w["q0_tree"], q_last_ia=w["q_last"],
                            bound=w["bound"], rule_met=(w["status"] == "pass")))
    subrows.append(dict(name="zero_extrapolation[IA+MEC]", status="recorded",
                        in_verdict=False, n_fits_with_root=ex["n_with_root"],
                        min_root=ex["min_root"], max_root=ex["max_root"],
                        **{key + "_root": f["root"] for key, f in ex["fits"].items()}))
    for k in ("ho", "vmc-ft"):
        c = m["ceiling"][k]
        subrows.append(dict(name="c0_ceiling[%s]" % k, status="recorded", in_verdict=False,
                            max_tree_f2_above_q_last=c["max_f2"], at_q=c["at_q"],
                            floor_f2=FLOOR_F2, above_floor=c["above_floor"]))
    for calc in ("IA+MEC", "IA"):
        s_ = m["slope"][calc]
        subrows.append(dict(name="slope_radius[%s]" % calc, status="recorded",
                            in_verdict=False, r_ws=s_["ws"]["r"],
                            r_ws_err_dig=s_["ws"]["r_err_dig"],
                            r_ws_variants_min=s_["ws"]["r_min_variants"],
                            r_ws_variants_max=s_["ws"]["r_max_variants"],
                            r_ho_same_fit=s_["ho"]["r_same_fit"],
                            r_vmcft_same_fit=s_["vmc-ft"]["r_same_fit"],
                            r_ho_richardson=s_["ho"]["r_richardson"],
                            r_vmcft_richardson=s_["vmc-ft"]["r_richardson"]))
    for (k, calc), sh in sorted(m["shape"].items()):
        subrows.append(dict(name="c0_shape[%s vs %s]" % (k, calc), status="recorded",
                            in_verdict=False, chi2=sh["chi2"], ndf=sh["ndf"], s=sh["s"],
                            ratio_at_q1=sh["ratio_q1"]["ratio"],
                            ratio_at_q2=sh["ratio_q2"]["ratio"],
                            ratio_at_q25=sh["ratio_q25"]["ratio"],
                            ratio_at_q_last=sh["ratio_qlast"]["ratio"]))
    subrows.append(dict(name="m1_shape[fm vs IA+MEC F_T^2]", status="recorded",
                        in_verdict=False, chi2=m1["chi2"], ndf=m1["ndf"], p=m1["p"],
                        s_hat=m1["s_hat"], s_err=m1["s_err"], chi2_s1=m1["chi2_s1"],
                        ndf_s1=m1["ndf_s1"], p_s1=m1["p_s1"],
                        ws_zero_symbol_q=m1["ws_zero_symbol_q"],
                        tree_qz=m["tree"]["options"]["fm_qz_fm"]))
    subrows.append(dict(name="c2[fq vs IA+MEC C2]", status="recorded", in_verdict=False,
                        ratio_low_q_mean=c2["ratio_low_q_mean"],
                        ratio_low_q_min=c2["ratio_low_q_min"],
                        ratio_low_q_max=c2["ratio_low_q_max"],
                        ws_zero_lo=c2["ws_zero_bracket"][0], ws_zero_hi=c2["ws_zero_bracket"][1],
                        tree_c2_zero=c2["tree_c2_zero"],
                        ws_second_lobe_q=c2["ws_second_lobe_q"],
                        ws_second_lobe_f2=c2["ws_second_lobe_f2"],
                        tree_f2_at_ws_lobe=c2["tree_f2_at_ws_lobe"]))
    for key in ("mu[IA+MEC]", "mu[IA]"):
        c = cl[key]
        subrows.append(dict(name="closure_%s" % key, status=_word(c["holds"]),
                            in_verdict=False, mu=c["mu"], printed=c["printed"], rel=c["rel"]))
    for calc in ("IA+MEC", "IA"):
        c = cl["mu_nominal_q[%s]" % calc]
        subrows.append(dict(name="closure_mu[%s, nominal q, post hoc]" % calc,
                            status=_word(c["holds"]), in_verdict=False, post_hoc=True,
                            mu=c["mu"], printed=c["printed"], rel=c["rel"],
                            max_offset_from_grid=c["max_offset_from_grid"]))
    mn = m["m1_nominal_q"]
    subrows.append(dict(name="m1_shape[nominal q, post hoc]", status="recorded",
                        in_verdict=False, post_hoc=True, chi2=mn["chi2"], ndf=mn["ndf"],
                        p=mn["p"], s_hat=mn["s_hat"]))
    c = cl["Q[IA+MEC]"]
    subrows.append(dict(name="closure_Q[IA+MEC]", status=_word(c["holds"]), in_verdict=False,
                        abs_q_fm2=c["abs_q_fm2"], printed_fm2=c["printed_fm2"], rel=c["rel"],
                        times_sqrt2_over_3=c["ratio_sqrt2_over_3"]))
    c = cl["lines_vs_symbols[IA+MEC]"]
    subrows.append(dict(name="closure_lines_vs_symbols[IA+MEC]", status=_word(c["holds"]),
                        in_verdict=False, n=c["n"], max_abs_pull=c["max_abs_pull"]))
    c = cl["rfy_table1[Stanford F_T^2]"]
    subrows.append(dict(name="closure_rfy_table1[q]", status=_word(c["q_holds"]),
                        in_verdict=False, max_abs_dq=c["max_abs_dq"], limit=RFY_Q_TOL_FM))
    subrows.append(dict(name="closure_rfy_table1[ratio]", status=_word(c["constant"]),
                        in_verdict=False, ratio_mean=c["ratio_mean"],
                        ratio_min=c["ratio_min"], ratio_max=c["ratio_max"],
                        spread_limit=RFY_RATIO_SPREAD))
    import lipolgen as lg  # noqa: WPS433
    config = dict(lipolgen=getattr(lg, "__version__", None), python=sys.version.split()[0],
                  data_file_sha256=DATA_SHA256, rfy_file_sha256=RFY_SHA256,
                  design_block_sha256=m["design_block_sha256"],
                  ws_pdf_sha256=WS_PDF_SHA256, ws_page=WS_PAGE,
                  slot_rows=m["slot"]["n"], slot_q_min=m["slot"]["q_min"],
                  slot_q_max=m["slot"]["q_max"], m_a_fm=m["tree"]["m_a_fm"])
    report = dict(name=NAME, generator=generator, reference=reference,
                  tolerance=tolerance, status=m["status"], reason=None,
                  subrows=subrows, config=config, q_last=ql, q0=q0,
                  verdict=v, verdict_ia=m["verdict_ia"], extrapolation=m["extrapolation"],
                  ceiling=m["ceiling"], slope=m["slope"],
                  shape={"%s vs %s" % k: x for k, x in m["shape"].items()},
                  m1=m1, m1_nominal_q=m["m1_nominal_q"], c2=c2, closures=cl,
                  tree_options=m["tree"]["options"],
                  provenance=m["tree"]["provenance"])
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        print("reference: %d record rows, %d slot rows (q = %.4f .. %.4f fm^-1); "
              "design block %s..." % (len(m["ref"]["record"]), m["slot"]["n"],
                                      m["slot"]["q_min"], m["slot"]["q_max"],
                                      m["design_block_sha256"][:12]))
        print("WS98 C0 drawn non-zero to q_last = %.6f (IA+MEC), %.6f (IA) fm^-1"
              % (ql["IA+MEC"], ql["IA"]))
        for k in ("ho", "vmc-ft"):
            print("  %-6s q0 = %s; rule q0 >= %.6f -> %s" % (
                k, qs(q0[k]), v[k]["bound"], v[k]["status"].upper()))
        print("(r1) zero extrapolations of the IA+MEC branch: %s" % ", ".join(
            "%s %s" % (key, "none<4" if f["root"] is None else "%.4f" % f["root"])
            for key, f in ex["fits"].items()))
        for k in ("ho", "vmc-ft"):
            c = m["ceiling"][k]
            print("(r2) %-6s max (fc/Z)^2 on (q_last, 4] = %.3e at q = %.3f (floor 1e-7)"
                  % (k, c["max_f2"], c["at_q"]))
        for calc in ("IA+MEC", "IA"):
            s_ = m["slope"][calc]
            print("(r3) slope radius %s: WS %.4f +- %.4f (dig) [variants %.4f .. %.4f]; "
                  "same fit ho %.4f vmc-ft %.4f; Richardson ho %.5f vmc-ft %.5f fm"
                  % (calc, s_["ws"]["r"], s_["ws"]["r_err_dig"], s_["ws"]["r_min_variants"],
                     s_["ws"]["r_max_variants"], s_["ho"]["r_same_fit"],
                     s_["vmc-ft"]["r_same_fit"], s_["ho"]["r_richardson"],
                     s_["vmc-ft"]["r_richardson"]))
        for (k, calc), sh in sorted(m["shape"].items()):
            print("(r4) shape %-6s vs %-6s chi2/ndf = %.4g/%d, s = %.4f; WS/(s tree) at "
                  "q=%.2f: %.3f, q=%.2f: %.3f, q=%.2f: %.3f, q_last=%.3f: %.3f"
                  % (k, calc, sh["chi2"], sh["ndf"], sh["s"], sh["ratio_q1"]["q"],
                     sh["ratio_q1"]["ratio"], sh["ratio_q2"]["q"], sh["ratio_q2"]["ratio"],
                     sh["ratio_q25"]["q"], sh["ratio_q25"]["ratio"], sh["ratio_qlast"]["q"],
                     sh["ratio_qlast"]["ratio"]))
        print("(r5) F_M1: chi2 = %.3f/%d (p = %.3g), s = %.4f +- %.4f; s = 1: %.3f/%d "
              "(p = %.3g); WS F_T zero bracket %s, tree q_z = %.2f"
              % (m1["chi2"], m1["ndf"], m1["p"], m1["s_hat"], m1["s_err"], m1["chi2_s1"],
                 m1["ndf_s1"], m1["p_s1"], "/".join("%.4f" % x for x in m1["ws_zero_bracket"]),
                 m["tree"]["options"]["fm_qz_fm"]))
        print("(r6) C2: |WS/tree| at q <= 0.5: mean %.3f (%.3f .. %.3f); WS zero in "
              "(%.4f, %.4f), tree %s; WS 2nd lobe %.3e at %.3f (tree %.3e)"
              % (c2["ratio_low_q_mean"], c2["ratio_low_q_min"], c2["ratio_low_q_max"],
                 c2["ws_zero_bracket"][0], c2["ws_zero_bracket"][1], qs(c2["tree_c2_zero"]),
                 c2["ws_second_lobe_f2"], c2["ws_second_lobe_q"], c2["tree_f2_at_ws_lobe"]))
        for key in ("mu[IA+MEC]", "mu[IA]"):
            c = cl[key]
            print("(r7) mu %s: %.4f vs printed %.3f (%+.2f %%) -> %s"
                  % (key, c["mu"], c["printed"], 100 * c["rel"], "HOLDS" if c["holds"] else "DOES NOT HOLD"))
        for key in ("mu_nominal_q[IA+MEC]", "mu_nominal_q[IA]"):
            c = cl[key]
            print("     POST HOC %s: %.4f vs printed %.3f (%+.2f %%) -> %s"
                  % (key, c["mu"], c["printed"], 100 * c["rel"],
                     "within the rule" if c["holds"] else "outside the rule"))
        mn = m["m1_nominal_q"]
        print("     POST HOC F_M1 at nominal q: chi2 = %.3f/%d (p = %.3g), s = %.4f"
              % (mn["chi2"], mn["ndf"], mn["p"], mn["s_hat"]))
        c = cl["Q[IA+MEC]"]
        print("(r7) Q: |Q| = %.4f fm^2 vs printed 0.23 (%+.1f %%; x sqrt2/3 = %.4f) -> %s"
              % (c["abs_q_fm2"], 100 * c["rel"], c["ratio_sqrt2_over_3"],
                 "HOLDS" if c["holds"] else "DOES NOT HOLD"))
        c = cl["lines_vs_symbols[IA+MEC]"]
        print("(r7) lines vs symbols: %d points, max |pull| %.2f -> %s"
              % (c["n"], c["max_abs_pull"], "HOLDS" if c["holds"] else "DOES NOT HOLD"))
        c = cl["rfy_table1[Stanford F_T^2]"]
        print("(r7) RFY Table I: max |dq| = %.4f fm^-1 -> %s; ratio %.4f (%.4f .. %.4f) -> %s"
              % (c["max_abs_dq"], "HOLDS" if c["q_holds"] else "DOES NOT HOLD",
                 c["ratio_mean"], c["ratio_min"], c["ratio_max"],
                 "constant" if c["constant"] else "NOT constant"))
        for r in subrows:
            vals = ", ".join("%s = %s" % (k, ("%.6g" % x) if isinstance(x, float) else x)
                             for k, x in r.items() if k not in ("name", "status", "in_verdict"))
            print("SUBROW | %s[%s] | %s | %s" % (NAME, r["name"], vals,
                                                 str(r["status"]).upper()))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, generator, reference, tolerance, m["status"].upper()))
    return report


# ---------------------------------------------------------------- figure

def plot(report, out_dir=None):
    """Two panels sharing q.  Top: |F_C0|^2 per nucleus on a log axis -- WS98's
    IA+MEC C0 line (ink, the reference; it ends where it meets the figure's
    1e-7 floor) with its C2 line for context (grey), against the two C0 edges
    of the tree, (fc/Z)^2, dashed (`ho` the shipped default's slot, `vmc-ft`
    the alternative's); the zeros the figure allows, q >= q_last - (delta +
    0.1), shaded.  Bottom: F_T^2 -- WS98's IA+MEC symbols with their Monte
    Carlo bars (ink) against the tree's (4/3) eta (1 + eta) fm^2/Z^2 times the
    record's free scale (dashed).  The curves are re-evaluated on a q grid for
    drawing only; every number the figure states is the report's."""
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    from _plots import BAND, INK, MUTED, mpl, save, status_title, style_series
    plt = mpl()
    tree = tree_side()
    ref = load_reference(DATA, tree["m_a_fm"], tree["z"])
    rec = ref["record"]
    edges = {k: Edge(v, tree) for k, v in tree["edges"].items()}
    fig, (ax, axt) = plt.subplots(2, 1, figsize=(6.4, 5.8), sharex=True,
                                  layout="constrained",
                                  gridspec_kw={"height_ratios": [3, 2]})
    fig.get_layout_engine().set(rect=(0.0, 0.03, 1.0, 0.905), hspace=0.02)
    v = report["verdict"]["ho"]
    q_last = report["q_last"]["IA+MEC"]
    ax.axvspan(v["bound"], 4.0, color=BAND, alpha=0.8, linewidth=0, zorder=1,
               label="zeros the figure allows: q >= q_last - (delta + 0.1) = %.4f"
                     % v["bound"])
    ax.axhline(FLOOR_F2, color=MUTED, linewidth=0.8, linestyle=":", zorder=2)
    ax.text(0.05, FLOOR_F2 * 1.3, "WS98 Fig. 1 floor (1e-7)", fontsize=7,
            color=MUTED, va="bottom")
    c0 = pick(rec, "C0", "IA+MEC")
    ax.plot([r["q"] for r in c0], [r["F2"] for r in c0], color=INK, linewidth=1.6,
            zorder=5, label="WS98 IA+MEC C0 (drawn to q_last = %.4f fm^-1)" % q_last)
    c2 = pick(rec, "C2", "IA+MEC")
    for br in (0, 1):
        pts = [r for r in c2 if r["branch"] == br]
        ax.plot([r["q"] for r in pts], [r["F2"] for r in pts], color=MUTED,
                linewidth=1.0, zorder=4, label="WS98 IA+MEC C2 (context)" if br == 0 else None)
    qs = [0.005 * i for i in range(1, 801)]
    for k, slot in (("ho", "default"), ("vmc-ft", "alt")):
        q0 = report["q0"][k]
        grid = sorted(set(qs) | ({q0} if q0 is not None else set()))
        kw = style_series(slot, label="HoSpin1FF %s: zero %s" % (
            k, "%.4f fm^-1" % q0 if q0 is not None else "none"))
        ax.plot(grid, [max(edges[k].c0(q) ** 2, 1e-12) for q in grid],
                linestyle="--", linewidth=1.6, **kw)
    ax.set_yscale("log")
    ax.set_ylim(1e-9, 2.0)
    ax.set_xlim(0.0, 4.0)
    ax.set_ylabel("|F_C0(q)|^2 per nucleus (F_C0(0) = 1)")
    ax.legend(loc="upper right", fontsize=7)
    sym = pick(rec, "FT2", "IA+MEC", "symbol")
    up = [0.0 if r["F2_hi"] is None else r["F2_hi"] - r["F2"] for r in sym]
    lo = [0.0 if r["F2_lo"] is None else min(r["F2"] - r["F2_lo"], 0.99 * r["F2"])
          for r in sym]
    axt.errorbar([r["q"] for r in sym], [r["F2"] for r in sym], yerr=[lo, up],
                 **style_series("reference", label="WS98 IA+MEC F_T^2 (MC errors)"))
    m1 = report["m1"]
    axt.plot(qs, [m1["s_hat"] * edges["ho"].ft2(q) for q in qs], linestyle="--",
             linewidth=1.6, **style_series("default", label=(
                 "HoSpin1FF F_T^2 (q_z = 1.30, b = 1.85; both edges) x %.3f"
                 % m1["s_hat"])))
    axt.set_yscale("log")
    axt.set_ylim(1e-10, 1e-3)
    axt.set_xlabel("q [fm^-1]")
    axt.set_ylabel("F_T^2 = B/Z^2")
    axt.legend(loc="lower left", fontsize=7)
    q0 = report["q0"]["ho"]
    status_title(fig, NAME, report["status"],
                 "C0 zero: ho %s vs WS98 non-zero to %.4f fm^-1; PASS iff q0 >= %.4f "
                 "(one-sided); F_T chi2 %.1f/%d"
                 % ("%.4f" % q0 if q0 is not None else "none", q_last, v["bound"],
                    m1["chi2"], m1["ndf"]))
    return save(fig, NAME, out_dir)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if "--redigitize" in argv:
        pdf = argv[argv.index("--redigitize") + 1]
        new = redigitize(pdf)
        old = body_of_file()
        if new == old:
            print("redigitize: %s re-read -> the vendored body, identical (%d lines)"
                  % (pdf, new.count("\n")))
            return 0
        a, b = old.splitlines(), new.splitlines()
        diff = [i for i in range(max(len(a), len(b)))
                if i >= len(a) or i >= len(b) or a[i] != b[i]]
        print("redigitize: %d line(s) differ; first at body line %d:\n  file: %s\n  pdf:  %s"
              % (len(diff), diff[0] + 1, a[diff[0]] if diff[0] < len(a) else "<none>",
                 b[diff[0]] if diff[0] < len(b) else "<none>"))
        return 1
    rep = run()
    if "--plot" in argv:
        print("figure: %s" % plot(rep))
    return 0


if __name__ == "__main__":
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    # `--plot` also writes the row's figure (plot()); the same rule holds.
    # `--redigitize PDF`: 0 identical, 1 different, 2 broken.
    try:
        sys.exit(main())
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        sys.exit(2)
