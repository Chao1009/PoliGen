#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""docs/benchmarking/02_data_nucleon_deuteron.md sec. 5.4 -- the d(e,e'p)n
reduced cross section of Boeglin et al. 2011 against the tree's deuteron
momentum distributions n(k) (tier T2).

New 2026-10-26 (run 2026-10-26 phase C5, PLAN.md "C5"; under
docs/benchmarking/08_channel_priority.md sec. 3.7, the tagged d-p control,
whose development order starts with "d(e,e'p)n momentum distributions
(Boeglin), with the paper's own 10-20 % FSI allowance declared", and sec. 4
ingredient G, the deuteron wave function default).

THE DOCUMENT.  W. U. Boeglin et al. (Hall A), "Probing the high momentum
component of the deuteron at high Q^2", PRL 107 (2011) 262501,
arXiv:1106.0275v3: 2H(e,e'p)n at Q^2 = 3.5 (GeV/c)^2.  Fig. 2 (PDF page 5)
is a VECTOR drawing: sigma_red = sigma_exp/(k sigma_cc1) [fm^3] against the
missing momentum p_m at theta_nq = 35, 45, 75 deg (+-5 deg), data with error
bars, and the MS model's CD-Bonn and Paris curves (PWIA and FSI), each a
polyline with one vertex at every data point's p_m.  All of it was read back
from the PDF (sha256 b9abd413...9c70f91, the author's corpus copy, NOT in this
repository) into data/boeglin2011_deuteron_pmiss.csv: 35 rows, every number
with the PDF coordinates it is read from.  The file's provenance block gives
the frames, the calibration and its check, the convention and the commands;
`--redigitize PDF` re-reads the PDF (dev time, PyMuPDF) and compares.  The
sha256 of every byte is DATA_SHA256 below, re-checked on every read; every
number is re-derived from its coordinates.

WHAT THE PAPER SAYS ABOUT THE WINDOW (p. 2, 4), stated before the design:
"all predict small FSI contributions (10 - 20%) for 35 <= theta_nq <= 45";
"for 0.2 <= p_m <= 0.5 GeV/c on average MS predicts between -18% and -5%
FSI, JM between -9% and +5% and JVO between -12% and -16% FSI.  This
indicates a kinematic region where the cross section is dominated by PWIA
and should therefore reflect directly the deuteron momentum distribution";
the PWIA calculations "differ at larger p_m due to differences in the wave
functions and in details of the off-shell electron nucleon interaction
used"; sigma_red divides data and theory alike by k sigma_cc1, the de Forest
CC1 off-shell cross section.  In PWIA, sigma_red(p_m) = n(p_m), the proton
momentum distribution in the deuteron (INT n d^3p = 1, fm^3).

THE TREE SIDE.  n(k) of the three deuteron waves this tree ships, at
k = p_m (the paper's own PWIA identification; instant form, which is also
how the generator places the spectator), in fm^3:
  hulthen  the analytic Hulthen pair of deuteron_channel() (beta = 0.30 GeV,
           P_D = 0.045) -- the tagged channel's DEFAULT;
  av18     deuteron_channel(source = ClusterWaveSource.VmcAV18);
  cdbonn   cdbonn_wave() (Machleidt 2001 App. D), n = (psi_s^2 + psi_d^2)/(4 pi)
           with INT (psi_s^2 + psi_d^2) p^2 dp = 1;
the first two as the M-averaged TaggedModel(...).n_of_kc(M, k, 0) (isotropic
-- t2_deeps_spectator_tail.py checks it) times (hbar c)^3 (t2_deeps's Tree).

WHAT IS COMPARED, AND THE TOLERANCE -- fixed in this header on 2026-10-27
(UTC) before this harness computed any comparison (PLAN.md ground rule 4).
The block between the two marker lines is hashed (DESIGN_BLOCK_SHA256
below); `design_block_sha256()` re-checks it, and the numerical rules
against its sentences, on every run; the harness exits 2 if either moved.
==== DESIGN BLOCK BEGIN ====
THE VERDICT, one per wave; the row's status is the DEFAULT wave's
(hulthen), av18 and cdbonn are sub-rows by the same rule.  The points: the
22 of panels (a) theta_nq = 35 deg (10) and (b) 45 deg (12), p_m = 0.055 ..
0.538 GeV/c -- the paper's 35-45 deg window up to 0.55 GeV/c; panel (c),
75 deg, is FSI-dominated by the paper's own account and is NOT compared.
d_i = sigma_red as plotted; sigma_bar_i = half the bar, (hi - lo)/2 (the
bars are symmetric in sigma_red).  The model m_i = n_wave(k = p_m_i) in
fm^3, times ONE free scale s common to both panels (the paper: "The PWIA
results are for all angles identical"); s absorbs the 4.5 % overall
normalisation and any constant FSI or off-shell factor.
THE ALLOWANCE BAND, declared before any comparison: A = 0.20 -- the upper
end of the paper's own "10 - 20%" FSI estimate for this window -- standing
for FSI AND the de Forest CC1 off-shell prescription together (the paper
gives no separate number for the off-shell spread, so none is added beyond
A), as a relative error on each point, uncorrelated, in quadrature:
    sigma_i^2 = sigma_bar_i^2 + (A d_i)^2 + (ln 10 x 0.000307 d_i)^2,
the last term the digitization error (0.01 pt).  s in closed form,
s = sum(d m/sigma^2)/sum(m^2/sigma^2); chi2 = sum((d - s m)/sigma)^2;
ndf = 22 - 1 = 21.
TOLERANCE: PASS iff p(chi2_min, ndf = 21) >= 0.01.  p is scipy's chi2.sf
when scipy is importable, else a regularized-gamma fallback.
RECORDED, NEVER IN THE VERDICT -- each rule fixed here, before computing:
 (r1) the absolute closure of the reading: the paper's CD-Bonn PWIA curve
      over the tree's cdbonn n(p_m) (fm^3, normalised to 1, no free scale)
      at the 22 points: mean, min, max.  Near 1 means sigma_red = n(p_m) in
      PWIA is read right; its p_m dependence is the MS model's off-shell
      cross section against CC1 (and the bins' kinematics).  And the
      paper's Paris PWIA over its CD-Bonn PWIA: mean, min, max (the
      wave-function spread the paper speaks of).
 (r2) theory against theory: each wave against the paper's CD-Bonn PWIA
      curve, one free scale fitted in ln (uniform weights): the min and max
      of curve/(s m) over the 22 points.
 (r3) the verdict without the allowance (A = 0): chi2, p, s per wave.
 (r4) per panel, with the allowance: (a) alone (ndf = 9) and (b) alone
      (ndf = 11), each with its own scale: chi2, p, s per wave.
 (r5) the MS model's own FSI in the window: cdbonn_fsi/cdbonn_pwia over the
      22 points, min and max, against the declared A.
 (r6) the ratio route: d_i over the paper's CD-Bonn PWIA, per panel mean,
      min, max; and the verdict's chi2 with the model replaced by
      s m_i (cdbonn_fsi_i/cdbonn_pwia_i) -- the MS FSI factor applied to the
      tree's PWIA -- and A = 0: chi2, p, s per wave.
 (r7) panel (c), 75 deg: its 13 points are vendored and drawn, never
      compared.
WHAT WAS SEEN BEFORE THIS BLOCK WAS WRITTEN (by the worker, 2026-10-27): the
page rendered at 150 dpi; every drawing operator, and the extraction's
printout of all 35 points with their bars and the four curves at each --
so the data/PWIA ratios were in view (at 35 deg about 0.65-0.85 of the
CD-Bonn PWIA curve over the whole range; at 45 deg about 0.6-0.9 up to
0.46 GeV/c and 1.2-1.3 at the last two points); on the tree side, the
hulthen and av18 n(k) at k = 0.05, 0.1, 0.2, 0.3, 0.4, 0.5 GeV (printed:
343, 39.6, 2.26, 0.309, 0.0627, 0.0152 and 376, 38.9, 1.71, 0.196, 0.0520,
0.0214 GeV^-3), their normalisation integrals (1.085, 1.088), the cdbonn
normalisation (1.0000) and (psi_s^2 + psi_d^2)(0.25 fm^-1) = 32.44; and the
committed Deeps row's finding that the hulthen tail falls too fast at
p_s = 0.30-0.46 GeV/c.  No chi2, scale or ratio of the tree to these data
was computed before this block.  The allowance and the tolerance are the
task's (the paper's 10-20 %, its upper end; p >= 0.01, this run's standard
for a shape chi2), fixed before the data were read.
==== DESIGN BLOCK END ====

POST-HOC ADDENDUM (2026-10-27, written after the first run of 00:45:08Z and
OUTSIDE the hashed block, which is unchanged).  The first figure showed the
hulthen and av18 lines stepping at low k.  That is the tree's model as the
block declared it: TaggedModel.n_of_kc looks up the cell AT OR BELOW k
(tagged.hpp: "the grid IS part of the model"; k grid 0.0001 .. 1.2 GeV, 280
nodes, 0.0043 GeV apart), and the table value of a cell is the density at
its lower node.  At the data's p_m the cell value is 1.0051 .. 1.1605
(hulthen) and 1.0032 .. 1.1808 (av18) times the same table read smoothly
(linear in ln n between node values) -- 1.16 / 1.18 at p_m = 0.055 GeV/c,
1.08 / 1.09 at 0.0925.  A DIAGNOSTIC, never a verdict: the verdict's chi2
with the smooth reading gives hulthen 23.42 / 21 (p = 0.32), av18 14.61
(p = 0.84), cdbonn unchanged (analytic) -- the verdict does not turn on
the grid.  The figure draws every wave through its node values.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE (unpolarized).
  2. isoscalar denominator -- NOT APPLICABLE; its analogue, the struck
     proton's off-shell e-p cross section, is divided out by the paper's own
     reduction (k sigma_cc1) on both the data and its theory.
  3. per nucleus / per nucleon -- sigma_red is PER DEUTERON (one proton per
     deuteron, n normalised to 1); the tree's n(k) is per deuteron.  Only
     the shape is in the verdict (one free scale); the absolute
     normalisation is record (r1).

WHAT THIS DOES NOT VALIDATE.  The tagged spectrum at EIC kinematics (a
fixed-target exclusive reaction at Q^2 = 3.5 GeV^2, the neutron
unobserved); the D-state separately (n sums S and D); the polarized or
tensor structure (A_zz^tag); the LC versus instant-form question at these
momenta (the instant form is the paper's own reduction); FSI itself (an
allowance, not a calculation -- the tree's GlauberFsiWeight is not used
here); and p_m above 0.55 GeV/c.

Run:  source env.sh && python3 validation/benchmarks/t2_boeglin_nk.py [--plot]
      (dev time: ... --redigitize PATH/TO/1106.0275.pdf)

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
DATA = os.path.join(HERE, "data", "boeglin2011_deuteron_pmiss.csv")
#: sha256 of EVERY byte of the vendored file, as written 2026-10-27T00:42:49Z;
#: re-checked on every read
DATA_SHA256 = "07a6fedecabeadb5124b2dcd3bdce0865e4986151448cad48b13066b9b653334"
DESIGN_BEGIN = "==== DESIGN BLOCK BEGIN ===="
DESIGN_END = "==== DESIGN BLOCK END ===="

NAME = "t2_boeglin_nk"
#: a library evaluation on a vendored table, well under 30 s, no event sample
EXPENSIVE = False

# ------------------------------------------------------------ the digitization
#: the PDF the numbers were read from (the author's corpus copy of
#: arXiv:1106.0275v3); NOT in this repository -- `--redigitize PDF` re-reads
#: it (dev time, PyMuPDF) and compares with the vendored body
BG_PDF_SHA256 = "b9abd4132e4736ac195d608d43354d3774c2ca233a587495adb1e31395c70f91"
BG_PAGE = 5                         # 1-based: Fig. 2
#: the three plot boxes (x0, y0, x1, y1) in PDF points, y down, as the four
#: 0.64-pt strokes of each frame draw them; the theta_nq of each panel
FRAMES = {"a": (114.973, 63.567, 251.737, 258.945),
          "b": (251.737, 63.567, 388.502, 258.945),
          "c": (388.502, 63.567, 525.267, 258.945)}
THETA_NQ = {"a": 35, "b": 45, "c": 75}
#: the axis values at the box edges: p_m 0 .. 0.7 GeV/c (ticks 0.1 .. 0.6
#: every 19.538 pt) and log10 sigma_red -5 .. +1 (fm^3; decade ticks every
#: 32.563 pt, minor ticks at log10 2 .. 9)
PM_RANGE = (0.0, 0.7)
LOG_RANGE = (-5.0, 1.0)
#: the coordinates are float32 output of a vector plotting program, exact to
#: ~1e-4 pt and written here to 0.001 pt; 0.01 pt is the declared
#: digitization error of every coordinate
DIG_PT = 0.01
DPM_GEV = DIG_PT * (PM_RANGE[1] - PM_RANGE[0]) / (FRAMES["a"][2] - FRAMES["a"][0])
DLOG10 = DIG_PT * (LOG_RANGE[1] - LOG_RANGE[0]) / (FRAMES["a"][3] - FRAMES["a"][1])

BLACK, BLUE, GREEN = (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.8, 0.0)
#: (colour, dash) -> curve; the caption: "CD-Bonn potential: dashed (blue)
#: lines PWIA, solid (blue) lines FSI.  Paris potential: dash-dot (green)
#: lines PWIA, long-dashed (green) lines FSI" (a [.63849 0] dash is solid)
CURVES = {(BLUE, "[ 3.19245 .63849 ] 0"): "cdbonn_pwia",
          (BLUE, "[ .63849 0 ] 0"): "cdbonn_fsi",
          (GREEN, "[ 3.19245 .63849 1.27698 .63849 ] 0"): "paris_pwia",
          (GREEN, "[ 6.3849 1.9154701 ] 0"): "paris_fsi"}
CURVE_KEYS = ("cdbonn_pwia", "cdbonn_fsi", "paris_pwia", "paris_fsi")
COLUMNS = (["panel", "theta_nq_deg", "point", "pm_gev", "sigma_red_fm3",
            "sigma_red_lo_fm3", "sigma_red_hi_fm3"] + [k + "_fm3" for k in CURVE_KEYS]
           + ["x_pt", "y_pt", "y_bar_lo_pt", "y_bar_hi_pt"]
           + ["y_%s_pt" % k for k in CURVE_KEYS] + ["dpm_gev", "dlog10"])
EXPECTED_POINTS = {"a": 10, "b": 12, "c": 13}


def to_pm(x_pt, panel):
    x0, _, x1, _ = FRAMES[panel]
    return PM_RANGE[0] + (x_pt - x0) / (x1 - x0) * (PM_RANGE[1] - PM_RANGE[0])


def to_log10(y_pt, panel):
    _, y0, _, y1 = FRAMES[panel]
    return LOG_RANGE[0] + (y1 - y_pt) / (y1 - y0) * (LOG_RANGE[1] - LOG_RANGE[0])


def to_sigma(y_pt, panel):
    return 10.0 ** to_log10(y_pt, panel)


def _r3(v):
    return round(float(v), 3)


def _col(c):
    return None if c is None else tuple(round(float(v), 2) for v in c)


def _shape(d):
    return "".join(it[0] for it in d["items"])


def _inside(x, y, frame):
    x0, y0, x1, y1 = frame
    return x0 < x < x1 and y0 < y < y1


def digitize(pdf):
    """Every data point (with its error bar) and the four MS curves' values
    at it, per panel (DEV TIME: PyMuPDF and the PDF)."""
    import pymupdf  # noqa: WPS433 (dev-time only)
    with open(pdf, "rb") as f:
        got = hashlib.sha256(f.read()).hexdigest()
    if got != BG_PDF_SHA256:
        raise RuntimeError("%s: sha256 %s is not the arXiv 1106.0275v3 PDF this "
                           "table was read from (%s)" % (pdf, got, BG_PDF_SHA256))
    ds = pymupdf.open(pdf)[BG_PAGE - 1].get_drawings()
    # the frames: four 0.64-pt strokes each
    for panel, (x0, y0, x1, y1) in FRAMES.items():
        edges = set()
        for d in ds:
            if d.get("type") == "s" and _col(d.get("color")) == BLACK and \
                    abs((d.get("width") or 0) - 0.64) < 0.01 and _shape(d) == "l":
                it = d["items"][0]
                edges.add((_r3(min(it[1].x, it[2].x)), _r3(min(it[1].y, it[2].y)),
                           _r3(max(it[1].x, it[2].x)), _r3(max(it[1].y, it[2].y))))
        need = {(x0, y0, x1, y0), (x0, y1, x1, y1), (x0, y0, x0, y1), (x1, y0, x1, y1)}
        if not need <= edges:
            raise RuntimeError("frame %s is not drawn at the declared box" % panel)
    rows = []
    for panel, frame in FRAMES.items():
        circles = []
        for i, d in enumerate(ds):
            r = d["rect"]
            cx, cy = 0.5 * (r.x0 + r.x1), 0.5 * (r.y0 + r.y1)
            if (d.get("type") == "fs" and _col(d.get("color")) == BLACK
                    and _col(d.get("fill")) == BLACK and _shape(d) == "cccccccc"
                    and _inside(cx, cy, frame)):
                circles.append((i, _r3(cx), _r3(cy)))
        bars = []
        for i, d in enumerate(ds):
            if (d.get("type") == "s" and _col(d.get("color")) == BLACK
                    and abs((d.get("width") or 0) - 1.28) < 0.01 and _shape(d) == "l"):
                it = d["items"][0]
                if abs(it[1].x - it[2].x) < 1e-6 and _inside(it[1].x, 0.5 * (it[1].y + it[2].y), frame):
                    bars.append((_r3(it[1].x), _r3(max(it[1].y, it[2].y)),
                                 _r3(min(it[1].y, it[2].y))))
        curves = {}
        for d in ds:
            key = CURVES.get((_col(d.get("color")), d.get("dashes")))
            if key is None or d.get("type") != "s" or len(d["items"]) < 3:
                continue                       # (the legend handles are 1 segment)
            pts = []
            for it in d["items"]:
                for p in (it[1], it[2]):
                    q = (_r3(p.x), _r3(p.y))
                    if not pts or pts[-1] != q:
                        pts.append(q)
            if all(_inside(x, y, frame) for x, y in pts):
                if key in curves:
                    raise RuntimeError("panel %s: two %s curves" % (panel, key))
                curves[key] = pts
        if set(curves) != set(CURVE_KEYS):
            raise RuntimeError("panel %s: curves %r" % (panel, sorted(curves)))
        circles.sort(key=lambda c: c[1])
        for n, (i, x, y) in enumerate(circles, 1):
            mb = [b for b in bars if abs(b[0] - x) < 0.01]
            if len(mb) != 1:
                raise RuntimeError("panel %s point %d: %d bars at x = %.3f" % (panel, n, len(mb), x))
            cy = {}
            for k, pts in curves.items():
                hit = [p for p in pts if abs(p[0] - x) < 0.01]
                if len(hit) != 1:
                    raise RuntimeError("panel %s point %d: curve %s has %d vertices "
                                       "at its x" % (panel, n, k, len(hit)))
                cy[k] = hit[0][1]
            rows.append(dict(panel=panel, point=n, drawing=i, x=x, y=y,
                             ylo=mb[0][1], yhi=mb[0][2], curves=cy))
        if any(len(pts) != len(circles) for pts in curves.values()):
            raise RuntimeError("panel %s: a curve has a vertex off the data's p_m" % panel)
    return rows


def body_text(rows):
    """The vendored file's body (the column line and every row), exactly as
    written."""
    buf = io.StringIO()
    buf.write(",".join(COLUMNS) + "\n")
    for r in rows:
        p = r["panel"]
        buf.write(",".join(
            [p, "%d" % THETA_NQ[p], "%d" % r["point"], "%.5f" % to_pm(r["x"], p),
             "%.6g" % to_sigma(r["y"], p), "%.6g" % to_sigma(r["ylo"], p),
             "%.6g" % to_sigma(r["yhi"], p)]
            + ["%.6g" % to_sigma(r["curves"][k], p) for k in CURVE_KEYS]
            + ["%.3f" % r["x"], "%.3f" % r["y"], "%.3f" % r["ylo"], "%.3f" % r["yhi"]]
            + ["%.3f" % r["curves"][k] for k in CURVE_KEYS]
            + ["%.6f" % DPM_GEV, "%.6f" % DLOG10]) + "\n")
    return buf.getvalue()


def redigitize(pdf):
    """DEV TIME: re-read the PDF and return the body this harness would vendor."""
    return body_text(digitize(pdf))


# ------------------------------------------------------------ the design block

#: sha256 of this header's DESIGN BLOCK (its two marker lines included), as
#: written at 2026-10-27T00:43:44Z, before this harness computed any
#: comparison.  `design_block_sha256()` re-checks it, and the rule constants
#: below against the block's own sentences, on every run.
DESIGN_BLOCK_SHA256 = (
    "437eab47375bd87d96e194c1b5706dd87a52f7b116101776e7adf4920089625d")
#: the verdict's rules (design block)
ALLOWANCE = 0.20
P_MIN = 0.01
NDF_VERDICT = 21
WINDOW_PANELS = ("a", "b")
DIG_LOG10 = 0.000307
WAVES = ("hulthen", "av18", "cdbonn")
DEFAULT_WAVE = "hulthen"
WAVE_LABEL = {"hulthen": "Hulthen pair, beta = 0.30 GeV, P_D = 0.045 (deuteron_channel() default)",
              "av18": "AV18 (deuteron_channel(source=VmcAV18))",
              "cdbonn": "CD-Bonn (cdbonn_wave(), Machleidt 2001 App. D)"}

_RULES = (
    (r"THE ALLOWANCE BAND, declared before any comparison: A = ([0-9.]+) --", "ALLOWANCE"),
    (r"PASS iff p\(chi2_min, ndf = (\d+)\) >= ([0-9.]+)\.", "TOL"),
    (r"the last term the digitization error \(0\.01 pt\)", None),
    (r"\(ln 10 x ([0-9.]+) d_i\)\^2", "DIG_LOG10"),
    (r"ndf = 22 - 1 = 21\.", None),
)


def design_block(source=None):
    if source is None:
        with open(os.path.abspath(__file__), encoding="utf-8") as f:
            source = f.read()
    i0 = source.index(DESIGN_BEGIN)
    i1 = source.index(DESIGN_END) + len(DESIGN_END)
    return source[i0:i1]


def rules_from_design_block(block):
    got = {}
    for rx, key in _RULES:
        m = re.search(rx, block)
        if m is None:
            raise RuntimeError("the design block no longer states %r" % rx)
        if key is not None:
            got[key] = m.groups()
    return dict(ALLOWANCE=float(got["ALLOWANCE"][0]), NDF_VERDICT=int(got["TOL"][0]),
                P_MIN=float(got["TOL"][1]), DIG_LOG10=float(got["DIG_LOG10"][0]))


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
    used = dict(ALLOWANCE=ALLOWANCE, NDF_VERDICT=NDF_VERDICT, P_MIN=P_MIN,
                DIG_LOG10=DIG_LOG10)
    for k, v in used.items():
        if v != stated[k]:
            raise RuntimeError("rule %s = %r differs from the hashed design "
                               "block's %r" % (k, v, stated[k]))
    if abs(DIG_LOG10 - DLOG10) > 5e-7:
        raise RuntimeError("the design block's digitization error is not 0.01 pt")
    return got


# ------------------------------------------------------------ the reference

def read_vendored(path=DATA, sha256=None):
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    want = DATA_SHA256 if sha256 is None else sha256
    if got != want:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s" % (path, got, want))
    return raw


def _close(a, b, rel):
    return abs(a - b) <= rel * max(abs(a), abs(b))


def load_reference(path=DATA):
    """The 35 rows, checked: sha256 of every byte, the columns, the points
    per panel, and every number re-derived from its PDF coordinates."""
    raw = read_vendored(path).decode("utf-8")
    lines = [l for l in raw.splitlines() if l.strip() and not l.startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames != COLUMNS:
        raise RuntimeError("%s: columns %r, expected %r" % (path, reader.fieldnames, COLUMNS))
    rows, counts = [], {}
    for r in reader:
        p = r["panel"]
        if p not in FRAMES or int(r["theta_nq_deg"]) != THETA_NQ[p]:
            raise RuntimeError("%s: panel %r / theta %r" % (path, p, r["theta_nq_deg"]))
        counts[p] = counts.get(p, 0) + 1
        row = dict(panel=p, theta_nq=THETA_NQ[p], point=int(r["point"]),
                   pm=float(r["pm_gev"]), d=float(r["sigma_red_fm3"]),
                   lo=float(r["sigma_red_lo_fm3"]), hi=float(r["sigma_red_hi_fm3"]),
                   x=float(r["x_pt"]), y=float(r["y_pt"]),
                   ylo=float(r["y_bar_lo_pt"]), yhi=float(r["y_bar_hi_pt"]))
        for k in CURVE_KEYS:
            row[k] = float(r[k + "_fm3"])
            if not _close(to_sigma(float(r["y_%s_pt" % k]), p), row[k], 6e-6):
                raise RuntimeError("%s: %s at %s/%d is not its coordinate" % (path, k, p, row["point"]))
        if abs(to_pm(row["x"], p) - row["pm"]) > 6e-6:
            raise RuntimeError("%s: p_m %s/%d is not x_pt calibrated" % (path, p, row["point"]))
        for yk, vk in (("y", "d"), ("ylo", "lo"), ("yhi", "hi")):
            if not _close(to_sigma(row[yk], p), row[vk], 6e-6):
                raise RuntimeError("%s: %s at %s/%d is not %s calibrated" % (path, vk, p, row["point"], yk))
        if abs(float(r["dpm_gev"]) - DPM_GEV) > 1e-6 or abs(float(r["dlog10"]) - DLOG10) > 1e-6:
            raise RuntimeError("%s: a row's digitization error is not 0.01 pt" % path)
        rows.append(row)
    if counts != EXPECTED_POINTS:
        raise RuntimeError("%s: points per panel %r, expected %r" % (path, counts, EXPECTED_POINTS))
    for p in FRAMES:
        pts = [r for r in rows if r["panel"] == p]
        if [r["point"] for r in pts] != list(range(1, len(pts) + 1)) or \
                any(b["pm"] <= a["pm"] for a, b in zip(pts, pts[1:])):
            raise RuntimeError("%s: panel %s not in increasing p_m" % (path, p))
    return rows


def body_of_file(path=DATA):
    raw = read_vendored(path).decode("utf-8")
    return "".join(l + "\n" for l in raw.splitlines() if not l.startswith("#"))


# ------------------------------------------------------------ the tree side

class Tree:
    """n(k) of the three waves, in fm^3 (k in GeV/c)."""

    def __init__(self):
        import lipolgen as lg  # noqa: WPS433 (env.sh puts it on the path)
        self.lg = lg
        self.hbarc = lg.HBARC_GEV_FM
        self.models = {"hulthen": lg.TaggedModel(lg.deuteron_channel()),
                       "av18": lg.TaggedModel(lg.deuteron_channel(
                           source=lg.ClusterWaveSource.VmcAV18))}
        self.cdb = lg.cdbonn_wave()

    def n(self, wave, k):
        if wave == "cdbonn":
            p = k / self.hbarc
            return (self.cdb.psi_s(p) ** 2 + self.cdb.psi_d(p) ** 2) / (4.0 * math.pi)
        t = self.models[wave]
        return sum(t.n_of_kc(m, k, 0.0) for m in (1.0, 0.0, -1.0)) / 3.0 * self.hbarc ** 3

    def nodes(self, wave="hulthen"):
        """The TaggedModel's k grid.  tagged.hpp: n_of_kc looks up the cell
        AT OR BELOW k (searchsorted - 1), so the table value T_j -- the
        density AT node k_j, (psihat_0^2 + psihat_2^2)/(4 pi) there -- holds
        on (k_j, k_j+1], and a query exactly at k_j returns T_j-1."""
        t = self.models.get(wave, self.models["hulthen"])     # cdbonn: the same grid
        return [float(x) for x in t.k]

    def n_node(self, wave, j):
        """The density at node k_j (a query inside its cell)."""
        ks = self.nodes(wave)
        if wave == "cdbonn":
            return self.n(wave, ks[j])
        return self.n(wave, ks[j] + 0.5 * (ks[1] - ks[0]))

    def n_smooth(self, wave, k):
        """POST HOC (header addendum): the same table read smoothly -- linear
        in ln n between the node values bracketing k.  CD-Bonn is analytic."""
        if wave == "cdbonn":
            return self.n(wave, k)
        ks = self.nodes(wave)
        j = max(0, min(len(ks) - 2, sum(1 for x in ks if x < k) - 1))
        a, b = self.n_node(wave, j), self.n_node(wave, j + 1)
        t = (k - ks[j]) / (ks[j + 1] - ks[j])
        return math.exp(math.log(a) + (math.log(b) - math.log(a)) * t)


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


def scaled_chi2(pts, model, allowance=ALLOWANCE):
    """The verdict's chi2 (design block): one free scale in closed form,
    sigma_i^2 = bar^2 + (A d)^2 + digitization^2."""
    d = [r["d"] for r in pts]
    m = [model(r) for r in pts]
    sg = [math.sqrt((0.5 * (r["hi"] - r["lo"])) ** 2 + (allowance * r["d"]) ** 2
                    + (math.log(10.0) * DIG_LOG10 * r["d"]) ** 2) for r in pts]
    sw = sum(mi * mi / s ** 2 for mi, s in zip(m, sg))
    s_hat = sum(di * mi / s ** 2 for di, mi, s in zip(d, m, sg)) / sw
    chi2 = sum(((di - s_hat * mi) / s) ** 2 for di, mi, s in zip(d, m, sg))
    ndf = len(d) - 1
    return dict(chi2=chi2, ndf=ndf, p=chi2_sf(chi2, ndf), s=s_hat,
                s_err=1.0 / math.sqrt(sw), n=len(d),
                pulls=[(di - s_hat * mi) / s for di, mi, s in zip(d, m, sg)],
                ratio=[di / (s_hat * mi) for di, mi in zip(d, m)])


# ------------------------------------------------------------ the comparisons

def measure(path=DATA):
    design = design_block_sha256()
    rows = load_reference(path)
    tree = Tree()
    win = [r for r in rows if r["panel"] in WINDOW_PANELS]
    if len(win) != NDF_VERDICT + 1:
        raise RuntimeError("the window holds %d points, the block says %d"
                           % (len(win), NDF_VERDICT + 1))
    nk = {w: (lambda r, w=w: tree.n(w, r["pm"])) for w in WAVES}
    verdict = {}
    for w in WAVES:
        f = scaled_chi2(win, nk[w])
        f["status"] = "pass" if f["p"] >= P_MIN else "fail"
        verdict[w] = f
    # (r1) the absolute closure, and the paper's own wave spread
    r1 = [r["cdbonn_pwia"] / tree.n("cdbonn", r["pm"]) for r in win]
    pp = [r["paris_pwia"] / r["cdbonn_pwia"] for r in win]
    # (r2) theory against theory, ln fit with uniform weights
    r2 = {}
    for w in WAVES:
        lr = [math.log(r["cdbonn_pwia"] / tree.n(w, r["pm"])) for r in win]
        c = sum(lr) / len(lr)
        rr = [math.exp(x - c) for x in lr]
        r2[w] = dict(s=math.exp(c), min=min(rr), max=max(rr))
    # (r3) no allowance; (r4) per panel
    r3 = {w: scaled_chi2(win, nk[w], allowance=0.0) for w in WAVES}
    r4 = {}
    for p in WINDOW_PANELS:
        pts = [r for r in win if r["panel"] == p]
        r4[p] = {w: scaled_chi2(pts, nk[w]) for w in WAVES}
    # (r5) the MS model's own FSI in the window
    fsi = [r["cdbonn_fsi"] / r["cdbonn_pwia"] for r in win]
    # (r6) the ratio route
    r6 = {}
    for p in WINDOW_PANELS:
        pts = [r for r in win if r["panel"] == p]
        rr = [r["d"] / r["cdbonn_pwia"] for r in pts]
        r6[p] = dict(mean=sum(rr) / len(rr), min=min(rr), max=max(rr))
    r6c = {w: scaled_chi2(win, (lambda r, w=w: tree.n(w, r["pm"]) * r["cdbonn_fsi"]
                                / r["cdbonn_pwia"]), allowance=0.0) for w in WAVES}
    smooth = {w: scaled_chi2(win, (lambda r, w=w: tree.n_smooth(w, r["pm"])))
              for w in WAVES}
    cell = {w: [tree.n(w, r["pm"]) / tree.n_smooth(w, r["pm"]) for r in win]
            for w in ("hulthen", "av18")}
    return dict(design_block_sha256=design, rows=rows, window=win, verdict=verdict,
                smooth_table=smooth,
                cell_bias={w: dict(min=min(c), max=max(c),
                                   at_lowest_two=[c[i] for i, r in enumerate(win)
                                                  if r["pm"] < 0.1])
                           for w, c in cell.items()},
                status=verdict[DEFAULT_WAVE]["status"],
                closure=dict(mean=sum(r1) / len(r1), min=min(r1), max=max(r1),
                             per_point=r1),
                paris_over_cdbonn=dict(mean=sum(pp) / len(pp), min=min(pp), max=max(pp)),
                theory=r2, no_allowance=r3, per_panel=r4,
                ms_fsi=dict(min=min(fsi), max=max(fsi), mean=sum(fsi) / len(fsi)),
                ratio_route=r6, fsi_corrected=r6c,
                tree_n={w: [tree.n(w, r["pm"]) for r in win] for w in WAVES})


def _word(ok):
    return "pass" if ok else "fail"


def run(path=DATA, verbose=True):
    t0 = time.perf_counter()
    m = measure(path)
    v = m["verdict"]
    win = m["window"]
    generator = ("n(k = p_m) x s, one free scale: hulthen (default) chi2 = %.2f/%d "
                 "(p = %.3g, s = %.3f); av18 %.2f (p = %.3g); cdbonn %.2f (p = %.3g)"
                 % (v["hulthen"]["chi2"], v["hulthen"]["ndf"], v["hulthen"]["p"],
                    v["hulthen"]["s"], v["av18"]["chi2"], v["av18"]["p"],
                    v["cdbonn"]["chi2"], v["cdbonn"]["p"]))
    reference = ("Boeglin et al. PRL 107 (2011) 262501 Fig. 2, d(e,e'p)n sigma_red "
                 "at Q^2 = 3.5 GeV^2, theta_nq = 35 deg (10) and 45 deg (12 points), "
                 "p_m = %.3f-%.3f GeV/c, digitized from the vector PDF"
                 % (min(r["pm"] for r in win), max(r["pm"] for r in win)))
    tolerance = ("p(chi2, %d) >= %g with one free scale, the bars (+) an allowance "
                 "A = %g of each point for FSI and the CC1 off-shell prescription "
                 "(declared in the hashed design block)" % (NDF_VERDICT, P_MIN, ALLOWANCE))
    subrows = []
    for w in WAVES:
        f = v[w]
        subrows.append(dict(name="nk[%s]" % w, status=f["status"],
                            in_verdict=(w == DEFAULT_WAVE), chi2=f["chi2"], ndf=f["ndf"],
                            p=f["p"], s=f["s"], s_err=f["s_err"],
                            max_abs_pull=max(abs(x) for x in f["pulls"])))
    c = m["closure"]
    subrows.append(dict(name="closure[paper CD-Bonn PWIA / tree cdbonn n]",
                        status="recorded", in_verdict=False, mean=c["mean"],
                        min=c["min"], max=c["max"]))
    c = m["paris_over_cdbonn"]
    subrows.append(dict(name="paper[Paris PWIA / CD-Bonn PWIA]", status="recorded",
                        in_verdict=False, mean=c["mean"], min=c["min"], max=c["max"]))
    for w in WAVES:
        t = m["theory"][w]
        subrows.append(dict(name="theory[paper CD-Bonn PWIA vs %s]" % w, status="recorded",
                            in_verdict=False, s=t["s"], min_ratio=t["min"], max_ratio=t["max"]))
    for w in WAVES:
        f = m["no_allowance"][w]
        subrows.append(dict(name="no_allowance[%s]" % w, status="recorded", in_verdict=False,
                            chi2=f["chi2"], ndf=f["ndf"], p=f["p"], s=f["s"]))
    for p in WINDOW_PANELS:
        for w in WAVES:
            f = m["per_panel"][p][w]
            subrows.append(dict(name="panel_%s[%d deg, %s]" % (p, THETA_NQ[p], w),
                                status="recorded", in_verdict=False, chi2=f["chi2"],
                                ndf=f["ndf"], p=f["p"], s=f["s"]))
    c = m["ms_fsi"]
    subrows.append(dict(name="ms_fsi[CD-Bonn FSI / PWIA]", status="recorded",
                        in_verdict=False, min=c["min"], max=c["max"], mean=c["mean"],
                        allowance=ALLOWANCE))
    for p in WINDOW_PANELS:
        c = m["ratio_route"][p]
        subrows.append(dict(name="ratio_route[data / paper CD-Bonn PWIA, %d deg]" % THETA_NQ[p],
                            status="recorded", in_verdict=False, mean=c["mean"],
                            min=c["min"], max=c["max"]))
    for w in WAVES:
        f = m["fsi_corrected"][w]
        subrows.append(dict(name="fsi_corrected[%s x MS FSI/PWIA, A = 0]" % w,
                            status="recorded", in_verdict=False, chi2=f["chi2"],
                            ndf=f["ndf"], p=f["p"], s=f["s"]))
    for w in WAVES:
        f = m["smooth_table"][w]
        subrows.append(dict(name="nk_smooth_table[%s, post hoc]" % w,
                            status="pass" if f["p"] >= P_MIN else "fail", in_verdict=False,
                            post_hoc=True, chi2=f["chi2"], ndf=f["ndf"], p=f["p"], s=f["s"]))
    for w, c in m["cell_bias"].items():
        subrows.append(dict(name="cell_bias[%s, post hoc]" % w, status="recorded",
                            in_verdict=False, post_hoc=True, min=c["min"], max=c["max"]))
    subrows.append(dict(name="panel_c[75 deg]", status="not compared", in_verdict=False,
                        n_points=sum(1 for r in m["rows"] if r["panel"] == "c"),
                        reason="FSI-dominated by the paper's own account; vendored and drawn"))
    import lipolgen as lg  # noqa: WPS433
    config = dict(lipolgen=getattr(lg, "__version__", None), python=sys.version.split()[0],
                  data_file_sha256=DATA_SHA256, design_block_sha256=m["design_block_sha256"],
                  bg_pdf_sha256=BG_PDF_SHA256, bg_page=BG_PAGE, allowance=ALLOWANCE,
                  k_identification="instant form, k = p_m (the paper's PWIA reduction)")
    points = [dict(panel=r["panel"], theta_nq=r["theta_nq"], pm=r["pm"], d=r["d"],
                   err=0.5 * (r["hi"] - r["lo"]), cdbonn_pwia=r["cdbonn_pwia"],
                   cdbonn_fsi=r["cdbonn_fsi"], paris_pwia=r["paris_pwia"],
                   paris_fsi=r["paris_fsi"],
                   **{"n_" + w: m["tree_n"][w][i] for w in WAVES})
              for i, r in enumerate(win)]
    report = dict(name=NAME, generator=generator, reference=reference, tolerance=tolerance,
                  status=m["status"], reason=None, subrows=subrows, config=config,
                  verdict={w: {k: x for k, x in f.items()} for w, f in v.items()},
                  closure=m["closure"], paris_over_cdbonn=m["paris_over_cdbonn"],
                  theory=m["theory"], no_allowance=m["no_allowance"],
                  per_panel=m["per_panel"], ms_fsi=m["ms_fsi"],
                  ratio_route=m["ratio_route"], fsi_corrected=m["fsi_corrected"],
                  smooth_table=m["smooth_table"], cell_bias=m["cell_bias"],
                  points=points,
                  panel_c=[dict(pm=r["pm"], d=r["d"], err=0.5 * (r["hi"] - r["lo"]))
                           for r in m["rows"] if r["panel"] == "c"])
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        print("reference: %d points (35 deg: %d, 45 deg: %d; 75 deg: %d not compared); "
              "design block %s..." % (len(m["rows"]), EXPECTED_POINTS["a"], EXPECTED_POINTS["b"],
                                      EXPECTED_POINTS["c"], m["design_block_sha256"][:12]))
        print("   panel  p_m     sigma_red    +-       CDB PWIA   CDB FSI    hulthen    av18       cdbonn  (tree n, fm^3)")
        for i, r in enumerate(win):
            print("   %s      %.4f  %.4e  %.2e  %.4e %.4e %.4e %.4e %.4e" % (
                r["panel"], r["pm"], r["d"], 0.5 * (r["hi"] - r["lo"]), r["cdbonn_pwia"],
                r["cdbonn_fsi"], m["tree_n"]["hulthen"][i], m["tree_n"]["av18"][i],
                m["tree_n"]["cdbonn"][i]))
        for w in WAVES:
            f = v[w]
            print("verdict %-8s chi2 = %.3f/%d, p = %.4g, s = %.4f +- %.4f -> %s%s"
                  % (w, f["chi2"], f["ndf"], f["p"], f["s"], f["s_err"], f["status"].upper(),
                     "  (the row's)" if w == DEFAULT_WAVE else ""))
            print("   ratios d/(s n): %s" % " ".join("%.2f" % x for x in f["ratio"]))
        c = m["closure"]
        print("(r1) paper CD-Bonn PWIA / tree cdbonn n: mean %.4f (%.4f .. %.4f); "
              "paper Paris/CD-Bonn PWIA %.4f (%.4f .. %.4f)" % (
                  c["mean"], c["min"], c["max"], m["paris_over_cdbonn"]["mean"],
                  m["paris_over_cdbonn"]["min"], m["paris_over_cdbonn"]["max"]))
        for w in WAVES:
            t = m["theory"][w]
            print("(r2) paper CD-Bonn PWIA vs %-8s s = %.4f, curve/(s n) %.3f .. %.3f"
                  % (w, t["s"], t["min"], t["max"]))
        for w in WAVES:
            f = m["no_allowance"][w]
            print("(r3) A = 0: %-8s chi2 = %.2f/%d, p = %.3g, s = %.4f" % (w, f["chi2"], f["ndf"], f["p"], f["s"]))
        for p in WINDOW_PANELS:
            for w in WAVES:
                f = m["per_panel"][p][w]
                print("(r4) %d deg alone: %-8s chi2 = %.2f/%d, p = %.3g, s = %.4f"
                      % (THETA_NQ[p], w, f["chi2"], f["ndf"], f["p"], f["s"]))
        c = m["ms_fsi"]
        print("(r5) MS CD-Bonn FSI/PWIA in the window: %.3f .. %.3f (mean %.3f); A = %g"
              % (c["min"], c["max"], c["mean"], ALLOWANCE))
        for p in WINDOW_PANELS:
            c = m["ratio_route"][p]
            print("(r6) data / paper CD-Bonn PWIA, %d deg: mean %.3f (%.3f .. %.3f)"
                  % (THETA_NQ[p], c["mean"], c["min"], c["max"]))
        for w in WAVES:
            f = m["fsi_corrected"][w]
            print("(r6) tree x MS FSI/PWIA, A = 0: %-8s chi2 = %.2f/%d, p = %.3g, s = %.4f"
                  % (w, f["chi2"], f["ndf"], f["p"], f["s"]))
        for w in WAVES:
            f = m["smooth_table"][w]
            print("POST HOC, the table read smoothly: %-8s chi2 = %.3f/%d, p = %.4g, s = %.4f"
                  % (w, f["chi2"], f["ndf"], f["p"], f["s"]))
        for w, c in m["cell_bias"].items():
            print("POST HOC, cell value / smooth reading at the 22 points, %s: %.4f .. %.4f "
                  "(p_m < 0.1: %s)" % (w, c["min"], c["max"],
                                       ", ".join("%.4f" % x for x in c["at_lowest_two"])))
        for r in subrows:
            vals = ", ".join("%s = %s" % (k, ("%.6g" % x) if isinstance(x, float) else x)
                             for k, x in r.items() if k not in ("name", "status", "in_verdict"))
            print("SUBROW | %s[%s] | %s | %s" % (NAME, r["name"], vals, str(r["status"]).upper()))
        print("REPORT | %s | %s | %s | %s | %s" % (NAME, generator, reference, tolerance,
                                                   m["status"].upper()))
    return report


# ---------------------------------------------------------------- figure

def plot(report, out_dir=None):
    """Top: sigma_red against p_m on a log axis -- the 35 and 45 deg points
    (ink; the reference) with their bars, the 75 deg points (grey, record
    (r7): drawn, never compared), each wave's n(k = p_m) times its
    own fitted scale (a line in its slot: the default `hulthen`, `av18` the
    alternative's, `cdbonn` the third's), and the declared allowance band
    s n (1 +- A) around the default wave (grey).  Bottom: the same over
    s n_hulthen -- the data with their bars (divided by the verdict's model,
    the TaggedModel cell value at each p_m), the band 1 +- A, and the other
    two waves' fitted shapes over the default's.  The lines go through the
    TaggedModel's k nodes (its density is piecewise constant between them),
    for display only; every number the figure states is the report's."""
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    from _plots import BAND, INK, MUTED, band, new_figure, save, status_title, style_series
    tree = Tree()
    v = report["verdict"]
    pts = report["points"]
    fig, (ax, axr) = new_figure(two_panel=True, height=4.2)
    lo, hi = 0.04, 0.56
    # the TaggedModel's own k nodes (its density is piecewise constant
    # between them, tagged.hpp): every curve drawn through its node values
    idx = [j for j, k in enumerate(tree.nodes()) if lo <= k <= hi]
    grid = [tree.nodes()[j] for j in idx]
    s0 = v[DEFAULT_WAVE]["s"]
    n0 = [s0 * tree.n_node(DEFAULT_WAVE, j) for j in idx]
    band(ax, grid, [x * (1 - ALLOWANCE) for x in n0], [x * (1 + ALLOWANCE) for x in n0],
         label="allowance A = %g around s n (hulthen)" % ALLOWANCE)
    for p, fmt, lab in (("a", "o", "Boeglin 2011, theta_nq = 35 deg"),
                        ("b", "s", "Boeglin 2011, theta_nq = 45 deg")):
        q = [r for r in pts if r["panel"] == p]
        kw = style_series("reference", label=lab, marker=fmt)
        if p == "b":
            kw.update(markerfacecolor="white")
        ax.errorbar([r["pm"] for r in q], [r["d"] for r in q], yerr=[r["err"] for r in q], **kw)
        rel = [r["d"] / (s0 * r["n_" + DEFAULT_WAVE]) for r in q]
        erel = [r["err"] / (s0 * r["n_" + DEFAULT_WAVE]) for r in q]
        kw = style_series("reference", label=None, marker=fmt)
        if p == "b":
            kw.update(markerfacecolor="white")
        axr.errorbar([r["pm"] for r in q], rel, yerr=erel, **kw)
    pc = report["panel_c"]                     # record (r7): drawn, never compared
    kw = style_series("reference", label="75 deg (FSI-dominated, not compared)", marker="^")
    kw.update(color=MUTED, markerfacecolor="none", markersize=3, zorder=3)
    ax.errorbar([r["pm"] for r in pc], [r["d"] for r in pc], yerr=[r["err"] for r in pc], **kw)
    for w, slot in (("hulthen", "default"), ("av18", "alt"), ("cdbonn", "alt2")):
        f = v[w]
        nw = [f["s"] * tree.n_node(w, j) for j in idx]
        ax.plot(grid, nw, linewidth=1.6,
                **style_series(slot, label="%s: chi2 %.1f/%d, p = %.2g, s = %.3f"
                               % (w, f["chi2"], f["ndf"], f["p"], f["s"])))
        if w != DEFAULT_WAVE:
            axr.plot(grid, [a / b for a, b in zip(nw, n0)], linewidth=1.4,
                     **style_series(slot, label=None))
    axr.axhspan(1 - ALLOWANCE, 1 + ALLOWANCE, color=BAND, alpha=0.8, linewidth=0, zorder=1)
    axr.axhline(1.0, color=INK, linewidth=0.8, zorder=2)
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylabel("sigma_red = sigma_exp/(k sigma_cc1) [fm^3]")
    ax.legend(loc="upper right", fontsize=7)
    axr.set_ylim(0.4, 2.5)
    axr.set_xlabel("p_m = k [GeV/c]")
    axr.set_ylabel("/ (s n_hulthen)")
    status_title(fig, NAME, report["status"],
                 "35 + 45 deg, 22 points, one free scale, bars (+) A = %g: hulthen p = %.2f "
                 "(PASS iff p >= %g); A = 0: all fail"
                 % (ALLOWANCE, v[DEFAULT_WAVE]["p"], P_MIN))
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
