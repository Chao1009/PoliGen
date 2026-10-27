#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 6 -- Rand-Frosch-Yearian: the MEASURED 6Li
magnetic form factor against the tree's F_m shape (tier T3); 7Li recorded.

Re-wired 2026-10-26 (run 2026-10-26 phase B2; PLAN.md phase B2, rewritten
after review; docs/benchmarking/08_channel_priority.md sec. 4, ingredient F,
and BENCHMARK_PLAN.md sec. 6, U-10).  Until then the row was BLOCKED on a
document and its header asked for a CSV of magnetization RADII.  The document
is now held and it prints no radius: it prints the measured F_M^2(q^2) points
(Tables I, II) and the parameters of the harmonic-well and cluster models it
fitted to them (Tables VII, VIII, XIV).  So the row now compares the POINTS.

THE DOCUMENT.  R. E. Rand, R. F. Frosch, M. R. Yearian, "Elastic Electron
Scattering from the Magnetic Multipole Distributions of Li6, Li7, Be9, B10,
B11, and N14", Phys. Rev. 144 (1966) 859-873, doi:10.1103/PhysRev.144.859
(docs/references/REFERENCES.md sec. 1 #8): the APS PDF, sha256 7cbb2d21...
d5bbe2, fetched 2026-10-26 through the institution's licence and kept OUTSIDE
the repository.  Vendored from it, read from the page images:
  data/rfy1966_li_magnetic_ff.csv     Tables I (6Li, 7 points) and II (7Li,
                                      24 points), p. 865: E0, q^2, F_M^2 +- err
  data/rfy1966_li_magnetization.csv   Tables VII, VIII, XIV (pp. 869, 870,
                                      873): a0 and the coupling parameters of
                                      every 6Li / 7Li fit, verbatim
Each file's provenance block states the publishing convention; the sha256 of
every byte of each is a constant below, re-checked on every read (an altered
file raises RuntimeError and the harness exits 2).

THE PHYSICS OF THE TWO SIDES.  RFY measure at 180 degrees, where Coulomb
scattering is suppressed, and quote F_M^2(q^2) of their Eqs. (5)-(6),
normalised to 1 at q^2 = 0 through the measured dipole moment mu0 in Eq. (9);
q is the FOUR-momentum transfer at theta = pi (p. 860).  For 6Li (J = 1)
lambda <= 2J allows M1 only: F_M^2 = F_M1^2.  For 7Li (J = 3/2) M3 adds, and
the unresolved 0.478-MeV level is inside the measured number (Table II's
caption).  The tree's side is rc.hpp's
    F_m(t) = (M_A/m_p)(mu/mu_N) F_mag(q) * tail_tensor_scale,
    F_mag(q) = (1 - q^2/q_z^2) exp(-q^2 b^2/4),   q^2 = t/(hbar c)^2,
with q_z = LI6_FF_FM_QZ_FM = 1.30 fm^-1 and b = LI6_FF_FM_B_FM = 1.85 fm --
"UNFITTED STARTING VALUES" in rc.hpp's own words (q_z the midpoint of a
design range, b the geometric mean of two single-landmark solutions against
Wiringa-Schiavilla's F_T) -- and NO nucleon folding on F_m.  So
F_m(t)/F_m(0) = F_mag(q) is a NORMALISED magnetic form factor, the same kind
of object as RFY's F_M, and the comparison is of F_mag^2 with F_M^2.  Its
zero q_z^2 = 1.69 fm^-2 lies INSIDE Table I's q^2 range (0.72-1.92 fm^-2),
which is why these seven points can test it.  The measured F_M includes the
nucleons' own magnetic size; the tree's F_mag is meant to be the whole F_m it
puts into POLRAD's elastic tail, so it is compared as it is used.

WHAT IS COMPARED, AND THE TOLERANCE -- fixed in this header on 2026-10-26
before this harness computed anything (PLAN.md ground rule 4).  The block
between the two marker lines is hashed (DESIGN_BLOCK_SHA256 below);
`design_block_sha256()` re-checks it, and the tolerance constants against
its sentences, on every run, and the harness exits 2 if either moved.
==== DESIGN BLOCK BEGIN ====
THE VERDICT.  The seven 6Li points of RFY Table I (q^2 = 0.72 ... 1.92 F^-2,
q = 0.85 ... 1.39 fm^-1), d_i = F_M^2 +- sigma_i as printed, against
    m_i = [F_m(t_i)/F_m(0)]^2 = F_mag(q_i)^2
from the SHIPPED HoSpin1FF.for_ion(li6()).fm, at t_i = q_i^2 (hbar c)^2 with
q_i^2 the PRINTED q^2 (RFY's four-momentum transfer squared at 180 degrees;
the tree's q^2 = t/(hbar c)^2 is the same invariant).  m_i is also checked
against the closed form (1 - q^2/q_z^2)^2 exp(-q^2 b^2/2) to 1e-12, so the
row reads the tree's shape.  Nothing is fitted except ONE free scale s on
the model:
    chi2(s) = sum_i [(d_i - s m_i)/sigma_i]^2,
    s_hat = sum(d m/sigma^2)/sum(m^2/sigma^2)   (closed form, unconstrained),
diagonal, with the printed errors, uncorrelated (the paper states no
correlation) and Gaussian (also for the 0.001 +- 0.001 point); all seven
points, none dropped (the three above 120 MeV, dominated by the 7Li
impurity, are in the paper's own fits too).  ndf = 7 - 1 = 6.
TOLERANCE: PASS iff p(chi2_min, ndf = 6) >= 0.01 -- the plan's "shape chi2
with one free scale, the paper's own errors".  p is scipy's chi2.sf when
scipy is importable, else the Numerical Recipes gammq fallback of
t1_nmc_f2d.py, re-implemented here.
RECORDED, NEVER IN THE VERDICT -- each rule fixed here, before computing:
 (r1) s = 1 fixed: both sides normalised to 1 at q^2 = 0 (the data through
      mu0 = 0.822 nm in RFY Eq. (9), the tree through F_m(0) = (M_A/m_p) mu):
      chi2 and p with ndf = 7.
 (r2) the per-point pulls (d_i - s_hat m_i)/sigma_i.
 (r3) the verdict's chi2 and s_hat with q_i^2 RECOMPUTED from the printed E0
      at theta = pi: q^2 = 4 E0 E'/(hbar c)^2, E' = E0/eta, eta = 1 + 2 E0/M_A
      (RFY Eq. (1) at theta = pi; the electron mass neglected, E0 >= 70 MeV),
      M_A and hbar c the tree's (lg.li6().mass(), lg.li7().mass(),
      lg.HBARC_GEV_FM).  And over all 31 rows of Tables I and II the largest
      |q^2_recomputed - q^2_printed|: the q^2 column CHECK holds iff it is
      <= 0.006 F^-2 (half a unit of the printed last digit plus 1e-3 for 1966
      constants); a failure flags the transcription, it moves no verdict.
 (r4) 7Li, Table II (24 points): NOT IMPLEMENTED IN THE TREE --
      HoSpin1FF.for_ion(li7()) throws (no measured-moment block in rc.hpp);
      the points are recorded and drawn, never compared; the sub-row's status
      is "not compared".
 (r5) TRANSCRIPTION CLOSURES: the paper's own fits re-evaluated at the
      vendored points with the vendored parameters and NO free parameter
      (normalised to 1 at q^2 = 0, p. 868), chi2 over that table's points
      with the printed errors:
        6Li, every harmonic-well row (A), (B), (C) of Table VII: Eq. (17)
        with R2 = 0 (J = 1), F_M = [<j0> + R1 <j2>] F_N, <j0> = (1 - 2x/3)
        e^-x (c.m.), <j2> = (2/3) x e^-x (c.m.), x = q^2 a0^2/4, (c.m.) =
        e^(x/A) with A = 6, F_N = (1 + q^2/7.1)^-1 (p. 868);
        6Li alpha+d (1s): F_M = F_d exp[-(q^2/(4 q0^2))(1 - 2/A)], F_d =
        exp(-q^2 (1.8 F)^2/6), the Gaussian of rms radius 1.8 F (p. 869);
        7Li, every harmonic-well row of Table VIII: Eq. (18), F_M^2 = [A
        <j0>^2 + B <j0><j2> + C <j2>^2] F_N^2 with the printed A, B, C and
        (c.m.) = e^(x/7).
      A closure HOLDS iff chi2_table <= chi2_paper + 1: the paper's chi2 also
      counts other authors' points (each >= 0), and rounding the printed
      parameters moves chi2 by far less than 1.  One that does not hold
      flags the transcription or this reading of Eqs. (17)/(18) (for 7Li,
      also whether Table II's errors are those the fit used: pp. 870-871 say
      the S1, S2 errors were "included in the final experimental error"); it
      is resolved by re-reading the page image, NEVER by editing a number,
      and it does not change the row's status.
 (r6) a0 VERBATIM (data/rfy1966_li_magnetization.csv) and, per a0 row,
      sqrt(5/2) a0 -- RFY p. 868's "<r^2> = (5/2) a0^2", the rms radius of
      the 1p-shell orbital of the fitted well, labelled as such and NEVER
      compared with the tree.
 (r7) the tree's slope radius of F_mag, r_mag = sqrt(6/q_z^2 + (3/2) b^2),
      UNCOMPARED (the BLOCKED row printed it: 2.9469 fm in the committed
      report).
THE 2026-09-23 RADIUS RULE (the tree's r_mag inside the envelope of
vendored r_mag rows) IS RETIRED -- the dated note below this block.
WHAT WAS SEEN BEFORE THIS BLOCK WAS WRITTEN (by the worker, 2026-10-26):
every row of Tables I, II, VII, VIII and XIV while transcribing them; the
committed report's tree-side r_mag(6Li) = 2.9469 fm; and, while reading
Table I, a BY-HAND evaluation of F_mag^2 at the seven points with a chi2
estimate -- about 17.5 for 7 points at s = 1 and about 13 for 6 dof with the
free scale (p ~ 0.04 by eye).  The tolerance above is the one PLAN.md phase
B2 and the task fixed before that estimate; no rule was changed after it.
Also by hand before this block: q^2 from E0 for six rows (r3); sqrt(5/2) a0
for the Table VII/VIII rows, which reproduces to the printed digit the
"magnetization radii" docs/benchmarking/06_critic.md sec. 4.6 transcribed
from ADNDT 14 (1974) Table V (3.42(6), 2.81, 2.69(13), 2.72 fm); and a slope
radius from Eq. (17) for the Table VII rows -- which is one of the two
reasons no radius sub-row is kept.
==== DESIGN BLOCK END ====

POST-HOC ADDENDUM (2026-10-26, written after the first run of 23:00:57Z and
OUTSIDE the hashed block, which is unchanged).  (a) That run found ONE
closure not holding: 7Li Table VIII (C) a0 = 1.66, chi2 over the 24 Table II
points 26.91 against the paper's 25.6 + 1 (the other 12 hold).  The page
image was re-read at 450 dpi: the row is as printed (A = 1.15, B = 0.25,
C = 1.09, chi2 = 25.6, 23 dof) and so are the 24 points.  The block's premise
that rounding the printed parameters "moves chi2 by far less than 1" is not
guaranteed for three strongly correlated parameters, so a rounding-free
check is ADDED as a diagnostic -- never a verdict; the closure's own result
stays recorded as it came out.  Eq. (18) is LINEAR in (A, B, C): at the
printed a0 their weighted least-squares minimum over the 24 table points is
computed in closed form.  The paper's chi2 is the minimum over a SUPERSET of
these points (other authors' points, each >= 0), so chi2_min(24) <=
chi2_paper must hold if the transcription, this reading of Eq. (18) and the
errors are the paper's -- exactly for the (B) and (C) rows (a0 fixed by the
fit's design) and, for the (A) row, up to the rounding of a0 (+-0.005 of its
+-0.01: at most 0.25 in chi2).  Rule: the REFIT HOLDS iff chi2_min <=
chi2_paper, + 0.25 for the (A) row.  (b) The (r3) q^2 check's 0.006 F^-2
did not allow for the rounding of the integer E0 (+-0.5 MeV moves q^2 by
+-0.0145 F^-2 at 147 MeV); it held all the same (max 0.0058 F^-2, Table II
row 17, 147 MeV), so nothing turns on it; the limit is not changed.
(c) Outcome of (a), the run of 23:02:41Z: the refit holds for all five rows
(chi2_min 21.36 / 28.82 / 24.53 / 20.33 / 20.20 against 23.3 / 30.8 / 25.6 /
23.5 / 30.0).  At a0 = 1.66 the 24-point optimum lies +0.50 / +0.02 / +1.41
printed errors from the printed A / B / C -- farther than rounding (at most
0.25 of an error in C) -- which is consistent with the paper's fit including
the two other-author points of Fig. 6 (26 points, 23 dof), and not with a
transcription error.  The closure sub-row stays FAIL as recorded.

THE RADIUS RULE, RETIRED (dated note, 2026-10-26).  The 2026-09-23 header
planned "the tree's r_mag inside the envelope of ALL the vendored 6Li rows,
each widened by its own quoted error", reading r_mag rows from RFY or ADNDT
14 Table V.  Retired, for two reasons.  (i) RFY print no magnetization
radius.  The only radius of their model is p. 868's <r^2> = (5/2) a0^2, the
1p-shell orbital rms of the fitted well; p. 873 distinguishes it from the
slope of F_M at q^2 = 0 ("does not imply <r^2> = 0, as may be seen by
inspection of (17)"), and the ADNDT 14 numbers 06_critic.md sec. 4.6 quotes
are exactly sqrt(5/2) a0 (record (r6)).  A slope radius would need R1 as
well as a0, the c.m. factor e^(x/A) and the nucleon factor F_N, and a0 is
only known to 1.6-2.2 F for 6Li (p. 869, Table XIV).  (ii) The task allowed
a slope-radius sub-row only if its formula, the c.m. and nucleon-size choices
and the RFY error entering it were written in this header BEFORE anything
was derived; the worker evaluated such a formula by hand before this header
existed, so that condition is not met and no radius is compared or recorded
as a comparison.  a0 is vendored verbatim instead (record (r6)); the tree's
r_mag stays a tree-side number (record (r7)).

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations / A-scaling -- NOT APPLICABLE: no b1, no A-scaling.
  2. isoscalar denominator -- NOT APPLICABLE: no structure function, no
     nuclear PDF grid.  (F_mag is not nucleon-folded; RFY's model folds F_N,
     the same for both nucleon kinds -- recorded in (r5), not asserted.)
  3. per nucleus / per nucleon -- both sides PER NUCLEUS and normalised: RFY's
     F_M^2 is the whole nucleus's, normalised to 1 at q^2 = 0 through mu0 of
     the nucleus (Eq. (9)); the tree's F_m(0) = (M_A/m_p) mu is the
     Rosenbluth normalisation of the same per-nucleus moment, and the row
     compares F_m(t)/F_m(0), so the normalisation cancels (asserted by the
     pytest: tail_tensor_scale and M_A drop out).

WHAT THIS DOES NOT VALIDATE.  F_m outside q = 0.85-1.39 fm^-1 -- in
particular the second lobe near q = 2 fm^-1 that rc.hpp's comment says is
"where the tail lives" (RFY's 6Li data stop at q^2 = 1.92 F^-2; Bergstrom,
Kowalski, Neuhausen, PRC 25 (1982) 1156, REFERENCES.md #22, reaches higher
and is not held); F_m(0), which is the measured moment by construction; C0
and C2 (rows 5 and the de Vries radius); the q -> 0 slope (the data start at
q = 0.85 fm^-1 and the free scale absorbs a normalisation); 7Li, which the
tree does not model; and whether the RC tail itself is right (row 7 and the
POLRAD rows).  A PASS means the unfitted (q_z, b) shape, scaled once, is
statistically compatible with seven 1966 points with errors of 30-100 %; it
does not mean those starting values are fitted or right.

Run:  source env.sh && python3 validation/benchmarks/t3_li_magnetization_rfy.py [--plot]

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, an edited design block, a missing dependency, any
exception).  A harness is a measurement, not a CI gate: the pytests are the
gate.
"""

import csv
import hashlib
import math
import os
import re
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
#: RFY Tables I and II (6Li, 7Li F_M^2 points) -- the PRIMARY comparison
DATA_FF = os.path.join(HERE, "data", "rfy1966_li_magnetic_ff.csv")
#: RFY Tables VII, VIII, XIV (a0 and the coupling parameters), verbatim
DATA_A0 = os.path.join(HERE, "data", "rfy1966_li_magnetization.csv")
#: sha256 of EVERY byte of each vendored file, as written 2026-10-26
#: (22:53:49Z and 22:54:51Z); re-checked on every read
FF_SHA256 = "a3ada6b099c50e8fd6c90453e6fd964f8436b97675d1bdb17c0a75bc4e593aa8"
A0_SHA256 = "8f6874262d2ce32cab4645c8095c177be5e35c2f01f2bf270d13c64ac5bed1e8"
#: sha256 of this header's DESIGN BLOCK (its two marker lines included), as
#: written at 2026-10-26T22:57:20Z, before this harness computed anything.
#: `design_block_sha256()` re-checks it, and the tolerance constants below
#: against the block's own sentences, on every run.
DESIGN_BLOCK_SHA256 = (
    "d00b1307a0e397dcad2830513abdcbe21da9fb7eef3d1931e5976a7a0b3a3589")
DESIGN_BEGIN = "==== DESIGN BLOCK BEGIN ===="
DESIGN_END = "==== DESIGN BLOCK END ===="

NAME = "t3_li_magnetization_rfy"
#: the tolerance (design block): PASS iff p(chi2_min, ndf = NDF_VERDICT) >= P_MIN
P_MIN = 0.01
NDF_VERDICT = 6
#: (r5): a transcription closure holds iff chi2_table <= chi2_paper + CLOSURE_MARGIN
CLOSURE_MARGIN = 1.0
#: (r3): the q^2 column check holds iff max |q2_recomputed - q2_printed| <= this
Q2_CHECK_FM2 = 0.006
#: the tree's F_mag read through `fm` must equal the closed form to this
#: (absolute, on F_mag <= 1)
SHAPE_CHECK_ABS = 1e-12
#: RFY p. 868: the nucleonic form factor F_N(q^2) = (1 + q^2/7.1)^-1, q^2 in F^-2
RFY_FN_SCALE_FM2 = 7.1
#: RFY p. 869: the alpha+d model's Gaussian deuteron form factor, rms 1.8 F
RFY_DEUTERON_RMS_FM = 1.8
MASS_NUMBER = {"6Li": 6, "7Li": 7}
#: the post-hoc refit diagnostic (header addendum (a)): the (A) row's slack
#: for the rounding of its fitted a0
REFIT_SLACK_A = 0.25
#: under 1 s (0.4-0.5 s wall measured 2026-10-26); a library evaluation, no
#: event sample
EXPENSIVE = False

FF_COLUMNS = ["nucleus", "table", "row", "page", "e0_mev", "q2_fm2", "fm2",
              "fm2_err"]
A0_COLUMNS = ["nucleus", "fit", "a0_fm", "a0_err_fm", "coupling_note",
              "q_range_fm", "table", "page"]


# ------------------------------------------------------------ the design block

_TOL_RE = re.compile(r"PASS iff p\(chi2_min, ndf = (\d+)\) >= ([0-9.]+)")
_CLOSURE_RE = re.compile(r"HOLDS iff chi2_table <= chi2_paper \+ ([0-9.]+)")
_Q2CHECK_RE = re.compile(r"<= ([0-9.]+) F\^-2 \(half a unit")


def design_block(source=None):
    """The DESIGN BLOCK of this file's header, both marker lines included."""
    if source is None:
        with open(os.path.abspath(__file__), encoding="utf-8") as f:
            source = f.read()
    i0 = source.index(DESIGN_BEGIN)
    i1 = source.index(DESIGN_END) + len(DESIGN_END)
    return source[i0:i1]


def tolerances_from_design_block(block):
    """(ndf, p_min, closure_margin, q2_check) as the block's sentences state
    them -- RuntimeError if a sentence is missing."""
    t, c, q = _TOL_RE.search(block), _CLOSURE_RE.search(block), _Q2CHECK_RE.search(block)
    if t is None or c is None or q is None:
        raise RuntimeError("the design block no longer states its tolerance, "
                           "closure margin or q^2 check")
    return (int(t.group(1)), float(t.group(2)), float(c.group(1).rstrip(".")),
            float(q.group(1)))


def design_block_sha256(source=None):
    """sha256 of the design block -- RuntimeError unless it is the recorded
    one AND the block's stated tolerances are the constants this code uses."""
    block = design_block(source)
    got = hashlib.sha256(block.encode("utf-8")).hexdigest()
    if got != DESIGN_BLOCK_SHA256:
        raise RuntimeError("the header's design block sha256 %s != recorded %s: "
                           "the design or the tolerance was edited"
                           % (got, DESIGN_BLOCK_SHA256))
    stated = tolerances_from_design_block(block)
    used = (NDF_VERDICT, P_MIN, CLOSURE_MARGIN, Q2_CHECK_FM2)
    if stated != used:
        raise RuntimeError("the tolerance constants %r differ from the hashed "
                           "design block's %r" % (used, stated))
    return got


# ------------------------------------------------------------ the references

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


def _csv(raw, columns, path):
    lines = [l for l in raw.decode("utf-8").splitlines()
             if l.strip() and not l.startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames != columns:
        raise RuntimeError("%s: columns %r, expected %r"
                           % (path, reader.fieldnames, columns))
    return list(reader)


def load_ff(path=DATA_FF):
    """(6Li points of Table I, 7Li points of Table II), each a dict with the
    printed numbers as floats and as text."""
    pts = []
    for r in _csv(read_vendored(path, FF_SHA256), FF_COLUMNS, path):
        pts.append(dict(
            nucleus=r["nucleus"], table=r["table"], row=int(r["row"]),
            page=int(r["page"]), e0_mev=float(r["e0_mev"]),
            q2_fm2=float(r["q2_fm2"]), fm2=float(r["fm2"]),
            fm2_err=float(r["fm2_err"]), q2_as_printed=r["q2_fm2"],
            as_printed="%s +- %s" % (r["fm2"], r["fm2_err"])))
    li6 = [p for p in pts if p["nucleus"] == "6Li" and p["table"] == "I"]
    li7 = [p for p in pts if p["nucleus"] == "7Li" and p["table"] == "II"]
    if (len(li6), len(li7), len(pts)) != (7, 24, 31):
        raise RuntimeError("%s: expected 7 Table I and 24 Table II rows, found "
                           "%d / %d of %d" % (path, len(li6), len(li7), len(pts)))
    if any(not p["fm2_err"] > 0.0 for p in pts):
        raise RuntimeError("%s: a non-positive error" % path)
    return li6, li7


_VALUE_RE = re.compile(r"^(-?[0-9]+(?:\.[0-9]+)?)(?:\+-([0-9]+(?:\.[0-9]+)?))?$")
_RANGE_RE = re.compile(r"^([0-9]+(?:\.[0-9]+)?)-([0-9]+(?:\.[0-9]+)?)$")


def parse_note(note):
    """The coupling_note's machine part -> {KEY: (value, err or None)}, and
    the paper's words after ' | '."""
    machine, _, text = note.partition(" | ")
    params = {}
    for item in machine.split(";"):
        item = item.strip()
        if not item:
            continue
        key, eq, val = item.partition("=")
        m = _VALUE_RE.match(val.strip())
        if not eq or m is None or not key.strip():
            raise RuntimeError("coupling_note item %r is not KEY=VALUE[+-ERR]" % item)
        params[key.strip()] = (float(m.group(1)),
                               None if m.group(2) is None else float(m.group(2)))
    return params, text.strip()


def load_a0(path=DATA_A0):
    """Every row of Tables VII, VIII and XIV, verbatim, parsed."""
    out = []
    for r in _csv(read_vendored(path, A0_SHA256), A0_COLUMNS, path):
        a0 = r["a0_fm"].strip()
        rng = _RANGE_RE.match(a0)
        params, text = parse_note(r["coupling_note"])
        out.append(dict(
            nucleus=r["nucleus"], fit=r["fit"], a0_as_printed=a0,
            a0_fm=(float(a0) if a0 and rng is None else None),
            a0_range_fm=((float(rng.group(1)), float(rng.group(2)))
                         if rng is not None else None),
            a0_err_fm=(float(r["a0_err_fm"]) if r["a0_err_fm"].strip() else None),
            params=params, note=text, q_range_fm=r["q_range_fm"],
            table=r["table"], page=int(r["page"])))
    n6 = sum(1 for r in out if r["nucleus"] == "6Li")
    n7 = sum(1 for r in out if r["nucleus"] == "7Li")
    if (n6, n7) != (13, 11):
        raise RuntimeError("%s: expected 13 6Li and 11 7Li rows, found %d / %d"
                           % (path, n6, n7))
    return out


# ------------------------------------------------------------ statistics

def chi2_sf_fallback(chi2, ndf):
    """Q(ndf/2, chi2/2) without scipy -- Numerical Recipes' gammq, as
    t1_nmc_f2d.py's `chi2_sf_fallback` (series for P below a + 1, the
    continued fraction for Q itself above), re-implemented here so this
    harness imports no other harness."""
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
    """Upper-tail probability of chi2 with ndf dof (scipy if present)."""
    try:
        from scipy.stats import chi2 as _c
    except ImportError:                             # pragma: no cover
        return chi2_sf_fallback(chi2, ndf)
    return float(_c.sf(float(chi2), int(ndf)))


def shape_fit(d, sig, m):
    """chi2(s) = sum [(d - s m)/sig]^2 minimised over ONE free scale s, in
    closed form: (s_hat, its 1-sigma, chi2_min, pulls)."""
    smm = sum((mi / si) ** 2 for mi, si in zip(m, sig))
    sdm = sum(di * mi / si ** 2 for di, mi, si in zip(d, m, sig))
    if not smm > 0.0:
        raise RuntimeError("the model is zero at every point")
    s_hat = sdm / smm
    pulls = [(di - s_hat * mi) / si for di, mi, si in zip(d, m, sig)]
    return dict(s_hat=s_hat, s_err=1.0 / math.sqrt(smm),
                chi2=sum(p * p for p in pulls), pulls=pulls)


def chi2_fixed(d, sig, m, s=1.0):
    return sum(((di - s * mi) / si) ** 2 for di, mi, si in zip(d, m, sig))


# ------------------------------------------------------------ the tree side

def tree_side():
    """The tree's F_m: its shape parameters, its slope radius (UNCOMPARED,
    record (r7)), and 7Li's absence.  The keys fm0, r_mag_li6,
    r_mag_li6_closed, fm_qz_fm, fm_b_fm and li7 are 2026-09-23's
    (python/tests/test_bench_nuclear.py reads them)."""
    import lipolgen as lg  # noqa: WPS433 (env.sh puts it on the path)
    hbarc = lg.HBARC_GEV_FM
    ff = lg.HoSpin1FF.for_ion(lg.li6())
    fm0 = ff.fm(0.0)
    eps = 1e-3
    r2 = -6.0 * (ff.fm((eps * hbarc) ** 2) / fm0 - 1.0) / eps ** 2
    o = ff.options
    qz, b = o.fm_qz_fm, o.fm_b_fm
    r2_closed = 6.0 / qz ** 2 + 1.5 * b ** 2
    try:
        lg.HoSpin1FF.for_ion(lg.li7())
        li7 = "HoSpin1FF.for_ion(li7()) constructed (unexpected: re-read rc.hpp)"
        li7_implemented = True
    except Exception as exc:
        li7 = "not implemented in the tree: HoSpin1FF.for_ion(li7()) throws (%s)" % (
            str(exc).splitlines()[0][:120])
        li7_implemented = False
    return dict(fm0=fm0, r_mag_li6=math.sqrt(r2), r_mag_li6_closed=math.sqrt(r2_closed),
                fm_qz_fm=qz, fm_b_fm=b, li7=li7, li7_implemented=li7_implemented,
                mu_n=o.mu_n, m_a_gev=o.m_a_gev,
                tail_tensor_scale=o.tail_tensor_scale, hbarc_gev_fm=hbarc,
                m_li6_gev=lg.li6().mass(), m_li7_gev=lg.li7().mass(),
                m_proton_gev=lg.PROTON_MASS, provenance=ff.provenance())


def tree_fmag(q2_fm2, ff=None):
    """[F_m(t)/F_m(0) at t = q^2 (hbar c)^2 for each q^2], read through the
    SHIPPED HoSpin1FF.for_ion(li6()).fm, and the largest |difference| from
    the closed form (1 - q^2/q_z^2) exp(-q^2 b^2/4)."""
    import lipolgen as lg  # noqa: WPS433
    ff = lg.HoSpin1FF.for_ion(lg.li6()) if ff is None else ff
    hc2 = lg.HBARC_GEV_FM ** 2
    o = ff.options
    fm0 = ff.fm(0.0)
    vals, worst = [], 0.0
    for q2 in q2_fm2:
        f = ff.fm(q2 * hc2) / fm0
        closed = (1.0 - q2 / o.fm_qz_fm ** 2) * math.exp(-0.25 * q2 * o.fm_b_fm ** 2)
        worst = max(worst, abs(f - closed))
        vals.append(f)
    return vals, worst


def q2_at_pi(e0_mev, m_a_gev, hbarc_gev_fm):
    """Record (r3): the four-momentum transfer squared [fm^-2] of elastic
    scattering at theta = pi, q^2 = 4 E0 E'/(hbar c)^2 with E' = E0/eta,
    eta = 1 + 2 E0/M_A (RFY Eq. (1) at theta = pi; m_e neglected)."""
    e0 = 1e-3 * e0_mev
    e1 = e0 / (1.0 + 2.0 * e0 / m_a_gev)
    return 4.0 * e0 * e1 / hbarc_gev_fm ** 2


# ------------------------------------------------------------ RFY's own models

def rfy_ho_j(q2, a0, mass_number):
    """RFY p. 868, harmonic well: (<j0>, <j2>) with x = q^2 a0^2/4 and the
    c.m. factor e^(x/A)."""
    x = q2 * a0 * a0 / 4.0
    g = math.exp(-x) * math.exp(x / mass_number)
    return (1.0 - 2.0 * x / 3.0) * g, (2.0 / 3.0) * x * g


def rfy_fn(q2):
    """RFY p. 868: F_N(q^2) = (1 + q^2/7.1)^-1."""
    return 1.0 / (1.0 + q2 / RFY_FN_SCALE_FM2)


def rfy_eq17_fm2(q2, a0, r1, mass_number, r2=0.0):
    """RFY Eq. (17): [(<j0> + R1 <j2>)^2 + R2^2 <j2>^2] F_N^2."""
    j0, j2 = rfy_ho_j(q2, a0, mass_number)
    return ((j0 + r1 * j2) ** 2 + (r2 * j2) ** 2) * rfy_fn(q2) ** 2


def rfy_eq18_fm2(q2, a0, a, b, c, mass_number):
    """RFY Eq. (18) (7Li, the 0.478-MeV level included):
    [A <j0>^2 + B <j0><j2> + C <j2>^2] F_N^2."""
    j0, j2 = rfy_ho_j(q2, a0, mass_number)
    return (a * j0 * j0 + b * j0 * j2 + c * j2 * j2) * rfy_fn(q2) ** 2


def rfy_alpha_d_fm2(q2, q0, mass_number, rd=RFY_DEUTERON_RMS_FM):
    """RFY p. 869, alpha+d (1s): F_M = F_d exp[-(q^2/4 q0^2)(1 - 2/A)], F_d
    the Gaussian of rms radius rd; squared."""
    f = math.exp(-q2 * rd * rd / 6.0) * math.exp(
        -(q2 / (4.0 * q0 * q0)) * (1.0 - 2.0 / mass_number))
    return f * f


def closures(li6, li7, a0_rows):
    """Record (r5): every fit of Tables VII and VIII the header names,
    re-evaluated at the vendored points with NO free parameter."""
    out = []
    for r in a0_rows:
        nuc, fit, pr = r["nucleus"], r["fit"], r["params"]
        pts = li6 if nuc == "6Li" else li7
        if r["table"] == "VII" and fit in ("A", "B", "C"):
            model = [rfy_eq17_fm2(p["q2_fm2"], r["a0_fm"], pr["R1"][0], MASS_NUMBER[nuc])
                     for p in pts]
            label = "%s VII (%s) a0=%s R1=%.2f" % (nuc, fit, r["a0_as_printed"], pr["R1"][0])
        elif r["table"] == "VII" and fit == "alpha+d(1s)":
            model = [rfy_alpha_d_fm2(p["q2_fm2"], pr["q0"][0], MASS_NUMBER[nuc])
                     for p in pts]
            label = "%s VII alpha+d(1s) q0=%g" % (nuc, pr["q0"][0])
        elif r["table"] == "VIII" and fit in ("A", "B", "C"):
            model = [rfy_eq18_fm2(p["q2_fm2"], r["a0_fm"], pr["A"][0], pr["B"][0],
                                  pr["C"][0], MASS_NUMBER[nuc]) for p in pts]
            label = "%s VIII (%s) a0=%s" % (nuc, fit, r["a0_as_printed"])
        else:
            continue
        chi2 = chi2_fixed([p["fm2"] for p in pts], [p["fm2_err"] for p in pts], model)
        paper = pr["chi2"][0]
        out.append(dict(name=label, nucleus=nuc, table=r["table"], fit=fit,
                        a0_fm=r["a0_fm"],
                        params={k: v[0] for k, v in pr.items()},
                        n_points=len(pts), chi2_table=chi2, chi2_paper=paper,
                        dof_paper=int(pr["dof"][0]),
                        holds=chi2 <= paper + CLOSURE_MARGIN, model=model))
    return out


def _solve(mat, vec):
    """Gaussian elimination with partial pivoting (a 3 x 3 normal system)."""
    n = len(vec)
    a = [list(row) + [v] for row, v in zip(mat, vec)]
    for i in range(n):
        k = max(range(i, n), key=lambda r: abs(a[r][i]))
        a[i], a[k] = a[k], a[i]
        if a[i][i] == 0.0:
            raise RuntimeError("singular normal equations")
        for r in range(i + 1, n):
            f = a[r][i] / a[i][i]
            for c in range(i, n + 1):
                a[r][c] -= f * a[i][c]
    x = [0.0] * n
    for i in reversed(range(n)):
        x[i] = (a[i][n] - sum(a[i][c] * x[c] for c in range(i + 1, n))) / a[i][i]
    return x


def closure_refit_7li(li7, a0_rows):
    """POST-HOC DIAGNOSTIC (header addendum (a)): for every harmonic-well row
    of Table VIII, the weighted least-squares minimum over (A, B, C) of Eq.
    (18) at the printed a0 on the 24 Table II points -- linear, closed form --
    against the paper's chi2 (a minimum over a superset of the points)."""
    out = []
    for r in a0_rows:
        if not (r["nucleus"] == "7Li" and r["table"] == "VIII"
                and r["fit"] in ("A", "B", "C")):
            continue
        a0, pr = r["a0_fm"], r["params"]
        basis = []
        for p in li7:
            j0, j2 = rfy_ho_j(p["q2_fm2"], a0, MASS_NUMBER["7Li"])
            fn2 = rfy_fn(p["q2_fm2"]) ** 2
            basis.append((j0 * j0 * fn2, j0 * j2 * fn2, j2 * j2 * fn2))
        w = [1.0 / p["fm2_err"] ** 2 for p in li7]
        mat = [[sum(wi * bi[i] * bi[j] for wi, bi in zip(w, basis)) for j in range(3)]
               for i in range(3)]
        vec = [sum(wi * bi[i] * p["fm2"] for wi, bi, p in zip(w, basis, li7))
               for i in range(3)]
        coef = _solve(mat, vec)
        chi2 = sum(((p["fm2"] - sum(c * b for c, b in zip(coef, bi))) / p["fm2_err"]) ** 2
                   for p, bi in zip(li7, basis))
        paper = pr["chi2"][0]
        slack = REFIT_SLACK_A if r["fit"] == "A" else 0.0
        printed = [pr[k] for k in ("A", "B", "C")]
        out.append(dict(
            name="7Li VIII (%s) a0=%s" % (r["fit"], r["a0_as_printed"]), fit=r["fit"],
            a0_fm=a0, chi2_min=chi2, chi2_paper=paper, slack=slack,
            holds=chi2 <= paper + slack, abc_fit=coef,
            abc_printed=[v for v, _ in printed],
            abc_shift_in_printed_errors=[(c - v) / e for c, (v, e) in zip(coef, printed)]))
    return out


def a0_record(a0_rows):
    """Record (r6): the a0 rows verbatim, each with sqrt(5/2) a0 -- RFY p.
    868's <r^2> = (5/2) a0^2, the 1p-shell orbital rms of the fitted well,
    NOT the slope of F_M and NEVER compared with the tree."""
    k = math.sqrt(2.5)
    out = []
    for r in a0_rows:
        row = dict(nucleus=r["nucleus"], fit=r["fit"], table=r["table"],
                   a0_as_printed=r["a0_as_printed"], a0_err_fm=r["a0_err_fm"])
        if r["a0_fm"] is not None:
            row["r_1p_orbital_rms_fm"] = k * r["a0_fm"]
            if r["a0_err_fm"] is not None:
                row["r_1p_orbital_rms_err_fm"] = k * r["a0_err_fm"]
        elif r["a0_range_fm"] is not None:
            row["r_1p_orbital_rms_range_fm"] = [k * v for v in r["a0_range_fm"]]
        out.append(row)
    return out


# ------------------------------------------------------------ the measurement

def measure(ff_path=DATA_FF, a0_path=DATA_A0):
    design = design_block_sha256()
    li6, li7 = load_ff(ff_path)
    a0_rows = load_a0(a0_path)
    tree = tree_side()

    d = [p["fm2"] for p in li6]
    sig = [p["fm2_err"] for p in li6]
    fmag, worst = tree_fmag([p["q2_fm2"] for p in li6])
    if not worst <= SHAPE_CHECK_ABS:
        raise RuntimeError("HoSpin1FF.fm / fm(0) differs from the closed-form "
                           "F_mag by %.3g: the harness no longer reads the "
                           "shape its header describes" % worst)
    m = [f * f for f in fmag]
    ndf = len(li6) - 1
    if ndf != NDF_VERDICT:
        raise RuntimeError("ndf %d != the design block's %d" % (ndf, NDF_VERDICT))
    fit = shape_fit(d, sig, m)
    p = chi2_sf(fit["chi2"], ndf)
    status = "pass" if p >= P_MIN else "fail"

    chi2_abs = chi2_fixed(d, sig, m, 1.0)
    p_abs = chi2_sf(chi2_abs, len(li6))

    # (r3) q^2 recomputed from E0 at theta = pi
    q2_re = {}
    for nuc, pts, mass in (("6Li", li6, tree["m_li6_gev"]), ("7Li", li7, tree["m_li7_gev"])):
        q2_re[nuc] = [q2_at_pi(pt["e0_mev"], mass, tree["hbarc_gev_fm"]) for pt in pts]
    dq = [abs(a - pt["q2_fm2"]) for nuc, pts in (("6Li", li6), ("7Li", li7))
          for a, pt in zip(q2_re[nuc], pts)]
    fmag_re, _ = tree_fmag(q2_re["6Li"])
    fit_re = shape_fit(d, sig, [f * f for f in fmag_re])
    p_re = chi2_sf(fit_re["chi2"], ndf)

    clos = closures(li6, li7, a0_rows)
    refit = closure_refit_7li(li7, a0_rows)
    return dict(design_block_sha256=design, li6=li6, li7=li7, a0_rows=a0_rows,
                tree=tree, fmag=fmag, model=m, shape_check_abs=worst,
                fit=fit, ndf=ndf, p=p, status=status,
                chi2_abs=chi2_abs, ndf_abs=len(li6), p_abs=p_abs,
                q2_recomputed=q2_re, q2_check_max=max(dq),
                q2_check_holds=max(dq) <= Q2_CHECK_FM2,
                fit_recomputed=fit_re, p_recomputed=p_re,
                closures=clos, closure_refit=refit,
                a0_record=a0_record(a0_rows))


def _status_word(ok):
    return "pass" if ok else "fail"


def run(verbose=True, ff_path=DATA_FF, a0_path=DATA_A0):
    t0 = time.perf_counter()
    m = measure(ff_path, a0_path)
    t, f = m["tree"], m["fit"]
    li6, li7 = m["li6"], m["li7"]
    c6 = [c for c in m["closures"] if c["nucleus"] == "6Li"]
    c7 = [c for c in m["closures"] if c["nucleus"] == "7Li"]

    generator = ("HoSpin1FF F_mag^2 (q_z = %.2f fm^-1, b = %.2f fm; unfitted "
                 "starting values) at RFY Table I's 7 points: chi2 = %.2f / %d "
                 "with one free scale s = %.3f +- %.3f (p = %.3g); s = 1: chi2 = "
                 "%.2f / %d (p = %.3g); 7Li: not implemented in the tree"
                 % (t["fm_qz_fm"], t["fm_b_fm"], f["chi2"], m["ndf"], f["s_hat"],
                    f["s_err"], m["p"], m["chi2_abs"], m["ndf_abs"], m["p_abs"]))
    reference = ("Rand-Frosch-Yearian, Phys. Rev. 144 (1966) 859, Table I: 6Li "
                 "F_M^2 at q^2 = %s-%s fm^-2 (q = %.2f-%.2f fm^-1), 7 points, "
                 "read from the page image; Table II (7Li, 24 points) recorded"
                 % (li6[0]["q2_as_printed"], li6[-1]["q2_as_printed"],
                    math.sqrt(li6[0]["q2_fm2"]), math.sqrt(li6[-1]["q2_fm2"])))
    tolerance = ("p(chi2, %d) >= %g with one free scale, the paper's errors, "
                 "diagonal (declared in the hashed design block)" % (NDF_VERDICT, P_MIN))

    subrows = [
        dict(name="li6_shape[free scale]", status=m["status"], in_verdict=True,
             chi2=f["chi2"], ndf=m["ndf"], p=m["p"], s_hat=f["s_hat"],
             s_err=f["s_err"], n_points=len(li6), pulls=f["pulls"]),
        dict(name="li6_absolute[s = 1]", status="recorded", in_verdict=False,
             chi2=m["chi2_abs"], ndf=m["ndf_abs"], p=m["p_abs"]),
        dict(name="li6_shape[q2 from E0]", status="recorded", in_verdict=False,
             chi2=m["fit_recomputed"]["chi2"], ndf=m["ndf"], p=m["p_recomputed"],
             s_hat=m["fit_recomputed"]["s_hat"]),
        dict(name="q2_column_check", status=_status_word(m["q2_check_holds"]),
             in_verdict=False, max_abs_dq2_fm2=m["q2_check_max"],
             limit_fm2=Q2_CHECK_FM2, n_rows=len(li6) + len(li7)),
        dict(name="li7[Table II]", status="not compared", in_verdict=False,
             reason=t["li7"], n_points=len(li7),
             q2_min_fm2=li7[0]["q2_fm2"], q2_max_fm2=li7[-1]["q2_fm2"]),
        dict(name="closure[6Li Table VII]",
             status=_status_word(all(c["holds"] for c in c6)), in_verdict=False,
             fits=[c["name"] for c in c6], chi2_table=[c["chi2_table"] for c in c6],
             chi2_paper=[c["chi2_paper"] for c in c6], n_points=len(li6),
             margin=CLOSURE_MARGIN),
        dict(name="closure[7Li Table VIII]",
             status=_status_word(all(c["holds"] for c in c7)), in_verdict=False,
             fits=[c["name"] for c in c7], chi2_table=[c["chi2_table"] for c in c7],
             chi2_paper=[c["chi2_paper"] for c in c7], n_points=len(li7),
             margin=CLOSURE_MARGIN),
        dict(name="closure_refit[7Li Table VIII, post-hoc]",
             status=_status_word(all(c["holds"] for c in m["closure_refit"])),
             in_verdict=False, fits=[c["name"] for c in m["closure_refit"]],
             chi2_min=[c["chi2_min"] for c in m["closure_refit"]],
             chi2_paper=[c["chi2_paper"] for c in m["closure_refit"]],
             slack_a_row=REFIT_SLACK_A, n_points=len(li7)),
    ]
    import lipolgen as lg  # noqa: WPS433
    config = dict(
        lipolgen=getattr(lg, "__version__", None), python=sys.version.split()[0],
        data_ff_sha256=FF_SHA256, data_a0_sha256=A0_SHA256,
        design_block_sha256=m["design_block_sha256"],
        fm_qz_fm=t["fm_qz_fm"], fm_b_fm=t["fm_b_fm"],
        tree_r_mag_li6_fm_uncompared=t["r_mag_li6"],
        radius_rule="retired 2026-10-26 (harness header)",
        q2_abscissa="printed (RFY Table I); recomputed recorded")
    points = dict(
        li6=[dict(row=p["row"], e0_mev=p["e0_mev"], q2_fm2=p["q2_fm2"],
                  fm2=p["fm2"], fm2_err=p["fm2_err"], tree_fmag2=mi,
                  pull=pu, q2_recomputed_fm2=qr)
             for p, mi, pu, qr in zip(li6, m["model"], f["pulls"],
                                       m["q2_recomputed"]["6Li"])],
        li7=[dict(row=p["row"], e0_mev=p["e0_mev"], q2_fm2=p["q2_fm2"],
                  fm2=p["fm2"], fm2_err=p["fm2_err"], q2_recomputed_fm2=qr)
             for p, qr in zip(li7, m["q2_recomputed"]["7Li"])])
    report = dict(name=NAME, generator=generator, reference=reference,
                  tolerance=tolerance, status=m["status"], subrows=subrows,
                  config=config, points=points, tree=t,
                  closures_detail=[{k: v for k, v in c.items() if k != "model"}
                                   for c in m["closures"]],
                  closure_refit_detail=m["closure_refit"],
                  a0_record=m["a0_record"], shape_check_abs=m["shape_check_abs"])
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        print("tree: %s" % t["provenance"][:160])
        print("tree: F_mag(q_z = %.4f fm^-1, b = %.4f fm), read through fm/fm(0); "
              "closed-form check %.1e; slope radius r_mag(6Li) = %.4f fm "
              "(UNCOMPARED, record (r7))" % (t["fm_qz_fm"], t["fm_b_fm"],
                                             m["shape_check_abs"], t["r_mag_li6"]))
        for p in points["li6"]:
            print("  6Li I.%d  E0 = %3.0f MeV  q^2 = %.2f fm^-2 (from E0: %.4f)  "
                  "F_M^2 = %.3f +- %.3f   tree F_mag^2 = %.3e   pull(s_hat) = %+.2f"
                  % (p["row"], p["e0_mev"], p["q2_fm2"], p["q2_recomputed_fm2"],
                     p["fm2"], p["fm2_err"], p["tree_fmag2"], p["pull"]))
        print("verdict: chi2 = %.4f / %d, s_hat = %.4f +- %.4f, p = %.4g -> %s "
              "(tolerance p >= %g)" % (f["chi2"], m["ndf"], f["s_hat"], f["s_err"],
                                       m["p"], m["status"].upper(), P_MIN))
        print("7Li: %s; %d Table II points recorded, not compared" % (t["li7"], len(li7)))
        for c in m["closures"]:
            print("  closure %-34s chi2(%2d table points) = %7.3f   paper %5.1f / "
                  "%d dof -> %s" % (c["name"], c["n_points"], c["chi2_table"],
                                    c["chi2_paper"], c["dof_paper"],
                                    "holds" if c["holds"] else "DOES NOT HOLD"))
        for c in m["closure_refit"]:
            print("  post-hoc refit %-22s (A, B, C) free at the printed a0: chi2_min(24) = "
                  "%7.3f   paper %5.1f%s -> %s; fitted A, B, C = %s (printed %s; "
                  "shift %s printed errors)" % (
                      c["name"], c["chi2_min"], c["chi2_paper"],
                      " + %.2f" % c["slack"] if c["slack"] else "",
                      "holds" if c["holds"] else "DOES NOT HOLD",
                      " / ".join("%.3f" % v for v in c["abc_fit"]),
                      " / ".join("%.2f" % v for v in c["abc_printed"]),
                      " / ".join("%+.2f" % v for v in c["abc_shift_in_printed_errors"])))
        for r in m["a0_record"]:
            if "r_1p_orbital_rms_fm" in r or "r_1p_orbital_rms_range_fm" in r:
                print("  a0 %s %-12s %-8s -> sqrt(5/2) a0 = %s fm (1p-orbital rms, "
                      "RFY p. 868; NOT compared)" % (
                          r["nucleus"], r["fit"], r["a0_as_printed"],
                          ("%.3f" % r["r_1p_orbital_rms_fm"]) if "r_1p_orbital_rms_fm" in r
                          else "%.3f-%.3f" % tuple(r["r_1p_orbital_rms_range_fm"])))
        for s in subrows:
            if s["name"] == "li6_shape[free scale]":
                g = "chi2 = %.4f / %d, s = %.4f +- %.4f" % (
                    s["chi2"], s["ndf"], s["s_hat"], s["s_err"])
                ref, tol = "RFY Table I, 7 points", "p >= %g" % P_MIN
                extra = " (p = %.4g)" % s["p"]
            elif s["name"] == "li6_absolute[s = 1]":
                g = "chi2 = %.4f / %d" % (s["chi2"], s["ndf"])
                ref, tol, extra = "RFY Table I, 7 points", "recorded", " (p = %.4g)" % s["p"]
            elif s["name"] == "li6_shape[q2 from E0]":
                g = "chi2 = %.4f / %d, s = %.4f" % (s["chi2"], s["ndf"], s["s_hat"])
                ref, tol = "q^2 = 4 E0 E'/(hbar c)^2 at theta = pi", "recorded"
                extra = " (p = %.4g)" % s["p"]
            elif s["name"] == "q2_column_check":
                g = "max |q2(E0) - q2 printed| = %.4f fm^-2 over %d rows" % (
                    s["max_abs_dq2_fm2"], s["n_rows"])
                ref, tol, extra = "RFY Tables I, II", "<= %g fm^-2" % Q2_CHECK_FM2, ""
            elif s["name"].startswith("closure_refit"):
                g = "chi2_min(24) = %s" % " / ".join("%.2f" % v for v in s["chi2_min"])
                ref = "paper chi2 %s" % " / ".join("%.1f" % v for v in s["chi2_paper"])
                tol, extra = ("chi2_min <= chi2(paper) (+ %g for (A)); post-hoc diagnostic"
                              % REFIT_SLACK_A), ""
            elif s["name"] == "li7[Table II]":
                g = "not implemented in the tree"
                ref, tol, extra = "RFY Table II, %d points" % s["n_points"], "n/a", ""
            else:
                g = "chi2(table) = %s" % " / ".join("%.2f" % v for v in s["chi2_table"])
                ref = "paper chi2 %s" % " / ".join("%.1f" % v for v in s["chi2_paper"])
                tol, extra = "chi2(table) <= chi2(paper) + %g each" % CLOSURE_MARGIN, ""
            print("SUBROW | %s[%s] | %s | %s | %s | %s%s" % (
                NAME, s["name"], g, ref, tol, s["status"].upper(), extra))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, generator, reference, tolerance, m["status"].upper()))
    return report


# ------------------------------------------------------------ the figure

def plot(report, out_dir=None):
    """F_M^2 against q on a log scale: RFY's 6Li points (ink) and 7Li points
    (open, recorded), the tree's F_mag^2 x s_hat (slot 'default') and at
    s = 1 (dotted), RFY's own fit (A) (published, ink dashed), the tree's zero
    q_z marked; under it the 6Li pulls against the scaled tree."""
    sys.path.insert(0, HERE)
    from _plots import INK_SECONDARY, MUTED, new_figure, save, status_title, style_series
    import numpy as np
    fig, (ax, axp) = new_figure(two_panel=True, height=4.0)
    t = report["tree"]
    li6, li7 = report["points"]["li6"], report["points"]["li7"]
    sub = {s["name"]: s for s in report["subrows"]}
    v = sub["li6_shape[free scale]"]
    s_hat = v["s_hat"]
    qz, b = t["fm_qz_fm"], t["fm_b_fm"]
    q = np.linspace(0.30, 2.05, 600)
    fmag2 = ((1.0 - q * q / qz ** 2) * np.exp(-0.25 * q * q * b * b)) ** 2
    fmag2 = np.where(fmag2 > 1e-7, fmag2, np.nan)      # the zero: break the log line
    ax.plot(q, s_hat * fmag2, **style_series(
        "default", label="tree F_mag^2 x s, s = %.2f" % s_hat))
    kw = style_series("default", label="tree F_mag^2, s = 1 (as generated)")
    kw.update(linestyle=":", linewidth=1.4)
    ax.plot(q, fmag2, **kw)
    a = [r for r in report["closures_detail"]
         if r["nucleus"] == "6Li" and r["table"] == "VII" and r["fit"] == "A"][0]
    a0, r1 = a["a0_fm"], a["params"]["R1"]
    rfy = [rfy_eq17_fm2(x * x, a0, r1, 6) for x in q]
    ax.plot(q, rfy, color=INK_SECONDARY, linestyle="--", linewidth=1.2, zorder=4,
            label="RFY fit (A), a0 = %.2f F (published)" % a0)
    ax.errorbar([math.sqrt(p["q2_fm2"]) for p in li6], [p["fm2"] for p in li6],
                yerr=[p["fm2_err"] for p in li6],
                **style_series("reference", label="6Li, RFY Table I"))
    k7 = style_series("reference", label="7Li, Table II (not in the tree)",
                      marker="s")
    k7.update(markerfacecolor="white", color=MUTED, zorder=4)
    ax.errorbar([math.sqrt(p["q2_fm2"]) for p in li7], [p["fm2"] for p in li7],
                yerr=[p["fm2_err"] for p in li7], **k7)
    ax.axvline(qz, color=MUTED, linestyle=(0, (1, 2)), linewidth=1.0, zorder=0,
               label="tree's zero q_z = %.2f fm^-1" % qz)
    axp.axvline(qz, color=MUTED, linestyle=(0, (1, 2)), linewidth=1.0, zorder=0)
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1.5)
    ax.set_xlim(0.30, 2.05)
    ax.set_ylabel("F_M^2(q), normalised to 1 at q = 0")
    ax.legend(loc="lower left", fontsize=7)
    axp.plot([math.sqrt(p["q2_fm2"]) for p in li6], [p["pull"] for p in li6],
             **style_series("default", marker="o"))
    axp.set_ylabel("6Li pull")
    axp.set_xlabel("q, four-momentum transfer at 180 deg  (fm^-1)")
    lim = max(3.0, max(abs(p["pull"]) for p in li6) * 1.2)
    axp.set_ylim(-lim, lim)
    status_title(fig, NAME, report["status"],
                 "6Li: chi2 = %.2f / %d with one free scale s = %.2f (p = %.3g); "
                 "tolerance p >= %g; 7Li not in the tree"
                 % (v["chi2"], v["ndf"], s_hat, v["p"], P_MIN))
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
