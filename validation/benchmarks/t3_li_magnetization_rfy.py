#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 6 -- Rand-Frosch-Yearian 6Li / 7Li
magnetization radii vs the tree's F_m (tier T3).  STATUS: BLOCKED.

WHY BLOCKED.  The reference is not on disk and is not open access.  Checked
2026-09-23: neither document is in PolarizedLithiumSim/refs/ (204 entries; no
file, README row or refs_dict.json key for either), and
docs/references/REFERENCES.md sec. 1 lists both as PLEASE DOWNLOAD:

  [#8] R. E. Rand, R. F. Frosch, M. R. Yearian, "Elastic Electron Scattering
       from the Magnetic Multipole Distributions of Li6, Li7, Be9, B10, B11
       and N14", Phys. Rev. 144 (1966) 859 -- https://doi.org/10.1103/PhysRev.144.859
  [#9] C. W. de Jager, H. de Vries, C. de Vries, "Nuclear charge- and
       magnetization-density-distribution parameters from elastic electron
       scattering", At. Data Nucl. Data Tables 14 (1974) 479 --
       https://doi.org/10.1016/S0092-640X(74)80002-1  (Table V, p. 502)

The UVa machine-readable archive carries NO magnetization rows (REFERENCES.md
#9; `FB_data.dat`, `SOG_data.dat`, `tabletr.txt` are charge-only).  The
survey (06_critic.md sec. 4.6) read ADNDT 14 Table V from a third-party
mirror; that copy is not a licensed open-access source and is NOT vendored
here, so the rows it transcribed are NOT used by this harness.

WHAT UNBLOCKS IT.  Download either document, transcribe the lithium rows of
ADNDT 14 Table V (or RFY66's own fit table) into
`validation/benchmarks/data/rfy1966_li_magnetization.csv` with a provenance
header (document, page, date, who transcribed) and columns
    nucleus,mu_nm,r_mag_fm,r_mag_err_fm,model_note,q_range_fm
one row per published fit -- BOTH 6Li rows (the free fit and the one with
a0 fixed from charge scattering; 06 sec. 4.6 reports they differ by ~22 %,
UNVERIFIED here) and every 7Li row.  The harness then runs unchanged.

THE COMPARISON IT WILL MAKE (written out so it is ready when the file lands):

  Tree side.  rc.hpp builds F_m(t) = (M_A/m_p)(mu/mu_N) F_mag(q),
  F_mag(q) = (1 - q^2/q_z^2) exp(-q^2 b^2/4), q_z = LI6_FF_FM_QZ_FM = 1.30
  fm^-1, b = LI6_FF_FM_B_FM = 1.85 fm (UNFITTED placeholders, rc.hpp), NO
  nucleon folding on F_m.  F_m(0) is anchored on the MEASURED moment mu
  (00_in_tree_checks.md C-3); the q -> 0 SLOPE is not anchored on anything:
      <r^2>_mag = -6 d[F_m/F_m(0)]/dq^2 |_(q=0) = 6/q_z^2 + (3/2) b^2 .
  This harness evaluates it from `HoSpin1FF.for_ion(li6()).fm` numerically
  (and checks it against the closed form).  For 7Li the tree HAS NO F_m:
  `HoSpin1FF.for_ion(li7())` throws (no measured-moment block), so the 7Li
  half of the comparison is "not implemented in the tree", reported as such.

  Reference side.  <r^2>^1/2_mag from the magnetization-density fits.
  Pass rule (6Li): the tree's r_mag lies inside the envelope of ALL the
  vendored 6Li rows, each widened by its own quoted error (the two model
  choices are a band, not a number -- 06 sec. 4.6).  7Li: recorded only.

  What it does NOT validate: the M1 shape beyond the fits' narrow q range
  (06 sec. 4.6 reports 0.85-1.39 fm^-1 for 6Li, UNVERIFIED here), C0, C2, or
  any tensor observable.  The q_z zero of F_mag is not tested by a radius.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: no b1, no A-scaling assumption.
  2. isoscalar denominator -- NOT APPLICABLE: no structure function.
  3. per nucleus / per nucleon -- both sides PER NUCLEUS: mu in mu_N of the
     whole nucleus, <r^2>^1/2_mag in fm of the nuclear magnetization density;
     the tree's F_m(0) = (M_A/m_p) mu is the POLRAD/Rosenbluth normalisation
     of the same per-nucleus moment.  Radii are compared as F_m/F_m(0), so
     the normalisation cancels.

Run:  source env.sh && python3 validation/benchmarks/t3_li_magnetization_rfy.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import csv
import math
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "rfy1966_li_magnetization.csv")

NAME = "t3_li_magnetization_rfy"
BLOCKED_MESSAGE = (
    "BLOCKED: reference not on disk and not open access. Download "
    "R. E. Rand, R. F. Frosch, M. R. Yearian, Phys. Rev. 144 (1966) 859, "
    "https://doi.org/10.1103/PhysRev.144.859 "
    "OR C. W. de Jager, H. de Vries, C. de Vries, At. Data Nucl. Data Tables "
    "14 (1974) 479 (Table V), https://doi.org/10.1016/S0092-640X(74)80002-1 "
    "(docs/references/REFERENCES.md sec. 1 items 8 and 9), then vendor the "
    "lithium rows as validation/benchmarks/data/rfy1966_li_magnetization.csv")


def tree_side():
    """The tree's q -> 0 magnetization radius, and 7Li's absence."""
    import lipolgen as lg  # noqa: WPS433
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
    except Exception as exc:
        li7 = "not implemented in the tree: HoSpin1FF.for_ion(li7()) throws (%s)" % (
            str(exc).splitlines()[0][:120])
    return dict(fm0=fm0, r_mag_li6=math.sqrt(r2), r_mag_li6_closed=math.sqrt(r2_closed),
                fm_qz_fm=qz, fm_b_fm=b, li7=li7)


def load_reference(path=DATA):
    rows = []
    with open(path) as f:
        lines = [l for l in f if l.strip() and not l.startswith("#")]
    for r in csv.DictReader(lines):
        rows.append(dict(nucleus=r["nucleus"].strip(), mu=float(r["mu_nm"] or "nan"),
                         r_mag=float(r["r_mag_fm"]),
                         err=float(r["r_mag_err_fm"] or 0.0),
                         note=r.get("model_note", ""), q_range=r.get("q_range_fm", "")))
    return rows


def run(path=DATA, verbose=True):
    tree = tree_side()
    if not os.path.isfile(path):
        report = dict(name=NAME,
                      generator="r_mag(6Li) = %.4f fm (tree F_m slope; UNCOMPARED)"
                                % tree["r_mag_li6"],
                      reference="none on disk", tolerance="n/a",
                      status="blocked", reason=BLOCKED_MESSAGE, tree=tree)
        if verbose:
            print(BLOCKED_MESSAGE)
            print("tree side, measured and NOT compared: r_mag(6Li) = %.4f fm "
                  "(closed form %.4f; q_z = %.2f fm^-1, b = %.2f fm, both unfitted); "
                  "7Li: %s" % (tree["r_mag_li6"], tree["r_mag_li6_closed"],
                               tree["fm_qz_fm"], tree["fm_b_fm"], tree["li7"]))
            print("REPORT | %s | %s | %s | %s | BLOCKED" % (
                NAME, report["generator"], report["reference"], report["tolerance"]))
        return report
    ref = load_reference(path)
    li6 = [r for r in ref if r["nucleus"].lower() in ("6li", "li6")]
    li7 = [r for r in ref if r["nucleus"].lower() in ("7li", "li7")]
    if not li6:
        raise RuntimeError("%s carries no 6Li row" % path)
    lo = min(r["r_mag"] - r["err"] for r in li6)
    hi = max(r["r_mag"] + r["err"] for r in li6)
    inside = lo <= tree["r_mag_li6"] <= hi
    status = "pass" if inside else "fail"
    report = dict(name=NAME,
                  generator="r_mag(6Li) = %.4f fm (tree F_m slope)" % tree["r_mag_li6"],
                  reference="6Li envelope [%.3f, %.3f] fm over %d rows" % (lo, hi, len(li6)),
                  tolerance="inside the envelope of all 6Li rows +- their errors",
                  status=status, tree=tree, li6_rows=li6, li7_rows=li7,
                  distance=(0.0 if inside else
                            (tree["r_mag_li6"] - hi if tree["r_mag_li6"] > hi
                             else tree["r_mag_li6"] - lo)))
    if verbose:
        for r in li6 + li7:
            print("%s  r_mag = %.3f +- %.3f fm  %s  q %s" % (
                r["nucleus"], r["r_mag"], r["err"], r["note"], r["q_range"]))
        print("7Li: %s" % tree["li7"])
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, report["generator"], report["reference"], report["tolerance"],
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
