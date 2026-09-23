#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 1 -- EPIOS Table II source modes vs the spin
bookkeeping (tier T5/T1).

WHAT IS COMPARED.  The eight (P_z, P_zz) configurations of the COSY/ANKE
polarized DEUTERON source (EPIOS, arXiv:2510.10794 Table II; vendored in
`data/epios2026_table2_deuteron_source_modes.json` with provenance) are fed to
LiPolGen's (P_z, P_zz) -> (n+, n0, n-) inversion `spin1_populations` -- the
function `--pzz-mode typed` (`HelicityFlipOptions::use_explicit_pzz = true`)
fills a spin-1 category with.  ASSERTED, per mode: the populations are
non-negative and sum to 1, and the moments recomputed from them by two
independent routes (`moments_along_axis`, and the density matrix
`rho_from_populations` -> `vector_polarization` / `tensor_polarization`)
return the table's (P_z, P_zz) to 1e-12.

RECORDED, NOT ASSERTED: the distance of the equal-spin-temperature fill
(`spin_temperature_pzz(1, P_z)`, the `--pzz-mode ladder` default) from each
mode, at the ideal P_z and at the measured P_z^LEP.  That is the measurement
behind the plan's "the ladder fill is the wrong model for an EIC ion beam".
The EST fill has P_zz >= 0 for every P_z, so it cannot represent the three
modes with P_zz < 0 at any P_z; it is ALSO undefined in the tree at |P_z| = 1
(`populations_maxent` throws "maxent populations need |pz| < 1"), where the
closed form's limit is P_zz = 1 = the table's value.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: no b1 on either side.
  2. isoscalar denominator -- NOT APPLICABLE: no structure function.
  3. per nucleus / per nucleon -- both sides are PER ION: the reference is the
     polarization of one spin-1 deuteron; LiPolGen's `spin1_populations` is the
     substate occupation of one spin-1 ion.  Convention on both sides:
     CARTESIAN, P_z = n+ - n-, P_zz = n+ + n- - 2 n0 in [-2, 1] (NOT the
     spherical t20 -- 06_critic.md sec. 5.1 trap T-8).  The reference's P_zz
     column is IDEAL (the RF-transition scheme's value), not measured; only
     P_z^LEP is measured.

Run:  source env.sh && python3 validation/benchmarks/t5_epios_source_modes.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import json
import math
import os
import sys
import traceback
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "epios2026_table2_deuteron_source_modes.json")

NAME = "t5_epios_source_modes"
TOL_ROUNDTRIP = 1e-12


def _lib():
    import lipolgen._lipolgen as _l  # noqa: WPS433 (import here: env.sh)
    return _l


def load_table(path=DATA):
    with open(path) as f:
        doc = json.load(f)
    rows = []
    for r in doc["rows"]:
        rows.append(dict(r, pz=Fraction(r["pz_ideal"]),
                         pzz=Fraction(r["pzz_ideal"])))
    return rows


def est_pzz(_l, pz):
    """(tree value or None, closed-form value, note)."""
    closed = 2.0 - math.sqrt(4.0 - 3.0 * float(pz) ** 2)
    try:
        return _l.spin_temperature_pzz(1.0, float(pz)), closed, ""
    except Exception as exc:  # the tree refuses |pz| = 1
        return None, closed, "tree throws: %s" % exc


def measure(path=DATA):
    _l = _lib()
    out = []
    for r in load_table(path):
        pz, pzz = float(r["pz"]), float(r["pzz"])
        pops = list(_l.spin1_populations(pz, pzz))
        m = _l.moments_along_axis(1.0, pops)
        rho = _l.rho_from_populations(1.0, pops)
        vz = _l.vector_polarization(rho, 1.0)[2]
        tz = _l.tensor_polarization(rho, 1.0)
        # exact rationals of the closed-form inversion, for the record only
        exact = [(2 + 3 * r["pz"] + r["pzz"]) / 6, (1 - r["pzz"]) / 3,
                 (2 - 3 * r["pz"] + r["pzz"]) / 6]
        est_tree, est_closed, est_note = est_pzz(_l, pz)
        est_lep_tree, est_lep_closed, _ = est_pzz(_l, r["pz_lep"])
        est_val = est_tree if est_tree is not None else est_closed
        if est_tree is not None:
            lad = list(_l.spin_temperature_ladder(1.0, pz).populations)
        else:  # |pz| = 1: the pure state is the only fill
            lad = [1.0, 0.0, 0.0] if pz > 0 else [0.0, 0.0, 1.0]
        out.append(dict(
            mode=r["mode"], pz=pz, pzz=pzz, pops=pops,
            pops_exact=[str(x) for x in exact],
            sum_err=abs(sum(pops) - 1.0),
            min_pop=min(pops),
            res_axis=max(abs(m.vector - pz), abs(m.tensor - pzz)),
            res_rho=max(abs(vz - pz), abs(tz - pzz)),
            pz_lep=r["pz_lep"], pz_lep_stat=r["pz_lep_stat"],
            realisation=(r["pz_lep"] / pz) if pz != 0 else None,
            est_pzz=est_val, est_from_tree=est_tree is not None,
            est_note=est_note,
            est_minus_ideal=est_val - pzz,
            est_pzz_at_lep=(est_lep_tree if est_lep_tree is not None
                            else est_lep_closed),
            ladder_pops=lad,
            ladder_linf=max(abs(a - b) for a, b in zip(lad, pops)),
        ))
    return out


def run(path=DATA, verbose=True):
    rows = measure(path)
    worst = max(max(x["res_axis"], x["res_rho"], x["sum_err"]) for x in rows)
    positive = all(x["min_pop"] >= 0.0 for x in rows)
    status = "pass" if (worst <= TOL_ROUNDTRIP and positive) else "fail"
    on_est = [x["mode"] for x in rows if abs(x["est_minus_ideal"]) < 1e-9]
    report = dict(
        name=NAME,
        generator="max moment round-trip residual %.3e over 8 modes "
                  "(min population %.3g)" % (worst, min(x["min_pop"] for x in rows)),
        reference="EPIOS Table II ideal (P_z, P_zz), 8 modes",
        tolerance=TOL_ROUNDTRIP,
        status=status,
        worst_residual=worst,
        rows=rows,
        est_modes_matched=on_est,
    )
    if verbose:
        print("mode  Pz      Pzz    (p+, p0, p-) exact      round-trip  "
              "EST Pzz   EST-ideal  EST@PzLEP  Linf(ladder)  PzLEP/Pz")
        for x in rows:
            print("%3d  %+7.4f %+5.2f  (%s)  %.1e  %s  %+8.5f   %8.5f   %8.5f    %s" % (
                x["mode"], x["pz"], x["pzz"], ", ".join(x["pops_exact"]),
                max(x["res_axis"], x["res_rho"]),
                ("%8.5f" % x["est_pzz"]) + ("" if x["est_from_tree"] else "*"),
                x["est_minus_ideal"], x["est_pzz_at_lep"], x["ladder_linf"],
                "n/a" if x["realisation"] is None else "%.3f" % x["realisation"]))
        print("  * |P_z| = 1: the tree's spin_temperature_pzz throws "
              "(populations_maxent needs |pz| < 1); the closed-form limit is shown.")
        print("  EST reproduces the ideal P_zz (to 1e-9) at modes %s only." % on_est)
        print("REPORT | %s | %s | %s | tol %.0e | %s" % (
            NAME, report["generator"], report["reference"], TOL_ROUNDTRIP,
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
