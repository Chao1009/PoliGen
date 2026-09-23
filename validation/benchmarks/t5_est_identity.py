#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 9 -- the equal-spin-temperature identity
(tier T5).

WHAT IS COMPARED.  `spin_temperature_pzz(1, P_z)` -- the rank-2 moment the
`--pzz-mode ladder` default records as `pzz_true` -- against the polarized-
target community's closed form for a spin-1 system in thermal equilibrium at
one spin temperature,

    P_zz = 2 - sqrt(4 - 3 P_z^2),

which follows from p_m ~ exp(beta m): with x = cosh(beta),
P_z = 2 sinh(beta)/(2x + 1) and P_zz = (2x - 2)/(2x + 1); eliminating beta
gives the line above.  The reference is ANALYTIC (vendored as a formula with
its sources in `data/est_identity_sources.json`); no measured number enters.

CITATIONS.
  * J. Koivuniemi et al., "Polarization Build Up in COMPASS 6LiD Target",
    SPIN 2004 proceedings
    (wwwcompass.cern.ch/compass/publications/technical/koivuniemi_spin2004.pdf):
    D and 6Li have J = 1 and 7Li J = 3/2; the nuclear polarizations follow the
    Brillouin function of the common spin temperature ("all the nuclei share
    the same spin temperature").  This is the construction of
    `spin_temperature_ladder(j, pz)`: p_m ~ t^m, P_z = B_J(J ln t).
  * D. Keller, D. Crabb, D. Day, "Enhanced Tensor Polarization in Solid-
    State Targets", NIM A 981 (2020) 164504, arXiv:2008.09515: the EST
    relation above is the
    equilibrium baseline that RF manipulation of the quadrupole-split line
    deliberately breaks.  Abstract read; the body's (P_z, P_zz) table is NOT
    read here.
  * The 7Li statement (J = 3/2): the same ladder with J = 3/2 is asserted here
    to be the Brillouin function B_{3/2} and geometric in m; its rank-2
    moment against the partition-function closed form.

THE COMPASS NUMBERS, CHECKED (recorded AND asserted at the printed digit;
`data/est_identity_sources.json`).  Koivuniemi's Fig. 1 prints the deuteron
populations 63.6 / 25.9 / 10.6 % at P = 53.0 % and "T = 1/2(1 - 3p0) = 11 %"
-- NOTE THE 1/2: their T is HALF the Cartesian P_zz = 1 - 3 p0 this tree uses.
From the tree's own ladder at P_z = 0.530 and the common spin temperature
(beta_i = beta_d (mu_i/J_i)/(mu_d/J_d)), the harness recomputes the
populations, T, T_S = h f0/(kB beta_d), and P(6Li), P(7Li) -- the latter two
by inverting `spin_temperature_ladder(j, .)` at the scaled ladder ratio.

THE LIMITATION, stated where the identity is: in the COMPASS 6LiD target the
unsplit NMR line measures p+ - p- only, so the tensor polarization there was
ESTIMATED from EST (Fig. 1), not measured; no measurement of the tensor
polarization of 6Li was found by the searches of 06_critic.md secs. 4.1-4.2
(a finding of absence, not a proof).  An identity is not a measurement.

WHAT WAS MEASURED (2026-09-23).  The critic measured |diff| <= 7.7e-15 at five
points.  On a 20001-point grid of [-1, 1] (19999 interior points) the
maximum is 1.138e-14 at P_z = -0.5996, and 496 of the 19999 points exceed
1e-14: the requested 1e-14 does NOT hold everywhere.  The residual is the bisection of `populations_maxent` (bracket
width < tol = 1e-13 in beta), bounded by 0.5 * 1e-13 * max|dP_zz/dbeta|
= 0.5e-13 * 0.42609 = 2.13e-14 (the maximum of 6 sinh b/(2 cosh b + 1)^2,
at b = 1.113, P_z = 0.621); the closed form's own double rounding is
<= 3.2e-16 (checked against 40-digit arithmetic).  The asserted tolerance is
that bound plus 8 x 2.2e-16 of rounding, TOL = 2.31e-14.  At |P_z| = 1 the tree THROWS
("maxent populations need |pz| < 1"); the closed-form limit there is 1.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE (no b1).
  2. isoscalar denominator -- NOT APPLICABLE (no structure function).
  3. per nucleus / per nucleon -- both sides are PER NUCLEUS (per ion):
     the moments of one spin-J nucleus, Cartesian P_zz = <3 Jz^2 - 2> at
     J = 1, P_zz in [-2, 1].

Run:  source env.sh && python3 validation/benchmarks/t5_est_identity.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import math
import os
import sys
import traceback

NAME = "t5_est_identity"
N_GRID = 20001                      # grid of [-1, 1]; the endpoints are tested apart
BISECTION_TOL_BETA = 1e-13          # populations_maxent's default `tol`
MAX_DPZZ_DBETA = 0.4260945789955075 # max_b 6 sinh b / (2 cosh b + 1)^2
TOL = 0.5 * BISECTION_TOL_BETA * MAX_DPZZ_DBETA + 8 * 2.2e-16  # = 2.31e-14
TOL_REQUESTED = 1e-14


def closed_form(pz):
    return 2.0 - math.sqrt(4.0 - 3.0 * pz * pz)


def brillouin(j, y):
    a = (2.0 * j + 1.0) / (2.0 * j)
    b = 1.0 / (2.0 * j)
    return a / math.tanh(a * y) - b / math.tanh(b * y)


def est_scan(n=N_GRID):
    import lipolgen._lipolgen as _l
    worst, at, n_over = 0.0, None, 0
    for i in range(1, n - 1):
        pz = -1.0 + 2.0 * i / (n - 1)
        d = abs(_l.spin_temperature_pzz(1.0, pz) - closed_form(pz))
        if d > TOL_REQUESTED:
            n_over += 1
        if d > worst:
            worst, at = d, pz
    endpoints = {}
    for pz in (-1.0, 1.0):
        try:
            endpoints[pz] = "returned %r" % _l.spin_temperature_pzz(1.0, pz)
        except Exception as exc:
            endpoints[pz] = "throws: %s" % exc
    return dict(worst=worst, at=at, n_points=n - 2, n_over_requested=n_over,
                endpoints=endpoints)


def ladder_j32(pzs=(-0.9, -0.5, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)):
    """7Li: the J = 3/2 ladder is geometric in m and IS the Brillouin
    function B_{3/2}; its rank-2 moment against the partition function."""
    import lipolgen._lipolgen as _l
    worst_b = worst_geo = worst_t = 0.0
    for pz in pzs:
        lad = _l.spin_temperature_ladder(1.5, pz)
        p, t = list(lad.populations), lad.t
        beta = math.log(t)
        worst_b = max(worst_b, abs(brillouin(1.5, 1.5 * beta) - pz))
        for a, b in zip(p, p[1:]):          # m = +3/2 ... -3/2
            worst_geo = max(worst_geo, abs(a / b - t) / t)
        ms = (1.5, 0.5, -0.5, -1.5)
        w = [math.exp(beta * m) for m in ms]
        z = sum(w)
        m2 = sum(m * m * x for m, x in zip(ms, w)) / z
        t_closed = (3.0 * m2 - 3.75) / 3.0   # T = <3Jz^2 - J(J+1)>/3
        worst_t = max(worst_t, abs(_l.spin_temperature_pzz(1.5, pz) - t_closed))
    return dict(brillouin=worst_b, geometric=worst_geo, rank2=worst_t)


DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                    "est_identity_sources.json")


def _pz_at_ratio(_l, j, t_target):
    """Invert the tree's ladder: the P_z whose `spin_temperature_ladder(j, .)`
    ratio is t_target (bisection; the ratio is monotone in P_z)."""
    lo, hi = -1.0 + 1e-12, 1.0 - 1e-12
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _l.spin_temperature_ladder(j, mid).t < t_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def compass_crosscheck(path=DATA):
    import json
    import lipolgen._lipolgen as _l
    with open(path) as f:
        d = json.load(f)
    k = d["koivuniemi2004"]
    mu = d["magnetic_moments_muN"]
    c = d["physical_constants"]
    fig = k["fig1_downstream"]
    pz = fig["P_percent"] / 100.0
    lad = _l.spin_temperature_ladder(1.0, pz)
    pops = [100.0 * x for x in lad.populations]            # m = +1, 0, -1
    printed = [fig["p_aligned_percent"], fig["p0_percent"],
               fig["p_anti_percent"]]
    beta_d = math.log(lad.t)
    g = lambda s: mu[s]["mu"] / mu[s]["J"]                 # noqa: E731
    t_s_mk = c["h_Js"] * k["f0_deuteron_MHz"] * 1e6 / (c["kB_JK"] * beta_d) * 1e3
    p6 = _pz_at_ratio(_l, 1.0, math.exp(beta_d * g("6Li") / g("deuteron")))
    p7 = _pz_at_ratio(_l, 1.5, math.exp(beta_d * g("7Li") / g("deuteron")))
    pzz = _l.spin_temperature_pzz(1.0, pz)
    return dict(
        pops=pops, printed=printed,
        pops_maxdiff=max(abs(a - b) for a, b in zip(pops, printed)),
        pzz_cartesian=pzz, t_half=0.5 * pzz, t_printed=k["T_half_convention_percent"] / 100.0,
        t_s_mk=t_s_mk, t_s_printed=k["T_S_mK"],
        p6=p6, p6_printed=k["P_6Li_percent"] / 100.0,
        p7=p7, p7_printed=k["P_7Li_percent"] / 100.0,
    )


def run(verbose=True, n=N_GRID):
    s = est_scan(n)
    j32 = ladder_j32()
    cx = compass_crosscheck()
    ok_end = all(v.startswith("throws") for v in s["endpoints"].values())
    status = "pass" if (s["worst"] <= TOL and ok_end) else "fail"
    report = dict(
        name=NAME,
        generator="max|spin_temperature_pzz(1,Pz) - closed form| = %.3e at "
                  "Pz = %.4f over %d interior points" % (s["worst"], s["at"],
                                                         s["n_points"]),
        reference="P_zz = 2 - sqrt(4 - 3 P_z^2) (EST, analytic)",
        tolerance=TOL, status=status, scan=s, j32=j32, compass=cx,
        requested_tolerance=TOL_REQUESTED,
        meets_requested=s["worst"] <= TOL_REQUESTED,
    )
    if verbose:
        print("EST identity J = 1: max |diff| = %.4e at P_z = %.6f over %d "
              "interior points; %d points exceed the requested %.0e" % (
                  s["worst"], s["at"], s["n_points"], s["n_over_requested"],
                  TOL_REQUESTED))
        print("  bisection bound 0.5*tol*max|dPzz/dbeta| = %.4e; asserted TOL "
              "= %.4e" % (0.5 * BISECTION_TOL_BETA * MAX_DPZZ_DBETA, TOL))
        for pz, v in sorted(s["endpoints"].items()):
            print("  P_z = %+.0f: %s (closed-form limit 1)" % (pz, v))
        print("7Li J = 3/2 ladder: |B_3/2(1.5 ln t) - P_z| <= %.2e, geometric "
              "ratio rel. err <= %.2e, rank-2 vs partition function <= %.2e" % (
                  j32["brillouin"], j32["geometric"], j32["rank2"]))
        print("COMPASS 6LiD (Koivuniemi SPIN 2004, Fig. 1): tree ladder at P = "
              "0.530 -> (%.2f, %.2f, %.2f) %% vs printed (%.1f, %.1f, %.1f) %% "
              "(max |diff| %.3f %%)" % (tuple(cx["pops"]) + tuple(cx["printed"])
                                         + (cx["pops_maxdiff"],)))
        print("  Cartesian P_zz = %.5f; their T = P_zz/2 = %.5f vs printed "
              "%.2f; T_S = %.4f mK vs %.2f; P(6Li) = %.4f vs %.3f; "
              "P(7Li) = %.4f vs %.3f" % (
                  cx["pzz_cartesian"], cx["t_half"], cx["t_printed"],
                  cx["t_s_mk"], cx["t_s_printed"], cx["p6"], cx["p6_printed"],
                  cx["p7"], cx["p7_printed"]))
        print("REPORT | %s | %s | %s | tol %.2e | %s" % (
            NAME, report["generator"], report["reference"], TOL, status.upper()))
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
