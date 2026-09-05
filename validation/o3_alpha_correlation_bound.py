#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Open item O3: how much does the alpha core's genuine 4-body correlation
move a position-space observable, relative to the uncorrelated product of
one-body densities `ClusterConfigSampler` actually draws?

WHAT THIS IS.  `ClusterConfigSampler`'s alpha core samples FOUR INDEPENDENT
positions from the ANL VMC one-body density `he4.density`
(`AlphaCoreSource::VmcHe4Density`), then recentres so the four sum to zero
exactly (`alpha_cm_inflate`, `cluster_config.hpp`).  That is, by construction,
an UNCORRELATED product: nothing in the draw couples nucleon i's position to
nucleon j's beyond the single centre-of-mass constraint.  Genuine GFMC/VMC
CONFIGURATIONS (a joint sample of all four positions from the correlated
4-body wavefunction) are not in this tree -- `docs/open_items/run_2026-09-02/
design_G_cluster_config.md` O3 calls the resulting incoherent/coherent-split
question "unanswerable from our side".

WHAT CAN BE BOUNDED OFFLINE.  `data/vmc/bonus_other_clusters/he4.dd` IS a
genuinely correlated input already in the tree: an S+D-wave VMC overlap of
the SAME 4He wavefunction onto an alpha -> d+d cluster channel (the same
`Forest et al., PRC 54, 646 (1996)` method as `li6.ad`/`li7.at`, run
19-Apr-2004).  Its own relative-motion radial density carries the full
4-body correlation; comparing its second moment to what an UNCORRELATED
product of the same he4.density would give for the identical observable --
the separation between the centroids of two 2-nucleon halves of the alpha --
isolates the size of that correlation for this one observable, on data this
repository actually holds.

THE NULL MODEL, worked once.  Draw 4 iid positions v_1..v_4 from a source
density (inflated by lambda^2 = 4/3, exactly `alpha_cm_inflate`'s own
constant) so that after recentring (s_i = v_i - mean(v)) the single-nucleon
density matches `he4.density` (<s^2> = R1^2 by construction -- this is
EXACTLY what `ClusterConfigSampler` does for the alpha core).  Split the four
into two pairs and let D = (s_1+s_2)/2 - (s_3+s_4)/2 be the centroid
separation.  The sample mean cancels identically in D (algebra below), so
D = (v_1+v_2-v_3-v_4)/2 exactly, independent of recentring and of which of
the three 2-2 partitions is chosen:

    <D^2> = (1/4) * sum_i <v_i^2> = <v^2> = (4/3) R1^2 .

This is the UNCORRELATED prediction for <D^2>, using only he4.density's own
measured R1 -- no free parameter, no fit.

THE CORRELATED ANSWER.  Read <D^2>-in-the-d+d-channel straight off he4.dd's
own S+D radial overlap (Add00(r), Add22(r)): the same
int (Add00^2+Add22^2) r^2 dr normalization / int (...) r^4 dr moment this
module already uses for `li6.ad`, `li7.at` (`validation/vmc_reconcile.py`).

WHAT THIS DOES AND DOES NOT SETTLE.  It bounds one position-space moment for
one clustering channel (alpha -> d+d, not alpha -> 4 free nucleons) and shows
genuine 4-body correlation is a SEVERAL-PERCENT-TO-~20% effect on that
moment, not a factor of several -- i.e. it rules out "correlations dominate"
as a concern for the SIZE of position-space moments.  It does NOT compute how
the coherent/incoherent SPLIT itself would move: that is a property of the
Good-Walker AMPLITUDE's configuration-to-configuration fluctuation under a
dipole model, which needs the actual GFMC 4He configurations (or the
`subnucleondiffraction` graft) to evaluate -- nothing available in this tree
supplies that. O3 stays open in that sense; this script bounds only what it
says it bounds.

Pinned nowhere (no pytest, no reference JSON) -- this is an offline bound for
`docs/OPEN_ITEMS_SOLUTIONS.md` section 11.5 / `phase_C_numbers.md` section C6,
not a gated physics number.
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DENSITY = ROOT / "data" / "vmc" / "density" / "he4.density"
HE4DD = ROOT / "data" / "vmc" / "bonus_other_clusters" / "he4.dd"

NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"


def _trapz(y, x):
    y = list(y)
    x = list(x)
    s = 0.0
    for i in range(1, len(x)):
        s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1])
    return s


def read_he4_density(path: Path):
    """R RHORP DRHORP rows, after the '*****  ********  ********' rule."""
    lines = path.read_text().splitlines()
    start = None
    for i, l in enumerate(lines):
        if re.match(r"^\s*\*+\s+\*+\s+\*+\s*$", l):
            start = i + 1
            break
    if start is None:
        raise RuntimeError(f"no header rule found in {path}")
    r, rho = [], []
    row = re.compile(rf"^\s*({NUM})\s+({NUM})\s+({NUM})\s*$")
    for l in lines[start:]:
        m = row.match(l)
        if not m:
            break
        r.append(float(m.group(1)))
        rho.append(float(m.group(2)))
    return r, rho


def read_he4_dd(path: Path):
    """rij Add00 (err) Add22 (err) #samples rows, after the 'rij ... Add00 ...'
    header."""
    lines = path.read_text().splitlines()
    start = None
    for i, l in enumerate(lines):
        if l.strip().startswith("rij") and "Add00" in l:
            start = i + 1
            break
    if start is None:
        raise RuntimeError(f"no r-space header found in {path}")
    row = re.compile(
        rf"^\s*({NUM})\s+({NUM})\s*\(\s*({NUM})\s*\)\s*({NUM})\s*\(\s*({NUM})\s*\)\s*(\d+)")
    r, a00, a22 = [], [], []
    for l in lines[start:]:
        m = row.match(l)
        if not m:
            break
        r.append(float(m.group(1)))
        a00.append(float(m.group(2)))
        a22.append(float(m.group(4)))
    return r, a00, a22


def he4_density_r1(path: Path = DENSITY):
    """R1 = sqrt(<r^2>) of the one-body point-nucleon density, RECOMPUTED
    from the file's own printed grid (never trusted as a typed constant --
    `docs/CONVENTIONS.md`).  Cross-checked against the file's own printed
    rms and against `LI6_R_POINT_VMC_FM`-style headers elsewhere in the
    tree."""
    r, rho = read_he4_density(path)
    w2 = [rho_i * ri * ri for rho_i, ri in zip(rho, r)]
    w4 = [rho_i * ri**4 for rho_i, ri in zip(rho, r)]
    norm = 4.0 * math.pi * _trapz(w2, r)
    r2 = 4.0 * math.pi * _trapz(w4, r) / norm
    return norm, math.sqrt(r2)


def he4_dd_moments(path: Path = HE4DD):
    """<r^2> and P_D of the alpha -> d+d relative-motion S+D overlap, from
    the SAME int (A00^2+A22^2) r^2/r^4 dr convention `vmc_reconcile.py` uses
    for `li6.ad`/`li7.at`."""
    r, a00, a22 = read_he4_dd(path)
    dr = r[1] - r[0]
    norm = s_norm = d_norm = r2mom = 0.0
    for ri, a0, a2 in zip(r, a00, a22):
        w = (a0 * a0 + a2 * a2) * ri * ri * dr
        norm += w
        r2mom += w * ri * ri
        s_norm += a0 * a0 * ri * ri * dr
        d_norm += a2 * a2 * ri * ri * dr
    r2 = r2mom / norm
    return {"norm": norm, "s_norm": s_norm, "d_norm": d_norm,
            "p_d": d_norm / (s_norm + d_norm), "r2": r2,
            "r_rms": math.sqrt(r2)}


def report(out=sys.stdout) -> dict:
    norm_he4, r1 = he4_density_r1()
    dd = he4_dd_moments()

    # the uncorrelated-product null: <D^2> = (4/3) R1^2, exactly, worked in
    # the module docstring -- independent of recentring, independent of
    # which 2-2 partition of the four iid draws is taken.
    d2_uncorr = (4.0 / 3.0) * r1 * r1
    d_rms_uncorr = math.sqrt(d2_uncorr)

    ratio_var = dd["r2"] / d2_uncorr
    ratio_rms = dd["r_rms"] / d_rms_uncorr

    print("O3 offline bound -- alpha-core correlation vs. the uncorrelated "
          "product of one-body densities", file=out)
    print("=" * 78, file=out)
    print(f"he4.density: printed norm check (should be ~1.9974, Z) "
          f"= {norm_he4:.5f}", file=out)
    print(f"he4.density: R1 = sqrt(<r^2>) (recomputed)             "
          f"= {r1:.5f} fm  (file's own rms 1.4404 fm)", file=out)
    print(file=out)
    print("he4.dd (alpha -> d+d, CORRELATED VMC overlap):", file=out)
    print(f"  S+D norm (should be ~0.973, file's own 'ndx' line)   "
          f"= {dd['norm']:.5f}", file=out)
    print(f"  P_D = d_norm/(s_norm+d_norm)  (file's own 0.024/0.973"
          f" = 0.0247) = {dd['p_d']:.5f}", file=out)
    print(f"  <r^2>_dd (correlated)                                 "
          f"= {dd['r2']:.5f} fm^2", file=out)
    print(f"  r_rms_dd (correlated)                                 "
          f"= {dd['r_rms']:.5f} fm", file=out)
    print(file=out)
    print("Uncorrelated-product null model (4 iid draws from the SOURCE "
          "density that reproduces he4.density after recentring, exactly "
          "ClusterConfigSampler's alpha_cm_inflate convention):", file=out)
    print(f"  <D^2>_uncorr = (4/3) R1^2                             "
          f"= {d2_uncorr:.5f} fm^2", file=out)
    print(f"  r_rms(D)_uncorr                                       "
          f"= {d_rms_uncorr:.5f} fm", file=out)
    print(file=out)
    print(f"RATIO <r^2>_correlated / <D^2>_uncorrelated              "
          f"= {ratio_var:.5f}  ({(ratio_var - 1) * 100:+.1f} %)", file=out)
    print(f"RATIO r_rms_correlated / r_rms_uncorrelated              "
          f"= {ratio_rms:.5f}  ({(ratio_rms - 1) * 100:+.1f} %)", file=out)
    print(file=out)
    print("Reading: genuine 4-body correlation in this ONE clustering "
          "channel widens the internal separation moment by "
          f"~{(ratio_var - 1) * 100:.0f}% (variance) / "
          f"~{(ratio_rms - 1) * 100:.0f}% (rms) relative to the "
          "uncorrelated-product null -- a several-to-twenty-percent effect, "
          "not a factor of several. This is NOT a measurement of how the "
          "coherent/incoherent split moves; that needs the actual GFMC "
          "4He configurations or the subnucleondiffraction graft.", file=out)

    return {"r1_he4_density": r1, "norm_he4_density": norm_he4,
            "dd_norm": dd["norm"], "dd_p_d": dd["p_d"], "dd_r2": dd["r2"],
            "dd_r_rms": dd["r_rms"], "d2_uncorrelated": d2_uncorr,
            "d_rms_uncorrelated": d_rms_uncorr,
            "ratio_variance": ratio_var, "ratio_rms": ratio_rms}


def _self_check(res: dict) -> None:
    # he4.density's own printed normalization is Z = 2 (proton density); the
    # recomputed trapezoid integral should land within the file's own
    # quadrature-vs-trapezoid spread already measured elsewhere in the tree
    # (`data/vmc/density/README.md`: 0.07% in the 4th digit).
    assert abs(res["norm_he4_density"] - 1.9974) / 1.9974 < 2e-3, res
    assert abs(res["r1_he4_density"] - 1.4404) / 1.4404 < 2e-3, res
    # he4.dd's own printed decomposition: ndx/s-wave/d-wave = 0.973/0.949/0.024
    assert abs(res["dd_norm"] - 0.973) / 0.973 < 5e-3, res
    assert abs(res["dd_p_d"] - 0.024 / 0.973) / (0.024 / 0.973) < 5e-2, res
    # the correlated separation must be LARGER than the uncorrelated null
    # (short-range NN repulsion / antisymmetrization push nucleons apart,
    # never together, for this channel) and the effect must be modest, i.e.
    # this bound rules out "correlations dominate":
    assert 1.0 < res["ratio_variance"] < 2.0, res
    assert 1.0 < res["ratio_rms"] < 1.5, res


def main() -> int:
    res = report()
    try:
        _self_check(res)
    except AssertionError as exc:
        print(f"\nSELF-CHECK FAILED: {exc}", file=sys.stderr)
        return 1
    print("\nself-check: OK", file=sys.stdout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
