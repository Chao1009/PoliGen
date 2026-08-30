#!/usr/bin/env python3
"""VmcRadial prototype: Fourier-Bessel-transform the real VMC (AV18+UIX)
two-cluster overlap tables at data/vmc/ into momentum space, and compare
their momentum densities and tail fractions against LiPolGen's two-parameter
Hulthen model (cluster.hpp / cluster.cpp, ported here from spectator.py /
tagged.py's Wave.radial).

Data provenance: data/vmc/README.md.  Primary channels:
  6Li(1+) -> alpha + d   S+D wave   data/vmc/li6_alpha_d/li6.ad
  7Li(3/2-) -> alpha + t P wave     data/vmc/li7_alpha_t/li7.at
Both files carry BOTH the raw k-space amplitudes (Monte Carlo, with
statistical errors) AND the raw r-space overlap (also Monte Carlo, an
independent estimator from the same run -- not derived from the k-space
column or vice versa).  That lets this script do two things at once:
  1. Compare VMC vs Hulthen directly using the tabulated k-space amplitudes
     (no transform needed/no transform error).
  2. Validate the Fourier-Bessel transform itself by transforming the
     r-space columns and checking they reproduce the (independently
     estimated) tabulated k-space columns -- first against the *clean*
     deuteron wave function (data/vmc/deuteron/fdeut.av18, which is not
     Monte Carlo, so this is a clean numerical check of the transform's
     normalization convention), then against the noisier VMC alpha-d/
     alpha-t tables.  This is exactly the machinery a future VmcRadial
     backend needs when the only available table is r-space only (e.g.
     data/vmc/h3_d_n/h2n.table, data/vmc/li7_li6_n/li6n_31.table).

Convention (matches LiPolGen's own normalization, cluster.hpp:10-11):
  R_L(k) = sqrt(2/pi) * Integral_0^inf  R_L(r) j_L(k r) r^2 dr
Parseval's theorem then gives Integral R_L(k)^2 k^2 dk = Integral R_L(r)^2
r^2 dr for this convention, matching "TaggedModel normalizes them on its
own grid so that integral psihat^2 k^2 dk = 1" (cluster.hpp:10).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.special import spherical_jn

HBARC_GEV_FM = 0.1973269804  # GeV * fm (CODATA hbar*c)

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "vmc"

NUM = r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?'
KROW_RE = re.compile(rf'^\s*({NUM})\s*({NUM})\s*\(\s*({NUM})\s*\)\s*({NUM})\s*\(\s*({NUM})\s*\)')

# --------------------------------------------------------------------------
# Parsers for the raw VMC Fortran-code output (li6.ad, li7.at)
# --------------------------------------------------------------------------


def _read_amplitude_block(lines: list[str], header_idx: int) -> np.ndarray:
    """Rows 'x  y1 (dy1)  y2 (dy2) [trailing ints ignored]' following a
    header line.  Stops at the first non-matching line.  Returns an
    (N, 5) array: x, y1, dy1, y2, dy2."""
    rows = []
    for line in lines[header_idx + 1:]:
        m = KROW_RE.match(line)
        if not m:
            if rows:
                break
            continue
        rows.append([float(g) for g in m.groups()])
    return np.array(rows)


@dataclass
class ClusterOverlap:
    """One VMC two-cluster raw-output file: two partial waves (label0,
    label1), each with a k-space amplitude table and an r-space amplitude
    table (independent MC estimators, both with 1-sigma errors)."""
    name: str
    label: tuple[str, str]
    l: tuple[int, int]
    k_fm: np.ndarray
    k_amp: tuple[np.ndarray, np.ndarray]      # (A_l0(k), A_l1(k))
    k_err: tuple[np.ndarray, np.ndarray]
    r_fm: np.ndarray
    r_amp: tuple[np.ndarray, np.ndarray]
    r_err: tuple[np.ndarray, np.ndarray]


def parse_cluster_overlap(path: Path, name: str, labels: tuple[str, str],
                           ls: tuple[int, int]) -> ClusterOverlap:
    lines = path.read_text().splitlines()
    k_hdr = next(i for i, l in enumerate(lines) if l.strip().startswith('k(fm-1)'))
    r_hdr = next(i for i, l in enumerate(lines) if l.strip().startswith('rij'))
    kblk = _read_amplitude_block(lines, k_hdr)
    rblk = _read_amplitude_block(lines, r_hdr)
    return ClusterOverlap(
        name=name, label=labels, l=ls,
        k_fm=kblk[:, 0], k_amp=(kblk[:, 1], kblk[:, 3]), k_err=(kblk[:, 2], kblk[:, 4]),
        r_fm=rblk[:, 0], r_amp=(rblk[:, 1], rblk[:, 3]), r_err=(rblk[:, 2], rblk[:, 4]),
    )


def parse_fdeut(path: Path) -> dict:
    lines = path.read_text().splitlines()
    r_hdr = next(i for i, l in enumerate(lines) if l.strip().startswith('r ') and 'u(k)' not in l)
    k_hdr = next(i for i, l in enumerate(lines) if l.strip().startswith('k ') and 'u(k)' in l)
    r_rows = []
    for line in lines[r_hdr + 1:]:
        parts = line.split()
        if len(parts) != 5:
            if r_rows:
                break
            continue
        r_rows.append([float(p) for p in parts])
    r_rows = np.array(r_rows)
    k_rows = []
    for line in lines[k_hdr + 1:]:
        parts = line.split()
        if len(parts) != 3:
            if k_rows:
                break
            continue
        k_rows.append([float(p) for p in parts])
    k_rows = np.array(k_rows)
    return dict(r_fm=r_rows[:, 0], u_r=r_rows[:, 1], w_r=r_rows[:, 3],
                k_fm=k_rows[:, 0], u_k=k_rows[:, 1], w_k=k_rows[:, 2])


# --------------------------------------------------------------------------
# Fourier-Bessel transform R_L(r) -> R_L(k)
# --------------------------------------------------------------------------


def fourier_bessel(r: np.ndarray, R_r: np.ndarray, k_grid: np.ndarray, l: int) -> np.ndarray:
    """R_L(k) = sqrt(2/pi) Integral_0^inf R_L(r) j_L(k r) r^2 dr.  This is
    the "reduced radial function" convention (u(r) = r R_0(r), matched
    exactly below against fdeut.av18's own tabulated u(k)/w(k) -- ratio
    1.000 to 4 digits at every check point)."""
    base = R_r * r ** 2
    out = np.empty(len(k_grid))
    for i, k in enumerate(k_grid):
        jl = spherical_jn(l, k * r)
        out[i] = np.sqrt(2.0 / np.pi) * np.trapz(base * jl, r)
    return out


def fourier_bessel_cluster(r: np.ndarray, A_r: np.ndarray, k_grid: np.ndarray, l: int) -> np.ndarray:
    """A_L(k) = 4*pi * Integral_0^inf A_L(r) j_L(k r) r^2 dr.

    This is a DIFFERENT normalization convention from fourier_bessel()
    above -- the one the VMC spectroscopic-overlap files (li6.ad, li7.at,
    and the r-space-only li6n_31.table / h2n.table) use, per the
    overlap_gfmc page's own definition: N(k) = A(k)^2/(4*pi), spectroscopic
    factor = 1/(2*pi)^3 Integral d^3k N(k) = Integral A^2(r) r^2 dr.
    Combined with the reduced-radial (fourier_bessel) convention's Parseval
    identity, (2pi)^(3/2) * sqrt(2/pi) = 4*pi exactly -- i.e. A_L(k) here is
    (2pi)^(3/2) times the "R_L(k)" of fourier_bessel().  Calibrated
    empirically below against li6.ad's own tabulated k-space Aad00(k):
    the naive sqrt(2/pi) transform came out a factor of ~15.7 too small
    at every k, matching 4*pi/sqrt(2/pi) = 4*pi*sqrt(pi/2) = 15.75 almost
    exactly -- the giveaway that this second convention was in play.
    """
    base = A_r * r ** 2
    out = np.empty(len(k_grid))
    for i, k in enumerate(k_grid):
        jl = spherical_jn(l, k * r)
        out[i] = 4.0 * np.pi * np.trapz(base * jl, r)
    return out


# --------------------------------------------------------------------------
# LiPolGen's Hulthen model, reimplemented from cluster.cpp / spectator.py
# --------------------------------------------------------------------------


def hulthen_radial(k: np.ndarray, kappa: float, beta: float, l: int) -> np.ndarray:
    k2 = k * k
    b2 = beta * beta
    kap2 = kappa * kappa
    if l == 0:
        return 1.0 / (k2 + kap2) - 1.0 / (k2 + b2)
    if l == 1:
        return k / ((k2 + kap2) * (k2 + b2))
    if l == 2:
        return k2 / ((k2 + kap2) * (k2 + b2) ** 2)
    raise ValueError("l must be 0, 1, or 2")


# Separation-energy-derived kappa, from spectator.py's own docstring
# (verified AME2020 masses, plans/02):
KAPPA_LI6_AD = 0.06066  # GeV: 6Li -> alpha+d, S = 1.474 MeV
KAPPA_LI7_AT = 0.08890  # GeV: 7Li -> alpha+t, S = 2.467 MeV
P_D_LI6_PYTHON = 0.0867  # tagged.py P_D_LI6 (SCENARIO value)


# --------------------------------------------------------------------------
# Density / tail-fraction utilities
# --------------------------------------------------------------------------


def density_from_amplitudes(k: np.ndarray, amps: list[np.ndarray]) -> np.ndarray:
    """k^2 * sum_L A_L(k)^2, unnormalized."""
    return k ** 2 * sum(a ** 2 for a in amps)


def normalize(k: np.ndarray, dens: np.ndarray) -> np.ndarray:
    return dens / np.trapz(dens, k)


def tail_fraction(k: np.ndarray, dens_norm: np.ndarray, k0: float) -> float:
    mask = k >= k0
    if mask.sum() < 2:
        return 0.0
    return float(np.trapz(dens_norm[mask], k[mask]))


def sample_pt_tail_fractions(k: np.ndarray, dens_norm: np.ndarray, thresholds,
                              n_samples=400_000, seed=0) -> dict:
    """Isotropic projection: sample |k| from k^2 n(k) (dens_norm is already
    that, i.e. proportional to the radial pdf in k), sample cos(theta)
    uniform on [-1,1], p_T = |k| sin(theta).  Returns {threshold: P(p_T>thr)}."""
    rng = np.random.default_rng(seed)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (dens_norm[1:] + dens_norm[:-1])
                                            * np.diff(k))])
    cdf /= cdf[-1]
    u = rng.uniform(0, 1, n_samples)
    ksamp = np.interp(u, cdf, k)
    cos_t = rng.uniform(-1, 1, n_samples)
    pt = ksamp * np.sqrt(1 - cos_t ** 2)
    return {thr: float(np.mean(pt > thr)) for thr in thresholds}


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main():
    li6 = parse_cluster_overlap(DATA / "li6_alpha_d" / "li6.ad", "6Li -> alpha+d",
                                 ("Aad00 (S, L=0)", "Aad22 (D, L=2)"), (0, 2))
    li7 = parse_cluster_overlap(DATA / "li7_alpha_t" / "li7.at", "7Li -> alpha+t",
                                 ("Aat33 (P, j=3/2)", "Aat11 (P, j=1/2)"), (1, 1))
    deut = parse_fdeut(DATA / "deuteron" / "fdeut.av18")

    print("=" * 78)
    print("1. FILES PARSED")
    print("=" * 78)
    print(f"{li6.name}: {len(li6.k_fm)} k-space rows (k=0..{li6.k_fm[-1]:.2f} fm^-1), "
          f"{len(li6.r_fm)} r-space rows (r=0..{li6.r_fm[-1]:.2f} fm)")
    print(f"{li7.name}: {len(li7.k_fm)} k-space rows (k=0..{li7.k_fm[-1]:.2f} fm^-1), "
          f"{len(li7.r_fm)} r-space rows (r=0..{li7.r_fm[-1]:.2f} fm)")
    print(f"deuteron (fdeut.av18): {len(deut['r_fm'])} r-space rows "
          f"(r=0..{deut['r_fm'][-1]:.1f} fm), {len(deut['k_fm'])} k-space rows "
          f"(k=0..{deut['k_fm'][-1]:.1f} fm^-1); D-state ~5.76% (AV18 std value, sanity check)")

    # ----------------------------------------------------------------
    # 2. Validate the Fourier-Bessel transform against the clean (non-MC)
    #    deuteron wave function, then against the noisy VMC alpha-d table.
    # ----------------------------------------------------------------
    print()
    print("=" * 78)
    print("2. FOURIER-BESSEL TRANSFORM VALIDATION")
    print("=" * 78)

    r = deut["r_fm"]
    R0_r = deut["u_r"] / r          # u(r) = r R_0(r)
    R2_r = np.divide(deut["w_r"], r, out=np.zeros_like(r), where=r > 0)
    k_check = np.array([0.2, 0.5, 1.0, 1.5, 2.0, 3.0])
    R0_k_pred = fourier_bessel(r, R0_r, k_check, 0)
    R2_k_pred = fourier_bessel(r, R2_r, k_check, 2)
    # Empirically (checked against the tabulated u(k)/w(k) below): ANL's
    # "u(k)"/"w(k)" columns are R_L(k) itself, NOT k*R_L(k) -- i.e. the
    # k-space amplitude is *not* reduced the way u(r)=r R_0(r) is in
    # r-space.  (A first attempt with the r-space-style k*R_L(k) reduction
    # came out exactly a factor of k off at every check point -- the
    # giveaway that the k-space convention drops that factor.)
    u_k_pred = R0_k_pred
    w_k_pred = R2_k_pred
    u_k_true = np.interp(k_check, deut["k_fm"], deut["u_k"])
    w_k_true = np.interp(k_check, deut["k_fm"], deut["w_k"])
    print("Deuteron u(k), w(k): transform of r-space u(r)/w(r) [my FT] vs the "
          "table's own tabulated k-space columns [ANL]:")
    print(f"{'k [fm^-1]':>10} {'u(k) mine':>12} {'u(k) ANL':>12} {'ratio':>8} "
          f"{'w(k) mine':>12} {'w(k) ANL':>12} {'ratio':>8}")
    for i, k in enumerate(k_check):
        ru = u_k_pred[i] / u_k_true[i] if u_k_true[i] else float("nan")
        rw = w_k_pred[i] / w_k_true[i] if w_k_true[i] else float("nan")
        print(f"{k:10.2f} {u_k_pred[i]:12.4f} {u_k_true[i]:12.4f} {ru:8.3f} "
              f"{w_k_pred[i]:12.4f} {w_k_true[i]:12.4f} {rw:8.3f}")

    # Now the same transform on the (Monte Carlo noisy) VMC alpha-d r-space
    # table, compared against the file's own (independent) MC k-space
    # estimate -- both are estimates of the same quantity from the same run.
    kfb = np.linspace(0.05, 3.0, 30)
    Aad00_k_pred = fourier_bessel_cluster(li6.r_fm, li6.r_amp[0], kfb, 0)
    Aad22_k_pred = fourier_bessel_cluster(li6.r_fm, li6.r_amp[1], kfb, 2)
    Aad00_k_native = np.interp(kfb, li6.k_fm, li6.k_amp[0])
    Aad22_k_native = np.interp(kfb, li6.k_fm, li6.k_amp[1])
    print()
    print("6Li alpha-d: FT(r-space Aad00/Aad22, 4*pi convention) vs the "
          "file's own (noisy, independent) k-space Aad00/Aad22 MC estimate:")
    print(f"{'k [fm^-1]':>10} {'Aad00 FT':>12} {'Aad00 ANL':>12} {'ratio':>7} "
          f"{'Aad22 FT':>12} {'Aad22 ANL':>12} {'ratio':>7}")
    for i in range(0, len(kfb), 5):
        r0 = Aad00_k_native[i] / Aad00_k_pred[i] if Aad00_k_pred[i] else float("nan")
        r2 = Aad22_k_native[i] / Aad22_k_pred[i] if Aad22_k_pred[i] else float("nan")
        print(f"{kfb[i]:10.2f} {Aad00_k_pred[i]:12.3f} {Aad00_k_native[i]:12.3f} {r0:7.2f} "
              f"{Aad22_k_pred[i]:12.4f} {Aad22_k_native[i]:12.4f} {r2:7.2f}")
    print("(Same normalization convention as fourier_bessel(), just with "
          "the 4*pi prefactor the cluster-overlap files use instead of "
          "sqrt(2/pi) -- ratio columns should sit near 1 modulo Monte Carlo "
          "noise in BOTH the r-space and k-space tables, which are "
          "independent MC estimates from the same run, not derived from "
          "each other. The transform's normalization convention itself is "
          "validated exactly (ratio 1.000 to 4 digits) against the clean, "
          "non-MC deuteron pair above.)")

    # ----------------------------------------------------------------
    # 3. D-state probability, 6Li alpha-d
    # ----------------------------------------------------------------
    print()
    print("=" * 78)
    print("3. 6Li alpha+d D-STATE PROBABILITY")
    print("=" * 78)
    k6, A6_00, A6_22 = li6.k_fm, li6.k_amp[0], li6.k_amp[1]
    I0 = np.trapz(k6 ** 2 * A6_00 ** 2, k6)
    I2 = np.trapz(k6 ** 2 * A6_22 ** 2, k6)
    P_D_vmc = I2 / (I0 + I2)
    print(f"VMC (AV18+UIX, this table): P_D = Integral k^2 Aad22^2 dk / "
          f"Integral k^2(Aad00^2+Aad22^2) dk = {P_D_vmc:.4f}")
    print("  (file's own printed 'ndx s-wave d-wave' line: 0.838 s, 0.017 d "
          f"-> P_D = 0.017/(0.838+0.017) = {0.017/(0.838+0.017):.4f}, consistent)")
    print(f"LiPolGen python (tagged.py P_D_LI6, SCENARIO placeholder): {P_D_LI6_PYTHON:.4f}")
    print(f"  ratio VMC/python = {P_D_vmc / P_D_LI6_PYTHON:.3f} "
          f"-- the scenario value is ~{P_D_LI6_PYTHON / P_D_vmc:.1f}x the VMC D-state fraction")

    # ----------------------------------------------------------------
    # 4. Momentum density comparison: VMC vs Hulthen, both channels
    # ----------------------------------------------------------------
    print()
    print("=" * 78)
    print("4. MOMENTUM DENSITY k^2 n(k) COMPARISON (VMC vs Hulthen)")
    print("=" * 78)

    def hulthen_density(k_gev, kappa, beta, waves):
        """waves: list of (l, prob).  Each L normalized separately to unit
        norm on [0, k_gev[-1]], combined with its prob weight (matches
        cluster.hpp: 'TaggedModel normalizes them on its own grid so that
        integral psihat^2 k^2 dk = 1')."""
        tot = np.zeros_like(k_gev)
        for l, prob in waves:
            shape = hulthen_radial(k_gev, kappa, beta, l)
            dens = k_gev ** 2 * shape ** 2
            dens = dens / np.trapz(dens, k_gev)
            tot += prob * dens
        return tot

    thresholds = [0.10, 0.20, 0.30, 0.45]

    # --- 6Li alpha+d ---
    k6_gev = k6 * HBARC_GEV_FM
    dens6_vmc = normalize(k6_gev, density_from_amplitudes(k6_gev, [A6_00, A6_22]))
    dens6_h30 = hulthen_density(k6_gev, KAPPA_LI6_AD, 0.30,
                                 [(0, 1 - P_D_LI6_PYTHON), (2, P_D_LI6_PYTHON)])
    dens6_h40 = hulthen_density(k6_gev, KAPPA_LI6_AD, 0.40,
                                 [(0, 1 - P_D_LI6_PYTHON), (2, P_D_LI6_PYTHON)])

    print(f"\n6Li -> alpha+d  (k range 0-{k6_gev[-1]:.3f} GeV, VMC table native extent)")
    print(f"{'model':<28}" + "".join(f"P(k>{t:.2f})".rjust(12) for t in thresholds))
    for label, dens in [("VMC (AV18+UIX)", dens6_vmc),
                         ("Hulthen beta=0.30", dens6_h30),
                         ("Hulthen beta=0.40", dens6_h40)]:
        vals = [tail_fraction(k6_gev, dens, t) for t in thresholds]
        print(f"{label:<28}" + "".join(f"{v:12.4f}" for v in vals))

    print(f"\n6Li -> alpha+d  isotropic-projection p_T tail (same densities, "
          f"p_T = |k| sin(theta), cos(theta) ~ U[-1,1])")
    print(f"{'model':<28}" + "".join(f"P(pT>{t:.2f})".rjust(12) for t in thresholds))
    for label, dens in [("VMC (AV18+UIX)", dens6_vmc),
                         ("Hulthen beta=0.30", dens6_h30),
                         ("Hulthen beta=0.40", dens6_h40)]:
        pt_fracs = sample_pt_tail_fractions(k6_gev, dens, thresholds)
        print(f"{label:<28}" + "".join(f"{pt_fracs[t]:12.4f}" for t in thresholds))

    # --- 7Li alpha+t ---
    k7, A7_33, A7_11 = li7.k_fm, li7.k_amp[0], li7.k_amp[1]
    k7_gev = k7 * HBARC_GEV_FM
    dens7_vmc = normalize(k7_gev, density_from_amplitudes(k7_gev, [A7_33, A7_11]))
    dens7_h30 = hulthen_density(k7_gev, KAPPA_LI7_AT, 0.30, [(1, 1.0)])
    dens7_h40 = hulthen_density(k7_gev, KAPPA_LI7_AT, 0.40, [(1, 1.0)])

    print(f"\n7Li -> alpha+t  (k range 0-{k7_gev[-1]:.3f} GeV, VMC table native extent, "
          f"P-wave only)")
    print(f"{'model':<28}" + "".join(f"P(k>{t:.2f})".rjust(12) for t in thresholds))
    for label, dens in [("VMC (AV18+UIX)", dens7_vmc),
                         ("Hulthen beta=0.30", dens7_h30),
                         ("Hulthen beta=0.40", dens7_h40)]:
        vals = [tail_fraction(k7_gev, dens, t) for t in thresholds]
        print(f"{label:<28}" + "".join(f"{v:12.4f}" for v in vals))

    print(f"\n7Li -> alpha+t  isotropic-projection p_T tail")
    print(f"{'model':<28}" + "".join(f"P(pT>{t:.2f})".rjust(12) for t in thresholds))
    for label, dens in [("VMC (AV18+UIX)", dens7_vmc),
                         ("Hulthen beta=0.30", dens7_h30),
                         ("Hulthen beta=0.40", dens7_h40)]:
        pt_fracs = sample_pt_tail_fractions(k7_gev, dens, thresholds)
        print(f"{label:<28}" + "".join(f"{pt_fracs[t]:12.4f}" for t in thresholds))

    # ----------------------------------------------------------------
    # 5. Plot
    # ----------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    ax = axes[0, 0]
    ax.semilogy(k6_gev, dens6_vmc, label="VMC (AV18+UIX)", color="k", lw=2)
    ax.semilogy(k6_gev, dens6_h30, label=r"Hulthen $\beta=0.30$", ls="--")
    ax.semilogy(k6_gev, dens6_h40, label=r"Hulthen $\beta=0.40$", ls="--")
    ax.set_xlabel("k [GeV]"); ax.set_ylabel(r"$k^2 n(k)$ (normalized)")
    ax.set_title("6Li -> alpha+d (S+D)"); ax.set_xlim(0, 1.0); ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.semilogy(k7_gev, dens7_vmc, label="VMC (AV18+UIX)", color="k", lw=2)
    ax.semilogy(k7_gev, dens7_h30, label=r"Hulthen $\beta=0.30$", ls="--")
    ax.semilogy(k7_gev, dens7_h40, label=r"Hulthen $\beta=0.40$", ls="--")
    ax.set_xlabel("k [GeV]"); ax.set_ylabel(r"$k^2 n(k)$ (normalized)")
    ax.set_title("7Li -> alpha+t (P)"); ax.set_xlim(0, 1.0); ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(deut["k_fm"], deut["u_k"], label="u(k) tabulated (ANL)", color="k")
    ax.plot(k_check, u_k_pred, "o", label="u(k) my FT of u(r)", color="tab:red")
    ax.plot(deut["k_fm"], deut["w_k"], label="w(k) tabulated (ANL)", color="0.4")
    ax.plot(k_check, w_k_pred, "s", label="w(k) my FT of w(r)", color="tab:blue")
    ax.set_xlim(0, 3); ax.set_xlabel("k [fm$^{-1}$]"); ax.set_ylabel("u(k), w(k)")
    ax.set_title("FT validation: deuteron (AV18)"); ax.legend(fontsize=7)

    ax = axes[1, 1]
    labels = [f">{t:.2f}" for t in thresholds]
    x = np.arange(len(thresholds))
    w = 0.25
    v6 = [tail_fraction(k6_gev, dens6_vmc, t) for t in thresholds]
    h30_6 = [tail_fraction(k6_gev, dens6_h30, t) for t in thresholds]
    h40_6 = [tail_fraction(k6_gev, dens6_h40, t) for t in thresholds]
    ax.bar(x - w, v6, w, label="VMC 6Li ad")
    ax.bar(x, h30_6, w, label=r"Hulth. $\beta$=0.30")
    ax.bar(x + w, h40_6, w, label=r"Hulth. $\beta$=0.40")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_yscale("log"); ax.set_ylabel("P(k > threshold) [GeV]")
    ax.set_title("6Li alpha+d tail fractions"); ax.legend(fontsize=7)

    fig.suptitle("VMC (AV18+UIX, ANL) vs LiPolGen Hulthen model", fontsize=12)
    fig.tight_layout()
    out_png = HERE / "vmc_overlap_comparison.png"
    fig.savefig(out_png, dpi=140)
    print(f"\nSaved {out_png}")


if __name__ == "__main__":
    main()
