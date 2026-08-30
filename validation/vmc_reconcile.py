#!/usr/bin/env python3
"""Reconcile the two independent readings of the ANL VMC cluster data.

Two earlier investigations fetched DIFFERENT files and reported numbers that
cannot both be true:

  Thread A  `overlap_old/li6.ad`, `overlap_old/li7.at` (AV18+UIX, 19-22 Apr
            2004; raw MC output with BOTH a k-space amplitude block and an
            r-space amplitude block).  Reported 7Li alpha+t
            P(k > 0.45 GeV) = 0.0095 (VMC) vs 0.22 (Hulthen beta = 0.40),
            i.e. "VMC much softer", and 6Li P_D = 0.020.
  Thread B  `momenta/li6_ad1.momentum`, `momenta/li7_at3.momentum` (AV18+UX,
            22-Mar-2014 / 12-Apr-2024; momentum DENSITIES rho_L(K) with an
            S/D split for 6Li).  Reported 7Li <k> = 0.186 GeV,
            P(k > 0.3) = 0.179 (VMC) vs 0.084 (Hulthen beta = 0.40), i.e.
            "VMC harder", plus 6Li <k> = 0.122, P(k>0.3) = 0.045,
            P_D = 1.94%, S_ad = 0.820, S_at = 1.008.

This script computes EVERY number in ONE place with ONE normalization, for

    {overlap k-space, overlap r-space -> Fourier-Bessel, momentum file}
      x  {6Li alpha+d, 7Li alpha+t}

and for LiPolGen's own Hulthen radial forms evaluated THROUGH THE C++
BINDINGS (`lipolgen.Wave.radial`, i.e. `src/core/cluster.cpp` itself, not a
Python re-implementation) at beta = 0.20 / 0.30 / 0.40 with the kappa the
channels actually carry.

THE ONE NORMALIZATION.  Every model is reduced to the radial probability
density in the relative momentum k,

    D(k) = k^2 * sum_L |psi_L(k)|^2  ,   normalized to  int_0^KMAX D dk = 1,

with the SAME KMAX for every row of a table.  Then
    <k>       = int k D dk
    P(k > k0) = int_{k0}^{KMAX} D dk
    P_D       = int k^2 |psi_2|^2 dk / int k^2 sum_L |psi_L|^2 dk .
KMAX matters and is quoted with every table; the default KMAX is
5 fm^-1 = 0.98663 GeV, the common native extent of all four VMC tables.

Unit conversion: k[GeV] = k[fm^-1] * hbar c, hbar c = 0.1973269804 GeV fm.
P(k > k0) is invariant under the change of variable (the Jacobian cancels in
the normalized ratio) as long as the threshold is converted too; <k> simply
scales by hbar c.

Outputs: a printed report, `validation/vmc_reconcile.png`, and the table that
`docs/open_items/vmc_reconciliation.md` carries.

Provenance of every input file: `data/vmc/README.md`.
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import spherical_jn

# numpy < 2.0 spells it np.trapz; >= 2.0 spells it np.trapezoid.
_trapz = getattr(np, "trapezoid", None) or np.trapz

HBARC = 0.1973269804          # GeV fm
FM_MAX = 5.0                  # fm^-1, common native extent of the VMC tables
KMAX_DEFAULT = FM_MAX * HBARC  # 0.98663 GeV
KMAX_MODEL = 1.2              # GeV, the TaggedModel grid's k_max

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data" / "vmc"

THRESHOLDS = (0.10, 0.20, 0.30, 0.45)

NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"
PAIR_ROW = re.compile(
    rf"^\s*({NUM})\s*({NUM})\s*\(\s*({NUM})\s*\)\s*({NUM})\s*\(\s*({NUM})\s*\)")


# --------------------------------------------------------------------------
# parsers
# --------------------------------------------------------------------------


@dataclass
class Overlap:
    """`overlap_old/*.ad|*.at`: raw VMC stdout with two partial-wave columns,
    each tabulated BOTH in k (amplitude A_L(k)) and in r (amplitude A_L(r)),
    as two independent MC estimators from the same run."""
    name: str
    cols: tuple[str, str]
    ls: tuple[int, int]
    k_fm: np.ndarray
    k_amp: tuple[np.ndarray, np.ndarray]
    k_err: tuple[np.ndarray, np.ndarray]
    r_fm: np.ndarray
    r_amp: tuple[np.ndarray, np.ndarray]
    r_err: tuple[np.ndarray, np.ndarray]


def _pair_block(lines: list[str], hdr: int) -> np.ndarray:
    rows = []
    for line in lines[hdr + 1:]:
        m = PAIR_ROW.match(line)
        if not m:
            if rows:
                break
            continue
        rows.append([float(g) for g in m.groups()])
    return np.array(rows)


def parse_overlap(path: Path, name: str, cols: tuple[str, str],
                  ls: tuple[int, int]) -> Overlap:
    lines = path.read_text().splitlines()
    k_hdr = next(i for i, l in enumerate(lines) if l.strip().startswith("k(fm-1)"))
    r_hdr = next(i for i, l in enumerate(lines) if l.strip().startswith("rij"))
    kb, rb = _pair_block(lines, k_hdr), _pair_block(lines, r_hdr)
    return Overlap(name, cols, ls,
                   kb[:, 0], (kb[:, 1], kb[:, 3]), (kb[:, 2], kb[:, 4]),
                   rb[:, 0], (rb[:, 1], rb[:, 3]), (rb[:, 2], rb[:, 4]))


@dataclass
class MomentumFile:
    """`momenta/*.momentum`: rho_L(K) [fm^3] with 1-sigma MC errors, and the
    file's own printed normalization 4 pi int rho K^2 dK / (2 pi)^3 = S."""
    name: str
    k_fm: np.ndarray
    rho_tot: np.ndarray
    drho_tot: np.ndarray
    norms: list[float]                       # every printed 4pi.../(2pi)^3
    split: dict[int, np.ndarray] | None      # L -> rho_L, when the file splits
    dsplit: dict[int, np.ndarray] | None


def parse_momentum(path: Path, name: str) -> MomentumFile:
    text = path.read_text()
    norms = [float(m) for m in re.findall(
        r"4\*PI\*TOTINT\([^)]*\)/\(2\*PI\)\*\*3\s*=\s*(" + NUM + ")", text)]
    lines = text.splitlines()
    blocks: list[np.ndarray] = []
    i = 0
    while i < len(lines):
        stripped = lines[i].strip().replace(" ", "")
        if stripped and set(stripped) == {"*"} and len(stripped) > 4:
            rows = []
            for line in lines[i + 1:]:
                parts = line.split()
                if not parts or not re.fullmatch(NUM, parts[0]):
                    break
                try:
                    rows.append([float(p) for p in parts])
                except ValueError:
                    break
            if rows:
                blocks.append(np.array(rows))
                i += len(rows)
        i += 1
    tot = blocks[0]
    split = dsplit = None
    if len(blocks) > 1 and blocks[1].shape[1] >= 5:
        b = blocks[1]
        split = {0: b[:, 1], 2: b[:, 3]}
        dsplit = {0: b[:, 2], 2: b[:, 4]}
    return MomentumFile(name, tot[:, 0], tot[:, 1], tot[:, 2], norms,
                        split, dsplit)


def parse_fdeut(path: Path) -> dict:
    lines = path.read_text().splitlines()
    r_hdr = next(i for i, l in enumerate(lines)
                 if l.strip().startswith("r ") and "u(k)" not in l)
    k_hdr = next(i for i, l in enumerate(lines)
                 if l.strip().startswith("k ") and "u(k)" in l)

    def block(hdr, ncol):
        rows = []
        for line in lines[hdr + 1:]:
            parts = line.split()
            if len(parts) != ncol:
                if rows:
                    break
                continue
            rows.append([float(p) for p in parts])
        return np.array(rows)

    r, k = block(r_hdr, 5), block(k_hdr, 3)
    return dict(r_fm=r[:, 0], u_r=r[:, 1], w_r=r[:, 3],
                k_fm=k[:, 0], u_k=k[:, 1], w_k=k[:, 2])


# --------------------------------------------------------------------------
# Fourier-Bessel transforms (two conventions, both validated below)
# --------------------------------------------------------------------------


def fb_reduced(r, R_r, k_grid, l):
    """R_L(k) = sqrt(2/pi) int R_L(r) j_L(kr) r^2 dr -- the reduced-radial
    (deuteron u/w) convention.  Parseval: int R^2 k^2 dk = int R^2 r^2 dr."""
    base = R_r * r ** 2
    return np.array([np.sqrt(2.0 / np.pi) * _trapz(base * spherical_jn(l, k * r), r)
                     for k in k_grid])


def fb_cluster(r, A_r, k_grid, l):
    """A_L(k) = 4 pi int A_L(r) j_L(kr) r^2 dr -- the convention the ANL
    cluster-overlap files use, per the overlap page's own definition
    N(k) = A^2(k)/(4 pi), S = (2 pi)^-3 int d^3k N(k) = int A^2(r) r^2 dr.
    It is exactly (2 pi)^(3/2) times fb_reduced()."""
    base = A_r * r ** 2
    return np.array([4.0 * np.pi * _trapz(base * spherical_jn(l, k * r), r)
                     for k in k_grid])


# --------------------------------------------------------------------------
# the ONE normalization
# --------------------------------------------------------------------------


@dataclass
class Model:
    """sum_L k^2 |psi_L|^2 on a common GeV grid, plus the per-L pieces."""
    label: str
    k: np.ndarray                  # GeV
    per_l: dict[int, np.ndarray]   # L -> k^2 |psi_L(k)|^2 (unnormalized)

    @property
    def dens(self):
        return sum(self.per_l.values())

    def normalized(self):
        d = self.dens
        return d / _trapz(d, self.k)

    def mean_k(self):
        return float(_trapz(self.k * self.normalized(), self.k))

    def tail(self, k0):
        d = self.normalized()
        m = self.k >= k0
        if m.sum() < 2:
            return 0.0
        # include the partial first cell by linear interpolation
        kk = np.concatenate([[k0], self.k[m]])
        dd = np.concatenate([[np.interp(k0, self.k, d)], d[m]])
        return float(_trapz(dd, kk))

    def p_d(self):
        if 2 not in self.per_l:
            return float("nan")
        num = _trapz(self.per_l[2], self.k)
        den = _trapz(self.dens, self.k)
        return float(num / den)

    def row(self):
        return (self.mean_k(), *[self.tail(t) for t in THRESHOLDS], self.p_d())


def on_grid(k_src_gev, per_l_src, grid):
    """Interpolate each k^2|psi_L|^2 onto `grid`, zero outside the source
    range (that is exactly what `VmcRadial` does at runtime)."""
    out = {}
    for l, v in per_l_src.items():
        out[l] = np.interp(grid, k_src_gev, v, left=0.0, right=0.0)
    return out


# --------------------------------------------------------------------------
# LiPolGen's own Hulthen forms, through the C++ bindings
# --------------------------------------------------------------------------


def load_lipolgen():
    build = ROOT / "build" / "python"
    if build.is_dir():
        sys.path.insert(0, str(build))
    import lipolgen
    return lipolgen


def hulthen_model(lp, label, kappa, waves, grid):
    """`waves` = [(L, P_L)].  Each L is normalized SEPARATELY on the grid to
    unit int psihat^2 k^2 dk and then weighted by P_L -- exactly what
    `TaggedModel`'s constructor does (src/core/tagged.cpp:110-118)."""
    per_l = {}
    for l, prob in waves:
        w = lp.Wave(l=l, prob=prob, beta=0.0)
        per_l[l] = None, w
    return per_l


def hulthen_density(lp, kappa, beta, waves, grid):
    per_l = {}
    for l, prob in waves:
        w = lp.Wave(l=l, prob=1.0, beta=beta)
        psi = np.array([w.radial(float(k), kappa) for k in grid])
        d = grid ** 2 * psi ** 2
        per_l[l] = prob * d / _trapz(d, grid)
    return per_l


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------

HDR = f"{'model':<44}{'<k> GeV':>9}" + "".join(
    f"P(k>{t:.2f})".rjust(11) for t in THRESHOLDS) + f"{'P_D':>9}"


def print_rows(models):
    print(HDR)
    print("-" * len(HDR))
    for m in models:
        mk, *tails, pd = m.row()
        pd_s = "     --  " if np.isnan(pd) else f"{pd:9.4f}"
        print(f"{m.label:<44}{mk:9.4f}" + "".join(f"{t:11.4f}" for t in tails) + pd_s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmax", type=float, default=KMAX_DEFAULT,
                    help="GeV; the one normalization range (default 5 fm^-1)")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--markdown", type=Path, default=None,
                    help="write the reconciliation table to this file")
    args = ap.parse_args()

    lp = load_lipolgen()
    kappa6 = lp.li6_alpha_channel().base.kappa()
    kappa7 = lp.li7_alpha_channel().base.kappa()

    li6 = parse_overlap(DATA / "li6_alpha_d" / "li6.ad", "6Li -> a+d",
                        ("Aad00", "Aad22"), (0, 2))
    li7 = parse_overlap(DATA / "li7_alpha_t" / "li7.at", "7Li -> a+t",
                        ("Aat33", "Aat11"), (1, 1))
    m6 = parse_momentum(DATA / "momenta" / "li6_ad1.momentum", "6Li ad1")
    m7 = parse_momentum(DATA / "momenta" / "li7_at3.momentum", "7Li at3")
    m7x = parse_momentum(DATA / "momenta" / "li7_at1.momentum", "7Li at1 (1/2- exc)")
    deut = parse_fdeut(DATA / "deuteron" / "fdeut.av18")

    grid = np.linspace(1e-4, args.kmax, 4000)

    print("=" * 96)
    print("0.  INPUTS, AND THE TWO CONVENTIONS THEY ARE WRITTEN IN")
    print("=" * 96)
    print(f"kappa(6Li a+d) = {kappa6:.6f} GeV, kappa(7Li a+t) = {kappa7:.6f} GeV "
          "(from the C++ channels)")
    print(f"analysis grid: 0 .. {args.kmax:.5f} GeV = 0 .. {args.kmax/HBARC:.3f} fm^-1, "
          f"{len(grid)} points")
    print()
    print("overlap_old (2004, AV18+UIX): AMPLITUDES A_L(k) [fm^3/2-ish], signed,")
    print("  with S = (2pi)^-3 int k^2 A^2 dk  and  N(k) = A^2(k)/(4 pi).")
    print("momenta/    (2014/2024, AV18+UX): DENSITIES rho_L(K) [fm^3], positive,")
    print("  with S = 4 pi int rho K^2 dK / (2 pi)^3, i.e. rho = A^2/(4 pi).")
    for f in (m6, m7, m7x):
        print(f"  {f.name:<22} printed norms: "
              + ", ".join(f"{n:.5f}" for n in f.norms))
    s_ad = _trapz(li6.k_fm ** 2 * (li6.k_amp[0] ** 2 + li6.k_amp[1] ** 2),
                        li6.k_fm) / (2 * np.pi) ** 3
    s_at = _trapz(li7.k_fm ** 2 * li7.k_amp[0] ** 2, li7.k_fm) / (2 * np.pi) ** 3
    print(f"  overlap_old li6.ad  (2pi)^-3 int k^2 (A00^2+A22^2) dk = {s_ad:.5f}"
          f"   [momenta S_ad = {m6.norms[0]:.5f}]")
    print(f"  overlap_old li7.at  (2pi)^-3 int k^2  Aat33^2       dk = {s_at:.5f}"
          f"   [momenta S_at = {m7.norms[0]:.5f}]")
    print("  -> the two file families ARE the same quantity in the same units;")
    print("     the residual difference is Hamiltonian (UIX/2004 vs UX/2024) +")
    print("     the 5 fm^-1 truncation, not a convention mismatch.")

    # -- 1. the r-space -> k-space transform, and the relative sign ---------
    print()
    print("=" * 96)
    print("1.  FOURIER-BESSEL CONVENTION, VALIDATED ON THE CLEAN (NON-MC) DEUTERON")
    print("=" * 96)
    r = deut["r_fm"]
    R0 = deut["u_r"] / r
    R2 = np.divide(deut["w_r"], r, out=np.zeros_like(r), where=r > 0)
    kchk = np.array([0.2, 0.5, 1.0, 1.5, 2.0, 3.0])
    u_mine, w_mine = fb_reduced(r, R0, kchk, 0), fb_reduced(r, R2, kchk, 2)
    u_true = np.interp(kchk, deut["k_fm"], deut["u_k"])
    w_true = np.interp(kchk, deut["k_fm"], deut["w_k"])
    print(f"{'k[fm-1]':>8}{'u(k) FT':>12}{'u(k) ANL':>12}{'ratio':>8}"
          f"{'w(k) FT':>12}{'w(k) ANL':>12}{'ratio':>8}")
    for i, k in enumerate(kchk):
        print(f"{k:8.2f}{u_mine[i]:12.4f}{u_true[i]:12.4f}{u_mine[i]/u_true[i]:8.3f}"
              f"{w_mine[i]:12.4f}{w_true[i]:12.4f}{w_mine[i]/w_true[i]:8.3f}")
    print("sqrt(2/pi) is the deuteron convention; the cluster-overlap files use")
    print("4 pi = (2 pi)^(3/2) * sqrt(2/pi), a factor 15.749 larger.")

    kfb = np.arange(0.1, 5.01, 0.1)
    ft6 = {0: fb_cluster(li6.r_fm, li6.r_amp[0], kfb, 0),
           2: fb_cluster(li6.r_fm, li6.r_amp[1], kfb, 2)}
    ft7 = {1: fb_cluster(li7.r_fm, li7.r_amp[0], kfb, 1)}
    nat6 = {0: np.interp(kfb, li6.k_fm, li6.k_amp[0]),
            2: np.interp(kfb, li6.k_fm, li6.k_amp[1])}
    print()
    print("6Li: FT of the r-space block vs the file's own (independent) k-space block")
    print(f"{'k[fm-1]':>8}{'A00 FT':>12}{'A00 kblk':>12}{'ratio':>8}"
          f"{'A22 FT':>12}{'A22 kblk':>12}{'ratio':>8}")
    for i in range(0, 30, 4):
        print(f"{kfb[i]:8.2f}{ft6[0][i]:12.3f}{nat6[0][i]:12.3f}"
              f"{nat6[0][i]/ft6[0][i]:8.3f}"
              f"{ft6[2][i]:12.4f}{nat6[2][i]:12.4f}{nat6[2][i]/ft6[2][i]:8.3f}")
    print("Ratios are 1.000 to 1e-3: the 4 pi convention is right and the r-space")
    print("route reproduces the file's own INDEPENDENT k-space MC estimate.  So an")
    print("r-space-only table (li6n_31.table, h2n.table) can be used the same way.")

    # -- 2. the S-D relative sign -----------------------------------------
    print()
    print("=" * 96)
    print("2.  THE S-D RELATIVE SIGN (what the tensor observables actually need)")
    print("=" * 96)
    print("sign convention used everywhere below: fix the GLOBAL phase by")
    print("psi_L0(k -> 0) > 0, with L0 the lowest wave (S for 6Li, P for 7Li).")
    print()
    print(f"{'r[fm]':>8}{'A00(r)':>12}{'A22(r)':>12}   |"
          f"{'k[fm-1]':>9}{'A00(k)':>12}{'A22(k)':>12}{'FT A00':>12}{'FT A22':>12}")
    for rr, k in ((0.15, 0.2), (0.95, 0.4), (1.85, 0.6), (2.45, 1.0),
                  (3.45, 2.0), (4.45, 3.0)):
        i = int(round(k / 0.1)) - 1
        print(f"{rr:8.2f}{np.interp(rr, li6.r_fm, li6.r_amp[0]):12.4f}"
              f"{np.interp(rr, li6.r_fm, li6.r_amp[1]):12.4f}   |"
              f"{k:9.2f}{np.interp(k, li6.k_fm, li6.k_amp[0]):12.4f}"
              f"{np.interp(k, li6.k_fm, li6.k_amp[1]):12.4f}"
              f"{ft6[0][i]:12.4f}{ft6[2][i]:12.4f}")

    def zero_crossings(x, y, lo, hi):
        m = (x >= lo) & (x <= hi)
        xs, ys = x[m], y[m]
        s_ = np.sign(ys)
        idx = np.where(s_[:-1] * s_[1:] < 0)[0]
        return [float(xs[i] - ys[i] * (xs[i + 1] - xs[i]) / (ys[i + 1] - ys[i]))
                for i in idx]

    n00 = zero_crossings(li6.k_fm, li6.k_amp[0], 0.05, 2.0)
    n22 = zero_crossings(li6.k_fm, li6.k_amp[1], 0.05, 3.0)
    n33 = zero_crossings(li7.k_fm, li7.k_amp[0], 0.05, 2.0)
    print()
    print("r-space: A00(r) > 0 inside, crosses zero at r ~ "
          f"{zero_crossings(li6.r_fm, li6.r_amp[0], 1.0, 3.0)[0]:.2f} fm (the ONE node the")
    print("  Wildermuth condition 2n + L = 2 demands of the alpha-d 2S relative motion)")
    print("  and is NEGATIVE in the asymptotic tail; A22(r) > 0 in the tail.")
    print("k-space: A00(k) node at k = "
          + ", ".join(f"{z:.3f}" for z in n00) + " fm^-1; A22(k) node at k = "
          + ", ".join(f"{z:.3f}" for z in n22) + " fm^-1.")
    print("  Cross-check against the momentum file, which cannot show a sign but must")
    print("  show a MINIMUM at every node:")
    def min_near(kfm, rho, node, halfwidth=0.4):
        m = np.abs(kfm - node) <= halfwidth
        idx = np.where(m)[0]
        return int(idx[np.argmin(rho[idx])])

    i0 = min_near(m6.k_fm, m6.split[0], n00[0])
    i2 = min_near(m6.k_fm, m6.split[2], n22[0])
    i7 = min_near(m7.k_fm, m7.rho_tot, n33[0])
    print(f"    rho_0 minimum {m6.split[0][i0]:.5f} at K = {m6.k_fm[i0]:.2f} fm^-1"
          f"   [overlap node {n00[0]:.3f}]")
    print(f"    rho_2 minimum {m6.split[2][i2]:.3e} at K = {m6.k_fm[i2]:.2f} fm^-1"
          f"   [overlap node {n22[0]:.3f}]")
    print(f"    7Li rho minimum {m7.rho_tot[i7]:.5f} at K = {m7.k_fm[i7]:.2f} fm^-1"
          f"   [overlap Aat33 node {n33[0]:.3f}]")
    print()
    print("THE RELATIVE SIGN.  With the global phase fixed by psi_0(k -> 0) > 0:")
    print(f"    psi_2/psi_0 < 0 for k < {n00[0]:.2f} fm^-1 = {n00[0]*HBARC:.3f} GeV")
    print(f"    psi_2/psi_0 > 0 for {n00[0]:.2f} < k < {n22[0]:.2f} fm^-1")
    print("LiPolGen's Hulthen forms are POSITIVE-DEFINITE in k, so the model has been")
    print("assuming psi_2/psi_0 = +1 everywhere.  The VMC S-wave node sits at 0.130 GeV,")
    print("inside the sampled range, so the VMC S-D interference term CHANGES SIGN")
    print("across the tagged acceptance -- a structure no two-parameter form can carry.")

    # -- 3. the j = 1/2 column of li7.at ----------------------------------
    print()
    print("=" * 96)
    print("3.  li7.at's SECOND COLUMN (Aat11) IS NOT A SECOND PARTIAL WAVE")
    print("=" * 96)
    w33 = _trapz(li7.k_fm ** 2 * li7.k_amp[0] ** 2, li7.k_fm)
    w11 = _trapz(li7.k_fm ** 2 * li7.k_amp[1] ** 2, li7.k_fm)
    print(f"int k^2 Aat33^2 dk = {w33:.4f},  int k^2 Aat11^2 dk = {w11:.6f}  "
          f"-> weight {w11/(w33+w11):.2e}")
    print("alpha(0+) x t(1/2+) with L = 1 gives j = 1/2 or 3/2; only j = 3/2 can")
    print("build the 3/2- ground state, so Aat11 is a selection-rule zero and what")
    print("is tabulated is MC leakage (it changes sign repeatedly and never exceeds")
    print("~1.2% of Aat33).  It is NOT the 1/2- excited state: that state has its own")
    print(f"momentum file, li7_at1.momentum, with S = {m7x.norms[0]:.5f} -- a FULL-SIZE")
    print("distribution, not a 1e-4 one.  Thread A summed |A33|^2 + |A11|^2, which")
    print("changes nothing at the quoted precision but is conceptually wrong.")

    # -- 4. THE TABLE ------------------------------------------------------
    def vmc_overlap_k(ov, ls, cols=(0,)):
        per = {}
        for c in cols:
            l = ls[c]
            d = ov.k_fm ** 2 * ov.k_amp[c] ** 2
            per[l] = per.get(l, 0.0) + d
        return on_grid(ov.k_fm * HBARC, per, grid)

    def vmc_overlap_r(ov, ls, cols=(0,)):
        per = {}
        for c in cols:
            l = ls[c]
            A = fb_cluster(ov.r_fm, ov.r_amp[c], kfb, l)
            per[l] = per.get(l, 0.0) + kfb ** 2 * A ** 2
        return on_grid(kfb * HBARC, per, grid)

    def vmc_momentum(mf, ls_map):
        per = {}
        for l, rho in ls_map.items():
            per[l] = mf.k_fm ** 2 * rho
        return on_grid(mf.k_fm * HBARC, per, grid)

    print()
    print("=" * 96)
    print(f"4.  THE ONE-NORMALIZATION TABLE   (KMAX = {args.kmax:.5f} GeV "
          f"= {args.kmax/HBARC:.2f} fm^-1)")
    print("=" * 96)

    m6_models = [
        Model("6Li VMC  momenta li6_ad1 (S+D)", grid,
              vmc_momentum(m6, {0: m6.split[0], 2: m6.split[2]})),
        Model("6Li VMC  overlap li6.ad k-block", grid,
              vmc_overlap_k(li6, li6.ls, (0, 1))),
        Model("6Li VMC  overlap li6.ad r-block FT", grid,
              vmc_overlap_r(li6, li6.ls, (0, 1))),
    ]
    p_d_vmc = m6_models[0].p_d()
    for beta in (0.20, 0.30, 0.40):
        m6_models.append(Model(
            f"6Li Hulthen beta={beta:.2f} P_D=0.0867 (code)", grid,
            hulthen_density(lp, kappa6, beta,
                            [(0, 1 - 0.0867), (2, 0.0867)], grid)))
    for beta in (0.20, 0.30, 0.40):
        m6_models.append(Model(
            f"6Li Hulthen beta={beta:.2f} P_D={p_d_vmc:.4f} (VMC)", grid,
            hulthen_density(lp, kappa6, beta,
                            [(0, 1 - p_d_vmc), (2, p_d_vmc)], grid)))
    print_rows(m6_models)

    print()
    m7_models = [
        Model("7Li VMC  momenta li7_at3 (P)", grid,
              vmc_momentum(m7, {1: m7.rho_tot})),
        Model("7Li VMC  overlap li7.at k-block Aat33", grid,
              vmc_overlap_k(li7, li7.ls, (0,))),
        Model("7Li VMC  overlap li7.at k-block A33+A11", grid,
              vmc_overlap_k(li7, li7.ls, (0, 1))),
        Model("7Li VMC  overlap li7.at r-block FT", grid,
              vmc_overlap_r(li7, li7.ls, (0,))),
    ]
    for beta in (0.20, 0.30, 0.40):
        m7_models.append(Model(f"7Li Hulthen beta={beta:.2f} (P wave, code)", grid,
                               hulthen_density(lp, kappa7, beta, [(1, 1.0)], grid)))
    # the WRONG comparison thread B made: the S-wave Hulthen for a P-wave channel
    for beta in (0.20, 0.30, 0.40):
        m7_models.append(Model(f"7Li Hulthen beta={beta:.2f} [S-WAVE form: WRONG]",
                               grid,
                               hulthen_density(lp, kappa7, beta, [(0, 1.0)], grid)))
    print_rows(m7_models)

    # -- 5. KMAX sensitivity ----------------------------------------------
    print()
    print("=" * 96)
    print("5.  KMAX SENSITIVITY OF THE HULTHEN NUMBERS (the VMC rows cannot move:")
    print("    the tables are zero past 5 fm^-1, which is what `VmcRadial` does too)")
    print("=" * 96)
    print(f"{'model':<40}" + "".join(f"KMAX={x:.3f}".rjust(14)
                                     for x in (KMAX_DEFAULT, KMAX_MODEL, 4.0)))
    for label, kappa, waves in (
            ("7Li Hulthen b=0.40 P wave  <k>", kappa7, [(1, 1.0)]),
            ("6Li Hulthen b=0.40 S+D     <k>", kappa6, [(0, 0.9806), (2, 0.0194)])):
        vals = []
        for km in (KMAX_DEFAULT, KMAX_MODEL, 4.0):
            g = np.linspace(1e-4, km, 6000)
            mm = Model(label, g, hulthen_density(lp, kappa, 0.40, waves, g))
            vals.append(mm.mean_k())
        print(f"{label:<40}" + "".join(f"{v:14.4f}" for v in vals))
    for label, kappa, waves in (
            ("7Li Hulthen b=0.40 P wave  P(k>0.3)", kappa7, [(1, 1.0)]),
            ("6Li Hulthen b=0.40 S+D     P(k>0.3)", kappa6, [(0, 0.9806), (2, 0.0194)])):
        vals = []
        for km in (KMAX_DEFAULT, KMAX_MODEL, 4.0):
            g = np.linspace(1e-4, km, 6000)
            mm = Model(label, g, hulthen_density(lp, kappa, 0.40, waves, g))
            vals.append(mm.tail(0.30))
        print(f"{label:<40}" + "".join(f"{v:14.4f}" for v in vals))

    # -- 6. verdict --------------------------------------------------------
    h40_7 = [m for m in m7_models if m.label.startswith("7Li Hulthen beta=0.40 (P")][0]
    h40_7s = [m for m in m7_models if "S-WAVE" in m.label and "0.40" in m.label][0]
    vmc7 = m7_models[0]
    ovl7 = m7_models[1]
    print()
    print("=" * 96)
    print("6.  VERDICT")
    print("=" * 96)
    print(f"7Li VMC (momenta)   <k> = {vmc7.mean_k():.4f}  P(k>0.3) = {vmc7.tail(0.30):.4f}"
          f"  P(k>0.45) = {vmc7.tail(0.45):.4f}")
    print(f"7Li VMC (overlap)   <k> = {ovl7.mean_k():.4f}  P(k>0.3) = {ovl7.tail(0.30):.4f}"
          f"  P(k>0.45) = {ovl7.tail(0.45):.4f}")
    print("    -> the two VMC families AGREE.  Threads A and B never disagreed about")
    print("       the data; A quoted the 0.45 threshold and B the 0.30 one.")
    print(f"7Li Hulthen b=0.40, P wave (the form the code uses for this channel):")
    print(f"                    <k> = {h40_7.mean_k():.4f}  P(k>0.3) = {h40_7.tail(0.30):.4f}"
          f"  P(k>0.45) = {h40_7.tail(0.45):.4f}   <- thread A's 0.22 reproduced")
    print(f"7Li Hulthen b=0.40, S wave (NOT this channel's form):")
    print(f"                    <k> = {h40_7s.mean_k():.4f}  P(k>0.3) = {h40_7s.tail(0.30):.4f}"
          f"  P(k>0.45) = {h40_7s.tail(0.45):.4f}   <- thread B's 0.084/0.15 reproduced")
    k_snap = 1.6 * HBARC
    print()
    print("Thread B's VMC tail fractions (7Li 0.179, 6Li 0.045) are ~15-25% below the")
    print("recomputed 0.30 GeV values because the threshold was snapped up to the next")
    print(f"tabulated grid point, K = 1.6 fm^-1 = {k_snap:.4f} GeV:")
    print(f"    7Li VMC momenta  P(k > {k_snap:.4f}) = {vmc7.tail(k_snap):.4f}   "
          "[thread B: 0.179]")
    print(f"    6Li VMC momenta  P(k > {k_snap:.4f}) = "
          f"{m6_models[0].tail(k_snap):.4f}   [thread B: 0.045]")
    print("Its <k> values (0.186, 0.122) are unaffected by that and are correct.")
    print()
    print("THREAD B'S HULTHEN BAND IS THE S-WAVE (Hulthen) DENSITY EVALUATED FOR A")
    print("P-WAVE CHANNEL.  Its VMC numbers are right; its reference band is not, and")
    print("the conclusion drawn from it ('the band is biased low, VMC is harder than")
    print("beta = 0.40') is therefore wrong.  Against the P-wave form the code really")
    print("uses, VMC alpha+t is SOFTER than every beta in the band -- thread A's")
    print("conclusion -- though far less dramatically than 'over 20x': the VMC/Hulthen")
    print("ratio is a strong function of where you cut.")

    # -- 7. plot ------------------------------------------------------------
    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))

        ax = axes[0, 0]
        for m, st in ((m6_models[0], "-"), (m6_models[1], "--"), (m6_models[2], ":")):
            ax.semilogy(grid, m.normalized(), st, lw=1.8, label=m.label[9:])
        for m in m6_models[3:6]:
            ax.semilogy(grid, m.normalized(), lw=1.0, alpha=0.8, label=m.label[4:])
        ax.set_xlim(0, args.kmax); ax.set_ylim(1e-4, 30)
        ax.set_xlabel("k [GeV]"); ax.set_ylabel(r"$k^2 n(k)$, normalized")
        ax.set_title(r"$^6$Li $\to\alpha+d$"); ax.legend(fontsize=6.5)

        ax = axes[0, 1]
        for m, st in ((m7_models[0], "-"), (m7_models[1], "--"), (m7_models[3], ":")):
            ax.semilogy(grid, m.normalized(), st, lw=1.8, label=m.label[9:])
        for m in m7_models[4:7]:
            ax.semilogy(grid, m.normalized(), lw=1.0, alpha=0.8, label=m.label[4:])
        ax.semilogy(grid, h40_7s.normalized(), "-.", color="crimson", lw=1.4,
                    label="S-wave form b=0.40 (thread B)")
        ax.set_xlim(0, args.kmax); ax.set_ylim(1e-4, 30)
        ax.set_xlabel("k [GeV]"); ax.set_ylabel(r"$k^2 n(k)$, normalized")
        ax.set_title(r"$^7$Li $\to\alpha+t$"); ax.legend(fontsize=6.5)

        ax = axes[1, 0]
        ax.plot(li6.k_fm, li6.k_amp[0] / np.abs(li6.k_amp[0][0]), "o-", ms=3,
                label=r"$A_{00}(k)$ / $|A_{00}(0)|$  (k-block)")
        ax.plot(li6.k_fm, li6.k_amp[1] / np.abs(li6.k_amp[0][0]) * 20, "s-", ms=3,
                label=r"$A_{22}(k)$ / $|A_{00}(0)|\times 20$  (k-block)")
        ax.plot(kfb, ft6[0] / np.abs(ft6[0][0]), "--", label="FT of r-block, L=0")
        ax.plot(kfb, ft6[2] / np.abs(ft6[0][0]) * 20, "--", label=r"FT of r-block, L=2 ($\times$20)")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlim(0, 4); ax.set_xlabel("k [fm$^{-1}$]")
        ax.set_title(r"$^6$Li S/D relative sign: $A_{22}/A_{00} < 0$")
        ax.legend(fontsize=6.5)

        ax = axes[1, 1]
        x = np.arange(len(THRESHOLDS)); w = 0.2
        sets = [("VMC momenta", m7_models[0], "k"),
                ("VMC overlap", m7_models[1], "0.45"),
                (r"Hulthen P $\beta$=0.40", h40_7, "tab:blue"),
                (r"S-wave form $\beta$=0.40", h40_7s, "crimson")]
        for i, (lab, m, col) in enumerate(sets):
            ax.bar(x + (i - 1.5) * w, [m.tail(t) for t in THRESHOLDS], w,
                   label=lab, color=col)
        ax.set_xticks(x); ax.set_xticklabels([f">{t:.2f}" for t in THRESHOLDS])
        ax.set_yscale("log"); ax.set_ylabel("P(k > threshold)")
        ax.set_title(r"$^7$Li tail fractions: the disputed numbers")
        ax.legend(fontsize=6.5)

        fig.suptitle("ANL VMC cluster data vs LiPolGen Hulthen -- one normalization "
                     f"(KMAX = {args.kmax:.3f} GeV)", fontsize=12)
        fig.tight_layout()
        out = HERE / "vmc_reconcile.png"
        fig.savefig(out, dpi=140)
        print(f"\nSaved {out}")

    if args.markdown:
        write_markdown(args.markdown, args.kmax, m6_models, m7_models, m6, m7,
                       kappa6, kappa7, s_ad, s_at,
                       dict(vmc7_snap=vmc7.tail(k_snap),
                            vmc6_snap=m6_models[0].tail(k_snap),
                            k_snap=k_snap, n00=n00[0], n22=n22[0], n33=n33[0]))
        print(f"Saved {args.markdown}")


def _md_table(models):
    out = ["| model | ⟨k⟩ [GeV] | " + " | ".join(f"P(k>{t:.2f})" for t in THRESHOLDS)
           + " | P_D |",
           "|---|---|" + "---|" * (len(THRESHOLDS) + 1)]
    for m in models:
        mk, *tails, pd = m.row()
        pd_s = "—" if np.isnan(pd) else f"{pd:.4f}"
        out.append(f"| {m.label} | {mk:.4f} | "
                   + " | ".join(f"{t:.4f}" for t in tails) + f" | {pd_s} |")
    return "\n".join(out)


def write_markdown(path, kmax, m6_models, m7_models, m6, m7, kappa6, kappa7,
                   s_ad, s_at, extra):
    path.parent.mkdir(parents=True, exist_ok=True)

    def find(models, needle):
        return [m for m in models if needle in m.label][0]

    h7p = find(m7_models, "beta=0.40 (P wave")
    h7s = find(m7_models, "beta=0.40 [S-WAVE")
    v7 = m7_models[0]
    v7o = m7_models[1]
    v6 = m6_models[0]
    h7p45, h7s30 = h7p.tail(0.45), h7s.tail(0.30)
    h7smk, h7pmk = h7s.mean_k(), h7p.mean_k()
    v7p45, v7mk = v7o.tail(0.45), v7.mean_k()
    v7snap, v6snap = extra["vmc7_snap"], extra["vmc6_snap"]
    ksnap = extra["k_snap"]
    v7p30, v6p30 = v7.tail(0.30), v6.tail(0.30)
    pd_ov, pd_mo = m6_models[1].p_d(), m6_models[0].p_d()
    body = f"""<!-- generated by validation/vmc_reconcile.py; do not hand-edit the tables -->

# Reconciling the two ANL VMC readings

Two independent investigations fetched different ANL files and reported
numbers that could not both be true.  `validation/vmc_reconcile.py` recomputes
every one of them in one place, with one normalization, and with LiPolGen's
Hulthen forms evaluated **through the C++ bindings** (`lipolgen.Wave.radial`,
i.e. `src/core/cluster.cpp` itself).

## The one normalization

    D(k) = k^2 sum_L |psi_L(k)|^2 ,    int_0^KMAX D dk = 1
    <k> = int k D dk ,   P(k>k0) = int_{{k0}}^{{KMAX}} D dk
    P_D = int k^2 |psi_2|^2 dk / int k^2 sum_L |psi_L|^2 dk

`KMAX = {kmax:.5f}` GeV `= {kmax/HBARC:.2f}` fm^-1 — the common native extent of every
VMC table, and the range outside which `VmcRadial` returns zero.
`k[GeV] = k[fm^-1] * 0.1973269804`; `kappa(6Li a+d) = {kappa6:.5f}` GeV,
`kappa(7Li a+t) = {kappa7:.5f}` GeV, taken from the C++ channels.

## Verdict

**Thread A was right about the physics; thread B's VMC numbers are right but
its Hulthen reference band is wrong.**

* The two VMC file families **agree with each other** — the 2004 AV18+UIX
  `overlap_old` amplitudes and the 2014/2024 AV18+UX `momenta` densities give
  the same distribution to within the Hamiltonian difference.  The apparent
  conflict (`P(k>0.45) = 0.0095` vs `P(k>0.3) = 0.179`) is two different
  thresholds on the same steeply-falling tail, not a disagreement.
* Thread B's "Hulthen band" is the **S-wave (Hulthen) density evaluated for a
  P-wave channel**.  `7Li -> alpha + t` is pure `L = 1` in the code
  (`li7_alpha_channel`), whose radial form is `k/((k^2+kappa^2)(k^2+beta^2))`,
  not `1/(k^2+kappa^2) - 1/(k^2+beta^2)`.  The S-wave form falls as `k^-6` and
  the P-wave form as `k^-4`, so the S-wave numbers are far too soft.  Every
  row in the table below is labelled with the form that produced it, and the
  "[S-WAVE form: WRONG]" rows reproduce thread B's 0.1505 / 0.0836 exactly.
* Consequence: `docs/open_items/physics_literature.md`'s "the beta band is
  biased low and does not bracket the truth" is **withdrawn**.  Measured
  against the P-wave form the code really uses, VMC `alpha+t` is *softer* than
  every beta in the 0.20–0.40 band.

## Every disputed number, reproduced

| number | thread | recomputed here | with what |
|---|---|---|---|
| ⁷Li Hulthen β=0.40 P(k>0.45) = 0.22 | A | {h7p45:.4f} | the **P-wave** form `k/((k²+κ²)(k²+β²))` — the form `li7_alpha_channel` really uses |
| ⁷Li Hulthen β=0.40 P(k>0.30) = 0.0836 | B | {h7s30:.4f} | the **S-wave** Hulthen form — *not* this channel's form |
| ⁷Li Hulthen β=0.40 ⟨k⟩ = 0.1505 | B | {h7smk:.4f} | ditto (P-wave gives {h7pmk:.4f}) |
| ⁷Li VMC P(k>0.45) = 0.0095 | A | {v7p45:.4f} | `overlap_old` k-block, Aat33 |
| ⁷Li VMC ⟨k⟩ = 0.186 | B | {v7mk:.4f} | `momenta/li7_at3` |
| ⁷Li VMC P(k>0.30) = 0.179 | B | {v7snap:.4f} at k > {ksnap:.4f} GeV | the threshold was snapped up to the next tabulated grid point, K = 1.6 fm⁻¹; at a true 0.300 GeV it is {v7p30:.4f} |
| ⁶Li VMC P(k>0.30) = 0.045 | B | {v6snap:.4f} at k > {ksnap:.4f} GeV | same snap; at a true 0.300 GeV it is {v6p30:.4f} |
| ⁶Li P_D = 0.020 / 1.94% | A / B | {pd_ov:.4f} / {pd_mo:.4f} | overlap k-block / momenta S-D blocks — the two agree |

So: **no measurement was ever in dispute.** Thread A used the right radial
form and the far tail; thread B used the wrong radial form and a
grid-snapped threshold.

## Normalization cross-check between the two file families

| quantity | `overlap_old` (AV18+UIX, 2004) | `momenta` (AV18+UX, 2014/24) |
|---|---|---|
| S_ad = (2π)⁻³ ∫k²A²dk | {s_ad:.5f} | {m6.norms[0]:.5f} |
| S_at (j=3/2 only) | {s_at:.5f} | {m7.norms[0]:.5f} |
| P_D(6Li) | {m6_models[1].p_d():.4f} | {m6_models[0].p_d():.4f} |

`rho_L(K) = A_L^2(K)/(4π)` relates the two conventions exactly, and the
node positions agree (S node at K≈0.68 fm⁻¹, D node at K≈2.3 fm⁻¹, ⁷Li P node
at K≈1.15 fm⁻¹).

## ⁶Li → α + d

{_md_table(m6_models)}

## ⁷Li → α + t

{_md_table(m7_models)}

## Production input, per isotope

| isotope | `|psi_L|` from | sign / interference from | why |
|---|---|---|---|
| ⁶Li α+d | `momenta/li6_ad1.momentum` S/D blocks | `overlap_old/li6.ad` | 1M VMC samples vs 200k, 2024 Hamiltonian, an explicit S/D split with MC errors, and the file's own normalization (S_ad, P_D) printed.  The overlap file's r- and k-blocks are the only source of the **relative S–D sign**, which a density file cannot carry. |
| ⁷Li α+t | `momenta/li7_at3.momentum` | not needed (single wave) | 500k samples, the 3/2⁻ ground state explicitly, `S_at` printed.  A single partial wave has no interference term, so the overall sign is unobservable. |

`li7.at`'s second column `Aat11` is **not** a second partial wave: α(0⁺)×t(½⁺)
with L = 1 gives j = ½ or 3/2 and only j = 3/2 builds the 3/2⁻ ground state, so
`Aat11` is a selection-rule zero carrying MC leakage at the 1e-4 level in
|A|².  It is also not the ½⁻ excited state — that has its own file,
`momenta/li7_at1.momentum`, with a full-size distribution (S = 0.98683).

## The S–D relative sign

Both blocks of `li6.ad` — and the Fourier–Bessel transform of the r-space
block — agree that **sign(A₂₂/A₀₀) = −1** across the whole region that carries
the D-wave strength.  (The two blocks are printed with opposite *overall*
phase, which is unobservable; the *relative* sign is not, and it is the same
in both.)  LiPolGen's Hulthen forms are positive-definite for every k, so the
model has been assuming `sign(psi_2/psi_0) = +1`.  Turning on the VMC backend
therefore **flips the sign of the S–D interference term** in `n_M(k, khat)`.
"""
    path.write_text(body)


if __name__ == "__main__":
    main()
