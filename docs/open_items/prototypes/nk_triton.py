#!/usr/bin/env python3
"""Ciofi degli Atti - Simula A=3 n(k) vs LiPolGen's Hulthen triton breakup.

SOURCE OF THE CS PARAMETERS (transcribed verbatim, no refit):
  /home/cpeng/Projects/polli/BeAGLE/source/src/dpmjet3.0-5F-new.f
  SUBROUTINE DT_KFERMI, lines 17216-17403, credited in the code to
  "Claudio Ciofi & S. Simula, Phys. Rev. C 53 (1996) 1689".

FUNCTIONAL FORM (k in fm^-1, n(k) in fm^3):

    n(k) = sum_{i=0,1,2}  A_i exp(-B_i k^2) / (1 + C_i k^2)^2

NORMALISATION.  The code INTEGRATES 4 pi k^2 n(k) dk (DT_KFERMI:17363) and
then renormalises numerically, so BeAGLE never notices what the A_i mean.
Measured here, the coded A_i satisfy

    int_0^inf n(k) k^2 dk  =  1.0031  (A=2)   0.6526  (A=3)   0.7996  (A=4)

i.e. the CS convention is int n k^2 dk = 1 WITHOUT the 4 pi, and the A=3 / A=4
sets are only the n_0(k) piece -- the part of n(k) in which the residual
(A-1) system is left in its GROUND state.  That is exactly what the source
comment says ("These are n0k parametrization not including n1k"), and it makes
the deficit a PHYSICS number, not a transcription error:

    A=3:  0.653 of the strength leaves a bound (A-1) = d      <- 2-body breakup
          0.347                        leaves the continuum   <- 3-body breakup
    A=4:  0.800 leaves a bound 3H/3He, 0.200 the continuum
    A=2:  1.003 -- the deuteron has no excitable remnant, as it must be.

So the CS parameterisation itself hands us the two-body/three-body split that
LiPolGen's triton branch is missing, and 0.653 sits right on the ~2/3 the
3He(e,e\'p)d literature quotes.

CONSEQUENCE FOR BeAGLE: because DT_KFERMI renormalises n_0 to unity, BeAGLE
samples the A=3 and A=4 Fermi momentum from the GROUND-STATE-REMNANT
distribution alone and drops the 35 % / 20 % correlated continuum -- its high-k
tail is too soft by construction.

COEFFICIENTS as coded:

  A=2 (deuteron): A0=157.4  B0=1.24  C0=18.3
                  A1=0.234  B1=1.27  C1=0
                  A2=0.00623 B2=0.220 C2=0
  A=3 (3He AND 3H - BeAGLE branches on the MASS NUMBER only, so the
       triton gets the 3He set; there is NO separate 3H set in BeAGLE):
                  A0=31.7   B0=1.32  C0=5.98
                  A1=0.00266 B1=0.365 C1=0
                  A2=0
  A=4 (alpha):    A0=4.33   B0=1.54  C0=0.419
                  A1=5.49   B1=4.90  C1=0
                  A2=0

WHAT LiPolGen DOES TODAY (src/core/breakup.cpp, the Triton branch):
  |k| ~ k^2 |psi_L0(k;kappa,beta)|^2 with the Hulthen form of cluster.hpp,
      psi_L0(k) = 1/(k^2+kappa^2) - 1/(k^2+beta^2),  beta = 0.30 GeV,
  kappa = sqrt(2 mu S) per channel:
      t* -> n + d      S = m_n + m_d   - m_t = 6.2572 MeV
      t* -> p + (nn)   S = m_p + 2 m_n - m_t = 8.4818 MeV  (mu with the nn
                       pair treated as a single fragment of mass 2 m_n)
  and the nn pair is then split at the virtual-state pole kappa_nn = 10.44 MeV.
"""

import numpy as np
from scipy.integrate import quad

HBARC = 0.1973269804        # GeV fm

CS = {
    2: [(157.4, 1.24, 18.3), (0.234, 1.27, 0.0), (0.00623, 0.220, 0.0)],
    3: [(31.7, 1.32, 5.98), (0.00266, 0.365, 0.0)],
    4: [(4.33, 1.54, 0.419), (5.49, 4.90, 0.0)],
}


def n_cs(k_fm, A):
    """CS n(k) [fm^3] at k [fm^-1]."""
    k2 = np.asarray(k_fm, float) ** 2
    tot = 0.0
    for (a, b, c) in CS[A]:
        tot = tot + a * np.exp(-b * k2) / (1.0 + c * k2) ** 2
    return tot


def norm_cs(A, kmax=50.0):
    """int_0^kmax n(k) k^2 dk in the CS convention (no 4 pi)."""
    f = lambda k: k * k * float(n_cs(k, A))
    v, _ = quad(f, 0.0, kmax, limit=800)
    return v


def p_above_cs(k_gev, A, kmax=10.0):
    """P(|k| > k_gev) for the CS distribution, k converted with hbar c."""
    k0 = k_gev / HBARC
    f = lambda k: 4.0 * np.pi * k * k * float(n_cs(k, A))
    num, _ = quad(f, k0, kmax, limit=400)
    den, _ = quad(f, 0.0, kmax, limit=400)
    return num / den


# ------------------------------------------------------------- LiPolGen side
M_P, M_N, M_D, M_T = 0.9382721, 0.9395654, 1.8756129, 2.8089210
BETA = 0.30


def kappa(mu, S):
    return np.sqrt(2.0 * mu * S)


def psi_hulthen(k, kap, beta=BETA):
    return 1.0 / (k**2 + kap**2) - 1.0 / (k**2 + beta**2)


def p_above_hulthen(k_gev, kap, beta=BETA, kmax=1.2):
    """P(|k|>k_gev) for k^2 |psi|^2, on the SAME k_max the C++ grid uses."""
    f = lambda k: k * k * psi_hulthen(k, kap, beta) ** 2
    num, _ = quad(f, k_gev, kmax, limit=400)
    den, _ = quad(f, 0.0, kmax, limit=400)
    return num / den


def mean_k_hulthen(kap, beta=BETA, kmax=1.2):
    f = lambda k: k * k * psi_hulthen(k, kap, beta) ** 2
    num, _ = quad(lambda k: k * f(k), 0.0, kmax, limit=400)
    den, _ = quad(f, 0.0, kmax, limit=400)
    return num / den


def mean_k_cs(A, kmax=10.0):
    f = lambda k: 4.0 * np.pi * k * k * float(n_cs(k, A))
    num, _ = quad(lambda k: k * f(k), 0.0, kmax, limit=400)
    den, _ = quad(f, 0.0, kmax, limit=400)
    return HBARC * num / den


def main():
    print("=" * 74)
    print("1. Normalisation of the CS sets:  int n(k) k^2 dk   (CS convention,")
    print("   no 4 pi).  The DEFICIT below 1 is the strength that leaves the")
    print("   residual (A-1) system EXCITED, i.e. the three-body-breakup share.")
    print("=" * 74)
    for A in (2, 3, 4):
        v = norm_cs(A)
        print(f"   A = {A}:  norm(n_0) = {v:7.4f}   -> P(2-body, bound "
              f"remnant) = {v*100:5.1f} %   P(3-body) = {(1-v)*100:5.1f} %")
    print()

    print("=" * 74)
    print("2. P(k > k_cut) tails.  CS = BeAGLE/Ciofi-Simula;  Hulthen = "
          "LiPolGen breakup.cpp")
    print("=" * 74)
    mu_nd = M_N * M_D / (M_N + M_D)
    s_nd = M_N + M_D - M_T
    kap_nd = kappa(mu_nd, s_nd)

    m_nn = 2.0 * M_N
    mu_pnn = M_P * m_nn / (M_P + m_nn)
    s_pnn = M_P + 2.0 * M_N - M_T
    kap_pnn = kappa(mu_pnn, s_pnn)

    mu_ad = 3.7273794 * M_D / (3.7273794 + M_D)
    kap_ad = kappa(mu_ad, 1.4743e-3)

    print(f"   kappa(t*->n+d)     = {kap_nd*1e3:7.2f} MeV   "
          f"(S = {s_nd*1e3:.4f} MeV)")
    print(f"   kappa(t*->p+(nn))  = {kap_pnn*1e3:7.2f} MeV   "
          f"(S = {s_pnn*1e3:.4f} MeV)")
    print(f"   kappa(6Li->a+d)    = {kap_ad*1e3:7.2f} MeV   (reference)")
    print()

    cuts = [0.10, 0.20, 0.30, 0.45]
    rows = [
        ("CS A=3 (3He/3H)", lambda c: p_above_cs(c, 3)),
        ("CS A=2 (deuteron)", lambda c: p_above_cs(c, 2)),
        ("CS A=4 (alpha)", lambda c: p_above_cs(c, 4)),
        ("Hulthen t*->n+d", lambda c: p_above_hulthen(c, kap_nd)),
        ("Hulthen t*->p+nn", lambda c: p_above_hulthen(c, kap_pnn)),
        ("Hulthen nn split", lambda c: p_above_hulthen(c, 0.0104399)),
    ]
    hdr = f"   {'distribution':<20}" + "".join(f"{'k>'+str(c):>12}" for c in cuts)
    print(hdr)
    print("   " + "-" * (20 + 12 * len(cuts)))
    vals = {}
    for name, fn in rows:
        v = [fn(c) for c in cuts]
        vals[name] = v
        print(f"   {name:<20}" + "".join(f"{x:12.4g}" for x in v))
    print()
    print("   ratio Hulthen(n+d) / CS(A=3):  " +
          "  ".join(f"{vals['Hulthen t*->n+d'][i]/vals['CS A=3 (3He/3H)'][i]:.3g}"
                    for i in range(len(cuts))))
    print()

    print("=" * 74)
    print("3. <k>")
    print("=" * 74)
    print(f"   CS A=3            : {mean_k_cs(3)*1e3:7.1f} MeV")
    print(f"   CS A=2            : {mean_k_cs(2)*1e3:7.1f} MeV")
    print(f"   CS A=4            : {mean_k_cs(4)*1e3:7.1f} MeV")
    print(f"   Hulthen t*->n+d   : {mean_k_hulthen(kap_nd)*1e3:7.1f} MeV")
    print(f"   Hulthen t*->p+nn  : {mean_k_hulthen(kap_pnn)*1e3:7.1f} MeV")
    print()
    print("   (the Hulthen numbers are truncated at the C++ grid ceiling "
          "k_max = 1.2 GeV;")
    print("    the CS numbers run to 10 fm^-1 = 1.97 GeV)")


if __name__ == "__main__":
    main()
