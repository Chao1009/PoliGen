# Phase G(i) — the numbers of `design_G_cluster_config.md` §8, reproduced

Measured on 2026-09-03 from the committed tables, through
`include/lipolgen/cluster_config.hpp` / `src/core/cluster_config.cpp` and the
`lipolgen-configs` CLI.  Every "measured" column below is what the code
returns; every "design" column is the value the design document states.
Where they differ the row is marked and the reason is given.

Quadrature convention, stated once because it explains every small
disagreement: the **analytic layer** is `np.trapz`-identical trapezoid
(`numerics.hpp`'s `trapezoid`) on **each table's own abscissa** — 0.1 fm for
`li6.ad` / `li6.adr.fit` / `he4.density`, 0.01 fm for `fdeut.av18` — with the
tables identically zero outside their own range.  That is the rule that
reproduces §8.  The sampler's own 512 × 96 cell grid is a separate,
midpoint-rule object and is compared to it in the "grid" rows below.

---

## 1. Inputs (the tables' own printed numbers vs the re-integration)

| quantity | file's own header | measured | note |
|---|---|---|---|
| ∫(u²+w²)dr | 1 | **0.9999981** | ✓ (5e-5 gate, T2) |
| P_D(d) `dstate` | 0.057599 | **0.0575986** | ✓ |
| Q_d = 𝒬/4 `qm` | 0.269673 | **0.2696703** fm² | ✓ 1.2e-5 rel |
| r_d = ½√⟨r_np²⟩ `rd` | 1.967364 | **1.9673386** fm | ✓ 1.3e-5 rel |
| ⟨r_np²⟩ | — | **15.48169** fm² | design 15.4817 ✓ |
| `li6.ad` S-wave norm | 0.838 | **0.837052** | ✓ (2e-3 gate, T3) |
| `li6.ad` D-wave norm | 0.017 | **0.017180** | ✓ |
| `li6.ad` `ndx` = S_αd | 0.856 | **0.854232** | design 0.85423 ✓ |
| `li6.adr.fit` ∫, P_D | — (no header) | **0.837301**, **0.025425** | design 0.83730 / 0.02542 ✓ |
| `FitRescaled` rescale factors | — | **s₀ = 1.01281, s₂ = 0.89835** | design 1.0128 / 0.8984 ✓ |
| `he4.density` 4π∫ρR²dR | 1.9974 | **1.99738** | ✓ |
| `he4.density` rms | 1.4404 | **1.44131** | design 1.4413 ✓ |
| ⟨s²⟩_α (recentred, λ² = 4/3) | — | **2.07739** fm² | design 2.0774 ✓ |
| `li6.density` 4π∫ρR²dR | 2.9991 | **2.99909** | ✓ |
| `li6.density` rms | **2.4433** | 2.44331 | = `LI6_R_POINT_VMC_FM` ✓ |
| α Gaussian ⟨s²⟩ = 3a² | — | **2.10019** fm² | design 2.1003 ✓ |

Both fetched files verify byte-for-byte against the design's md5s:
`he4.density` 4751 B `5d717cf0ed642dc7779aef3873e06317`, `li6.density`
5189 B `96c6b15f045785db0593c00e9a32f9b1`.

---

## 2. The §8 table, all three α–d sources

| quantity | design `FitRescaled` | **measured** | design `OverlapRaw` | **measured** | design `FitRaw` | **measured** |
|---|---|---|---|---|---|---|
| P_D(α–d) | 0.02011 | **0.02011** | 0.02011 | **0.02011** | 0.02542 | **0.02542** |
| ⟨R²⟩_αd [fm²] | 16.963 | **16.9633** | 17.548 | **17.5482** | 16.956 | **16.9564** |
| rms R [fm] | 4.119 | **4.1187** | 4.189 | **4.1891** | 4.118 | **4.1178** |
| 𝒬[R₀,R₂] [fm²] | −1.3204 | **−1.3204** | −1.4009 | **−1.4009** | −1.4895 | **−1.4895** |
| 𝒬 interference | −1.2570 | **−1.2573** | −1.3324 | **−1.3325** | −1.4094 | **−1.4098** |
| 𝒬 pure-D | −0.0631 | **−0.0631** | −0.0684 | **−0.0684** | −0.0798 | **−0.0798** |
| D_T | 0.98190 | **0.98190** | 0.98190 | **0.98190** | 0.97712 | **0.97712** |
| ⟨r²⟩ [fm²] | 6.4447 | **6.4447** | 6.5747 | **6.5747** | 6.4432 | **6.4432** |
| r_rms [fm] | 2.5386 | **2.5386** | 2.5641 | **2.5641** | 2.5383 | **2.5383** |
| Q_matter(±1) [fm²] | −1.2309 | **−1.2309** | −1.3382 | **−1.3382** | −1.4590 | **−1.4590** |
| Q_matter(0) [fm²] | +2.4618 | **+2.4618** | +2.6764 | **+2.6765** | +2.9180 | **+2.9181** |
| Q_charge(+1) [fm²] | −0.6154 | **−0.6154** | −0.6691 | **−0.6691** | −0.7295 | **−0.7295** |
| δ⊥ per nucleon [fm²] | −0.1026 | **−0.1026** | −0.1115 | **−0.1115** | −0.1216 | **−0.1216** |
| eps_b0 equiv (B = 52.0) | −0.0506 | **−0.0506** | −0.0550 | **−0.0550** | −0.0600 | **−0.0600** |
| a₂(±1) at \|t\| = 0.3 | +0.197 | **+0.1976** | +0.215 | **+0.2148** | +0.234 | **+0.2342** |
| α–d share of \|Q\| | 76.9 % | **76.9 %** | 77.9 % | **77.9 %** | 79.0 % | **79.0 %** |
| … of which S–D interference | 95.2 % | **95.2 %** | 95.1 % | **95.1 %** | 94.6 % | **94.6 %** |
| (G9) root s for Q_ch = −0.0818 | 0.4032 | **0.4032** | 0.3806 | **0.3806** | 0.3577 | **0.3576** |
| … the P_D′ it implies | 3.3e-3 | **3.325e-3** | 3.0e-3 | **2.965e-3** | 3.3e-3 | **3.325e-3** |
| per-config sd of Σ(3z²−r²) | 28.9 | **28.78** | ≈ 29 | **29.28** | ≈ 29 | **28.67** |
| per-config sd of ⟨r²⟩ | 3.66 | **3.642** | ≈ 3.7 | **3.677** | ≈ 3.7 | **3.642** |
| asymptotic η (R = 6–8 fm) | −0.042…−0.053 | **−0.0482** | −0.043…−0.072 | **−0.0538** | "as FitRescaled" | **−0.0544 ✗** |
| `match_li6_radius()` scale | — | **0.93488** | — | **0.91917** | — | **0.93507** |

**One row disagrees.** The design's η column says `FitRaw` is "as
`FitRescaled`"; it is not, and cannot be. `FitRescaled` multiplies the fit's
D wave by s₂ = 0.8984 and its S wave by s₀ = 1.0128, so every R₂/R₀ ratio —
and therefore η — is smaller by s₂/s₀ = 0.8870. Measured: −0.0544 × 0.8870 =
−0.0483, which is the `FitRescaled` value to 4 digits. The design's `FitRaw`
column is the only §8 entry that does not survive; the physics conclusion
(η ≈ −0.05, ≈ 2× the measured −0.025 ± 0.006 ± 0.010, not 5–15×) is unchanged
for the default source.

### 2.1 The Whittaker division, term by term

κ, η_c and μ are **derived** from `LI6_ALPHA_TAG()` (`spectator.hpp`), not
retyped: κ = ch.kappa()/`HBARC_GEV_FM`, η_c = 2 α_EM μ / κ[GeV].

| quantity | design | measured |
|---|---|---|
| E_sep(⁶Li → α+d) | 1.4743 MeV | **1.474300 MeV** |
| μ | 1247.7 MeV | **1247.75 MeV** |
| κ | 0.3074 fm⁻¹ | **0.30738 fm⁻¹** |
| η_c (Sommerfeld) | 0.3002 | **0.30023** |
| W₂/W₀ at R = 6 / 7 / 8 fm | 3.34 / 2.92 / 2.63 | **3.3434 / 2.9233 / 2.6265** |
| R₂/R₀ (`FitRescaled`) at 6/7/8 | −0.139 / −0.148 / −0.139 | **−0.1396 / −0.1490 / −0.1364** |
| η (`FitRescaled`) at 6/7/8 | −0.042 / −0.051 / −0.053 | **−0.0418 / −0.0510 / −0.0519** |
| R₂/R₀ (`OverlapRaw`) at 6/7/8 | −0.145 / −0.139 / −0.189 | **−0.1307 / −0.1527 / −0.1837** |
| η (`OverlapRaw`) at 6/7/8 | −0.043 / −0.048 / −0.072 | **−0.0391 / −0.0522 / −0.0699** |

The `OverlapRaw` per-point values differ from the design's in the third digit —
that column is the **noisy raw block** read at exactly R = 6, 7, 8 fm by linear
interpolation between the 5.95/6.05 nodes, and the design labels it
"indicative". The mean, −0.0538, sits inside the design's own quoted
−0.043…−0.072 band. `tricomi_u` was validated against SciPy's `hyperu`
(U(1.3002, 2, 3.6888): 0.16875783 vs 0.16875538).

### 2.2 The a₂ map against the only published polarized coherent calculation

`a2_from_quadrupole(2 Q_d, A = 2, |t|, m)`, zero free parameters:

| \|t\| [GeV²] | a₂(±1) pred | digitized | dev | a₂(0) pred | digitized | dev |
|---|---|---|---|---|---|---|
| 0.05 | **−0.0433** | −0.04 | 8 % | **+0.0866** | +0.08 | 8 % |
| 0.10 | **−0.0866** | −0.08 | 8 % | **+0.1731** | +0.15 | 15 % |
| 0.20 | **−0.1731** | −0.17 | 2 % | **+0.3463** | +0.30 | 15 % |
| 0.30 | **−0.2597** | −0.28 | 7 % | **+0.5194** | +0.43 | 21 % |

Design: "≤ 8 % at m = ±1, 8–21 % at m = 0" — reproduced exactly.
B(⁶Li) = `gaussian_slope(√LI6_R2_POINT_FM2)` = **52.0368 GeV⁻²** (design 52.04).

---

## 3. The grid quadrature vs the analytic layer (T7a / T8a)

The design asks for **1e-6**. That gate is not reachable and the reason is
arithmetic, not physics:

| quantity | design gate | measured (n_r = 512, n_c = 96) | at n_r = 4096 | at n_c = 3072 |
|---|---|---|---|---|
| ⟨r²⟩ grid / analytic − 1 | 1e-6 | **+6.46e-4** | +6.59e-4 | +6.46e-4 |
| Q_matter grid / analytic − 1 | 1e-6 | **+2.98e-3** | +2.97e-3 | **+8.7e-4** |
| ⟨P₂ \| M=+1, m_S=−1⟩ + 2/7 | 1e-6 | **+4.4e-8** ✓ | +4.4e-8 | +4.1e-14 |
| ⟨P₂ \| M=+1, m_S=0⟩ − 1/7 | 1e-6 | **+2.3e-4** | +2.3e-4 | **+2.3e-7** |

Two independent, purely numerical residues:

* **The radial floor, +6.5e-4, which does not shrink with `n_r`.** The
  analytic layer is a trapezoid on the table's own h = 0.1 fm nodes; the cell
  sum converges to the **exact** integral of the piecewise-linear interpolant.
  For the square of a piecewise-linear f the node trapezoid is **high** on
  ∫f² x² dx by ≈ (h²/6)∫(f′)² x² dx; it inflates the r⁴ moment *less* than the
  norm, so the analytic ⟨R²⟩ comes out **low** by ≈ 5e-4 relative to the exact
  integral of the interpolant — which is what the cell sum converges to.
  (Measured on the committed `FitRescaled` table: N_trap/N_exact − 1 =
  **+4.4e-4**, M_trap/M_exact − 1 = **−3.8e-5**, ⟨R²⟩ ratio **+4.8e-4**;
  +6.5e-4 on (G6) once the α and p–n pieces are in.)  An h² effect at the
  *table's* spacing, which no refinement of the *sampler's* grid can remove.
* **The cos θ midpoint residue, O(h_c²).** `|Θ₂¹|² ∝ (1−c²)c²` has a *simple*
  zero at c = ±1, so the midpoint rule leaves an O(h_c²) term; `|Θ₂²|² ∝
  (1−c²)²` has a double zero and leaves none. Refining n_c 96 → 3072 (32×,
  h_c² down 1024×) drives the m_S = 0 residue 2.3e-4 → 2.3e-7, which is the
  proof that it is quadrature.

**Gates as implemented**: T7a 2e-3, T8a 5e-3 (plus a 1.5e-3 gate at
n_c = 3072), T21 1e-6 for m_S = −1 and 1e-3 for m_S = 0 **plus** a 1e-6 gate on
the n_c = 3072 refinement. Every MC gate is unchanged at 5σ of the
sampler's own per-configuration sd, as the design requires.

---

## 4. The Monte Carlo closure (200 000 configurations, m = +1, seed 424242)

| source | sampled ⟨r²⟩ | (G6) | 5σ | sampled Q_matter | (G5) | 5σ |
|---|---|---|---|---|---|---|
| `FitRescaled` | **6.4483** | 6.4447 | 0.0407 | **−1.2216** | −1.2309 | 0.3217 |
| `OverlapRaw` | **6.5771** | 6.5747 | 0.0411 | **−1.3308** | −1.3382 | 0.3274 |
| `FitRaw` | **6.4468** | 6.4432 | 0.0407 | **−1.4492** | −1.4590 | 0.3206 |

At 100 000 configurations (the doctest gates): m = +1 Q = −1.2424 (5σ 0.454),
m = 0 Q = +2.5661 (5σ 0.477), m = −1 Q = −1.2325 (5σ 0.454); the equal-thirds
mixture is 0 within 5σ; ⟨r²⟩ = 6.4599 (5σ 0.0579); δ⊥ = −0.1035.
Σ_i r⃗_i = 0 to **2.2e-15 fm** worst component over 20 000 configurations
(gate 1e-12).

### 4.1 The m_S conditioning of R̂ (T21) — the error no moment test can see

| ⟨P₂(cos θ_R)⟩, M = +1 | exact | grid | MC (2 × 10⁵) |
|---|---|---|---|
| m_S = −1 (pure \|Y₂²\|²) | −2/7 = −0.285714 | **−0.2857142** | **−0.2783 ± 0.0052** |
| m_S = 0 (pure \|Y₂¹\|²) | +1/7 = +0.142857 | **+0.1430894** | **+0.1424 ± 0.0102** |
| m_S = +1 (S-dominated) | — | — | **−0.03575 ± 0.00098** |

An m_S-**marginal** draw would give the same, m_S-**summed** ⟨P₂⟩ in *every*
branch: Σ_{m_S} P(m_S|+1)⟨P₂⟩_{m_S} = **−0.0370** for the committed tables
(the design's own earlier estimate of that number is −0.0355).  The m_S = +1
branch alone is −0.0351, which is close to it only because that branch carries
98.2 % of the weight; the other two branches sit **46σ** (m_S = −1:
0.2428/0.00522) and **17σ** (m_S = 0: 0.1779/0.01022) away from the marginal.
Those σ are the T21 `MESSAGE`'s own per-branch errors. The conditioning is
doing what it must.

### 4.2 The α core (T5)

| gate | design | measured |
|---|---|---|
| ⟨s²⟩ | 2.0774 ± 0.015 fm² | **2.0774** (analytic), MC 2.0740 at n = 25 000 (gate now 5σ = **0.0348**, not the design's typed 0.015 — which is 2.2σ) |
| χ²/ndf vs the **closed-form recentred** density | < 2 | **1.21** on 29 bins (10⁵ core entries) |
| χ²/ndf vs `he4.density` itself | recorded ≈ 170 on 40 bins / 10⁶ entries | **21.6** on 31 bins / 10⁵ entries — χ² scales with the entry count, so ≈ **216** per 10⁶ |
| Gaussian core vs the analytic Gaussian | < 1.5 | **1.28** on 25 bins |

The closed form ρ̃_s(q) = f̃(3q/4) f̃(q/4)³ was independently validated against a
4 × 10⁵-configuration NumPy Monte Carlo of the recentring (χ²/ndf = 0.95).

**With the hard core on** (`--min-nn-separation 0.9`, T5h): rejection
correlates the four s⃗_i, so the closed form ⟨s²⟩ = (3/4)⟨v²⟩ does not apply
and the constructor measures it instead — ⟨s²⟩ = **2.2691 fm²** (free 2.0774),
(G6) ⟨r²⟩ = **6.5725 fm²** against the **6.5859** sampled at n = 5 × 10⁴,
seed 7 (5σ = 0.0813, i.e. 0.8σ).  The free closed form's 6.4447 would have
missed that same set by **8.7σ**, which is why the analytic layer follows the
sampler here rather than the other way round.  `rho_alpha_recentred()` and
T5(b) remain independent-draw statements and are not gated at 0.9 fm.

---

## 5. Timing (design §8's "to be measured and written back")

| measurement | value |
|---|---|
| `sample_set(100000)` in-process (doctest T20) | **0.54 µs/config** |
| `lipolgen-configs --n 100000` sampling loop | **1.31 µs/config**, 0.131 s |
| `lipolgen-configs --n 100000` total wall (build + sample + write 37 MB + md5) | **1.06 s** |
| target | < 5 µs/config ✓ |

Constructor cost (tables, 9 × (R, c) CDFs, the recentred-core transform) is
≈ 0.25 s and is paid once.

---

## 6. Reachability and the dial

| quantity | design | measured |
|---|---|---|
| Q_charge at s = 0 (the floor, = +Q_d) | +0.270 fm² | **+0.26967 fm²** |
| Q_charge at s = 1 (`FitRescaled`) | −0.615 fm² | **−0.61545 fm²** |
| root for Q_charge = `LI6_QUADRUPOLE_FM2` | s = 0.4032 | **0.40320**, closes to **1e-9** |
| the √(target/model) rule's answer | −0.048 fm² | **−0.0477 fm²** at s = 0.3646 |
| `validate()` throws above +0.270 / below −0.615 | yes | **yes** |

---

## 7. Files

* `data/vmc/density/{he4.density,li6.density,README.md}` — fetched from the
  Wayback Machine exactly as `data/vmc/README.md` documents; both md5s match.
* `docs/open_items/run_2026-09-02/example_li6_m1_configs.dat` — **100
  configurations, m = +1, seed 1**, 37 655 B (well under the 200 kB ceiling,
  so the 100-configuration form was used), plus its `.meta.json` sidecar with
  `md5` and `git` filled in by `lipolgen-configs`.
