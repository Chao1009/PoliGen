# Phase D — the §7 number placeholders, filled

Design: `docs/open_items/run_2026-09-02/design_D_b1_li6.md` §7.
Code: `include/lipolgen/b1_nuclear.hpp`, `src/core/b1_nuclear.cpp`.
Tests: `tests/test_b1_nuclear.cpp`.
Author: Agent A of design D §8, 2026-09-03; re-measured the same day after
review (the y-grid default and the b₁ᵈ-support fix below both move columns).

> **READ THE GATE FIRST.** `phase_D_gate.md` records that the A = 2 validation
> gate **fails G3b** at the gate's default settings: a factor **2.27** below
> the digitized CDKS Fig. 4 peak at CDKS Eq. (21)'s δ-function (3.68 at
> Eq. (17)'s κ = 1 form), with the nucleon PDF the dominant remaining item
> (×1.67 by the CT18NLO proxy, which brings it to a factor 1.39 — inside the
> window). Under design §5.4's *Escalation* clause **none of the ⁶Li numbers
> below may be published**: they are recorded here so the phase is
> reproducible and so a regression is catchable, not because they are results.
> Every one of them also carries the mandatory ±100 % band.

All values: **x·b₁ per nucleon, ⁶Li, Q² = 2.5 GeV²**, default options
(`Li6ConvolutionB1()`: raw digitized CDKS theory-1 b₁ᵈ, ToyF2 F₁ with
**`r_sigma_lt`** — the kernel's own toy R, *not* the `r1998` the A = 2 gate
runs with, see `docs/USAGE.md` §2a — κ = 1, `norm_target` =
`VMC_N_ALPHA_D_LI6`, w_CG = 1/10, w_αd = 1). `LI6_B1_PER_NUCLEON` = 2/6 on the
three deuteron terms and 4/6 on the struck-α term.

> **Two review fixes of 2026-09-03 moved numbers in this file.**
> (a) **The y grid.** `LightConeDensities::Options` now defaults to
> 2400 / 2400 / 3200 rather than the design's 600 / 2400 / 800. The struck
> deuteron's z-width scales with M_α/M_d = 1.987, so on the coarser outer
> segments term **(2d)** came out **2.3 % high at x = 0.10** and 0.8 % high at
> x = 0.30 (checked against an independent 2-D Gauss–Legendre evaluation of
> Eq. (21) with no y table at all, and against a ×4 refinement, which agree).
> On the new default (2d) is within 1.2e−3 of the ×4 value and **(2α)/(2d)
> lands on the analytic 0.5064** instead of 0.489–0.497. T2 now has a ⁶Li
> clause that would catch a regression.
> (b) **The support of b₁ᵈ.** `LightConeDensities::convolve` cut the y
> integral at x/y = 1 while its own header said it must not; the digitized
> column runs to x = 1.59 and the deuteron's per-nucleon support genuinely
> exceeds 1, so terms (1) and (3) were losing real strength near the top of
> the window — 0 % below x = 0.8, 1.6 % at x = 0.9, **5.1 % at x = 0.95**
> (T5b). Every number below x = 0.8 in this file is unaffected; ∫b₁ dx and any
> x ≥ 0.9 row are not.

## The four terms

| x | (1) embedded d, S | (3) CG D-wave | (2d) α–d D-wave, struck d | (2α) α–d D-wave, struck α | **total** | `Li6B1(MillerB1)` | `Li6B1(CdksB1)` |
|---|---|---|---|---|---|---|---|
| 0.05 | +1.5686e−6 | +3.0979e−9 | +1.3838e−6 | +6.9528e−7 | **+3.6507e−6** | +3.6503e−4 | +8.9916e−7 |
| 0.10 | −4.6926e−6 | −9.2965e−9 | +1.2673e−6 | +6.3805e−7 | **−2.7966e−6** | +4.1222e−4 | −2.6566e−6 |
| 0.20 | −2.3664e−5 | −4.6736e−8 | +1.6669e−6 | +8.3886e−7 | **−2.1205e−5** | +2.6700e−4 | −1.3491e−5 |
| 0.30 | −4.4967e−5 | −8.8342e−8 | +2.2133e−6 | +1.1092e−6 | **−4.1732e−5** | +6.7167e−6 | −2.5964e−5 |
| 0.50 | +4.2448e−5 | +8.5642e−8 | +6.1632e−6 | +3.0422e−6 | **+5.1739e−5** | −1.8650e−4 | +2.2792e−5 |

Reading:

* **Term (3) is 0.197 % of term (1)** at every x — the CG algebra's
  0.1 × P_D^{αd} = 0.1935 %. It is kept because it is free and because it is
  the analytic bridge to `LI6_B1_RANK2_TRANSFER`, not because it matters.
* **Terms (2d)+(2α) are 7–133 % of term (1)** (+1.326 / −0.406 / −0.106 /
  −0.074 / +0.217 as a signed ratio at the five x) and have the **opposite
  sign** over the whole window, because the α–d S–D relative sign is opposite
  to the deuteron's below the S node (66 % of the density). They very nearly
  cancel term (1) at x ≈ 0.1 and dominate below. *This is the ⁶Li quadrupole
  puzzle showing up in b₁, and it is the quantitative reason the 100 % band is
  mandatory.*
* **The struck-α term (2α) is half of that**: (2α)/(2d) = **0.5024 / 0.5035 /
  0.5033 / 0.5011 / 0.4936 / 0.5029** at x = 0.05/0.10/0.20/0.30/0.50/0.70,
  against the analytic 2(M_d/M_α)² = **0.5064**. (On the design's coarser y
  grid this column read 0.489–0.497; the 1–2 % offset was the grid, not the
  physics — see the review note above.) A regression that silently dropped the
  term would move b₁ by ~30 % (T15).
* Against the **production** default `Li6B1(MillerB1)` the new model differs by
  two orders of magnitude and by sign, because Miller (HERMES-like) and CDKS
  (convolution) are different *camps* for b₁ᵈ. **The default does not change.**

## SD / DD split of the orbital term (x·b₁, (2d)+(2α) together)

| x | SD | DD |
|---|---|---|
| 0.10 | +1.7537e−6 | +1.5162e−7 |
| 0.30 | +3.0574e−6 | +2.6510e−7 |
| 0.50 | +8.4690e−6 | +7.3642e−7 |

SD dominates DD by 11.50–11.57× throughout, as it must: DD is O(W²) and SD is
O(UW).

## Close–Kumano integrals — REPORTED, not enforced

| model | ∫b₁ dx |
|---|---|
| `Li6ConvolutionB1`, [0.01, 1.2], 241 points | **+1.45880e−4** |
| `Li6B1(MillerB1)`, same grid | +9.12429e−4 |
| `Li6B1(CdksB1)`, same grid | +6.92980e−5 |
| `close_kumano_integral(true)` (digitized CDKS, its own grid) | +4.59200e−4 |
| `close_kumano_integral(false)` (Miller, its own grid) | +5.91474e−3 |
| the A = 2 kernel of this phase, [0.01, 1.59], at its default κ | +2.15370e−4 |
| the same at Eq. (17)'s κ = 1 | +1.07915e−4 |

None is zero; none is enforced (`sf.hpp`: "Reported, not enforced").

## Band rows — `--b1-band-scale 0 / 1 / 2` (x·b₁ at x = 0.30)

| scale | x·b₁ |
|---|---|
| 0 | 0 (exactly) |
| 1 | −4.1732e−5 |
| 2 | −8.3465e−5 (bit-for-bit 2× the scale-1 row) |

**Never quote a single row.** Every published number from this backend is
{0, 1, 2} × b₁.

## Term-(2) knob rows — `w_alpha_d_dwave` 0 / 1 / 2

| x | w = 0 | w = 1 | w = 2 |
|---|---|---|---|
| 0.10 | −4.7019e−6 | −2.7966e−6 | −8.9130e−7 |
| 0.30 | −4.5055e−5 | −4.1732e−5 | −3.8410e−5 |
| 0.50 | +4.2534e−5 | +5.1739e−5 | +6.0945e−5 |

The knob is linear by construction (it scales (2d) and (2α) together, which is
right: they are one physical effect). At x = 0.10 the orbital term is 70 % of
the total's magnitude, so w = 2 nearly zeroes b₁ there.

## `finite_q_delta` on / off (⁶Li, x·b₁)

| x | κ = 1 (the ⁶Li default) | κ = √(1+γ²) | ratio |
|---|---|---|---|
| 0.10 | −2.7966e−6 | −2.7700e−6 | 0.9905 |
| 0.30 | −4.1732e−5 | −4.1256e−5 | 0.9886 |
| 0.50 | +5.1739e−5 | +5.5531e−5 | 1.0733 |

Far smaller than on the A = 2 gate (where it is 1.57–1.69 at x ≥ 0.5), because
in ⁶Li the term κ multiplies is the *orbital* one, which is a small fraction of
the total, and because M_d and M_α are large so the k_z κ shift matters less.
**This 1 % is why the two objects default differently** (design A10): the gate
is quoted at CDKS Eq. (21), which is their exact definition and worth a factor
1.62 there, while `Li6ConvolutionOptions::finite_q_delta` stays `false` and
keeps f(y) a function of y alone — one cached table reused at every x, against
a 0.24 s rebuild per x.

## Densities

| quantity | struck d | struck α |
|---|---|---|
| `norm()`, renormalised | 0.8194810 (= `VMC_N_ALPHA_D_LI6`, by construction) | 0.8194810 |
| `norm()`, raw (`renormalize = false`) | **0.8176766** | **0.8179966** |
| identity N_αd × ⟨E⟩/M_struck | 0.817421 (3.1e−4 apart, T3's own tolerance) | — |
| `mean_y()` = ⟨z²⟩/⟨z⟩ | **0.9997871** | **0.9984638** |
| closed form 1 − ε/M − ⟨k²⟩/(2 M_recoil M) | 0.9974859 | — |
| `p_d()` (y-weighted) | **0.0193299** | — |
| `p_d_momentum()` | **0.0193516** | — |
| `VMC_P_D_LI6` | 0.019354933 | — |
| y_max | 1.99286 | 1.25120 |

⟨z⟩ ≈ 1 for both, **not** 1/3: the fair-share normalisation of design §1.6,
which is what makes the ⁶Li convolution literally CDKS Eq. (16) one level up.

`p_d_momentum()` is the table's own trapezoid D/(S+D) and sits **1.7e−4 below**
`VMC_P_D_LI6`, which is built from the file's *printed* norms — that is the
file's own quadrature spread (design §2.1), not an error. Likewise the table
trapezoid S+D = 0.8194274 against the printed 0.8194810 (6.5e−5); part of that
6.5e−5 is the 1.5e−5-per-power inconsistency between `HBARC_GEV_FM` = 0.19733
and `cluster.cpp`'s file-local 0.1973269804 (three powers → 4.6e−5), recorded
in `read_fdeut_k`'s header comment.

## T14 — the angular-average truncation, MEASURED (not asserted)

Ratio of the dropped P₂ and P₄ remainders to term (1). The design's first
revision asserted "≲ 0.2 %" from wrong reasoning; these are measurements.

| x | P₂ remainder / (1) | P₄ remainder / (1) |
|---|---|---|
| 0.05 | −2.3128e−5 | +5.4028e−7 |
| 0.10 | −8.7307e−6 | −9.6625e−6 |
| 0.20 | −1.7030e−6 | −1.0482e−6 |
| 0.30 | +1.4109e−5 | −4.7182e−7 |
| 0.50 | −5.1057e−5 | −1.3610e−5 |

**Both are below 6e−5 of term (1) everywhere in the window**, i.e. four orders
of magnitude inside the 100 % band — but the reason is *not* "P_D × 1/10": the
P₄ coefficient is 15× the isotropic one that term (3) keeps. It is small
because the P₂- and P₄-weighted densities integrate to zero, so they enter only
through the curvature of b₁ᵈ(x/z), and because they multiply **b₁ᵈ, not F₁ᵈ**.
The remainders are themselves cancellations and move by up to 2 % between the
default grid and a 4×-refined one, which is why T14 pins them at rtol 3e−2
rather than the design's 1e−6.

The other dropped piece — the S–D interference in the b₁ᵈ sector — is term
(2d)'s SD structure with coefficient exactly 1/3 of it (G0f, exact) and b₁ᵈ in
place of F₁ᵈ, i.e. suppressed by b₁ᵈ/F₁ᵈ ~ 1e−3 on top: ~3e−4 of term (2d).

## The ±5 % N_αd systematic (design §2.1) and the Q4 knob

Three tabulations of the α–d spectroscopic factor span 5 %: 0.819481 (2014
`li6_ad1.momentum`, **used**), 0.856 (2004 `li6.ad`), 0.863 (Wiringa et al.,
PRC 89 (2014) 024305 §III: 0.846 + 0.017). The spread is **not explained** —
entries 1 and 3 are the same year and the same Hamiltonian family — and is
carried as a stated systematic through `norm_target`:

| N_αd | x·b₁ at x = 0.30 | relative |
|---|---|---|
| 0.95 × 0.819481 | −3.9646e−5 | −5 % |
| 0.819481 (default) | −4.1732e−5 | — |
| 1.05 × 0.819481 | −4.3819e−5 | +5 % |

b₁ is exactly linear in `norm_target`, so the systematic is ±5 % flat.

**Q4, `use_spectroscopic_factor = false`** (the "all of ⁶Li is α+d" variant):
x·b₁(0.30) = **−5.0925e−5**, i.e. **×1/N_αd = 1.2203** on the whole answer.
Both readings are defensible; the default is the conservative one (the 18 % of
⁶Li that is not α+d gets b₁ = 0) and the difference is well inside the band.

## The default's decomposition (design §1.7 table), confirmed

| | default `Li6B1` | `Li6ConvolutionB1` |
|---|---|---|
| P_D^{αd} | 0.0867 (Hulthén) → weight 0.921947 = `LI6_B1_RANK2_TRANSFER` | 0.0193549 (VMC) → 0.982581 |
| smearing in z | none | VMC α–d density (⟨z⟩ ≈ 1, so ≈ no change) |
| N_αd | 1 (implicit) | 0.819481 |
| α–d orbital terms | absent | (2d) **and** (2α) |

T6 pins the no-smearing limit: with `LightConeDensities::delta_limit(p_d, 1)`
the model reduces to `b1_li6_from_deuteron(b1_d, 1 − 0.9 p_d, 2/6)` to 1e-12,
and at P_D = `P_D_LI6` the weight 1 − 0.9 P_D equals `LI6_B1_RANK2_TRANSFER` to
1e-4 — the analytic bridge between this backend and the constant the library
default uses today.

## Alternative deuteron camp, for the plot legend

With `deuteron_b1 = MillerB1` instead of the raw CDKS column (legal, but it
must be said out loud — they are different camps):

| x | x·b₁ |
|---|---|
| 0.10 | +3.6104e−4 |
| 0.30 | +6.2303e−6 |
| 0.50 | −1.4956e−4 |

(`deuteron_b1_x_max` defaults to 1.0 for an injected camp, which is right for
Miller's analytic form; the digitized CDKS column's own 1.59 is used only when
`deuteron_b1` is null.)
