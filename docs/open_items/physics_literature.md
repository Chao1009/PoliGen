I have enough across all seven items. Writing up the report.

# LiPolGen theory-input research report
**Date:** 2026-08-29 · Scope: `plans/04` #6, #8, #9, #10, #14, #15, #16, #18 and `plans/08` D1/D2 · No code changed.

**Headline:** Item 6 is *closed* — three independent sources agree and the code's sign is wrong. Item 1 is much better than "not classified anywhere": a complete arbitrary-spin basis has existed since 1989 and the spin-3/2 parton-model rank-2/rank-3 functions were written down explicitly in 2022. Item 5 is a convention mismatch, not a physics disagreement. Items 2 and 7 are genuine first-mover gaps. Item 4 has two publishable precedents I had not expected to find.

---

## 1. Spin-3/2 inclusive structure-function basis (rank-2 and rank-3)

**Verdict: a complete basis exists — `plans/04` #14 should be downgraded from "not classified anywhere we can adopt" to "classified, needs transcription into our frame."**

| Ref | What it gives |
|---|---|
| **R. L. Jaffe, A. Manohar, "Deep Inelastic Scattering from Arbitrary Spin Targets", Nucl. Phys. B 321 (1989) 343** (no arXiv — INSPIRE recid 266865) | **The complete basis, for arbitrary J.** Lepton scattering off polarized targets of arbitrary spin J in terms of Lorentz-invariant structure functions; in the Bjorken limit the target is described by **2J+1 quark distributions per flavour**, evolving logarithmically. Parton model *and* OPE analysis, plus a convolution model for nuclear targets. Only ~45 forward citations — this is why it reads as missing. **Not on arXiv; needs a library copy.** |
| **D. Fu, B.-D. Sun, Y. Dong, "Generalized parton distributions in spin-3/2 particles", [arXiv:2209.12161](https://arxiv.org/abs/2209.12161), PRD 106 (2022) 116012** | **The explicit J=3/2 rank-2 and rank-3 functions.** §3.1: *"in the parton model for the spin-3/2 sector, there are four independent structure functions in DIS at leading twist and leading order in α_s."* Eqs. (19a)–(19d) give them as quark densities q↑^λ(x) over the four helicity states, citing JM89: F₁ (rank-0), **b₁ = ½Σe_q²{[q^{3/2}+q^{−3/2}] − [q^{1/2}+q^{−1/2}]}/2** (rank-2), g₁ (rank-1), and a fourth they call **g₂ = ½Σe_q²{[q^{3/2}−q^{−3/2}] − 3[q^{1/2}−q^{−1/2}]}/√20** — *"the new structure function as the spin goes from 1 up to 3/2"*, i.e. **the rank-3 octupole partner of g₁**. Eq. (21): sum rules ∫b₁ = ∫g₂ = 0 for an untensor-polarized sea. **Warning: their "g₂" collides with the twist-3 nucleon g₂ — rename it in our code (e.g. `g1_rank3`).** |
| **D. Fu, Y. Dong, S. Kumano, J.-J. Xie, "Generalizing the Soffer Bound: Positivity Constraints on PDFs of Spin-3/2 Particles", [arXiv:2602.11587](https://arxiv.org/abs/2602.11587), PRD (2026)** | Leading-twist PDF classification: **6 quark** (f₁, g₁, h₁, f₁LL, g₁LLL, h₁LLT) + **5 gluon** (f₁, g₁, f₁LL, h₁TT, g₁LLL) PDFs; spin-3/2 density matrix Eq. (21) ρ = ¼(1 + ⅘SⁱΣⁱ + ⅔TⁱʲΣⁱʲ + 8/9 RⁱʲᵏΣⁱʲᵏ). **Complete positivity bounds** — free sanity checks for any scenario shapes we invent. |
| **J. Zhao, Z. Zhang, Z.-t. Liang, T. Liu, Y.-j. Zhou**, [arXiv:2206.11742](https://arxiv.org/abs/2206.11742) (PRD 106 (2022) 094006) and [arXiv:2401.10031](https://arxiv.org/abs/2401.10031) (PRD 109 (2024) 074017) | The covariant spin-3/2 polarization parameterization: S^μ, T^{μν}, R^{μνρ} orthogonal to P (Eq. 9), decomposed in Eqs. (11)–(13) into **S_L, S_T, S_LL, S_LT, S_TT, S_LLL, S_LLT, S_LTT, S_TTT** — the 3+5+7 = 15 parameters. Written for the *produced* hadron, but the algebra is identical for a target. This is the geometry `spin.py` needs. |
| **W. Cosyn, C. Weiss, [arXiv:2603.23699](https://arxiv.org/abs/2603.23699) — Appendix D** (local: `refs/2603.23699.pdf`) | **The on-ramp, already in our refs directory.** "Extension to higher-spin targets": a spin-j density matrix decomposes into l = 0…2j multipoles; spin-j adds a new l = 2j multipole with 2(2j)+1 parameters, a symmetric traceless rank-2j Lorentz tensor (**pseudo-tensor for half-integer j**). Counting Eq. (D2) ⌈(3²/2)(2j+1)²⌉ → SIDIS series **5, 18, 41, 72, 113** for j = 0, ½, 1, 3/2, 2, i.e. **72 structure functions for spin-3/2 SIDIS**. Rule (i) fixes the parity: **odd-l multipoles give cos-terms only with a polarized lepton, and sin-terms (T-odd) with an unpolarized one.** |

**The operationally decisive consequence.** Combining Cosyn–Weiss rule (i) with time-reversal at Born level: in **inclusive** DIS with an **unpolarized electron beam**, only **even-l** multipoles survive. The l=1 vector and l=3 octupole sectors would have to appear as sin(nφ_S) single-spin terms, which vanish in the one-photon-exchange approximation. Therefore:

- The generator's **rank-≤2 truncation is exact, not an approximation**, for A_zz and the cos 2φ_TT observables that carry every published figure. `plans/04` #14's "Default: rank ≤ 2 truncation" can be restated as a theorem with a citation rather than a caveat.
- The rank-2 geometry (T_LL = Q_NN·P₂(cosΘ), T_TT = (3/2)Q_NN·sin²Θ — Cosyn [arXiv:2410.12764](https://arxiv.org/abs/2410.12764) Eq. 9) is **spin-independent** and carries to J=3/2 verbatim, confirming the "isotope-generic kernel" note of 2026-08-28. What changes for ⁷Li is the *value* of the alignment and the structure functions, not the machinery.
- The rank-3 sector only ever reaches a **beam-helicity-dependent** observable. It is not on the critical path for any published LiPolGen figure.

**Recommended path.** What is genuinely absent is the finite-γ *inclusive* decomposition for J=3/2 — the analogue of Cosyn 2410.12764 Eqs. (10) and (17a–e), i.e. how many b-like functions exist at finite Q² for spin-3/2 and their kinematic ε/γ prefactors. A short note would: (i) write ρ^{3/2} in the covariant form of 2401.10031 Eqs. (9)–(13); (ii) contract with Cosyn's basis tensors A_i^{μν} following Cosyn–Weiss Appendix D's iteration; (iii) name and normalize the resulting rank-2 block against a stated alignment tensor and give the P_zz → `spin`'s T map; (iv) identify the rank-3 block and show it multiplies λ_e only. Expect the spin-1 set unchanged in form plus a rank-3 block. **Effort: 5–10 person-days** for someone comfortable with the algebra, and Cosyn/Weiss wrote Appendix D precisely as this on-ramp — this is the cleanest co-authorship ask in the whole register.

---

## 2. b₁ for A > 2, and ⁶Li in particular

**Verdict: no b₁ calculation exists for any A > 2 nucleus. Confirmed by title, fulltext and citation searches. LiPolGen would be implementing the first one — a genuine first-mover result, and simultaneously a result with nothing to validate against.**

| Ref | What it gives |
|---|---|
| **W. Cosyn, Y.-B. Dong, S. Kumano, M. Sargsian, [arXiv:1702.05337](https://arxiv.org/abs/1702.05337), PRD 95 (2017) 074036** (local) | **The kernel to build on.** Two complete convolution formalisms, every step written out. Eq. (16) b₁ = ∫(dy/y) δ_T f(y) F₁ᴺ(x/y); Eq. (17) light-cone density; Eqs. (19)–(20) the S+D decomposition with CG coefficients ⟨2m_L:1m_S\|1H⟩; **Eq. (21)** the closed-form tensor density δ_T f(y) ∝ ∫d³p y[−(3/4√2π)φ₀φ₂ + (3/16π)\|φ₂\|²](3cos²θ−1)δ(…); Eqs. (40)–(46) the light-front virtual-nucleon variant and the A_zz relation. Figures 4/5 already digitized in our repo. |
| **M. Hirai, S. Kumano, K. Saito, T. Watanabe, [arXiv:1008.1313](https://arxiv.org/abs/1008.1313), PRC 83 (2011) 035202** | **The closest A > 2 analogue that exists.** Not b₁, but the same group convolving nucleonic SFs with **cluster momentum distributions** in light nuclei (⁹Be) from AMD wave functions. The template for writing a cluster light-cone density and folding it into a DIS convolution. |
| **H. Khan, P. Hoodbhoy, PRC 44 (1991) 1219**; Khan & Hoodbhoy, PLB 298 (1993) 181 | The original convolution b₁, as an explicit few-parameter parametrization of F₂ᴰ, g₁ᴰ, b₁ᴰ in terms of a handful of wave-function moments — the form one would re-fit for an α–d system. **Different normalization from Cosyn; do not mix.** |
| **A. Yu. Umnikov, [arXiv:hep-ph/9605291](https://arxiv.org/abs/hep-ph/9605291), PLB 391 (1997) 177** | The health warning: non-relativistic convolution gets small-x wrong and violates the Close–Kumano sum rule (PRD 42 (1990) 2377). Sets the validity floor at **x ≳ 0.1**. |

**"Embedded deuteron": supported at the spin-bookkeeping level, unsupported at the DIS level.** The α+d cluster picture is standard (Kubo & Hirata, PTP 45 (1971) 1786), and the polarized-target community routinely writes ⁶LiD as "an unpolarized α bound to a polarized deuteron" (COMPASS ⁶LiD targets, NIM A 498 (2003) 101). But nobody has convolved a deuteron b₁ with an α–d density, for ⁶Li or anything else.

**The caveat that must go in the systematics.** Q(⁶Li) = −0.0806 fm² is *negative and 3.5× smaller in magnitude* than Q_d = +0.286 fm²: the α–d relative D-wave nearly cancels the free deuteron's intrinsic D-state contribution in the closest measured observable. Since b₁ is driven entirely by D-state admixture, a naive P_zz-scaled b₁(d) ignores a second D-wave source known to enter with the opposite sign. The two are different operators (charge quadrupole vs light-cone momentum alignment), so the cancellation need not carry over — but **treat b₁(⁶Li) as uncertain at the 100% level in both sign and magnitude, and expose the α–d D-wave term as a knob.**

**Implementation.** Recouple \|⁶Li; 1H⟩ in the α–d basis exactly as Cosyn Eq. (19) does for the deuteron, with L = α–d relative orbital angular momentum and S = deuteron spin. Three terms result: (1) *embedded deuteron* — Eq. (16) with F₁ᴺ → b₁^d, convolved with the α–d light-cone density f_{d/⁶Li}(z) from Eq. (17) (peaked near z = 1/3); (2) *α–d D-wave* — Eq. (21) unchanged in form with φ₀,φ₂ → φ₀^{αd}, φ₂^{αd} and F₁ᴺ → F₁^d; (3) *CG depolarization* of the deuteron inside the L=2 component, an analytic P₂-weighted Clebsch–Gordan sum. Inputs: an α–d relative wave function reproducing S_αd = 1.474 MeV and the ⁶Li rms radius; **Wiringa 2014 ([arXiv:1309.3794](https://arxiv.org/abs/1309.3794)) §II gives the α–d spectroscopic factor N_αd = 0.86 = 0.846 (S) + 0.017 (D), i.e. an α–d D-state of ≈ 2%** — use it to normalize. Validate by reproducing Cosyn Figs. 4/5 exactly with A=2 inputs before trusting A=6, and check ∫b₁ against Close–Kumano. **Effort: 10–15 person-days** for the defensible three-term version; **4–6 days** for a term-1-only v0 with a documented 100% band.

---

## 3. Tensor-sector radiative corrections

**Verdict: better than the register assumes. Closed-form formulas and released code both exist — POLRAD 2.0 carries the full spin-1 tensor sector. What does *not* exist is any RC treatment for the φ-dependent (cos φ_TL, cos 2φ_TT) observables, or anything for spin-3/2.**

| Ref | What it gives |
|---|---|
| **I. Akushevich, A. Ilyichev, N. Shumeiko, A. Soroko, A. Tolkachev, POLRAD 2.0, [arXiv:hep-ph/9706516](https://arxiv.org/abs/hep-ph/9706516), CPC 104 (1997) 201** | **The tensor sector is in there.** Born Eqs. (9)/(10): longitudinal target (F₁ − (Q_N/3)b₁), transverse (F₁ + (Q_N/6)b₁) — note the −1/3 : +1/6 ratio is exactly P₂(cosΘ), our own rank-2 geometry. Eq. (A.1)/(A.3): eight generalized SFs ℑ₁…ℑ₈, with ℑ₅–ℑ₈ the quadrupole ones in terms of b₁–b₄. **Eq. (A.4): the tensor *elastic* radiative tail in closed form via the quadrupole form factor F_q.** Eq. (33)/(34): universal L₁…L₅ kernels stated to apply for F = F₁,₂, g₁,₂, **b₁₋₄**. Eqs. (110)/(111): iterative b₁ unfolding. Code: subroutine `B14SF`, fit file `BB1FIT.DAT`, `ITQUAD.DAT`, quadrupolarization-degree input line. **Limits: fixed-target, form factors for d/³He/C/O only, no φ-dependent observables.** |
| **G. I. Gakh, O. Shekhovtsova, [arXiv:hep-ph/0403262](https://arxiv.org/abs/hep-ph/0403262), JETP 99 (2004) 898** | The **only** dedicated analytic tensor-DIS RC calculation — model-independent, unpolarized beam on a tensor-polarized deuteron (exactly our b₁ configuration), including the elastic tail; Appendix A tabulates all hard-photon coefficients; numerics for HERMES kinematics. **Headline: at x ~ 10⁻³–10⁻², RC on the tensor part is 10–30% of Born and shifts the b₁ zero-crossing to smaller x.** **Honest flag: zero citations on INSPIRE — never independently checked or used.** |
| **HERMES, [arXiv:hep-ex/0506018](https://arxiv.org/abs/hep-ex/0506018), PRL 95 (2005) 242001** | **What was actually done.** RADGEN-based MC with the *unpolarized* photon spectrum; joint detector+radiative unfolding; coherent and quasi-elastic tails from deuteron form factors; **polarized QE tail explicitly neglected** because "there is no net tensor effect by inclusive scattering on weakly-bound spin-1/2 objects" (Zhou et al., PRL 82 (1999) 687). Quantitatively: radiative background ≈ 50% of statistics in the lowest-x bin; its subtraction inflates stat+syst by ~2× at low x; **residual RC systematic on A_zz ≈ 2×10⁻³ at low x, negligible at high x** (against \|A_zz\| ≤ 0.02). |
| **JLab E12-13-011 (PR12-13-011, Slifer et al.); Poudel, Bacchetta, Chen, Santiesteban, [arXiv:2506.04506](https://arxiv.org/abs/2506.04506), EPJ A 61 (2025)** | The minimal-treatment precedent: *"no polarized radiative corrections at the lepton vertex, and the unpolarized corrections are known to better than 1.5%"*; **RC budgeted at 1.5% of a 9.2% total A_zz systematic** (polarimetry 8%, dilution 4%). |

**Recommended minimal treatment** (which is also what both existing experiments did): (1) collinear ISR/FSR with the *unpolarized* photon spectrum applied identically to F₁, F₂ and b₁–b₄ — the lepton-vertex correction is spin- and tensor-blind, and because A_zz is a ratio the bulk cancels, leaving only shape migration, which the (x,Q²) shift captures correctly; (2) **do the tensor elastic tail explicitly** — this is the one piece that does not cancel, since it is a background carrying genuine tensor dependence through the quadrupole form factor: use POLRAD Eq. (A.4) with ⁶Li C0/C2 form factors from **Wiringa & Schiavilla, [arXiv:nucl-th/9807037](https://arxiv.org/abs/nucl-th/9807037), PRL 81 (1998) 4317** (VMC AV18+UIX ⁶Li longitudinal and transverse form factors); (3) ~~neglect the polarized QE tail, citing Zhou et al. — the argument applies to ⁶Li's six nucleons exactly as it did to the deuteron's two~~ — **the transfer clause is withdrawn (2026-09-04, task B3): Zhou et al. is a *deuteron* measurement and is not in this repository, nothing in it is an A = 6 statement, and once the per-nucleon factor and POLRAD's Pauli suppression were fixed the quasi-elastic term became 22 % / 73 % / 99.9 % of the tail at x = 0.01 / 0.10 / 0.30, so “neglect” is now the largest unpriced piece of it rather than a rounding choice. The tail is still not computed — no polarized quasi-elastic radiative-tail calculation exists for an A = 6 spin-1 nucleus — but it is now PRICED, opt-in and off by default, by `RcOptions::qe_tensor_scale`, which lends it the elastic tail's own tensor fraction as a *borrowed magnitude and not a derived bound* (`include/lipolgen/rc.hpp`; `run_2026-09-03/phase_B_numbers.md` §B3)**; (4) **quote 1.5% on A_zz for x ≳ 0.05** (the E12-13-011 published precedent) and escalate to a **10–30% band on the tensor part for x ≲ 0.01** citing Gakh–Shekhovtsova, flagged in the docs as a single uncited calculation. This replaces `plans/04` #10's "no band" with a defensible one for the vector-like part while keeping the honest "no band" for the φ-dependent terms, which nobody has computed. **Effort: 5–8 person-days** (+5 if POLRAD's full ℑ-based inelastic-tail integrals are transcribed rather than a leading-log radiator).

Also relevant and easy to miss: **Cosyn & Sargsian, [arXiv:1407.1653](https://arxiv.org/abs/1407.1653), J.Phys.Conf.Ser. 543 (2014) 012006** compute A_zz *with FSI* — a different correction to the same observable, and the only paper that touches the tensor asymmetry with any non-Born effect.

---

## 4. FSI for tagged spectator fragments

**Verdict: the nucleon-spectator deuteron case is solved twice over; the *cluster*-spectator case exists for ³He→d (including polarized) but has never been done for α, t, or any A > 3 cluster; and as of March 2026 even the polarized *deuteron* tagged case is still impulse-approximation only.**

| Ref | What it gives |
|---|---|
| **W. Cosyn, M. Sargsian, [arXiv:1704.06117](https://arxiv.org/abs/1704.06117), Int. J. Mod. Phys. E 26 (2017) 1730004** | **The parameterized rescattering amplitude, with fitted values.** Eq. (11): f_XN(t,Q²,x) = **σ_tot(Q²,M_X)·(i + ε)·e^{B t/2}**. Fitted to JLab "Deeps": **σ_tot rises 30 → 70 mb over W = 1.2 → 2.4 GeV**, *smaller* at Q² = 2.8 than 1.8 GeV² (a color-transparency-like trend), **ε fixed at −0.5**, B fitted. Eq. (19): the distorted momentum distribution S^D(p_s) = IA + FSI loop. §4 Eqs. (22)–(23): **pole extrapolation and the loop theorem** — only the IA term has the pole at t → M_B², so multiplying by I(α_s, p_s⊥, t) recovers the on-shell structure function. Companion: [arXiv:1311.3550](https://arxiv.org/abs/1311.3550) (PRC 89, 014612), GEA FSI in *inclusive* deuteron DIS. |
| **M. Strikman, C. Weiss, [arXiv:1706.02244](https://arxiv.org/abs/1706.02244), PRC 97 (2018) 035209** | The complementary intermediate-x (0.1 < x < 0.5) model: FSI arise **predominantly from the spectator interacting with *slow* hadrons (rest-frame momenta ≲ 1 GeV) in the target-fragmentation region**, built from measured DIS hadron distributions and low-energy hadron scattering amplitudes. Studies the recoil-momentum and angular dependence, the analytic properties, and the effect on on-shell extrapolation. **Note the framework disagreement with Ciofi–Kopeliovich below on whether NN or πN dominates — span both in the band.** |
| **C. Ciofi degli Atti, L. P. Kaptari, [arXiv:1011.5960](https://arxiv.org/abs/1011.5960), PRC 83 (2011) 044602** — and the polarized version **L. P. Kaptari, A. Del Dotto, E. Pace, G. Salmè, S. Scopetta, [arXiv:1307.2848](https://arxiv.org/abs/1307.2848), PRC 89 (2014) 035206** | **The cluster-spectator precedent, which I did not expect to find.** The first computes **³He(e,e′d)X — a *bound deuteron cluster* spectator** — with Glauber FSI: Eq. (6) distorted momentum distribution, **Eq. (7) survival factor S_FSI = Π_{i≥2}[1 − θ(z_i−z₁)Γ(b₁−b_i, z₁−z_i)]**, **Eq. (8) profile Γ = [(1−iα)σ_eff/(4πb₀²)]exp(−b²/2b₀²)**. The second is the **polarized ³He with a detected deuteron**, giving a *spin-dependent distorted spectral function*. **Structurally: the cluster is bound in the wave function, but the debris rescatters off its individual constituent nucleons — the cluster survives only if no constituent is struck.** |
| **W. Cosyn, C. Weiss, [arXiv:2006.03033](https://arxiv.org/abs/2006.03033) (PRC 102 (2020) 065204) and [arXiv:2603.23700](https://arxiv.org/abs/2603.23700) (2026)** (latter local) | The polarized tagged framework, IA only. **2603.23700 §VI C is the honest statement of the gap**, verbatim: *"The present treatment … is limited to the IA and should be extended to include FSI. FSI has a large effect on the unpolarized tagged DIS differential cross section at spectator momenta ≳ 300 MeV … Interesting questions are: (a) Can a possible spin dependence of the rescattering enable interference of S- and D-wave? (b) What is the effect of FSI on tensor-polarized asymmetries…? (c) What is the size of T-odd structure functions, which are zero in the IA…?"* Our exact question, open in the newest paper in the field. |

**Minimal implementable model.** Apply a multiplicative weight on the sampled spectator momentum. Primary (citable) version — the per-nucleon Glauber product of Ciofi Eqs. (7)–(8), folded with our α+d overlap:

> W_FSI(p_spec) = |∫d³r e^{i p·r} φ_{αd}(r) · Π_{i∈spectator}[1 − Γ(**b**₁−**b**_i, z₁−z_i)]|² / |φ̃_{αd}(p)|², with Γ = [(1−iα)σ_eff/(4πb₀²)]e^{−b²/2b₀²}

Starting parameters: **σ_eff = 30–70 mb** rising with W (Cosyn–Sargsian Deeps fit) or the time-dependent σ_eff(t) = σ_NN + σ_πN·n_M(t) with σ_NN = 40 mb, σ_πN = 20–30 mb, Δt ≈ 1 fm rising to ~100–160 mb; **α = ε = −0.5** (Cosyn–Sargsian) or −0.35 (JLab-tuned); **b₀ = 0.35–0.6 fm**, which converts to the diffractive slope as **B = b₀²** → 6.4 GeV⁻² at b₀ = 0.5 fm, consistent with Cosyn–Sargsian's B ≈ 6 GeV⁻². A coherent-α-amplitude variant (σ_tot(Xα) ≈ 2.5–3 × σ_XN after Glauber shadowing — *not* 4×; B_α ≈ 31 GeV⁻² from p–⁴He elastic, stable across p_lab = 5–400 GeV/c) makes a good systematic cross-check but has no published precedent for DIS debris.

**Two structural points specific to LiPolGen.** (a) FSI matters above p_spec ≈ 300 MeV and is strongly angular — precisely where the tagged tensor asymmetries are O(1), so this is a *leading*, not a residual, systematic on the tagged channels. (b) **The pole-extrapolation escape route is unusually favourable here**: the α–d separation energy is **1.474 MeV** and α–t is **2.467 MeV** (both already in `coherent.py` from TUNL A=6/A=7), *comparable to or smaller than* the deuteron's 2.2 MeV — so the extrapolation length into the unphysical region is similarly short for the α-tagged ⁶Li channel. That is a genuinely strong, citable argument (Cosyn–Sargsian §4, loop theorem) and worth making in any tagged-channel writeup.

**Honest limits:** nothing constrains the *spin dependence* of the rescattering, so the weight should be applied spin-independently and quoted as an unpolarized-shape systematic, never as a correction to A_zz. And **BeAGLE has no FSI for light nuclei at all** — [arXiv:2204.11998](https://arxiv.org/abs/2204.11998), PRD 106 (2022) 012007, §II: *"No Final-State Interactions (FSI) are present in the BeAGLE generator for lepton-deuteron collisions"*; its FLUKA Fermi break-up has no memory of ⁶Li's α–d structure. **Effort: 5–8 person-days** for the weight plus parameter band; +3 for a pole-extrapolation diagnostic.

---

## 5. ⁶Li and ⁷Li effective nucleon polarizations

**Verdict: 1/3 vs 0.81 was never a physics disagreement — it is a convention mismatch, and both ends are wrong. The ab-initio answer is 0.85.**

Every source in this chain uses the **whole-nucleus sum**, Cloët–Bentz–Thomas Eq. (24): P_α = ⟨J,H| Σ_{i∈α} σ_z(i) |J,H⟩ = N_{↑α} − N_{↓α}, **with no division by Z or N**, feeding g_{1A} = P_p g_{1p} + P_n g_{1n}. (**³He is NOT in that convention.** Bissey et al. Eq. (2) reads, verbatim,
`g1He(x, Q2) = Pn g1n(x, Q2) + 2Pp g1p(x, Q2)` with `Pn = 0.86±0.02` and
`Pp = −0.028±0.004` — the proton term carries an explicit **2**, so P_p =
−0.028 is **per proton** and the whole-nucleus proton sum is 2 P_p = −0.056.
That is exactly how `HE3()` stores it, and `TRITON()` mirrors it. This
parenthetical said "same convention" until 2026-09-16 and
`docs/PHYSICS_CHANNELS.md` recorded the disagreement as "an unresolved
convention conflict"; re-fetched from the source PDF that day, it is this line
that was wrong, not the code. A reader who had "corrected" `HE3()`/`TRITON()`
to half their proton values would have flipped the sign of g₁(³He): on the toy
backends at x = 0.3, Q² = 5 the shipped code gives −2.7445e−03 and the
whole-nucleus reading gives +8.3171e−04.)

| Ref | What it gives |
|---|---|
| **R. B. Wiringa, R. Schiavilla, S. C. Pieper, J. Carlson, [arXiv:1309.3794](https://arxiv.org/abs/1309.3794), PRC 89 (2014) 024305** | **Table I is the provenance of 0.866/−0.037.** VMC AV18+UX, spin-up/down proton/neutron counts at M_J = J: ⁶Li 1.924/1.076 both species → **P_p = P_n = 0.848**; ⁷Li p 1.934/1.066 → **0.868**, n 1.981/2.019 → **−0.038**. JLab PR12-14-001 cites exactly this for "P_p = 0.866, P_n = −0.037" (extra digit from the higher-precision Argonne online tables). The coincidence with the Cohen–Kurath 13/15 = 0.8667 is just that — the same shell model gives P_n = +2/15 = +0.133, *opposite in sign* to ab initio. |
| **M. Piarulli, S. Pastore, R. B. Wiringa et al., [arXiv:2210.02421](https://arxiv.org/abs/2210.02421), PRC 107 (2023) 014314** | The uncertainty band. VMC with five chiral **Norfolk NV2+3** Hamiltonians: ⁶Li → 0.86; ⁷Li → 0.88 / −0.02. Caption: *"variation among the different interactions … is less than 0.01."* **Swapping AV18 for chiral moves P by ≤ 0.02.** |
| **B. S. Pudliner, V. R. Pandharipande, J. Carlson, S. C. Pieper, R. B. Wiringa, [arXiv:nucl-th/9705009](https://arxiv.org/abs/nucl-th/9705009), PRC 56 (1997) 1720** | §VIII, the original with the physics discussion: ⁶Li integrated neutron densities **1.93 up / 1.07 down, "net polarization of 29%"** (whole-nucleus 0.86, per-neutron 0.287); ⁷Li P(p↑) = 1.94, neutrons 1.98/2.02. Table XIV: μ(⁶Li) = 0.828(1) vs exp 0.822 (**0.7% — validates P**), but **Q(⁶Li) = −0.33(18) vs exp −0.083** (does *not* validate the tensor structure). |
| **I. C. Cloët, W. Bentz, A. W. Thomas, [arXiv:nucl-th/0605061](https://arxiv.org/abs/nucl-th/0605061), PLB 642 (2006) 210** (local) | Where the polarizations enter the polarized EMC ratio: Eq. (23) R_As, Eq. (24) the definition, Eq. (25) P scales linearly with M_J, **Eqs. (26)/(27) the K=1 reduced-matrix-element factor √[(2J+1)(2J+2)/6J] = 1.491 for ⁷Li, 1.414 for ⁶Li** — mandatory if working in the multipole basis and easy to drop. **Trap: their published ⁷Li *curves* used shell-model P_n = +2/15, not −0.037**; the QMC value appears only in prose. |

**Resolution of the 1/3 vs 0.81:**

| Statement | Whole-nucleus P | Per-nucleon | What it actually is |
|---|---|---|---|
| Slides' "1/3" | **1.00** | 0.333 | Naive α+d, no D-state, no cluster loss. **16% too big.** |
| Slides' "0.87 × 0.93 = 0.81" | **0.805** | 0.268 | Cluster estimate, ~5% low — and probably a mild double-count, since SLAC E155's 87% traces to nucleon-level Faddeev/GFMC numbers that *already* contain the deuteron D-state. |
| **Ab-initio VMC** | **0.85** | **0.283** | Everything included. |

**Recommendation (whole-nucleus convention, no Z or N factor):**

| | ⁶Li (J=1, M_J=+1) | ⁷Li (J=3/2, M_J=+3/2) |
|---|---|---|
| P_p | **+0.85 ± 0.03** | **+0.87 ± 0.03** |
| P_n | **+0.85 ± 0.03** | **−0.03 ± 0.03** |
| per-nucleon (÷Z, ÷N) | +0.283 / +0.283 | +0.290 / −0.0075 |

Band construction: MC + inter-model spread ≤ 0.01, AV18-vs-chiral ≈ 0.012–0.018 observed, plus ~0.02 for the **unquantified** VMC→GFMC shift (no GFMC-quality polarization numbers exist for either isotope — every published value is VMC). **⁷Li's 0.866/−0.037 is confirmed** — keep it. **For ⁶Li, replace 1/3 with 0.85** (or 0.283 if the slot is per-nucleon — state which, in a comment beside the constant). Two caveats for the code: μ(⁶Li) is reproduced to 0.7% so the *vector* polarization is solid, but Q(⁶Li) is off by 4× — **do not derive a ⁶Li tensor polarization or b₁ input from these wave functions**; and ⁷Li's μ is 10% low in impulse approximation from missing T=½ MEC, which affects the moment, not ⟨σ_z⟩. **Effort: 0.5–1 person-day** to adopt the constants with an M_J-scaling docstring and a test reproducing R_pol = 0.866 − 0.037 g₁ⁿ/g₁ᵖ; +0.5–1 day for the K=1 factor.

---

## 6. Tensor sign convention (`TENSOR_LL_SIGN`) — **decided**

**Verdict: `TENSOR_LL_SIGN` must flip from `+1` to `−1`. Three independent sources agree; there is no residual ambiguity.** This closes `plans/08` D1.

**(a) Cosyn, Roldan Tomei, Sosa, Zec, [arXiv:2410.12764](https://arxiv.org/abs/2410.12764), EPJ A 61 (2025) 83** (local `refs/2410.12764v1.pdf`). Footnote 1 states their A_T is "elsewhere denoted as A_zz". Chain: Eq. (6) Q = n₊ + n₋ − 2n₀ (= P_zz); Eq. (14a) T_LL = t_zz[RF]; Eq. (20a) T_LL = Q/3 for polarization along **N**_q; Eq. (19) A_T = (1/Q)·2·[T_LL(F_{UT_LL,T} + εF_{UT_LL,L}) + …]/(F_{UU,T} + εF_{UU,L}); Bjorken limit Eqs. (25)/(26) F_{UT_LL,T} = −2b₁, rest zero; Eq. (16) F_{UU,T} = 2F₁. Result — **Eq. (27):**

> **A_T = −(2/3) b₁/F₁**, *"the relation that was used in the b₁ extraction of the HERMES result."*

**(b) HERMES, [arXiv:hep-ex/0506018](https://arxiv.org/abs/hep-ex/0506018), PRL 95 (2005) 242001.** Eq. (2): d²σ_P = d²σ[1 − P_z P_B D A₁ + **½ P_zz A_zz**]. Eq. (6):

> **b₁ᵈ = −(3/2) A_zzᵈ F₁ᵈ**, with F₁ᵈ = (1+γ²)F₂ᵈ/[2x(1+R)]

— i.e. **A_zz = −(2/3) b₁/F₁**, identical. Their normalizations also match exactly: substituting T_LL = Q/3 into Cosyn Eq. (10) gives σ_P/σ_U = 1 + (Q/2)A_T, which is HERMES Eq. (2) term-for-term with Q = P_zz. **The measured data carry the same sign**: A_zz negative at low x while b₁ is positive there.

**(c) Hoodbhoy–Jaffe–Manohar, Nucl. Phys. B 312 (1989) 571** — the origin of b₁. HJM/HERMES define b₁ = ½(q⁰ − q¹) with q¹ = q^{+1} = q^{−1}. HJM does not itself write A_zz (an experimental construct), but the derivation is two lines and independent: a pure m=0 state has P_zz = −2, so σ⁰/σ_U = 1 − A_zz; parton model gives σ⁰/σ_U = 3q⁰/(2q¹+q⁰), so −A_zz = 2(q⁰−q¹)/(2q¹+q⁰) = 4b₁/(6F₁), i.e. **A_zz = −(2/3) b₁/F₁**. ✓

**(d) Independent fourth check, POLRAD 2.0** ([arXiv:hep-ph/9706516](https://arxiv.org/abs/hep-ph/9706516)), which I checked because its appendix at first appeared to disagree. Born Eq. (9), target polarized **longitudinally**: σ ∝ (F₁ − (Q_N/3)b₁)xy² + …; Eq. (10), **transverse**: σ ∝ (F₁ + (Q_N/6)b₁)xy² + …. With Q_N = P_zz, Eq. (9) gives σ/σ_U = 1 − (P_zz/3)b₁/F₁ = 1 + (P_zz/2)A_zz ⟹ **A_zz = −(2/3) b₁/F₁**. ✓ And the −1/3 : +1/6 ratio between the two is exactly P₂(cos 0°) : P₂(cos 90°) = 1 : −½ — **the same rank-2 geometry the generator already implements**, which is a bonus validation of `_tensor_moments`.

**What the program should adopt.**
- Set `asymmetries.TENSOR_LL_SIGN = -1.0` and update `evgen/tests/test_tensor_convention.py`, whose assertion `A_T(θ_S=0)·(1 + ε(y)R) == TENSOR_LL_SIGN·(2/3)b₁/F₁` is already written in the right F₂/R-free form and whose guard message names this decision. B2 already verified the flip is a one-line change with only that test failing.
- **The flip also flips κ, and therefore the O(γ²) subtraction of D2 — which is why D2 was correctly gated on this.** D2 can now proceed.
- Two convention traps to record beside the constant: **POLRAD Eqs. (110)/(111) use a *third* asymmetry, A_q ≡ b₁/F₁**, neither A_zz nor A_T; and the widely quoted form that divides by [1+2(1−y)/y²] *and* multiplies by [1+2(1−y)(1+R)/y²] double-counts R (those brackets are identically 1+εR) and misses by 1.17 — already noted in `plans/08` B2, worth keeping.

---

## 7. Coherent diffraction on polarized light nuclei

**Verdict: `plans/04` #18's "lightest published is Ca" is now out of date at the rate/acceptance level — ⁷Li coherent J/ψ has been simulated, and the Sartre table bottleneck was removed three months ago. But the physics claim stands: no coherent-diffraction calculation exists for any *polarized* nucleus beyond the deuteron, and none for ⁶Li or any α-cluster nucleus at all.**

| Ref | What it gives |
|---|---|
| **H. Mäntysaari, F. Salazar, B. Schenke, C. Shen, W. Zhao, [arXiv:2408.13213](https://arxiv.org/abs/2408.13213), PLB 858 (2024) 139053** (local) | **The template, already digitized as `coherent.MANTYSAARI_A2_DEUTERON`.** Eq. (9) dσ/dΦd\|t\| ∝ 1 + 2Σa_n e^{inΦ} (cos 2Φ coefficient = **2a₂**); Fig. 4 a₂, a₄ per m-state; polarizations defined in the γ*d c.m. frame. **The only polarized coherent calculation in existence.** |
| **H. Mäntysaari, B. Schenke, C. Shen, W. Zhao, [arXiv:2303.04866](https://arxiv.org/abs/2303.04866), PRL 131 (2023) 062301**, "Multi-scale Imaging of Nuclear Deformation at the EIC"; and **Mäntysaari, Salazar, Schenke, [arXiv:2207.03712](https://arxiv.org/abs/2207.03712)**, "Nuclear geometry at high energy from exclusive VM production" | The deformation-in-coherent-diffraction machinery (IP-Glasma/Good–Walker, coherent vs incoherent separation). 2207.03712 is Pb-only; 2303.04866 is the deformed-nucleus paper. Neither treats light nuclei. |
| **W. Chang et al., [arXiv:2511.05638](https://arxiv.org/abs/2511.05638), PRD 113 (2026) 032018** (local) | **Coherent J/ψ, φ, ρ on ²H, ³He, ⁴He, ⁷Li, ⁹Be, ¹²C, ¹⁶O at the EIC IR-8**, 10M events each, 0.1 < Q² < 100 GeV², through the EIC afterburner, with \|t\| distributions and far-forward detection efficiencies. **Uses eSTARlight** (github.com/eic/estarlight), *not* Sartre. Unpolarized, Woods-Saxon/Gaussian densities, no cluster structure. |
| **T. Toll, T. Ullrich, [arXiv:1211.3048](https://arxiv.org/abs/1211.3048), CPC 185 (2014) 1835 (Sartre)**; and **T. Toll, D. Ghosh, A. Srivastav, [arXiv:2606.14633](https://arxiv.org/abs/2606.14633) (June 2026)** | **The decisive answer on Sartre.** Sartre needs a precomputed amplitude lookup table per (nucleus, vector meson); historically this was the blocker — *"a few CPU-year for each combination"*. 2606.14633 improves table production by **3–4 orders of magnitude**, so a new table is now **"produced in a few hours."** Sartre is a C++ class library, so adding a nucleus means supplying its density and regenerating the table — **now genuinely feasible.** |
| **A. Mondal, A. Kumar, D. Sarkar, [arXiv:2608.23445](https://arxiv.org/abs/2608.23445) (Aug 2026)** | Shell structure in the coherent \|t\|-differential cross section from self-consistent QMC-model densities, **"enhancing the secondary diffractive lobes for light nuclei"** — the closest thing to "realistic light-nucleus density → coherent \|t\| structure." Mean-field, not clusters; unpolarized. |

**What does *not* exist, checked several ways:** no coherent diffraction, exclusive VM, or DVCS calculation for ⁶Li; nothing for α-cluster nuclei in coherent diffraction (the α-clustering EIC work — Magdy et al., [arXiv:2405.07844](https://arxiv.org/abs/2405.07844), EPJ A 60 (2024) — is BeAGLE forward-multiplicity observables on ⁹Be/¹²C/¹⁶O, not diffraction); and no tensor cos 2φ of a coherent yield for any A > 2. `plans/06` §6.4b's ε_B0 ∈ −(0.04–0.13) scaling remains an extrapolation with no external support, and the note that it cannot be rescaled to ⁷Li at all stands.

**Recommended path, in increasing cost.** (1) **Unpolarized exclusive-VM rate scale, ~1–2 person-days** (done 2026-09-02; 85–97 % of eSTARlight's rate sits at Q² < 0.1 GeV², below `coherent.hpp`'s q2_min = 0.7, so it is not a coherent-channel baseline — reworded 2026-09-23): add ⁶Li to eSTARlight — Z = 3 puts it in the Z ≤ 6 branch that already uses a Gaussian mass distribution, and ⁷Li is already supported, so this is a table entry plus a validation run against 2511.05638's ⁷Li numbers. This replaces the hand-tuned f₀ = 0.04 ×2÷2 with something citable. (2) **Sartre with a VMC ⁶Li density, ~10–15 person-days:** supply the Argonne VMC one-body density (phy.anl.gov/theory/research — note the site is behind Cloudflare and needs a real browser) and regenerate the amplitude table using the 2606.14633 method; gives a dipole-model \|t\| distribution with realistic light-nucleus structure, still unpolarized. (3) **The tensor cos 2φ — no shortcut.** The amplitude must be built from an **α+d configuration-space density with m-dependence**, Good–Walker-averaged over orientations, which is exactly the Mäntysaari–Schenke deuteron setup with the pn density replaced by α–d. That is the ask in `plans/04` #18 and it is correctly scoped as a collaboration, not an implementation: **~1–2 person-months of their code, days of ours to consume the output.** ⁶Li remains an attractive proposal because Q(⁶Li) = −0.0806 fm² makes it a near-null test against the deuteron's +0.286. *(Annotation, 2026-09-04 — that sentence is now priced, and the price is high but payable. `OPEN_ITEMS_SOLUTIONS.md` §11.3 / `phase_C_numbers.md` §C2: a near-null target gives a near-null signal, and near-null turns out to be expensive rather than unmeasurable. Feeding the measured Q through the closed-form quadrupole → a₂ map gives a₂(±1) = +0.026 at |t| = 0.3 GeV², but the coherent sample lives at |t| ≈ 1/B ≈ 0.026 GeV², so the information-weighted modulation is **0.25 %** — which on 4.2 × 10⁵ events is 2.6 σ. With eSTARlight's own **σ = 11.971 nb over the whole Q² range** (`estarlight_li6.md` §2f), 10 fb⁻¹/u, **both J/ψ lepton channels** and arXiv:2511.05638's 17.75 % efficiency that is **2.62 σ**, and 3 σ needs **13.1 fb⁻¹/u** — inside the {1, 10, 100} band. **Read that as a BAND and never as that point: 2.63 σ at the band's low edge and 2.84 … 3.29 σ at its top, 3 σ at 8.3 … 13.0 fb⁻¹/u** (`OPEN_ITEMS_SOLUTIONS.md` §11.3b–c): the 17.75 % is a ⁷Li number at ⁷Li's own top energy 18 × 118 applied at 10 × 99.5 — ×1.12–1.16 conservative on that paper's own ³He energy scan — while the chain has no decay-lepton acceptance or reconstruction efficiency in it at all, which can only cost. Those two nearly cancel. **Whether the band's top crosses 3 σ is *not* established** (§11.3c): the ⁷Li → ⁶Li efficiency substitution is undetermined in direction, ×0.99–1.33 read off the four species entries that share a beam energy, and that is the whole of the top. **Path (3) is therefore worth proposing as a J/ψ measurement, in photoproduction** (Q² < 0.1 is 85 % of the rate, and it is also where [Mant24]'s own Fig. 4 lives). The light mesons have more statistics — ρ⁰ 293 σ, φ 39 σ — but both are below LiPolGen's own M_X ≥ 1.2 GeV floor and are large-dipole channels the published deuteron calculation does not cover, so J/ψ being measurable is what keeps the theory ask the easy one. Four things are still **not** established and are stated in the draft: no detection efficiency exists below Q² = 0.1 anywhere in this tree; the 17.75 % is transferred across a beam-energy change (conservative, and sized: ×1.12–1.16) **and a species change whose direction is not established at all** (×0.99–1.33); no decay-lepton reconstruction efficiency exists anywhere in this tree; and the far-forward working point is unchosen — that last being the one correction that on its own restores the NO. **A first version of this annotation said 0.75 σ / 160 fb⁻¹/u and "not worth proposing"; that came from one lepton channel in one Q² window and is withdrawn — `OPEN_ITEMS_SOLUTIONS.md` §11.3a.**)*

---

## Summary

| # | Item | Verdict | Effort |
|---|---|---|---|
| 6 | Tensor sign | **Closed. `TENSOR_LL_SIGN = −1`**, four-way confirmed (Cosyn Eq. 27, HERMES Eq. 6, HJM derivation, POLRAD Eqs. 9/10). Unblocks D2. | **1 line** + test update |
| 5 | ⁶Li/⁷Li polarization | **Closed. 1/3 vs 0.81 was a convention mismatch; adopt 0.85 ± 0.03 whole-nucleus for ⁶Li**, keep 0.87/−0.03 for ⁷Li. | 0.5–1 d |
| 1 | Spin-3/2 basis | **Exists** (Jaffe–Manohar NPB 321; explicit J=3/2 functions in 2209.12161 Eqs. 19a–d). **Rank-≤2 truncation is exact for our unpolarized-beam observables**, not an approximation. Gap = finite-γ inclusive decomposition. | 5–10 d note |
| 3 | Tensor RC | **Formulas and code exist** (POLRAD tensor sector; Gakh–Shekhovtsova). Adopt ISR shift + tensor elastic tail; quote 1.5% (x ≳ 0.05), 10–30% band (x ≲ 0.01). | 5–8 d |
| 4 | Cluster-spectator FSI | Deuteron solved; **³He→d cluster spectator done, including polarized** (1307.2848); α/t spectators never done. Pole extrapolation favourable (S_αd = 1.474 MeV). Apply spin-independently. | 5–8 d |
| 2 | b₁ for A > 2 | **Nothing exists for any A > 2.** First-mover. Three-term α–d convolution off the Cosyn kernel; **100% band mandatory** (⁶Li quadrupole cancellation). | 10–15 d (4–6 d v0) |
| 7 | Coherent ⁶Li | Unpolarized ⁷Li now simulated (eSTARlight); **Sartre tables went CPU-year → hours** (2606.14633). Polarized/tensor case still nonexistent beyond the deuteron. | 1–2 d / 10–15 d / collab |

Two caveats on this report. First, item 7's Sartre and light-nuclei findings are mine directly from the sources; a parallel search on that item had not returned when I wrote up, so there may be additional deformed-nucleus cos 2φ references I have not seen — the conclusion that no ⁶Li or polarized-A>2 coherent calculation exists was cross-checked independently and I am confident in it. Second, Gakh–Shekhovtsova (item 3) has **zero INSPIRE citations**: its 10–30% number is a plausibility band from a single unchecked calculation, and should be documented as such rather than applied as a correction.

---

Both agents are in, and I verified their key claims directly. Here are the revisions to the report I gave above — including one correction to my own item 7, and one open question that turns out to be **closed**.

---

# Revisions

## NEW — `plans/04` #15 (VMC two-cluster overlaps) is **closed for both isotopes**, and the β band is biased low

> **CORRECTION, 2026-08-29 (`docs/open_items/vmc_reconciliation.md`).** The
> data findings in this section stand and the files are now in `data/vmc/` and
> wired into the code (`ClusterWaveSource::VmcAV18`).  **The comparison table
> below and the "the β band is biased low" conclusion are WITHDRAWN.** The
> β = 0.20/0.30/0.40 numbers quoted for ⁷Li (⟨k⟩ 0.1113/0.1333/0.1505,
> P(k>0.3) 0.0231/0.0522/0.0836) are the **S-wave Hulthen** density; ⁷Li → α+t
> is a pure **P wave** in the code, whose form is
> `k/((k²+κ²)(k²+β²))`, and whose real numbers are ⟨k⟩ 0.2250/0.2813/0.3279
> and P(k>0.3) 0.2133/0.3485/0.4581.  Against the right form, VMC α+t is
> **softer** than every β in the band, not harder.  The VMC tail fractions
> quoted here were also read at a grid-snapped 1.6 fm⁻¹ = 0.3157 GeV rather
> than 0.300 GeV (at a true 0.300 the values are 0.212 for ⁷Li and 0.059 for
> ⁶Li).  The ⟨k⟩ values, S_αd, S_αt and P_D = 1.94 % are all correct.

The Argonne QMC momentum-distribution files exist for exactly what LiPolGen needs, with S/D decomposition and Monte-Carlo error bars. The ANL site now 403s curl behind Cloudflare, but the **Wayback mirror serves the raw ASCII** (I downloaded and integrated both):

- **`momenta/li6_ad1.momentum`** — ⁶Li(1⁺) α–deuteron relative momentum distribution, VMC AV18+UX, 1M samples, K = 0–5 fm⁻¹ in 0.1 steps. Normalization 4π∫ρ_K K²dK/(2π)³ = **0.81971** (= S_αd). **It has a second block splitting S and D waves**: RHOKA0 with **0.80362** and RHOKA2 with **0.015861**, each with 1σ error bars → **α–d D-state fraction = 1.94%**. The S/D node at K ≈ 0.7 fm⁻¹ is clearly resolved.
- **`momenta/li7_at3.momentum`** — ⁷Li(3/2⁻) **ground state** α–triton distribution, VMC AV18+UX, 500k samples, dated 12-Apr-2024, S_αt = **1.0084**. (`li7_at1` is the 1/2⁻ excited state, 0.98683.) The P-wave character is explicit — ρ(K→0) = 0.0012 versus ⁶Li's 1041.5.

Retrieve with the Wayback availability API then the `id_` raw form, e.g.
`http://web.archive.org/web/20250606203105id_/https://www.phy.anl.gov/theory/research/momenta/li6_ad1.momentum`

**I integrated both and compared against the project's two-parameter band** (converged by K_max = 3 fm⁻¹; numbers stable to 4 digits):

| | ⟨k⟩ (GeV/c) | P(k > 0.3 GeV/c) |
|---|---|---|
| ⁶Li α–d, VMC | **0.122** | **0.045** |
| ⁷Li α–t, VMC | **0.186** | **0.179** |
| ⁷Li, `plans/04` #15 band β = 0.20 / 0.30 / 0.40 | 0.1113 / 0.1333 / 0.1505 | 0.0231 / 0.0522 / 0.0836 |

**The VMC answer lies outside the assumed band on both moments** — 24% above the β = 0.40 edge on ⟨k⟩, and **a factor 2.1 above it on the tail probability**. The band is not merely wide, it is biased low and does not bracket the truth, and it is low precisely in the k > 0.3 GeV/c region where `plans/04` #15 notes the tagged tensor asymmetries are O(1). This should be treated as a correction to existing tagged acceptances, not just a tightening. **Effort: ~1–2 person-days** to replace the radial forms with interpolated VMC tables and re-run the tagged studies. This also feeds item 2 directly — Cosyn Eq. (21)'s δ_T f(y) needs exactly the φ₀φ₂ and |φ₂|² structures, and the S/D split is the momentum-space form of them.

Caveat: these are VMC, not GFMC, and AV18+UX; no GFMC cluster momentum distributions are published. There is still **no tabulated coordinate-space ⟨⁴He+d|⁶Li⟩ overlap** anywhere — every published ANL overlap (Brida–Pieper–Wiringa, [arXiv:1106.3121](https://arxiv.org/abs/1106.3121), PRC 84, 024319) is single-nucleon. The momentum-space file is the substitute.

## Item 4 — revised

Corrections and additions to what I wrote:

- **Released code exists: `github.com/wcosyn/physics-code`** (Cosyn), with `DIS/DeuteronCross` and `Glauber/FastParticle` carrying a self-contained energy-dependent NN parameterization closed by β² = σ_tot²(1+ε²)/(16πσ_el). That is a far better starting point than transcribing formulas. Add **Cosyn & Sargsian [arXiv:1012.0293](https://arxiv.org/abs/1012.0293), PRC 84, 014601 (2011)** to the GEA citations.
- **Strikman–Weiss structure** is more usable than the abstract suggests: a distorted spectral function S[dist] = S[IA] + S[FSI] + S[FSI²] with a **factorized approximation** (§VI.C) shown numerically excellent; S[FSI] < 0 (absorption), S[FSI²] > 0 (refraction); unitarity sum rule ∫dΓ S[FSI+FSI²] = 0; slow-hadron parameters Σc_h = 0.6–0.8, ζ₀ ≈ 0.2, B_h = 6–8 GeV⁻². **FSI vanishes at the pole — pole extrapolation is FSI-proof**, which strengthens the α–d argument I made (S_αd = 1.474 MeV). The no-loop theorem is **Sargsian–Strikman, [arXiv:hep-ph/0511054](https://arxiv.org/abs/hep-ph/0511054), PLB 639, 223**.
- **Measured α–p cross sections, which I did not have:** **σ_tot(αp) = 121.5 ± 2.9 mb, σ_el = 31.4 ± 2.8 mb, B = 31 ± 1 GeV⁻²** at T_p = 620 MeV in the α rest frame — Blinov et al., [arXiv:nucl-ex/9910012](https://arxiv.org/abs/nucl-ex/9910012), Phys. At. Nucl. 64, 907 (2001) (5 GeV/c ⁴He beam on a hydrogen bubble chamber; I verified the paper and kinematics, not the numbers themselves). Bujak et al., PRD 23, 1895 (1981) gives b = 24 + 1.13 ln s → 31.5 GeV⁻², so **B_αN ≈ 31 GeV⁻² is energy-stable**. Glauber A_eff ≈ 2.9 of 4 — **do not use 4×**.
- **The physics point I had missed:** σ_el/σ_tot ≈ 0.26 means **~74% of α rescatterings destroy the tag**. So in the weight, use σ_tot in the linear (absorptive) term but **σ_el in the quadratic (gain) term** — the α-tagged distortion is nearly purely absorptive, and ∫w dΓ then equals the tagged-α survival probability rather than 1. That deviation is physical; log it as a diagnostic.
- **GEA's distinctive prediction: FSI peaks at α_spec = 1, not θ = 90°.** Effect is zero at the pole, largest at p_T ≈ 100–250 MeV (|t′| ≈ 0.02–0.15 GeV²), with **20–50% depletion by p_T ≈ 100 MeV**, roughly 3× the nucleon-spectator case.
- **Preferred primary implementation: the per-nucleon Glauber product** (Ciofi Eqs. 36–37 folded with the α+d overlap), because it is citable; keep the coherent-cluster amplitude as the systematic variant. Effort confirmed at **5–8 person-days** (2–3 for a Hankel-transform weight grid, +3–5 for W/Q²-dependent σ and derived cluster composition); ~1 day for a crude ±30% envelope.

Unchanged: no eA-DIS FSI treatment with a *coherent* cluster amplitude exists; the σ_tot/B table for d/t/α spectators is extrapolation, not measurement.

## Item 7 — **I was wrong about Sartre; correcting**

My earlier claim that the Sartre route is "now genuinely feasible" was based on the abstract of [arXiv:2606.14633](https://arxiv.org/abs/2606.14633). Verification against the actual v1.39 source overturns it:

- **Sartre has no user-defined nuclei.** `src/Nucleus.cpp::init()` is a hard-coded `switch (mA)` with `default:` → `exit(1)`. Supported: A = 208, 197, 110, 90, 63, 40, 27, 16, 2, 1. **No ³He, ⁴He, ¹²C, or Li.** Precomputed tables exist only for **A = 1, 40, 197, 208** — and `tables/README.txt` says the Pb tables "should not be used" for EIC. **That table list is the true origin of the project's "lightest published is Ca" note: it is a Sartre table-availability fact, not a theory-prediction fact.**
- The add-a-nucleus and table-generator manual chapters were never written. The Vogt radius formula gives R = 1.28 fm for ⁶Li versus r_ch ≈ 2.55 fm — off by 2× — so `Nucleus::rho()` needs its own branch. Nucleon sampling is spherically symmetric by construction, so **Sartre cannot carry a tensor-polarization axis at all.** Total ≈ 1 person-month plus a cluster allocation, for an unpolarized spherical ⁶Li. **Not worth it.** 2606.14633's 2000–7100× speedup is real but the code is **unreleased** (contact Tobias Toll).
- **Use `github.com/hejajama/subnucleondiffraction` instead** — the actual code behind [arXiv:2408.13213](https://arxiv.org/abs/2408.13213), which already has an A-generic `Nucleons` amplitude, a `VMC` wave-function branch, and `SetHeId(i)` Good–Walker averaging over an *external configuration table* (`he3.dat`). That is exactly the architecture a ⁶Li α+d sampler needs. Author asks to be contacted before use.

Two project claims are now **outdated** (`plans/06:328`, `plans/04:240`, `docs/note_cos2phi_coherent_6Li.md:326`, `docs/note_7li_theory_questions.md:92`):
- **[arXiv:2605.00454](https://arxiv.org/abs/2605.00454), PRD 114, 014068 (2026)** (Mäntysaari, Roch, Schenke, Shen, Zhao) publishes coherent *and* incoherent J/ψ down to **A = 3 and A = 4**, comparing VMC/NLEFT/PGCM configurations and explicitly assessing **α clustering** in ¹⁶O/²⁰Ne. The Ca floor is gone. Also **[arXiv:2202.12200](https://arxiv.org/abs/2202.12200)** (Guzey et al.), coherent J/ψ on ³He and ⁴He.
- **[arXiv:2511.05638](https://arxiv.org/abs/2511.05638)** (already local) Fig. 6 gives generated *and* detected coherent J/ψ **|t| distributions for e+⁷Li** with tagging efficiency vs |t| — a directly comparable curve for `coherent_optics_scan.py`. It is an acceptance study with no saturation physics, so the f₀ scenario band survives, but "no ⁷Li coherent number anywhere" does not.

Still true, and checked exhaustively: **2408.13213 remains the only polarization-dependent coherent-diffraction calculation for any nucleus.** All 8 INSPIRE forward citations were pulled; none extends it to A > 2 and none is polarized. No ⁶Li diffractive work exists at all.

**One physics warning worth adding to `plans/06`:** there are **two distinct cos 2φ mechanisms and they must not be conflated** — (i) the geometric/tensor one of 2408.13213 (density anisotropy about the polarization axis), and (ii) the **photon-polarization** one from linearly polarized quasi-real photons (STAR [arXiv:2204.01625](https://arxiv.org/abs/2204.01625), Sci. Adv. 9, eabq3903 (2023); [arXiv:2006.06206](https://arxiv.org/abs/2006.06206)). At Q² > 0 in eA the second is a genuine background to the tensor a₂ and needs an explicit control — e.g. the modulation measured w.r.t. the lepton plane versus w.r.t. the spin axis. *(Annotation, 2026-09-04 — the control is now worked out and it is **cheaper than this warning implies**, `OPEN_ITEMS_SOLUTIONS.md` §11.3 / `phase_C_numbers.md` §C2.5. Two independent handles: (i) the two modulations are about **different axes** — the spin axis, fixed in the lab, versus the photon's polarisation direction, uniform about the beam — so the second averages to zero in a spin-referenced histogram under a φ_γ-uniform acceptance; and (ii) the tensor term is **odd** in P_zz and the photon term **even**, so `tensor_flip_plan`'s two-fill difference removes it at first order and costs 5.4 % against the optimal use of the same fills while **gaining 1.50×** on a single-fill run. The residual — the two fills' φ-averaged |t| spectra differ only at second order in ε = 0.0067 — matters only for a background amplitude above 1.9 (J/ψ), 0.13 (φ), 0.017 (ρ⁰) on the whole-Q² samples (7.1 / 0.60 / 0.10 on the smaller 0.1 < Q² < 100 ones), and for ρ⁰ it makes |t| binning mandatory. **Both handles survive Q² → 0**, which is what lets the photoproduction region — 85 % of the coherent rate — be used at all: below Q² = 0.1 the electron is undetected and φ_γ is unknown event by event, but (i) is a statement about a *spin-axis* histogram and (ii) never references the lepton plane. Note also that "it dies as Q² rises" is **not** a property this note should be read as asserting: the virtual photon's linear-polarisation degree is a function of y, not a decreasing function of Q²; what Q² ≳ 0.7 GeV² buys is that the scattered electron is detected, so the photon-plane azimuth is known event by event — measurability, not suppression, and it is a control this argument never needed.)*

**Revised path:** ~1 person-day to add ⁶Li to eSTARlight (Z ≤ 6 already uses a Gaussian, so it runs today) for tagging/acceptance/rate work — do this first, it is nearly free and replaces a scenario number with a simulated one. Then **≈ 10–15 person-days** to build a ⁶Li configuration sampler (α core from ANL `he4.density`/TrNucGen; p–n pair from `deuteron/fdeut.av18` u(r), w(r) oriented by the polarization axis — *this is where the tensor signal enters*; α–d separation from the Fourier transform of `li6_ad1.momentum` weighted by S_αd) and graft it into `subnucleondiffraction`. **Or ~5 days** if the Mäntysaari group will run their existing polarized-deuteron setup with an α+d density — still the highest-leverage ask on the register.

## Small corrections to items 2 and 5

- **Better ⁶Li quadrupole number.** I quoted Pudliner's Q(⁶Li) = −0.33(18); the current value is **GFMC AV18+IL7 Q = −0.20(6) e·fm² versus experiment −0.082(2)** — Pastore, Pieper, Schiavilla, Wiringa, [arXiv:1212.3375](https://arxiv.org/abs/1212.3375), PRC 87, 035503. Still a factor ~2.5 discrepancy (the long-standing ⁶Li quadrupole puzzle), so the conclusion is unchanged and if anything sharper: **any b₁(⁶Li) or a₂(⁶Li) built on these wave functions inherits that as an irreducible systematic.** The 100% band on item 2 stands.
- **α–d spectroscopic factor:** use **S_αd = 0.820** (AV18+UX, from `li6_ad1` normalization) rather than the 0.86 I quoted from Wiringa 2014 — different interaction and a different quantity; do not conflate them. The D-state fraction 1.94% supersedes my "≈ 2%" estimate with an exact number.
- **⁶Li charge density:** de Vries–de Jager–de Vries (ADNDT 36, 495) has **no model parameters for ⁶Li** — only rms radii 2.54–2.57 fm. The UVa archive's `FB_data.dat` does carry a ⁶Li Fourier–Bessel row (R = 6.0 fm, a₁…a₇), provenance probably Li–Yearian–Sick (1971) — **verify against the primary source**, but if it holds it is the otherwise-missing analytic ⁶Li charge density.

---

**Net effect on the register.** Item 6 closes (`TENSOR_LL_SIGN = −1`, four-way confirmed). Item 5 closes (0.85 ± 0.03 for ⁶Li; the 1/3-vs-0.81 dispute was a convention mismatch). **#15 closes for both isotopes, and reveals the current β band is biased low by a factor ~2 in the tail that matters.** Item 1 downgrades from "unclassified" to "classified since 1989, needs transcription" — with the bonus that the rank-≤2 truncation is provably exact for every published observable. Items 2 and 7 remain genuine first-mover gaps, and item 7's Sartre route is now firmly ruled out in favour of `subnucleondiffraction`.