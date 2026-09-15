# The Mäntysaari-group ask — reconciled, single draft (C6)

> **THIS IS A DRAFT ON DISK. NOTHING IN THIS FILE HAS BEEN SENT, POSTED, OR
> COMMUNICATED TO ANYONE.** No email has gone out, no issue has been opened,
> nobody has been contacted. Whether to send anything below, in what form,
> and when, is **the author's decision alone** — this document exists so that
> decision can be made with the numbers in front of it, not so that it can be
> skipped.

## 0. Why this file exists

Two earlier drafts of a collaboration request to the Mäntysaari group (H.
Mäntysaari, F. Salazar, B. Schenke, C. Shen, W. Zhao and collaborators — the
authors of [arXiv:2408.13213](https://arxiv.org/abs/2408.13213), PLB 858
(2024) 139053, the only published polarized coherent-diffraction calculation,
and of `github.com/hejajama/subnucleondiffraction`, the code behind it) were
committed independently: `docs/OPEN_ITEMS_SOLUTIONS.md` §11.2 and
`docs/open_items/run_2026-09-02/design_G_cluster_config.md` §10. They
disagreed on one number (whether the α+d model's ⁶Li point-radius agreement
is quoted as "3 %" or "4 %" — resolved in §3 below: both are correct, against
different reference radii, and the letter now says so instead of picking
one silently) and, more importantly, they were drafted **before** open item
O5 was answered. **This file is the single reconciled draft. Both committed
copies now point here instead of restating the ask; their own text is left
in place with the editor's notes already on it, as the record of what
changed and when — this file is not a rewrite of history, it is what
replaces the letter itself.**

The reconciliation was informed by task C2 (open item O5), completed
2026-09-04: **the answer is MARGINAL for the channel the drafted letter
actually asks about — a band whose low edge is 2.6 σ in one EIC year and
whose top is 2.8–3.3 σ, with 3 σ at 8.3–13.0 fb⁻¹/u.** See §6.

> **This file was written twice on 2026-09-04, and the two versions reached
> opposite conclusions.** The first said the answer was **NO**, that the
> measurement "does not exist at any luminosity the EIC is quoted at", and
> recommended not asking for the calculation. That was built on a J/ψ rate
> restricted to **one lepton channel** and to **0.1 < Q² < 100 GeV²**, a
> window `estarlight_li6.md` §1 records as *"arXiv:2511.05638's range
> verbatim"* — an acceptance study's kinematic range, not a physics one, and
> one that omits 85 % of the coherent J/ψ rate. Correcting both multiplies
> the rate by 13.5 and the significance by 3.5. **The NO is withdrawn**;
> `OPEN_ITEMS_SOLUTIONS.md` §11.3a is the record of what moved and by how
> much. Nothing was sent between the two versions, because nothing has been
> sent at all.

## 1. What this repository can supply

* **An α+d configuration sampler** (`include/lipolgen/cluster_config.hpp`,
  `src/core/cluster_config.cpp`, CLI `lipolgen-configs`) that emits
  nucleon-position tables in exactly the layout
  `Nucleons::InitializeTarget` reads for ³He in `subnucleondiffraction`: six
  positions per configuration, in fm, ion rest frame, centre of mass at the
  origin (Σr⃗ᵢ = 0 to 2.2 × 10⁻¹⁵ fm worst component over 2 × 10⁴
  configurations), one configuration per line, with the polarization axis
  and the magnetic substate m ∈ {+1, 0, −1} (plus an unpolarized mix)
  applied and recorded in a sidecar alongside every input table's md5 and
  every sampler option. `docs/USAGE.md` §9 has the CLI, the Python API and
  the byte-level output format.
* **The physics behind that sampler**: the α core from the ANL VMC ⁴He
  one-body point-nucleon density (`he4.density`), the α–d separation and its
  L = 0/2 orientation from the ANL VMC α–d overlap (`li6.ad` /
  `li6.adr.fit`), and the p–n pair inside the deuteron from AV18 u(r), w(r)
  (`fdeut.av18`) — all three correlated to one drawn deuteron spin
  projection m_S per configuration, reproducing the exact m-state moments
  ⟨P₂(cos θ_R)⟩ = −2/7 (m_S = −1) / +1/7 (m_S = 0) to 1e-6/2.3e-4.
* **An unpolarized coherent-rate baseline that IS citable**: eSTARlight
  (`github.com/eic/estarlight`, the same code arXiv:2511.05638 cites) run at
  this repository's own e 10 GeV × ⁶Li 99.5 GeV/u point with no code change
  (Z = 3 falls into the generic light-nucleus Gaussian branch). Cross
  sections at R_G = 1.2·A^{1/3} = 2.1805 fm, **over the whole Q² range** —
  J/ψ **11.971 nb**, φ 654.3 nb, ρ⁰ 17.8 μb — reproduce the analytic Gaussian
  slope to 3.7–5.2 % and the ⁷Li/⁶Li ratio to within 2 % of A^{4/3}, i.e.
  this baseline has passed its own internal cross-checks and is now in code
  as `estarlight_li6_coherent()` and `estarlight_li6_q2_floors()` rather than
  living only in prose. (Restricted to 0.1 < Q² < 100 GeV² the same runs give
  1.773 / 30.16 / 506.4 nb. **That window is arXiv:2511.05638's
  acceptance-study kinematic range, not a physics one**; 85 % of the coherent
  J/ψ rate is below Q² = 0.1 and the |t| slope is the same there to 0.4 %.) It is unpolarized and has no diffractive
  minimum at any A — a rate-and-slope baseline, not an imaging one — but it
  is a genuine external calculation, not a scenario number.
* **A closed-form tensor estimate** built from the sampler's own point-matter
  quadrupole, `a2_from_quadrupole` (`cluster_config.hpp`), which reproduces
  every digitized row of [arXiv:2408.13213]'s deuteron a₂(m = ±1) to 8 %
  (m = 0: 8–21 %) with **zero free parameters** — a genuine, if narrow,
  validation of the map, not a tuned fit.

## 2. What this repository cannot supply, stated plainly

* **`a2_from_quadrupole` is a closed form, not a dipole-model amplitude, and
  is NOT a Good–Walker amplitude.** It carries the target's point-matter
  quadrupole moment through [Mant24]'s own published |t| dependence for the
  deuteron and stops there: no Good–Walker average over configurations, no
  actual scattering amplitude, no saturation model, none of
  [Mant24]'s own statistical or theoretical uncertainties, and the *matter*
  quadrupole stands in for the transverse **gluon** anisotropy the real
  observable is sensitive to. Every number this repository quotes about ⁶Li's
  coherent a₂ inherits this and only this.
* **The α+d truncation overshoots the measured quadrupole by a factor ≈ 7.5**,
  and that factor is now decomposed rather than asserted (open item O1,
  `phase_C_numbers.md` §C1): it is **two** factors, not three —
  **3.3165** (the model → an η-matched dial, anchored on the *measured*
  asymptotic α–d D/S ratio η = −0.025 ± 0.006 ± 0.010, George & Knutson, PRC
  59, 598 (1999)) × **2.2686** (that dial → the measured quadrupole) =
  **7.524**, identically. The apparent 15 % "missing component" some earlier
  notes counted as a third factor is **not** one — 1/S_αd = 1.171 is already
  inside the model's −0.615 fm², because both waves are divided by √S_αd
  before any moment is taken. **The band, not the central value, is the
  physics**: GK's own 1σ on η maps onto model Q from −0.4005 fm² to
  **+0.0298 fm², through a sign change at 0.86σ** — so the 3.32× leg alone
  spans 1.54× to a flipped sign, and η cannot distinguish the measured
  quadrupole (+0.48σ) from GFMC AV18+IL7's −0.20(6) fm² (Pastore *et al.*,
  PRC 87, 035503 (2013), −0.07σ). **The rule this repository holds itself
  to: no tensor number is quoted from these wave functions without this
  band attached.**
* **The α core is an uncorrelated product of one-body densities**, not
  genuine correlated GFMC ⁴He configurations (open item O3). What can be
  bounded from data already in this tree, and what cannot, is in §5.
* **No Good–Walker loop, no `subnucleondiffraction` integration, no dipole
  amplitude of any kind exists in this repository.** The graft this letter
  used to ask for has not been attempted here in any form.

## 3. The one number the two earlier drafts disagreed on

Both drafts' closing paragraph said the α+d model "reproduces the ⁶Li point
radius to N %, but overshoots Q(⁶Li) by a factor ≈ 7.5" — `OPEN_ITEMS_SOLUTIONS.md`
said **3 %**, `design_G_cluster_config.md` said **4 %**, and the first draft
was introduced as quoting the second "verbatim", which it did not. Both
numbers are real, and the disagreement was only ever about which reference
radius: the sampler's r_rms = 2.5386 fm sits **3.0 %** above the measured
point radius `LI6_R2_POINT_FM2` = 2.4655 fm (from the Angeli–Marinova charge
radius, corrected for nucleon size) and **3.9 %** above the *ab initio* VMC
`li6.density` value `LI6_R_POINT_VMC_FM` = 2.4433 fm — the design's own T6
gate quotes exactly this second number ("the model is 4 % too large" against
`li6.density`, ratio 1.039 ± 0.005). **The letter below states both,
attributed, rather than picking one.**

## 4. THE DECIDING NUMBER (open item O5, task C2)

Coherent J/ψ off tensor-polarised ⁶Li, at this repository's own programme
luminosity `Scenario::lumi_fb_per_nucleon` = 10 fb⁻¹/u ("one EIC year" — hence
10/6 = 1.6667 fb⁻¹ of e+⁶Li, since the luminosity is quoted per nucleon), with
**σ_coh = 11.971 nb over the whole Q² range** (eSTARlight, §1), **J/ψ → e⁺e⁻
and μ⁺μ⁻ (branching 0.11932, eSTARlight's own `JpsiBree` + `JpsiBrmumu`)**,
arXiv:2511.05638's global 17.75 % coherent-J/ψ efficiency (a ⁷Li number, used
as the only stand-in that exists), and `tensor_flip_plan(0.6)` evaluated on
the **background-immune** two-fill difference (⟨P_zz²⟩ = 0.81, the estimator
the separation argument below actually requires, rather than the optimal 0.90):

> **N = 4.22 × 10⁵ reconstructed events → δa₂(|t| = 0.3 GeV²) = 0.0100
> against a predicted a₂ = +0.0263. S = 2.62 σ. Three sigma needs
> 13 fb⁻¹/u — inside the {1, 10, 100} fb⁻¹/u band this repository uses, and
> close to one EIC year.** With a perfect detector: 6.21 σ. On the measured
> charge radius R = 2.589 fm instead of eSTARlight's default 2.18 fm — the
> other defensible density, and the leading rate systematic — 1.58 σ and
> 36 fb⁻¹/u, still inside the band.

**QUOTE THAT AS A BAND, NEVER AS THAT POINT.** The 17.75 % above is a ⁷Li
number at ⁷Li's own *top* energy (18 × 117.9 GeV/u) used on a ⁶Li sample at
10 × 99.5, which the same paper's ³He energy scan makes **×1.12–1.16
conservative**; and the chain carries **no decay-lepton acceptance or
reconstruction efficiency at all**, which can only cost. Correcting both:

> **S = 2.63 σ at the band's low edge, 3 σ at 13.0 fb⁻¹/u — inside the band,
> and open below** (the per-lepton reconstruction efficiency is unbounded in
> this repository). The 2.62 σ point lands 0.3 % under that low edge, because
> the two omissions cancel to 0.7 %.

**AND THE BAND'S TOP IS A SPAN, NOT AN EDGE.** The same 17.75 % is a **⁷Li**
number used for **⁶Li**, and the size *and direction* of that substitution are
**not established**. Read off the four entries of arXiv:2511.05638's species
list that share a beam energy (²D, ⁴He, ¹²C, ¹⁶O, all at 137 GeV/u, bracketing
⁶Li's own fixed-rigidity 137.5), the substitution is **×0.99–1.33 — it
straddles 1**. So:

> **the band's TOP is 2.84 … 3.29 σ, 3 σ at 8.3 … 11.1 fb⁻¹/u**, and taken
> with the low edge **3 σ sits at 8.3 … 13.0 fb⁻¹/u on every form — inside
> the {1, 10, 100} band at both ends**. **Whether the top crosses 3 σ is NOT
> established.**

*(An earlier version of this file quoted the top as a single 3.15 σ /
9.1 fb⁻¹/u, read off the same species list treated as a pure species lever
because every entry sits at the same magnetic rigidity. That is a non
sequitur — at fixed rigidity E/u and Z both still vary, by ×1.55 and 1 → 8 —
and it is withdrawn; `OPEN_ITEMS_SOLUTIONS.md` §11.3c records it.)*

**The 0.25 % modulation is real, and it is not a verdict.** a₂ is *linear* in
|t|, quoted at |t| = 0.3 GeV², but the coherent sample lives at e^{−B|t|} with
B = 39–55 GeV⁻², so the information-weighted modulation is κ√⟨t²⟩ = **0.25 %**
at B = 50 — a factor **10.6** below the a₂(0.3) = 2.6 % everyone quotes, and
that factor follows from the coherent form factor alone. But a modulation is
not a significance: on the sample that would actually be taken (B = 38.8,
ā₂ = κ√⟨t²⟩ = 0.32 %) the reach is ā₂√(2⟨P_zz²⟩N) = 0.0032 × √(2 × 0.81 ×
4.22 × 10⁵) = **2.6 σ**. The |t| window is *not* where anything is lost
either: `COHERENT_T_MAX_DEFAULT` = 0.2 GeV² keeps 99.995 % of the rate and
99.72 % of the Fisher information at B = 50.

**This document's first version got this wrong, in the direction of
pessimism.** It quoted **0.749 σ and 160 fb⁻¹/u** and concluded NO. That row
used one lepton channel and the 0.1 < Q² < 100 GeV² window. One correction at
a time, at 10 fb⁻¹/u:

| | S | 3 σ at |
|---|---|---|
| Q² > 0.1, e⁺e⁻ only, on the optimal ⟨P_zz²⟩ = 0.90 — *the first version* | 0.749 σ | 160 fb⁻¹/u |
| + both lepton channels | 1.005 σ | 89 fb⁻¹/u |
| + a deliberately conservative Q² > 0.01 floor instead | 1.371 σ | 48 fb⁻¹/u |
| + both conservatisms, and the measured-radius density | 0.828 σ | 131 fb⁻¹/u |
| no Q² floor, both leptons, measured radius | 1.576 σ | 36 fb⁻¹/u |
| **no Q² floor, both leptons, default density** | **2.618 σ** | **13 fb⁻¹/u** |
| … and with the de-squeezed far-forward tagging optics on top | **0.731 σ** | **168 fb⁻¹/u** |
| corrected ε_det, BAND LOW (beam-energy leg × lepton pair stand-in) | **2.63 σ** | **13.0 fb⁻¹/u** |
| corrected ε_det, BAND TOP, ⁷Li → ⁶Li read low (×0.99) | **2.84 σ** | **11.1 fb⁻¹/u** |
| corrected ε_det, BAND TOP, ⁷Li → ⁶Li read high (×1.33) | **3.29 σ** | **8.3 fb⁻¹/u** |

Every conservatism on that ladder **except one** leaves 3 σ inside the band on
its own; it takes two of them stacked to put it back out. **The exception is
the optics row**, which does it alone — and it is a *choice* rather than a
measurement. That is what makes this **MARGINAL** rather than YES — and it is
why "does not exist at any luminosity the EIC is quoted at" was not a
defensible sentence.

> *(An earlier version of this file said "every conservatism on its own leaves
> 3 σ inside the band" **and left the optics row out of the table entirely**,
> so the counterexample was invisible here while `validation/o5_a2_reach.py`
> printed the same sentence three lines below its own bullet stating it.
> `OPEN_ITEMS_SOLUTIONS.md` §11.3b records the repair.)*

**And it is this repository's own factor-7.5 quadrupole overshoot (§2) that
sets the scale.** The identical J/ψ sample gives 19.7 σ fed with this α+d
geometry's own Q = −0.6154 fm², 6.40 σ fed with GFMC's −0.20, and **2.62 σ**
fed with the measured −0.0818 — the last of these being **the uncorrected
ε_det = 0.1775 chain's middle, kept so that the comparison is between
quadrupoles only, and not quotable on its own: the band is 2.63 σ at its low
edge and 2.84 … 3.29 σ at its top, 3 σ at 8.3 … 13.0 fb⁻¹/u** (the two boxes of §4 above). ⁶Li is attractive as a near-null test against
the deuteron's own +0.286 fm² — and this is the exact price of a near-null
target: near-null is expensive. What it is not, on these numbers, is
unmeasurable.

**Separation from the photon-polarisation cos 2φ background is free, was never
the obstacle, and survives Q² → 0.** Two independent handles make it so.
(a) The tensor modulation references the spin axis (fixed in the lab); the
photon-polarisation modulation references the photon's linear-polarisation
direction (uniform about the beam), so it averages to zero in a
spin-referenced histogram. (b) The tensor term is *odd* in P_zz — exactly,
because a₂(0) = −2 a₂(±1) makes Σ_m p_m a₂(m) = P_zz a₂(±1) identically — and
the photon term is *even*, so `tensor_flip_plan`'s two-fill difference is
background-free at first order and self-normalises against relative-luminosity
offsets. **Below Q² = 0.1 the scattered electron is not detected and φ_γ is
unknown event by event; neither handle needs it** — (a) is a statement about a
spin-axis histogram under a φ_γ-uniform acceptance, and (b) never references
the lepton plane. That is what makes the photoproduction region usable. The
difference estimator's cost is *negative*: 5.4 % worse in δ than the optimal
use of the same two fills (and the headline above is quoted on the difference,
not the optimum), but **1.50× better** than putting the whole luminosity into
a single +0.6 fill, because the m = 0-rich fill carries |P_zz| = 1.2. The
residual (the two fills' φ-averaged |t| spectra differ only at second order
in ε = |eps_b0_equivalent()| = 0.006728, i.e. 1.13 × 10⁻³ at |t| = 0.2 and
2.26 × 10⁻⁵ averaged) matches √(2/N) only above a background amplitude A_γ of
1.9 (J/ψ), 0.13 (φ), 0.017 (ρ⁰) on these larger samples — harmless for J/ψ by
a wide margin; for ρ⁰, |t|-binning is mandatory, not optional, because no
magnitude for A_γ exists anywhere in this tree to rule that out.

**The light mesons have more statistics and it no longer matters.** φ: 39 σ
with no Q² floor, 8.3 σ over 0.1 < Q² < 100, 2.1 σ inside this repository's
own Q² > 0.7 window. ρ⁰: 293 σ / 49 σ / 8.6 σ. Both sit **below**
`COHERENT_MX_MIN_DEFAULT` = 1.2 GeV — the shipped generator carries no rate
there at all — ρ⁰ is the channel where the photon-polarisation cos 2φ has
actually been *observed* (STAR's ultraperipheral measurement,
arXiv:2204.01625 — this tree records the mechanism and the reference, not the
channel attribution or a magnitude, so treat that attribution as this
document's, not as sourced), and both are large-dipole, strongly-absorbed
amplitudes at ⟨W⟩ = 14–18 GeV where the closed-form map is **least**
defensible. **J/ψ being measurable is what makes all of that irrelevant to the
ask**: the channel we want is the channel the published calculation already
covers.

**THE LIMITATION RIDES WITH THIS NUMBER.** `a2_from_quadrupole` is a closed
form and not a Good–Walker dipole-model amplitude (§2). A genuine dipole-model
run could move the signal by a factor: it would have to move it **up by 1.15**
to reach 3 σ at 10 fb⁻¹/u, **down by 1.31** to fall back to 2 σ, or up by 1.91
to reach 5 σ. *(The first version of this paragraph said 4.0 / 6.7 / ×16 in
luminosity; those belonged to the restricted row.)* One thing the wider Q²
window *improves*: the digitised deuteron a₂ the map is validated against is a
**photoproduction** calculation ([Mant24] Fig. 4), so Q² < 0.1 is where the
map is anchored and 0.1 < Q² < 100 was the extrapolation, not the reverse.

**FIVE THINGS ARE NOT ESTABLISHED, and they are not details.**

1. **No detection efficiency exists below Q² = 0.1 GeV² anywhere in this
   tree.** arXiv:2511.05638's 17.75 % was measured on a sample generated with
   0.1 < Q² < 100 GeV² (its p. 3), so no point of it was ever evaluated in the
   region that now carries 85 % of the rate. Applying it there is an
   extrapolation. Two things bound it and neither is a measurement: the
   efficiency is a *recoil-nucleus* far-forward tagging acceptance ("tagging
   efficiency × acceptance", their Figs. 2 and 5) and does **not** require the
   scattered electron, which at Q² < 0.1 goes down the beam pipe; and the
   recoil |t| spectrum is unchanged to 0.4 % between the two windows, with ⟨W⟩
   moving in the direction their own Fig. 2 says *raises* the efficiency. If
   the efficiency collapsed to zero below Q² = 0.1 the headline would return
   to 1.005 σ and 89 fb⁻¹/u.
2. **The efficiency is a ⁷Li number at ⁷Li's own TOP energy, used at
   10 × 99.5 — and that is CONSERVATIVE, by a stated amount.** 17.75 % was
   measured at 18 × 117.9 GeV/u, where this repository's own
   configuration-identical eSTARlight run has ⟨W⟩ = 43.2 GeV; the sample
   priced above is at 10 × 99.5 with ⟨W⟩ = 30.2. arXiv:2511.05638 measures
   that dependence in its §V.B on ³He — **32.23 % at 18 × 183, 54.38 % at
   10 × 100, 99.77 % at 5 × 41** — i.e. d ln ε/d ln E_ion = −0.87 and −0.68,
   so the step this transfer needs is worth **×1.12–1.16 UP**. We state it
   because it is the direction that flatters us and we would rather you heard
   it from us. It is *not* their 183 → 100 ratio of 1.687: that is a step
   3.6× larger in ln E than the one we need. *(An earlier version of this
   bullet added "×1.16–1.22 for ⁷Li → ⁶Li at fixed rigidity". That is
   withdrawn — see item 3.)*
3. **And it is a ⁷Li number used for ⁶Li, which we cannot even sign.** The
   17.75 % is the ⁷Li entry of arXiv:2511.05638's p. 4 species list. Every
   entry of that list is at the same magnetic rigidity, A/Z × E = 275 GeV/e —
   but that fixes neither the per-nucleon energy (E/u = R·Z/A, which runs 118
   to 183 GeV/u down the list) nor the charge (p_T,cut ∝ Z·R, Z = 1 to 8), so
   its A-ordering is **not** a species lever on its own, and the ×1.55 spread
   in E/u alone is worth ×1.35–1.46 in ε on the scan of item 2. Read off the
   **four entries that do share a beam energy** — ²D, ⁴He, ¹²C, ¹⁶O at
   137 GeV/u, which bracket ⁶Li's own fixed-rigidity 137.5 — the ⁷Li → ⁶Li
   substitution is **×0.99–1.33: it straddles 1**. That is the whole of the
   width of the band's top, and it is why we do not claim the top reaches 3 σ.
4. **There is no decay-lepton acceptance or reconstruction efficiency in the
   chain at all, and we cannot supply one.** The 17.75 % is the fraction of
   scattered **nuclei** inside the far-forward acceptance, under a simulation
   which by your co-authors' own words *"only accounts for the acceptance
   effect and does not incorporate the efficiencies of the detector.
   Additionally, we did not account for the efficiency and acceptance of the
   reconstructed distribution"*. The **geometry** we can bound, and it is not
   the problem: at 10 × 99.5 a J/ψ at ⟨W⟩ = 30.2 sits at y = −0.39, nearly at
   rest in the lab, so both 1.548 GeV decay leptons are inside |η| < 3.5 with
   acceptance **0.994**, and never below 0.89 anywhere the beams can reach.
   The per-lepton reconstruction efficiency is what we have nothing for — and
   on 2026-09-15 we surveyed our own detector-simulation repository
   (`../PolarizedLithiumSim`: all of `fastsim/polli_fastsim/`, `tools/fullsim/`,
   `tools/analysis/`, `plans/03`, `plans/09` and all 54 `refs/` entries) and it
   has none either, only a per-track stand-in it labels as one and a
   constructed scattered-electron ID profile, so this stays **unbounded rather
   than bounded** (`run_2026-09-06/phase_C_numbers.md` §C1). **It
   can only cost**, and it is why the number above is a band that is open
   below: a pair acceptance × efficiency of 0.52 would cost the "above 2 σ"
   half of MARGINAL, and 0.12 would put 3 σ outside the {1, 10, 100} band.
   *(Items 2 and 4 pull opposite ways and nearly cancel — ×1.12 against
   ×0.90 — which is why we quote neither on its own.)*
5. **The far-forward working point is unchosen.** The number above assumes the
   Yellow Report optics, `Optics::lumi_fraction` = 1. LiPolGen's own
   de-squeezed ⁶Li tagging point delivers 0.0781 of the machine luminosity at
   10 × 100, and `estarlight_li6.md` §2e instructs to multiply by it. This
   estimate does not, on the stated ground that the two are alternatives
   rather than multipliers — the efficiency above is arXiv:2511.05638's
   **IR-8 secondary-focus** number, and that paper's own case for IR-8 (its
   p. 4) is that the secondary focus buys the low-p_T acceptance *without* the
   β*_x de-squeeze that costs luminosity in IR-6. **If that is wrong, every
   significance here multiplies by √0.0781 = 0.279: the whole band becomes
   0.73–0.92 σ with 3 σ at 106–167 fb⁻¹/u, back outside the band at both
   ends.** This is the single largest unmade choice in the estimate, and
   **the only correction on the ladder above that restores a NO on its own.**

Full derivation, every assumption listed and priced, and the self-checking
script: `docs/open_items/run_2026-09-03/phase_C_numbers.md` §C2;
`docs/open_items/run_2026-09-02/estarlight_li6.md` §2f (the Q² scan);
`validation/o5_a2_reach.py`; pinned in `python/tests/test_o5_reach.py` and
`tests/test_coherent.cpp` T10a / T10c.

## 5. Open item O3 — what is bounded offline, and what genuinely needs the configurations

O3 asks: the α core is sampled as an uncorrelated product of one-body
densities; how much would the incoherent/coherent split move if genuine
correlated GFMC ⁴He configurations (as arXiv:2605.00454 uses) were used
instead? No GFMC configuration table for any nucleus is in this repository,
so the split itself — a property of a Good–Walker amplitude's
configuration-to-configuration fluctuation — cannot be computed here at all.
**That part of O3 stays open and genuinely needs the configurations
requested in §6(a) below.**

What *can* be bounded with data already in the tree: `data/vmc/
bonus_other_clusters/he4.dd` is a genuinely correlated VMC overlap of the
same ⁴He wavefunction onto an α → d+d cluster channel (same method as
`li6.ad`/`li7.at`, Forest *et al.*, PRC 54, 646 (1996)). Comparing its own
relative-motion second moment to what an **uncorrelated** product of
`he4.density` would give for the identical observable — the separation
between the centroids of two 2-nucleon halves of the α, worked out exactly
in `validation/o3_alpha_correlation_bound.py` (⟨D²⟩ = (4/3)R₁², independent
of recentring and of which 2-2 partition is taken) — measures:

| quantity | value |
|---|---|
| r_rms(d–d separation), correlated (`he4.dd`) | **1.809 fm** |
| r_rms(D), uncorrelated product of `he4.density` | **1.664 fm** |
| ratio, variance | **1.182 (+18.2 %)** |
| ratio, rms | **1.087 (+8.7 %)** |

**Reading**: genuine 4-body correlation inflates this one position-space
moment by a several-to-twenty-percent effect, not a factor of several — small
next to the factor-7.5 quadrupole gap (§2) and comparable in size to the
wave-function-choice systematic already measured elsewhere in this run
(6.6 % on the tagged tensor dilution, `phase_C_numbers.md` §C5.2). **This
does not answer O3** — it is a different observable (α → d+d, not α → 4 free
nucleons; a position-space moment, not the coherent/incoherent split of a
dipole amplitude) standing in for the general size of correlation effects in
this nucleus, on the only genuinely correlated data this repository holds.
What genuinely needs the actual configurations: any statement about how much
the incoherent cross section itself — which depends on the full
configuration-to-configuration fluctuation of a Good–Walker amplitude, not on
a single second moment — would move. That is exactly request (a) below.

## 6. The reconciled ask

Both earlier drafts asked, as their request (c), for the Mäntysaari group to
run their existing polarized-J/ψ machinery on this repository's ⁶Li
configuration tables. **§4 says that request is worth making**: coherent J/ψ
is a 2.6 σ measurement at the band's low edge and a 2.8–3.3 σ one at its top,
with 3 σ at 8.3–13.0 fb⁻¹/u, it is the
only channel inside our own M_X floor, and it is the only channel their
published calculation already covers — so it is also the *cheapest* version of
the ask. What §4 changes is the **kinematics** of it: the measurement lives in
photoproduction, Q² < 0.1 GeV², not in the 0.1 < Q² < 100 electroproduction
window the first version of this file priced it in — which is also where their
own Fig. 4 lives.

> *(An intermediate version of this file said "§4 shows that is a request to
> spend person-months measuring a channel that is already known, in-tree, to
> be blind at any EIC luminosity" and rewrote (c) into a question about
> whether to bother. That is withdrawn: the blindness was an artefact of one
> lepton channel and one Q² window. (c) is a request for a calculation
> again — with the four things that are not established (no detection efficiency below Q² = 0.1, where 85 % of the coherent rate sits; no decay-lepton reconstruction efficiency anywhere in this tree, which is what leaves the band OPEN BELOW; the ⁷Li → ⁶Li efficiency substitution, which straddles 1 and makes the band's top a span; and the far-forward working point, unchosen — `cluster_config.hpp:525-534`) stated inside
> it.)*

**What is offered, in one line**: the configuration sampler and its output
format (§1), the eSTARlight rate/slope baseline including the Q² scan (§1), and
— if (b) is accepted — a patch, in return for (a) the GFMC configurations and
(c) the polarized ⁶Li calculation.

> We have built an α+d configuration sampler for polarized ⁶Li that emits
> nucleon-position tables in exactly the format `Nucleons::InitializeTarget`
> reads for ³He (positions in fm, c.m. at the origin, one configuration per
> line), with the polarization axis and the substate m = +1, 0, −1 applied
> and recorded. The α core comes from your ANL VMC/GFMC ⁴He one-body
> density, the p–n pair from AV18 u(r), w(r) with the full m_S angular
> correlation (your Eqs. (7)–(8), which we reproduce exactly), and the α–d
> separation from the ANL VMC α–d overlap with the L = 2 orientation
> correlated to m through the CG recoupling.
>
> Three things, in increasing size.
>
> **(a)** Would you share the ⁴He GFMC nucleon configurations used in
> arXiv:2605.00454? Our α core is currently an uncorrelated product of
> one-body densities, which is the weakest part of the sampler and the one
> part you already have solved. We can bound the general size of what
> correlations do to a position-space moment using your own published α → d+d
> overlap against our own one-body density (it comes out to a
> several-to-twenty-percent effect, not a factor of several), but we cannot
> compute how the incoherent/coherent split itself would move without the
> actual configurations.
>
> **(b)** Would you accept a ~30-line generalization of the `A == 3` branch
> to an `-configfile/-configid` option for any A? We can send the patch; it
> makes the code A-agnostic and removes the 13698-configuration bound.
>
> **(c)** Would you run your polarized-deuteron J/ψ machinery on our ⁶Li
> tables — **as photoproduction, Q² < 0.1 GeV²**? We have priced the
> measurement as carefully as we can from outside a dipole model, and we would
> rather show you the number and its holes than ask you to take it on faith.
> Feeding our sampler's own point-matter quadrupole through a closed-form
> quadrupole → a₂ map that reproduces your published deuteron a₂(m = ±1) to
> 8 % (with zero free parameters — it is not a fit) predicts ⁶Li's tensor a₂
> at the *measured* quadrupole is about 10× smaller than the deuteron's and
> of the opposite sign. Combining that with an eSTARlight ⁶Li coherent rate
> (11.97 nb at 10 × 99.5 GeV/u with no Q² cut), both J/ψ lepton channels, and
> arXiv:2511.05638's 17.75 % far-forward tagging efficiency, one EIC year at
> 10 fb⁻¹/nucleon gives a band whose low edge is **2.6 σ** and whose top is
> **2.8–3.3 σ**, with 3 σ at roughly **8–13 fb⁻¹/nucleon**. On the measured
> charge radius rather than eSTARlight's default it is 1.6 σ and
> 36 fb⁻¹/nucleon. So it looks like a
> real measurement rather than a hopeless one — but a marginal one, quoted as
> a band on purpose, and four of the largest factors in it are things we
> cannot settle:
>
> * **the detection efficiency below Q² = 0.1 GeV².** That is 85 % of the
>   coherent rate and it is where the measurement would be made, but
>   arXiv:2511.05638 generated 0.1 < Q² < 100, so their 17.75 % was never
>   evaluated there. We are extrapolating it. The recoil |t| spectrum is the
>   same to 0.4 % between the two windows and the efficiency is a
>   recoil-nucleus tagging acceptance rather than an electron-arm one, which
>   is why we think the extrapolation is defensible — but it is not measured.
> * **the decay leptons, which that 17.75 % does not cover.** It is the
>   fraction of scattered *nuclei* in the far-forward acceptance — the paper
>   says in as many words that it carries no detector efficiency and no
>   reconstructed-distribution acceptance — so our chain has no central-
>   detector acceptance or reconstruction efficiency for the e⁺e⁻/μ⁺μ⁻ pair at
>   all. The geometry we can bound and it is benign (the J/ψ is nearly at rest
>   at 10 × 99.5, so 0.99 of pairs are inside |η| < 3.5); the reconstruction
>   efficiency we simply do not have — we searched our own detector-simulation
>   repository for one on 2026-09-15 and it has none either, only a per-track
>   stand-in it labels as one — and it can only cost. That is the open
>   bottom of the band above.
>   Pulling the other way, and we would rather tell you than have you find it:
>   the 17.75 % is a ⁷Li number at ⁷Li's own *top* energy, 18 × 118, and we
>   use it at 10 × 99.5 — your co-authors' own ³He scan (32.23 / 54.38 /
>   99.77 % at 18 × 183, 10 × 100, 5 × 41) makes that ×1.12–1.16 conservative.
>   **That and the missing lepton factor nearly cancel**, which is why the
>   point estimate barely moved when we found them and why the honest form is
>   a band.
> * **the ⁷Li → ⁶Li substitution itself, which we cannot sign.** The 17.75 %
>   is the ⁷Li entry of your co-authors' species list, and we need ⁶Li. Every
>   entry of that list is at the same magnetic rigidity, but that fixes
>   neither the per-nucleon energy (118 to 183 GeV/u across the list) nor the
>   charge (Z = 1 to 8), so the A-ordering in it is not a species lever on its
>   own. Read off the four entries that *do* share a beam energy — ²D, ⁴He,
>   ¹²C, ¹⁶O at 137 GeV/u, which bracket ⁶Li's own fixed-rigidity 137.5 — the
>   ⁷Li → ⁶Li substitution is ×0.99–1.33: it straddles 1, and we cannot even
>   give you its sign. That is the whole of the width of the band's top, and
>   it is why we do not claim the top reaches 3 σ.
> * **the far-forward working point.** Ours costs luminosity (a de-squeezed
>   β*_x for spectator tagging, 0.078 of the machine luminosity at 10 × 100);
>   the IR-8 secondary focus that the 17.75 % comes from does not. We have
>   assumed the latter. If the polarized ⁶Li run has to share the former,
>   everything above multiplies by 0.28 and the measurement goes away. This is
>   the one item on our list that does that on its own.
>
> Everything on our side is a closed form rather than an amplitude — no
> Good–Walker average, no saturation, and the *matter* quadrupole standing in
> for the *gluon* anisotropy — which is exactly the gap your code closes. Two
> honest numbers about our model while you decide: it reproduces ⁶Li's point
> radius to 3 % against the measured value (2.4655 fm) and 4 % against your
> own VMC density (2.4433 fm), but overshoots the *measured* quadrupole by a
> factor ≈ 7.5 — two factors, not one (3.32× from the model to an η-matched
> dial anchored on the measured asymptotic D/S ratio, × 2.27× from there to
> the measurement), with a 1σ band on the first factor that runs from 1.54× to
> a sign change in Q. If that overshoot is the more interesting problem than
> the measurement itself, we would rather hear that from you first.

## 7. Recommendation

**The letter is now a request for a calculation again, and the calculation is
worth requesting.** Request (a) is a low-cost, channel-independent ask that
stands on its own merit (it is also the only route to closing O3, §5);
request (b) is general infrastructure with no physics content; and request (c)
has gone back to asking for the polarized J/ψ run — in photoproduction, with
the two unestablished factors stated inside the ask rather than discovered by
the recipient.

**What changed three times.** The first version of this file asked for (c) as
a calculation. The second withdrew it, on a 0.75 σ estimate. The third
restored it, on 2.6 σ, because the 0.75 σ came from one lepton channel in one
Q² window. **The fourth (2026-09-04) turns that point into a band, 2.6–3.2 σ**,
because the efficiency chain under it had two defects of comparable size and
opposite sign — ε_det applied at ⁷Li's *top* beam energy, and no decay-lepton
acceptance at all — and neither was in the assumption list that claimed to be
complete. The lesson worth carrying into the letter — and it is why the letter
states its own holes — is that **the estimate was wrong in the direction of
the conclusion its author had already reached**, three times over: the Q²
window was inherited without being questioned; the μ⁺μ⁻ gain was computed,
written down and then not applied; and the one efficiency in the chain was
carried across a beam-energy change and a species change without either being
named, while the factor that could only cost was never introduced at all. The
fourth pass is the first that moved the number in *both* directions.

Whether to send the letter at all, to send only (a) and (b), or to wait for
further work, remains **the author's call** — this document's job was to make
sure that call is made with the deciding number in front of it, and with an
honest account of how that number moved.

## 8. Provenance

* Open item O5 (task C2): `docs/open_items/run_2026-09-03/phase_C_numbers.md`
  §C2; `docs/open_items/run_2026-09-02/estarlight_li6.md` §2f (the
  photoproduction Q² scan, 2026-09-04); `validation/o5_a2_reach.py`;
  `python/tests/test_o5_reach.py`; `tests/test_coherent.cpp` T10a and T10c.
* What the first two versions of this file got wrong and by how much:
  `docs/OPEN_ITEMS_SOLUTIONS.md` §11.3a.
* Open item O1 (task C1): `phase_C_numbers.md` §C1; `tests/test_cluster_config.cpp`
  T22b.
* Open item O3 (this file, §5): `validation/o3_alpha_correlation_bound.py`.
* The two earlier drafts, left in place with their editor's notes:
  `docs/OPEN_ITEMS_SOLUTIONS.md` §11.2; `docs/open_items/run_2026-09-02/
  design_G_cluster_config.md` §10.
* The design and every sampler number: `docs/open_items/run_2026-09-02/
  design_G_cluster_config.md`; `docs/open_items/run_2026-09-02/phase_G_numbers.md`.
