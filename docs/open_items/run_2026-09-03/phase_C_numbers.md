# Phase C — the ⁶Li quadrupole budget (O1), the EIC reach of its a₂ (O5), the α–d overlap default (O2), the ΔB convention (O4), the VMC inputs, O3's offline bound and the collaboration ask (C6)

Measured 2026-09-03 (§C1, §C3) and 2026-09-04 (§C2, §C4, §C5, §C6) from the committed tables
at commit `31ed181`, through `include/lipolgen/cluster_config.hpp` /
`src/core/cluster_config.cpp` / `include/lipolgen/coherent.hpp` and the
`build/python` bindings. Every number below is reproducible with
`PYTHONPATH=build/python`; the ones that matter are pinned in
`tests/test_cluster_config.cpp` **T22b** (§C1) and **T22c** (§C3),
`tests/test_coherent.cpp` **T10a** (§C2), and mirrored in
`python/tests/test_cluster_config.py` and `python/tests/test_o5_reach.py`.
§C2's whole arithmetic is `validation/o5_a2_reach.py`, which prints every
table in that section and self-checks. §C4 is pinned in
`tests/test_coherent.cpp` **T10b** and `tests/test_cluster_config.cpp`
**T23**; §C5 in `tests/test_tagged.cpp` **T24**, **T25** and **T26**, with
`python/tests/test_module.py` mirroring all five.

Nothing here changes a shipped default. The default α–d source is still
`FitRescaled`, `CoherentScenario`'s **values** are untouched, and every
previously pinned number is bit-for-bit unchanged; §C2 adds one new inert table
to `coherent.hpp` (`estarlight_li6_coherent()`, the single code home of the
eSTARlight cross sections that had lived only in prose), §C4 adds two derived
accessors and one inverse map, and §C5 adds four opt-in knobs whose zero rows
are today bit for bit.

---

## C1. Open item O1 — the factor 7.5 in Q(⁶Li) is **two** factors, not three

### C1.0 The mechanism no document had written down

`src/core/cluster_config.cpp:398-435` bisects the (G9) root for
`quadrupole_target_fm2`, then **mutates the D wave in place**:

```cpp
const double n = 1.0 - (1.0 - dial_s_ * dial_s_) * pd;
const double inv = 1.0 / std::sqrt(n);
for (double& v : ad_f2_) v *= dial_s_ * inv;
for (double& v : ad_f0_) v *= inv;
```

`asymptotic_ds_ratio()` (`cluster_config.cpp:771`) reads `ad_f2_ / ad_f0_`
**afterwards**, and the common `1/√n(s)` cancels in that ratio. Therefore

> **η(s) = s · η(1), exactly** — verified to 1e−13 in T22b.

Measured: η(1) = **−0.0482160903696**; at the anchor
`quadrupole_target_fm2 = −0.18557` the dial closes on
s = **0.520019565622** and η = **−0.02507331037**, and
0.520019565622 × (−0.0482160903696) = −0.02507331037 to the last digit.

Two consequences, and they point in opposite directions:

* **It buys a budget leg anchored on a measurement.** The first leg of the
  quadrupole gap can now be derived from George & Knutson's *measured*
  η = −0.025 ± 0.006 ± 0.010 (PRC 59, 598 (1999)) instead of from the GFMC
  number, which is what §2.7 of design G had to use.
* **It is not an independent check.** Dialling to a target Q and reading η back
  recovers the dial, not the wave function. η and `quadrupole_dial_s()` are the
  same number in different units.

`LI6_ETA_DS_GK` / `_STAT` / `_SYST` are new in `cluster_config.hpp`, bound in
`python/bindings.cpp`, and are the **single home** of that measurement; the
literal `-0.025` that T22 carried now reads through the constant.

### C1.1 The budget

Anchor: `quadrupole_target_fm2 = −0.18557 fm²`, the target that puts η on the
GK central value to 0.3 %. (The exact η = −0.025 root is
**−0.1842160147 fm²**; −0.18557 is the **regression anchor** the tests pin, not
a derived physics value.)

| leg | ratio | measured |
|---|---|---|
| model → η-matched dial | −0.615448256134 / −0.18557 | **3.3165288362** |
| η-matched → measured Q | −0.18557 / −0.0818 (`LI6_QUADRUPOLE_FM2`) | **2.26858190709** |
| **product** | −0.615448256134 / −0.0818 | **7.52381731215** |

3.3165288362 × 2.26858190709 = 7.52381731215 **identically** — the factorisation
telescopes, which is the whole point: any split of the 7.5 into a "D-wave" leg
and a "everything else" leg is fixed once the intermediate Q is named.

Reading of the two legs:

* **3.317× — "too much D wave."** The α–d relative motion carries a D
  amplitude 3.3× larger, in the Q-relevant sense, than the one the measured
  asymptotic D/S ratio implies. 𝒬 is *linear* in the D amplitude through the
  S–D interference term (which is 95.2 % of q_int), so a factor in the wave is
  very nearly a factor in Q.
* **2.269× — everything else.** Core/cluster polarization, the r⁴ weight of
  (G4) acting on a free deuteron's shape inside a compressing nucleus, and the
  missing 15 % non-α+d component. This leg is **not** attributed here; it is
  what is left over.

### C1.2 1/S_αd = 1.171 is **already inside** the −0.615 — it is a ceiling, not a factor

`cluster_config.cpp:369-374` divides **both** `a0` and `a2` by
`sqrt(s_alpha_d_)` **before** `ad_m_ = radial_moments(...)` is recomputed. Both
`q_int` (bilinear in a₀a₂) and `q_dd` (quadratic in a₂) therefore already carry
1/S_αd = 1/0.8542323487 = **1.1706416896**.

Writing it as a third multiplicative factor gives

    3.3165288362 × 1.1706416896 × 2.26858190709 = 8.80769421057

which **overshoots the measured ratio 7.52381731215 by 17 %**. T22b pins that
overshoot so the three-factor form cannot come back.

The honest statement is: **1.17× is a CEILING on what any coherent
missing-component model could add on top of the model Q**, since the model Q
already assumes the missing 15 % contributes nothing coherent. It never
multiplies the budget.

### C1.3 THE BAND IS THE PHYSICS. The central value is not.

George & Knutson: η = −0.025 ± 0.006 (stat) ± 0.010 (syst);
σ_comb = √(0.006² + 0.010²) = **0.0116619037897**.

Because η is exactly linear in the dial, that error bar maps *directly* onto a
model-Q interval. The dial scan (each row an actual `ClusterConfigSampler`):

| target Q [fm²] | `quadrupole_dial_s()` | η | `eps_b0_equivalent()` | `a2_from_geometry(0.3, 1)` |
|---|---|---|---|---|
| −0.09 | 0.4124242 | −0.0198855 | −0.0074028 | +0.0288912 |
| −0.13 | 0.4575291 | −0.0220603 | −0.0106929 | +0.0417318 |
| **−0.18557** | **0.5200196** | **−0.0250733** | **−0.0152637** | **+0.0595705** |
| −0.24 | 0.5810629 | −0.0280166 | −0.0197407 | +0.0770433 |
| −0.30 | 0.6481993 | −0.0312536 | −0.0246759 | +0.0963042 |

and the ±1 σ_comb edges, inverted through the same line:

| η | dial s | model Q_charge [fm²] | model/measured |
|---|---|---|---|
| −0.025 − σ_comb = **−0.0366619038** | 0.7603665811 | **−0.4004752268** | 1.5367948251× |
| −0.025 (central) | 0.5184991112 | −0.1842160147 | 3.341× |
| **−0.0149706307** | 0.3104903485 | **0.0000** | **0** (sign change) |
| −0.025 + σ_comb = **−0.0133380962** | 0.2766316412 | **+0.0297575059** | wrong sign |

> **The 1 σ band on η contains a sign change in Q.** Q_charge passes through
> zero at η = −0.0149706307, only **0.86 σ_comb** from the GK central value,
> because at s = 0.31 the α–d D wave has only just cancelled the deuteron's own
> +0.2697 fm². Below that the model quadrupole has the **opposite sign to the
> measurement**.

So the "3.317× leg" is a point estimate on a line whose 1 σ range runs from
**1.54×** at one edge to **an unbounded ratio through zero** at the other. It
is *not* "~1.5× to ~12×"; the near edge does not merely weaken, it flips.
Report the leg as **3.32× (1 σ: 1.54× … sign change)**, or report the band and
skip the leg.

### C1.4 …and therefore η does not discriminate between the two literature Q's

Running the dial the other way, onto each published quadrupole:

| target Q | dial s | implied η | pull vs GK |
|---|---|---|---|
| measured **−0.0818** (`LI6_QUADRUPOLE_FM2`) | 0.4031636152 | **−0.0194389733** | **+0.48 σ_comb** |
| GFMC AV18+IL7 **−0.20** (`LI6_QUADRUPOLE_GFMC_FM2`) | 0.5362175172 | **−0.0258543123** | **−0.07 σ_comb** |
| GFMC −0.14 (+1σ) | 0.4687884 | −0.0226031 | +0.21 σ_comb |
| GFMC −0.26 (−1σ) | 0.6034578 | −0.0290964 | −0.35 σ_comb |

**Both literature values sit inside GK's 1 σ.** The measured η is consistent
with a ⁶Li quadrupole anywhere from −0.0818 to −0.26 fm² and with zero. It
cannot be used to prefer one over the other, and it cannot be used to argue
that the model's −0.615 is "3.3× too big" to better than a factor ~2 either
way. What it *can* do — and this is the actual content of O1 — is establish
that the D-wave excess is **moderate and measured**, so the leading candidate
for the residual 2.3× is cluster polarization plus the missing component, not
an ANL normalization or phase-convention error.

### C1.5 The anchor point, pinned (T22b, and the pytest mirror)

At `quadrupole_target_fm2 = −0.18557`:

| quantity | value |
|---|---|
| `quadrupole_dial_s()` | 0.520019565622 |
| `asymptotic_ds_ratio()` | −0.02507331037 |
| `quadrupole_band_fm2()[2]` | −0.18557 (to 1e−12) |
| `eps_b0_equivalent()` | −0.0152637033625 |
| `a2_from_geometry(0.3, 1)` | +0.0595705416512 |
| `p_d_alpha_d()` | 0.00551970830165 |

Baseline (no dial), for the ratios: η = −0.0482160903696,
Q_charge = −0.615448256134, s = 1, S_αd = 0.8542323487.

**These are regression anchors, not measurements.** −0.18557 fm² is not a
prediction for Q(⁶Li); it is the point on the dial where η equals a number
whose own error bar spans a sign change (§C1.3).

---

## C2. Open item O5 — does ⁶Li's tensor a₂ survive EIC statistics?

**MARGINAL — and a BAND, not a point.** Coherent J/ψ off tensor-polarised ⁶Li,
over the **whole** Q² range and with **both** lepton channels reconstructed, at
the repository's own programme luminosity of 10 fb⁻¹/nucleon and its own
`tensor_flip_plan(0.6)`, is

> **S = 2.63 σ at the band's LOW EDGE and 2.84 … 3.29 σ at its TOP, with
> 3 σ at 8.3 … 13.0 fb⁻¹/u — inside the {1, 10, 100} fb⁻¹/u band at BOTH ends
> on every form, and close to one EIC year. Whether the TOP crosses 3 σ is
> NOT established.**

It is a band because the single largest multiplier in the chain, ε_det, was
being applied **at the wrong beam energy** and with **no decay-lepton factor at
all**, and the two pull opposite ways (§C2.3c); and its **top is a span rather
than an edge** because the ⁷Li → ⁶Li substitution inside ε_det straddles 1 once
it is read off entries that share a beam energy (§C2.3c(b), fourth pass). The uncorrected point estimate
— **2.62 σ, 13.1 fb⁻¹/u**, δa₂(|t| = 0.3) = 0.0100 against a predicted
a₂ = +0.0263 — lands 0.3 % *under* that band's low edge, because the two
omissions cancel to 0.7 % on the levers that can be quantified. **That cancellation is the reason
neither may be quoted without the other**: a reader told only about the beam
energy moves the verdict up, a reader told only about the leptons moves it
down, and both together barely move it.

> **This section previously said NO, and the NO was wrong.** It said
> *"0.75 σ … 3 σ needs 160 fb⁻¹/u … the measurement does not exist at any
> luminosity the EIC is quoted at."* That number came from **one lepton
> channel** in **one Q² window** — 0.1 < Q² < 100 GeV² — and both restrictions
> were artefacts, not physics. Correcting them multiplies the **rate** by
> 2.00 × 6.75 = **13.5**, and the significance goes from 0.749 σ to 2.62 σ —
> a factor 3.5, which is √13.5 = 3.67 less the 5.4 % this pass also pays by
> quoting the background-immune ⟨P_zz²⟩ = 0.81 instead of 0.90 (§C2.2). The
> old row is still arithmetically right *for that row* (§C2.4a); it was never
> the measurement.

Neither correction is a modelling choice. §C2.4b is the second lepton channel,
which §C2.4's own prose already granted and its arithmetic never used.
§C2.3a is the photoproduction region, measured with a fresh eSTARlight run of
the same build, beams and seed (`estarlight_li6.md` §2f, in code as
`estarlight_li6_q2_floors()`).

> **And this section then shipped a MARGINAL verdict whose efficiency chain
> had two defects, in opposite directions, and neither of them in §C2.8's
> "every assumption, in one list".** (i) ε_det was applied at the wrong beam
> energy — `COHERENT_JPSI_EFF_IR8_LI7` = 0.1775 is a ⁷Li number at **18 × 117.9
> GeV/u, that nucleus's own top energy** (⟨W⟩ = 43.2 in this tree's own row for
> those beams), used on a ⁶Li sample at **10 × 99.5** (⟨W⟩ = 30.2). §C2.8 item
> 3 listed ε_det as assumed flat in |t|, in φ_Δ, in Q² and across mesons;
> **beam energy was not in the list**, and `coherent.hpp` recorded "none at
> 10 × 99.5" without the direction or the size *the same paper supplies*.
> (ii) There was **no decay-lepton acceptance or reconstruction efficiency
> anywhere in the chain** — arXiv:2511.05638's number counts scattered
> *nuclei*, and its own text says the simulation carries no detector
> efficiency and no reconstructed-distribution acceptance. Both are now in the
> list, with their directions, in §C2.3c and §C2.8 items 3a and 3b. **They
> partly cancel, which is exactly why leaving both out was worse than leaving
> out either.**

The separation from the linearly-polarised-photon cos 2φ background is *not*
the problem — it is free, it survives Q² → 0, and the P_zz flip is a **1.5×
gain**, not a cost (§C2.5). What is still open is **not** statistics: it is
that **no detection efficiency exists anywhere in this tree below
Q² = 0.1 GeV²**, where 85 % of the rate now sits, that **no decay-lepton
reconstruction efficiency exists anywhere in this tree at all** — which is
what leaves the band open below — and that the far-forward working point has
not been chosen. That last one is **the single correction that on its own
restores the original NO** (§C2.8 items 3, 3a, 3b and 8).

The full estimate is `validation/o5_a2_reach.py` (run it; it prints every
table below and self-checks), pinned in `python/tests/test_o5_reach.py` and
`tests/test_coherent.cpp` **T10a**/**T10c**. Every number it multiplies is
read from code: `a2_from_quadrupole` at the **measured** `LI6_QUADRUPOLE_FM2`;
`estarlight_li6_coherent()` **and** `estarlight_li6_q2_floors()` in
`coherent.hpp`; `COHERENT_JPSI_EFF_IR8_LI7`; `Scenario` / `tensor_flip_plan`;
and `tagging_optics` / `yr_optics` for the luminosity fraction.

> *(That sentence used to read "Nothing is retyped and nothing outside those
> three enters", and it was false: `EPS_DET_CHANG = 0.1775` — the single
> largest multiplier in the chain — and the measured-radius slopes
> `{55.0, 54.8, 54.1}` were both typed into the script, the second with no
> code home and no pin at all, so a re-fit of `estarlight_li6.md` §2b could
> not have reached it. Both now live in `coherent.hpp` and are pinned in T10a
> and in pytest.)*

### C2.0 What this estimate is, and the limitation that rides with every number

`a2_from_quadrupole` is a **closed form, not a dipole-model amplitude**. It
carries the target's quadrupole moment through the deuteron's own published
|t| dependence ([Mant24] Eq. (9)) and reproduces that published a₂(m = ±1) to
8 % with zero free parameters — but it contains no Good–Walker average, no
amplitude, no saturation model and none of their uncertainties, and it treats
the *matter* quadrupole as a proxy for the transverse **gluon** anisotropy.
Everything below inherits that. Read it as *"what a measured quadrupole of
this size buys, if the map holds"*, never as a prediction of the coherent
cross section's azimuthal structure. A dipole-model run could move the
signal by a factor; it would now have to move it **down by 1.31** to fall
back to 2 σ at 10 fb⁻¹/u, or **up by 1.15** to reach 3 σ there — a much
narrower lever than the factor 4.0 the first pass needed.

One thing the photoproduction correction *improves*: `mantysaari_a2_deuteron()`
is digitised from a **photoproduction** calculation
([arXiv:2408.13213](https://arxiv.org/abs/2408.13213) Fig. 4, x_P = 1.7e−3, and
`coherent.hpp` says "photoproduction" in as many words). So Q² < 0.1 is where
the map is *anchored*; the 0.1 < Q² < 100 window this item used to stop at was
the extrapolation, not the other way round. (The x_P is still not ours: at
⟨W⟩ ≈ 30 GeV the eSTARlight sample sits near x_P ≈ 1e−2, an order of magnitude
above the digitised curve. The map has no x dependence, so this does not enter
the arithmetic — it enters the 8 % validation, and it is part of §C2.0's
caveat, not separate from it.)

### C2.1 The signal — and the |t| slope, which is what actually kills it

a₂(m = ±1, |t|) = κ|t| exactly, with κ = `a2_from_quadrupole(2Q, 6, 1, 1)`
linear in Q:

| Q_charge [fm²] | κ [GeV⁻²] | a₂(\|t\| = 0.3) |
|---|---|---|
| **measured** `LI6_QUADRUPOLE_FM2` = −0.0818 | **0.08752978** | **+0.026259** |
| GFMC AV18+IL7 `LI6_QUADRUPOLE_GFMC_FM2` = −0.20 | 0.21400924 | +0.064203 |
| this α+d geometry, −0.615448 | 0.65855806 | +0.197567 |

The +0.026 that §11.2 and design G quote is **a₂ at |t| = 0.3 GeV², and the
coherent sample does not live there.** dN/d|t| ∝ e^{−B|t|} with B ≈ 39–55, so
the information-weighted modulation is κ√⟨t²⟩, not κ × 0.3:

| B [GeV⁻²] | √⟨t²⟩ [GeV²] | **ā₂ = κ√⟨t²⟩** |
|---|---|---|
| 38.9 (eSTARlight's own R_G = 2.1805 fm) | 0.036065 | **0.00316** |
| 50.0 (`CoherentScenario::slope_b`) | 0.028246 | **0.00247** |
| 55.0 (measured R = 2.589 fm) | 0.025698 | **0.00225** |

**The observable is a 0.25 % azimuthal modulation, not a 2.6 % one** — a
factor **10.6** below the quoted a₂(0.3) at B = 50. That factor is real and it
is not a modelling uncertainty: it follows from the coherent form factor
alone. **What it is not is a verdict.** A modulation is not a significance:
S = ā₂·√(2⟨P_zz²⟩N), so on the sample that would actually be taken — the
no-floor J/ψ row, B = 38.8 GeV⁻², ā₂ = 0.00316, N = 4.22 × 10⁵,
⟨P_zz²⟩ = 0.81 — that is 0.00316 × 827.2 = **2.62 σ**. The factor 10.6 was
correctly identified in the first pass and then paired with a rate that was
13.5× too small.

Truncating at `COHERENT_T_MAX_DEFAULT` = 0.2 GeV² costs nothing — it keeps
**99.995 %** of the rate and **99.72 %** of the Fisher information at B = 50 —
so the |t| window is not where anything is lost, and the truncation keeps the
whole sample inside the |t| ≤ 0.30 range over which [Mant24]'s a₂ was
digitised and the linear form validated. (⟨|t|⟩ = 1/B = 0.0257 GeV² at
B = 38.9 reproduces eSTARlight's own generated ⟨|t|⟩ = 0.0254, so the
exponential is the generator's actual distribution, not an approximation of
it — a Gaussian form factor has no diffractive minimum, `estarlight_li6.md`
§2c.)

### C2.2 The luminosity, the run plan, and which ⟨P_zz²⟩ the headline is

`Scenario::lumi_fb_per_nucleon` = **10.0 fb⁻¹/u**, "one EIC year", band
{1, 10, 100} — the repository's own default (`sampler.hpp`), and the only
luminosity it has. It is **not** a Li number: *"Li luminosity. Confirmed gap
— no Li number exists in any document (EPIOS included)"*
(`docs/surveys/needs_survey.md` §3.8, quoting `plans/04:114-117` #4). It is
quoted **per nucleon**, so the e+⁶Li luminosity is 10/6 = **1.6667 fb⁻¹**, and
that division is where a factor 6 lives that is easy to drop.

`tensor_flip_plan(0.6)` (`bookkeeping.hpp`) is the run plan: m = ±1-rich
bunches at **P_zz = +0.6** and m = 0-rich bunches at **P_zz = −1.2**, equal
shares, axis transverse (θ_S = 90°), electrons unpolarised.

**Two ⟨P_zz²⟩, and the headline uses the smaller one.**

| estimator | ⟨P_zz²⟩ | what it assumes |
|---|---|---|
| optimal combination of the two fills | **0.90** = 2.5 × 0.6² | that the cos 2φ_γ background is absent, or already removed |
| **background-immune two-fill difference** | **0.81** = 0.90/(10/9) | nothing — §C2.5's parity argument, 5.4 % worse in δ |

The first pass quoted 0.90 without saying which it was, while §C2.5 computed
the 5.4 % cost of the estimator its own separation argument requires. **Every
headline in this section is now on 0.81**; the tables print both, and the
0.749 σ of §C2.4a is retained on 0.90 because that is the number the item
shipped with (on 0.81 it is 0.711 σ, 178 fb⁻¹/u).

**The far-forward working point is a choice, and it is not applied.**
`estarlight_li6.md` §2e instructs "multiply by `Optics::lumi_fraction`", which
for ⁶Li at 10 × 100 is **0.0781** (0.1467 at 5 × 41; 1 at the Yellow Report
envelopes). It is not applied here, for a stated reason: the efficiency this
chain uses is arXiv:2511.05638's **IR-8 secondary-focus** number, and that
paper's own argument for IR-8 (its p. 4) is that the secondary focus buys the
low-p_T far-forward acceptance **without** the β*_x de-squeeze that costs
luminosity in IR-6 — which is exactly what `lumi_fraction` prices. They are
alternatives, not multipliers. **If the tensor run must instead share
LiPolGen's own de-squeezed ⁶Li tagging optics, every significance below
multiplies by √0.0781 = 0.279 and every 3 σ luminosity by 12.8: the headline
becomes 0.73 σ and 168 fb⁻¹/u, back outside the band.** That is the single
largest open choice in this item.

### C2.3 The rate and the reach, per channel

Chain: σ × L(e+⁶Li) → produced; × branching → reconstructible; × ε_det.
σ and B always travel **together** — a bigger nucleus means both a steeper
slope and a smaller cross section, so the measured-radius rows use the
measured-radius slope (`EstarlightLi6Row::b_rmeas`, in code since 2026-09-04).

#### C2.3a The Q² window was an acceptance study's, not a physics window

`estarlight_li6_coherent()`'s σ stops at Q² > 0.1 GeV², and `estarlight_li6.md`
§1 says where that came from: *"this is arXiv:2511.05638's range verbatim"* —
copied so the ⁷Li cross-check would be configuration-identical. **It is not
eSTARlight's Q² reach**, and most of the rate is below it. The scan
(`estarlight_li6.md` §2f; same build, beams and seed, only `MIN_GAMMA_Q2`
changed; the 0.1 rows reproduce §2a/§2b to every printed digit):

| VM | Q² floor | σ_coh | ×(0.1) | B | σ_coh (R = 2.589) | ×(0.1) | B |
|---|---|---|---|---|---|---|---|
| J/ψ | 0.1 | 1.773 nb | 1.000 | 38.9 | 1.255 nb | 1.000 | 55.0 |
| J/ψ | 0.01 | 3.348 nb | 1.888 | 39.2 | 2.371 nb | 1.889 | 55.0 |
| J/ψ | **none** | **11.971 nb** | **6.752** | 38.8 | **8.458 nb** | **6.739** | 54.6 |
| φ | 0.1 | 30.16 nb | 1.000 | 39.2 | 22.08 nb | 1.000 | 54.8 |
| φ | 0.01 | 103.18 nb | 3.421 | 39.3 | 76.12 nb | 3.448 | 55.1 |
| φ | **none** | **654.34 nb** | **21.696** | 38.9 | **482.95 nb** | **21.874** | 54.4 |
| ρ⁰ | 0.1 | 506.4 nb | 1.000 | 38.6 | 379.5 nb | 1.000 | 54.1 |
| ρ⁰ | 0.01 | 2285 nb | 4.512 | 38.9 | 1750 nb | 4.611 | 54.1 |
| ρ⁰ | **none** | **17823 nb** | **35.190** | 38.7 | **13736 nb** | **36.194** | 54.0 |

"none" is `MIN_GAMMA_Q2 = 0`, i.e. eSTARlight's own kinematic limit
Q²_min = (m_e E_γ)²/(E_e(E_e − E_γ)) ≈ 10⁻⁹ GeV². **85 % of the coherent J/ψ
rate is below Q² = 0.1.** The multiplier is a *flux* effect and not a density
one — the two densities give 6.752 and 6.739 for J/ψ, agreeing to 0.2 % — and
it is bounded by the logarithm: at 10⁻⁹ GeV² the ∫dQ²/Q² has run out of
decades. (An adversarial Weizsäcker–Williams estimate calibrated to ⟨Q²⟩ =
0.906 predicted ×7.6; the generator says ×6.75, 13 % lower. The direction and
the order were right and the number is now measured, not estimated.)

**The |t| slope does not move**: J/ψ 38.94 → 38.78 GeV⁻² (−0.4 %),
⟨|t|⟩ 0.0254 → 0.0255, and 54.98 → 54.59 at the measured radius. So the recoil
p_T spectrum the far-forward acceptance cuts on is the **same sample** — which
is the one fact that makes an efficiency measured above Q² = 0.1 arguably
transferable below it. ⟨W⟩ falls, 32.2 → 30.2 GeV, and by arXiv:2511.05638's
own Fig. 2 the tagging efficiency *rises* as W falls, so the transfer is if
anything conservative on the W axis. **Neither of those is a measurement**;
see §C2.8 item 3.

#### C2.3b The reach

At 10 fb⁻¹/u, ⟨P_zz²⟩ = **0.81** (background-immune), |t| < 0.2 GeV²,
**ε_det = 0.1775** = `COHERENT_JPSI_EFF_IR8_LI7` (arXiv:2511.05638's *global*
far-forward tagging efficiency × acceptance for coherent J/ψ off e+⁷Li at
18 × 118 through the IR-8 secondary focus, its p. 4 — there is no ⁶Li number,
this is a stand-in applied to all three mesons), J/ψ with **both** lepton
channels (`branching_all` = 0.11932 = eSTARlight's own JpsiBree + JpsiBrmumu):

| VM | window | σ [nb] | produced | × BR | × ε_det | **δa₂(0.3)** | **S** | 3 σ at |
|---|---|---|---|---|---|---|---|---|
| J/ψ | Q² > 0.7 | 0.605 | 1.01e6 | 1.20e5 | 2.14e4 | 0.0447 | 0.587 | 261 |
| J/ψ | Q² > 0.1 | 1.773 | 2.96e6 | 3.52e5 | 6.26e4 | 0.0261 | 1.005 | 89.1 |
| J/ψ | Q² > 0.1, R = 2.589 | 1.255 | 2.09e6 | 2.50e5 | 4.43e4 | 0.0436 | 0.603 | 248 |
| J/ψ | Q² > 0.01 | 3.348 | 5.58e6 | 6.66e5 | 1.18e5 | 0.0192 | 1.371 | 47.9 |
| J/ψ | Q² > 0.01, R = 2.589 | 2.371 | 3.95e6 | 4.72e5 | 8.37e4 | 0.0317 | 0.828 | 131 |
| **J/ψ** | **no Q² floor** | **11.971** | **2.00e7** | **2.38e6** | **4.22e5** | **0.0100** | **2.618** | **13.1** |
| J/ψ | no floor, R = 2.589 | 8.458 | 1.41e7 | 1.68e6 | 2.99e5 | 0.0167 | 1.576 | 36.2 |
| φ | Q² > 0.7 | 1.897 | 3.16e6 | 1.55e6 | 2.75e5 | 0.0126 | 2.09 | 20.6 |
| φ | Q² > 0.1 | 30.16 | 5.03e7 | 2.46e7 | 4.37e6 | 0.0032 | 8.34 | 1.29 |
| φ | no Q² floor | 654.3 | 1.09e9 | 5.34e8 | 9.48e7 | 0.0007 | 39.1 | 0.06 |
| ρ⁰ | Q² > 0.7 | 15.32 | 2.55e7 | 2.55e7 | 4.53e6 | 0.0031 | 8.62 | 1.21 |
| ρ⁰ | Q² > 0.1 | 506.4 | 8.44e8 | 8.44e8 | 1.50e8 | 0.0005 | 49.5 | 0.04 |
| ρ⁰ | no Q² floor | 17823 | 2.97e10 | 2.97e10 | 5.27e9 | 0.00009 | 293 | 0.00 |

(The script prints all 21 rows, at three efficiencies and both ⟨P_zz²⟩. The
φ/ρ⁰ rows are unchanged by `branching_all` — neither has a second
reconstructible channel.)

With **ε_det = 1** — a perfect detector, an upper bound and not a projection —
the J/ψ no-floor row goes to **6.21 σ**, and Q² > 0.1 to 2.39 σ.

Events needed at ⟨P_zz²⟩ = 0.81: **3 σ at N = 5.6e5 (B = 38.9), 9.1e5
(B = 50), 1.10e6 (B = 55)**; 5 σ at 1.5e6 / 2.5e6 / 3.1e6.

#### C2.3c The efficiency chain itself — two defects, opposite directions

*(2026-09-04, third pass. Arithmetic in `validation/o5_a2_reach.py` §3b, the
two source tables in `coherent.hpp` as `chang26_he3_energy_scan()` and
`chang26_species_efficiency()`, pinned in **T10d** and in pytest.)*

Until this pass the chain was **σ × BR(ℓ⁺ℓ⁻) × ε_det and nothing else**, with
ε_det = 0.1775 taken flat in |t|, in φ_Δ, in Q², across mesons and — unstated —
**across beam energy**, and with **no factor at all for the decay leptons**.

**(a) THE BEAM ENERGY. UP, and it is measured.** `COHERENT_JPSI_EFF_IR8_LI7`
is arXiv:2511.05638's ⁷Li entry at **18 × 117.9 GeV/u** — ⁷Li's own *top*
energy, Z/A × 275, which the paper writes "18 × 118". This tree's
configuration-identical eSTARlight row for those beams has ⟨W⟩ = **43.2 GeV**
(`estarlight_li6.md` §2a). It is applied to a ⁶Li sample at **10 × 99.5**,
⟨W⟩ = **30.2**. The same paper measures that dependence, in its §V.B, on ³He:

| e × ³He | efficiency |
|---|---|
| 18 × 183 GeV/u | **32.23 %** |
| 10 × 100 GeV/u | **54.38 %** |
| 5 × 41 GeV/u | **99.77 %** |

*(Verified against the PDF: p. 5–6, "At the top collision energy, 32.23 % of
the scattered ³He nuclei occur within a safe distance from the beam … At the
energy of 10 × 100 GeV², the total detection efficiency is 54.38 % … For the
lowest collision energy, 5 GeV electron beams on 41 GeV ³He beams, the
detection efficiency is 99.77 %." The lower two ion energies are the paper's
own text, not Z/A-scaled the way its 183 is.)*

d ln ε / d ln E_ion = **−0.866** (183 → 100) and **−0.681** (100 → 41); the two
differ because ε is capped at 1 and must flatten, so both are carried as a band
rather than fitted to one power law. **For the step this transfer actually
requires, 117.9 → 99.5 GeV/u, that is ×1.122–1.158 UP: ε_det = 0.199–0.206,
not 0.1775.**

**It is NOT the paper's own 1.687.** 54.38/32.23 is their 183 → 100 ratio, a
step **3.56× larger in ln E** than the one needed here. Carrying it whole gives
S = 3.40 σ and 3 σ at 7.8 fb⁻¹/u — which is scaling by the wrong lever arm, not
a more honest number. (It would also have broken the old `2.0 < S < 3.0`
self-check pin, which is a symptom, not an argument.)

The W axis says the same thing and is **not independent of it**: ⟨W⟩ 43.2 →
30.2, and by the same paper's Fig. 2 lower W means higher efficiency. Not
applied on top — it is the same effect seen on the other axis.

**(b) THE SPECIES. DIRECTION UNDETERMINED — and this is the one factor in the
chain that was presented as quantified without ever being quantified.**
§C2.8 item 3 flagged the ⁷Li → ⁶Li substitution with no direction and no size;
the third pass gave it both, **wrongly**, and the fourth pass takes them back.

> ~~Every entry of the paper's p. 4 species list is at its own top energy
> Z/A × 275, so `A/Z × E_ion` = 274.0–275.3 GeV/e for all seven — the same
> magnetic rigidity. The A-ordering in that list is therefore not a rigidity
> effect and reads as a species lever on its own.~~ **RETRACTED, 2026-09-04,
> fourth pass.** The rigidity statement is true and checked in T10d; the
> conclusion drawn from it is a **non sequitur**. Eliminating rigidity
> eliminates rigidity and nothing else: at fixed R = A E/Z the per-nucleon
> energy is **E/u = R·Z/A** and the total beam momentum is **p_z = Z·R**, and
> both still vary down the list — E/u from **118 GeV/u (⁷Li) to 183 (³He)**,
> Z from **1 to 8**, p_z from **274 to 2192 GeV**. And the confound is the
> *size of the effect being read off it*: leg (a)'s own d ln ε/d ln E =
> −0.68 … −0.87 makes the list's **×1.55 spread in E/u worth ×1.35–1.46** in
> ε — as large as the whole ×1.16–1.22 that was being attributed to species.
> The list is a **joint (A, Z, E/u) lever**, and this chain asserted the
> decomposition without ever testing it.

**The confounded reading, kept only so the confound is visible.** Four closed
forms interpolating ⁴He (29.42 %, **at 137 GeV/u**) onto ⁷Li (17.75 %, **at
118 GeV/u**) at A = 6:

| construction | ε(⁶Li) | factor |
|---|---|---|
| log-linear in the Gaussian radius R_G = 1.2 A^{1/3} | 0.2060 | ×1.160 |
| log-linear in A | 0.2101 | ×1.183 |
| this repository's own `tag_acceptance` = exp(−B p_T,cut²) | 0.2101 | ×1.184 |
| linear in A | 0.2164 | ×1.219 |

→ ×1.16–1.22 — but its two anchors are **19 GeV/u apart in E/u**, which on
leg (a)'s own slopes is worth ×1.11–1.14 by itself.

**The reading with no energy step in it, and the one the band now uses.**
**Four** of the seven entries share a beam energy — ²D, ⁴He, ¹²C, ¹⁶O, all at
**137 GeV/u** — and ⁶Li's own fixed-rigidity energy is Z/A × 275 = **137.5**,
so those four bracket A = 6 with **no energy step at all**. ⁴He → ¹²C is the
adjacent pair:

| construction | ε(⁶Li, 137 GeV/u) | factor |
|---|---|---|
| `tag_acceptance`, p_T,cut inverted from ¹²C | 0.1763 | **×0.993** |
| log-linear in R_G | 0.1782 | ×1.004 |
| log-linear in A | 0.2006 | ×1.130 |
| `tag_acceptance`, p_T,cut inverted from ⁴He | 0.2012 | ×1.134 |
| linear in A | 0.2366 | **×1.333** |

**→ ×0.99–1.33. IT STRADDLES 1.** Neither the size nor the **direction** of
the ⁷Li → ⁶Li substitution is established, and that is the whole of the band's
top (§C2.4a, §C2.7). Chained with ⁶Li's own 137.5 → 99.5 GeV/u energy step
(×1.246–1.323) the two legs give **×1.238–1.763** in total, against the
×1.446–1.613 the third pass quoted. *(Linear-in-A is the crude end of that —
a straight chord across a 4.6× fall between A = 4 and A = 12 — but it is one
of the same four forms the confounded reading was quoted from, and dropping it
now would be choosing the answer.)*

**The one form that can be tested, read three ways — and the Z that was wrong
in it.** The form is ε = exp(−B(A)·p_T,cut²) with B = `gaussian_slope(1.2
A^{1/3})` and p_T,cut = **0.1958 GeV** inverted from the ⁷Li row. The criterion
arXiv:2511.05638 states is *"within a safe distance **from the beam**"* — a cut
on the **angle** — so **p_T,cut = θ·p_z = θ·Z·R: it carries the charge.** ³He
and ⁴He are Z = 2; ⁷Li is Z = 3.

| reading | ³He/⁴He (measured 0.3223/0.2942 = **1.0955**) | the list's absolute values, pred/meas over ²D … ⁹Be |
|---|---|---|
| **p_T,cut ∝ Z** — the criterion taken literally, (2/3)² on the pair | **1.0967 — 0.1 %** | poor: 1.95 (²D), 2.00 (³He), 2.00 (⁴He), 0.21 (⁹Be) |
| **p_T,cut Z-independent** | 1.2309 — *this is the shipped one*, "overstates by 12 %" | **good**: **1.003** (²D), **1.034** (⁴He), **1.047** (⁹Be); ³He off by **16 %** |

Under the second reading the only misfit among the light entries is ³He, and
³He is also **the one nucleus in the list at E/u = 183** rather than the
family's 137 — for which leg (a)'s own scan predicts an **18–22 % shortfall**,
which is what the 16 % is. *(And say what that reading does **not** do: the
form degrades with A, missing ¹²C by 32 % and ¹⁶O by a factor 3.1. It is a
**local** interpolation over the light end, not a description of the list —
A = 6 sits inside the range where it works, which is the most that can be
claimed for it.)* **Three readings of one table, and nothing in the table picks
between them.** The third pass took the one that yields the directional instruction
*"Read the low end"* — printed at three sites — and that instruction is now
**deleted**. `tests/test_coherent.cpp` T10d carries the Z-corrected prediction
beside the wrong-Z one so the choice stays visible.

**(c) THE DECAY LEPTONS. DOWN — and this factor was not in the chain at all.**
arXiv:2511.05638's number is the fraction of scattered **nuclei** in the
far-forward acceptance ("32.23 % of the scattered ³He nuclei occur within a
safe distance from the beam"), and its own p. 4 says the simulation *"only
accounts for the acceptance effect and does not incorporate the efficiencies
of the detector. Additionally, we did not account for the efficiency and
acceptance of the reconstructed distribution"*. So σ × BR(ℓ⁺ℓ⁻) × ε_recoil
carries **no central-detector acceptance and no reconstruction efficiency for
the e⁺e⁻/μ⁺μ⁻ pair**. It is ≤ 1 by construction.

**The geometry is bounded here, and it is not the problem.** At 10 × 99.5 the
photon that makes a J/ψ at ⟨W⟩ = 30.2 has only **E_γ = 2.29 GeV** in the lab,
so the J/ψ sits at **y = −0.39** — very nearly at rest — and its two 1.548 GeV
decay leptons are central. Both are inside |η| < `Scenario::eta_max` = 3.5
(this repository's *own* "crude central-detector acceptance", documented there
for the **scattered electron**, because no decay-lepton acceptance exists in
this tree at all) exactly when |cos θ*| < tanh(3.5 − |y_ψ|):

| W [GeV] | E_γ,lab | y(J/ψ) | A(η only) | A(+ p_T > 0.2) |
|---|---|---|---|---|
| 20.0 | 1.003 | +0.434 | 0.9935 | 0.9875 |
| **30.2** (⟨W⟩) | 2.289 | −0.391 | **0.9940** | 0.9875 |
| 32.2 | 2.603 | −0.519 | 0.9923 | 0.9875 |
| 40.0 | 4.018 | −0.954 | 0.9818 | 0.9818 |
| 50.0 | 6.279 | −1.400 | 0.9563 | 0.9563 |
| 63.1 (E_γ = E_e, the kinematic top) | 10.000 | −1.865 | 0.8940 | 0.8940 |

**A_geom = 0.994** at this sample's own ⟨W⟩, and **never below 0.89** anywhere
the beams can reach. SCHC transverse decay throughout — (3/8)(1 + cos²θ*), the
*pessimistic* weighting, which pushes leptons toward the beam; isotropic gives
0.9960. p_T(J/ψ) ≈ 0.16 GeV is neglected (it helps one lepton and hurts the
other). The p_T > 0.2 GeV column is the sibling `../PolarizedLithiumSim`'s own
`HfsModel(pt_min_track=0.2)` stand-in, not an in-tree number.

**The reconstruction efficiency is UNBOUNDED here, and that is the result.**
No per-lepton tracking or PID efficiency exists anywhere in this tree. The
sibling assumes `eff_track` = 0.95 and **labels it a stand-in in as many
words**; the pair would then cost 0.9025, and A_geom × that is **0.897**. That
0.897 is used below only as a *stand-in*, and it is what makes the band's low
end a band edge rather than a bound.

**(d) THEY PARTLY CANCEL.** ×1.122 (energy leg) × 0.897 (geometry × the
stand-in pair factor) = **1.007 — unity to 0.7 %.** With the species leg on top
the product is 1.447. That is why both had to be stated: **an unquantified
factor named as unquantified is a result; one left out of a list that claims
completeness is a defect**, and here the two defects were of the same size and
opposite sign, so a reader given either one alone would move the verdict in one
direction.

### C2.4 The deciding number

Coherent **J/ψ** is the channel the case is built on: the only one inside
`COHERENT_MX_MIN_DEFAULT` = 1.2 GeV, and the only one [Mant24] actually
computed. Over the whole Q² range, with both lepton channels, at one EIC year
it delivers **4.22 × 10⁵ reconstructed events**, hence

> **δa₂(\|t\| = 0.3) = 0.0100 against a predicted a₂ = +0.0263 — S = 2.62 σ.
> Three sigma needs 13.1 fb⁻¹/u, inside the {1, 10, 100} fb⁻¹/u band.**
> (On the optimal ⟨P_zz²⟩ = 0.90: 2.76 σ, 11.8 fb⁻¹/u.)

**That row is the chain with ε_det = 0.1775 flat, and §C2.3c is why it is the
middle of a band rather than the answer.** With the beam-energy leg and the
decay-lepton factor both applied the band's low edge is **2.63 σ, 3 σ at
13.0 fb⁻¹/u** — and the uncorrected row lands 0.3 % under it, because the two
omissions cancel to 0.7 %. Its **top is 2.84 … 3.29 σ, 3 σ at
8.3 … 11.1 fb⁻¹/u**, a span rather than an edge because the species leg
straddles 1 (§C2.3c(b)). Quote the band; never the point, and never one edge
of the band alone.

#### C2.4a The row this item shipped with, and why it was not the measurement

| ladder rung | S | 3 σ at |
|---|---|---|
| Q² > 0.1, e⁺e⁻ only, ⟨P_zz²⟩ = 0.90 — **the number this item shipped with** | **0.749 σ** | **160 fb⁻¹/u** |
| the same, on the background-immune 0.81 | 0.711 σ | 178 fb⁻¹/u |
| + both lepton channels | 1.005 σ | 89.1 fb⁻¹/u |
| + the conservative Q² > 0.01 floor instead | 1.371 σ | 47.9 fb⁻¹/u |
| + both, on the measured radius R = 2.589 fm | 0.828 σ | 131 fb⁻¹/u |
| no Q² floor, both leptons, measured radius | 1.576 σ | 36.2 fb⁻¹/u |
| **no Q² floor, both leptons, default density** | **2.618 σ** | **13.1 fb⁻¹/u** |
| … and with the de-squeezed tagging optics on top of that | **0.731 σ** | **168 fb⁻¹/u** |

**Every conservatism on this ladder EXCEPT ONE leaves 3 σ inside the band on
its own; it takes two of them stacked — the Q² > 0.01 floor *and* the
measured-radius density — to put it back out. THE ONE EXCEPTION IS THE LAST
ROW**: the de-squeezed far-forward optics do it alone, and they are a *choice*
rather than a measurement. That is what makes the verdict MARGINAL rather than
YES, and it is why "the measurement does not exist at any luminosity the EIC is
quoted at" was not a defensible sentence: it was a statement about one rung of
that ladder.

> *(2026-09-04, third pass — **this paragraph used to be a self-contradiction,
> at two of its three sites.** It read "Every single conservatism on its own
> leaves 3 σ inside the band" with the escape clause "or the far-forward
> working point on its own" tacked to the end of the same sentence, which is
> the same claim and its counterexample in one breath; `o5_a2_reach.py`
> printed it with **no** escape clause, three lines below its own bullet
> saying the optics row is "back outside the band"; and
> `mantysaari_collaboration_draft.md` repeated the version with no escape
> clause **and omitted the optics row from its ladder entirely**. All three
> now carry the exception as an exception, and the draft's ladder has the
> row.)*

**And the same ladder on the corrected ε_det of §C2.3c**, which is the honest
form. Every rung is the headline row with ε_det = 0.1775 multiplied by the
stated factor:

| rung | factor on ε_det | S | 3 σ at |
|---|---|---|---|
| as shipped — no beam-energy leg, no lepton factor | 1.000 | 2.618 σ | 13.1 fb⁻¹/u |
| + the beam-energy leg only | ×1.122–1.158 | 2.77 σ | 11.7 fb⁻¹/u |
| + the decay-lepton **geometry** only | ×0.9940 | 2.610 σ | 13.2 fb⁻¹/u |
| **BAND LOW** = beam lo × geometry × the pair stand-in | ×1.007 | **2.63 σ** | **13.0 fb⁻¹/u** |
| **BAND TOP**, species read low | ×1.179 | **2.84 σ** | **11.1 fb⁻¹/u** |
| **BAND TOP**, species read high | ×1.582 | **3.29 σ** | **8.3 fb⁻¹/u** |
| … de-squeezed optics on the BAND LOW rung | ×0.0787 | 0.734 σ | 167 fb⁻¹/u |
| … de-squeezed optics on the BAND TOP rung, species low | ×0.0921 | 0.794 σ | 143 fb⁻¹/u |
| … de-squeezed optics on the BAND TOP rung, species high | ×0.124 | 0.920 σ | 106 fb⁻¹/u |

**The TOP is a span because the species leg is.** The five closed forms of
§C2.3c(b) put it at **2.84, 2.86, 3.03, 3.04 and 3.29 σ** — two below 3 σ,
three above — so *"the band's top just crosses 3 σ"*, which the third pass
printed at six sites, is **not established** and is withdrawn. What *is*
established is that 3 σ sits at **8.3 … 13.0 fb⁻¹/u on every one of them**,
inside {1, 10, 100} at both ends.

**The band is open BELOW its low rung**, because the per-lepton reconstruction
efficiency is unbounded in this tree. What it takes to leave: a pair
acceptance × efficiency of **0.52** costs the "above 2 σ" half of MARGINAL, and
**0.117** — a detector reconstructing about one J/ψ in eight — puts 3 σ
outside {1, 10, 100}. So *"3 σ inside the band"* is robust to that unknown in a
way *"MARGINAL rather than NO"* is not, and both statements should be made
with that asymmetry attached.

#### C2.4b The second lepton channel: the write-up refuted itself

§C2.4 said, in the first pass, *"The row uses J/ψ → e⁺e⁻ only (branching
0.0597 …); adding μ⁺μ⁻ gains √2 and reaches 1.06 σ"* — three lines under a
headline of "3 σ needs 160 fb⁻¹/u". **160/2 = 80 fb⁻¹/u, which is inside the
{1, 10, 100} band**, so the document's own two numbers contradicted its
conclusion before any new run was made. The gain is √(0.11932/0.0597) =
**1.41374**, i.e. √2 to 0.05 %. Both branchings are eSTARlight's own
(`starlightconstants.h`: `JpsiBree` = 0.05971, `JpsiBrmumu` = 0.05961), and
`PROD_PID = 443013` reproduces the second one from the generator (106.801 pb
generated against 1.792 nb total = 0.05960).

The old self-check made this permanent rather than visible:
`assert res["lumi_for_3sigma_jpsi"] > 100.0` was computed from the e⁺e⁻-only,
Q² > 0.1 row, so it would have held forever no matter what the photoproduction
region or the second lepton did. It is now pinned in **both** directions —
that row still exceeds 100, and the headline row must be between 2 σ and 3 σ
with its 3 σ luminosity inside the band.

#### C2.4c The quadrupole is still the whole story

Run the *same* J/ψ sample (no floor, both leptons) at the other two
quadrupoles:

| Q_charge fed to the map | S at 10 fb⁻¹/u | 3 σ at |
|---|---|---|
| **measured −0.0818** | **2.62 σ** | **13.1 fb⁻¹/u** |
| GFMC AV18+IL7 −0.20 | 6.40 σ | 2.2 fb⁻¹/u |
| this α+d geometry −0.6154 | 19.7 σ | 0.2 fb⁻¹/u |

*(All three rows are the **uncorrected** ε_det = 0.1775 chain, so that the
comparison is between quadrupoles and nothing else; the first row is the
band's uncorrected middle and may not be quoted on its own — the band is
2.63 σ at its low edge and 2.84 … 3.29 σ at its top, §C2.4a.)*

**The factor 7.5 by which the geometry overshoots Q(⁶Li) — §C1's whole
subject — is still the difference between a comfortable measurement and a
marginal one.** ⁶Li is attractive as a *near-null* test against the deuteron's
+0.286 fm² (`physics_literature.md` §11 says so), and this is the price of
that: near-null is expensive. What it is *not*, on these numbers, is
unmeasurable.

### C2.5 The separation from the photon-polarisation cos 2φ — free, and not the problem

Two independent handles, and **neither is the bottleneck**.

**(a) The azimuthal reference.** The tensor modulation is in
cos 2(φ_Δ − φ_S), about the **spin axis**; the photon-polarisation one is in
cos 2(φ_Δ − φ_γ), about the **photon's linear-polarisation direction**, which
at Q² > 0 is the lepton plane. φ_S is fixed in the lab (vertical) and φ_γ is
uniform about the beam by rotational invariance, so
⟨cos 2(φ_γ − φ_S)⟩ = 0 and the background **averages to zero in a
spin-axis-referenced histogram** — provided the acceptance is uniform in φ_γ.
That much needs no run plan at all.

**(b) The P_zz parity, which is the real control.** The tensor term is
**odd** in P_zz (c₂ = 2κ|t|P_zz — exactly, because a₂(0) = −2a₂(±1) makes the
population average Σ_m p_m a₂(m) = P_zz·a₂(±1) identically, `a2_m_state`);
the photon-polarisation term is **even** in it, being a property of the photon
and of the target's unpolarised transverse structure. The difference of the
two fills' fitted amplitudes, divided by P_zz(1) − P_zz(2) = 1.8, is therefore
background-free at first order, and it cancels a relative-luminosity offset
exactly because each fill's amplitude is a self-normalised ratio (the
tolerance is on the **acceptance**, hence `tensor_flip_plan`'s bunch-by-bunch
alternation, not on the luminosities).

**BOTH SURVIVE Q² → 0, which is what lets §C2.3a be used at all.** Below
Q² = 0.1 the scattered electron is not detected and φ_γ is unknown event by
event. Neither handle needs it: (a) is a statement about a *spin-axis*
histogram under a φ_γ-uniform acceptance and does not require φ_γ to be known,
and (b) never references the lepton plane. What is lost at low Q² is the
ability to **project the background out directly** (§C2.6) — which was never
the control being used here.

**Its cost is negative.** Against the optimal use of the same two fills the
background-immune difference loses **1.111 in variance, 1.054 in δ** — 5.4 %,
and this section is the reason the headline is quoted on ⟨P_zz²⟩ = 0.81 and
not 0.90. Against putting the whole luminosity into one +0.6 fill it **gains**:
0.444 in variance, **0.667 in δ, a 1.50× improvement**, because the m = 0-rich
fill carries |P_zz| = 1.2. Asking "can we afford to flip P_zz?" has the answer
"you cannot afford not to."

**The residual, computed.** What breaks the parity is that the two fills'
|t| spectra are not identical: the same deformation that makes a₂ modulates
the slope, B → B(1 + ε cos 2φ) with ε = |`eps_b0_equivalent()`| = **0.006728**
at the measured Q. The **φ-averaged** spectra therefore differ only at
*second* order, ⟨e^{−εBt cos 2φ}⟩_φ = I₀(εBt) = 1 + (εBt)²/4 + …, i.e. by
**1.13 × 10⁻³ at |t| = 0.2 GeV²** and **2.26 × 10⁻⁵** averaged over the
window. A background of amplitude A_γ leaks at most that fraction of itself,
which matches the per-fill statistical error √(2/N) only above (no Q² floor,
both leptons, ε_det = 0.1775):

| VM (no Q² floor, ε_det = 0.1775) | N | √(2/N) | A_γ at which the residual matters |
|---|---|---|---|
| J/ψ | 4.22e5 | 2.2e−3 | **1.9** |
| φ | 9.48e7 | 1.5e−4 | **0.13** |
| ρ⁰ | 5.27e9 | 1.9e−5 | **0.017** |

So the flip is still exact enough for J/ψ by a wide margin, and the thresholds
for φ and ρ⁰ have tightened by exactly the √rate the photoproduction region
added: for **ρ⁰ it now bites at A_γ ≈ 0.02**, and there binning in |t| (which
removes the residual by construction) is mandatory rather than optional.
**No magnitude for A_γ exists anywhere in this tree** — `physics_literature.md`
records the mechanism (STAR
[arXiv:2204.01625](https://arxiv.org/abs/2204.01625)), not a number — which is
why the argument above is deliberately magnitude-free and rests on parity.

### C2.6 One correction to the framing this item was posed with

The item was posed with "the linearly-polarised-photon cos 2φ … is a real
photon effect that dies as Q² rises". **Nothing in this tree supports that,
and the standard virtual-photon result argues against it**: the degree of
linear polarisation of a virtual photon about the lepton plane is ε, a
function of **y** (→ 1 as y → 0), not a decreasing function of Q². What
changes with Q² is *measurability*, not magnitude — above Q² ≈ 0.7 GeV² the
scattered electron is in the central detector, so φ_γ is known event by event
and the background can be projected out directly, and below it it cannot. The
separation argument of §C2.5 does not use the premise and does not need it:
the P_zz parity does all the work, and §C2.5 now says so explicitly for
Q² → 0, because that is where the rate turned out to be. The second half of
the premise — that the background does not depend on the target's tensor
polarisation — is the half that is right, and it is the whole control.

Second, smaller: `estarlight_li6.md` §2a and its two quoting sites say the
fitted slopes reproduce the analytic B = R_G²/(3ħ²c²) = 40.7026 "to 4 %". The
measured deviations are **−3.69 % (φ), −4.43 % (J/ψ), −5.17 % (ρ)** — one-sided,
as the |t_min| floor and the electroproduction flux require, but spanning
3.7–5.2 %, not 4 %. Corrected at all three sites and pinned in T10a.

### C2.7 The verdict, and what it changes about the ask

**MARGINAL for coherent J/ψ — as a BAND — and no longer an argument against
the ask.**

* **J/ψ: MARGINAL.** **2.63 σ at the band's low edge and 2.84 … 3.29 σ at its
  top** in an EIC year over the whole Q² range with both leptons, **3 σ at
  8.3 … 13.0 fb⁻¹/u**, inside the {1, 10, 100} band at both ends. **Whether
  the top crosses 3 σ is not established** (§C2.3c(b)). 6.2 σ with a perfect detector *and* a perfect lepton
  pair, which is an upper bound and not a projection. This is the channel
  inside LiPolGen's own M_X floor and the only one [Mant24] computed. On the
  measured-radius density it is 1.6 σ and 36 fb⁻¹/u — still inside the band.
  It falls back below 1 σ only if two conservatisms are stacked, or if the
  de-squeezed tagging optics are imposed (§C2.2) — **that last one is the
  single correction that on its own restores the original NO**, and it does so
  at both ends of the band (0.73–0.92 σ, 106–167 fb⁻¹/u).
* **The band is a band, not an error bar.** Its top comes from one *measured*
  up-lever on ε_det (beam energy) times a species substitution whose direction
  is **not established** (×0.99–1.33, §C2.3c(b)); its bottom from a *stand-in*
  for the decay-lepton reconstruction efficiency that no source in this tree
  supplies. Below that bottom it is **open**: 0.52 on the
  pair costs "above 2 σ", 0.117 costs "3 σ inside the band" (§C2.3c, §C2.4a).
* **φ: YES on statistics.** 39 σ with no Q² floor, 8.3 σ over 0.1 < Q² < 100,
  and **2.1 σ** inside LiPolGen's own Q² > 0.7 window.
* **ρ⁰: YES on statistics alone.** 293 σ / 49 σ / 8.6 σ. But ρ⁰ sits **below**
  `COHERENT_MX_MIN_DEFAULT` = 1.2 GeV, so the shipped generator has no rate
  there at all; it is the channel in which the photon-polarisation cos 2φ has
  actually been *observed* (STAR's ρ⁰ ultraperipheral measurement — this tree
  records the reference, not its channel and not its magnitude, so treat the
  channel attribution as mine and not as sourced here) and where the
  P_zz-flip residual first bites, now at A_γ ≈ 0.02 (§C2.5); and it is
  where the closed-form map is **least** defensible — a large, strongly
  absorbed dipole at ⟨W⟩ = 14–18 GeV, where black-disc absorption saturates
  the very anisotropy the map is linear in, validated only against a **J/ψ**
  calculation.

**What that changes.** The drafted ask in §11.2 requests the Mäntysaari
group's polarised **J/ψ** setup on our ⁶Li configuration tables. The first
pass of this section concluded that this was the one channel that could not
see the effect and that the ask therefore argued against itself. **That
conclusion is withdrawn.** The honest recommendation is now:

1. **The ⁶Li tensor a₂ in coherent J/ψ is a marginal-to-good measurement at
   EIC luminosities**, and the ask is worth making on that channel — which is
   also the only channel the published calculation covers, so the theory ask
   stays the *easy* one.
2. **Propose it as a photoproduction measurement**, Q² < 0.1 included and
   dominant, not as the 0.1 < Q² < 100 electroproduction sample this item was
   first priced on. That is also where [Mant24]'s own Fig. 4 lives (§C2.0), so
   the request matches the published calculation's kinematics better than the
   first draft did.
3. **Say what is missing**, because it is not statistics: no detection
   efficiency exists below Q² = 0.1 in any source this tree has seen; **no
   decay-lepton acceptance or reconstruction efficiency exists in the chain at
   all**, and only its geometric half can be bounded here; and the far-forward
   working point is unchosen. All three are the collaboration's and the
   detector groups' to supply, and all three are stated in the draft. Say also
   that the ε_det actually used is a **top-energy ⁷Li** number applied to a
   ⁶Li sample at 10 × 99.5, i.e. conservative by ×1.12–1.16 on the beam-energy
   axis alone — because the recipient will otherwise find that themselves and
   read the rest of the estimate in that light.
4. The **null test** is still the cheap thing worth doing: ⁶Li's a₂ is
   predicted to be **positive** where the deuteron's is negative. A sign is a
   one-bit measurement and needs far less than a 3 σ magnitude — and a 2.6 σ
   sample does deliver that bit.

### C2.8 Every assumption, in one list

Stated because each of them is optimistic, and the answer is only MARGINAL.
Items 3a, 3b, 3c and 10 were added on 2026-09-04 (third pass): the list
claimed completeness while two factors of comparable size and opposite
sign were missing from it entirely.

1. **Transverse (vertical) tensor-polarised ⁶Li**, θ_S = 90°. The signal
   carries sin²θ_S (`delta_perp_analytic_fm2`); a longitudinal axis gives
   **exactly zero**. `tensor_flip_plan`'s default axis is transverse.
2. **Uniform acceptance in φ_Δ**, and a |t| resolution small enough not to
   smear it. **No in-tree coherent-VM p_T resolution exists**, so no smearing
   dilution is applied; the formula is D = exp(−2σ_φ²) with σ_φ ≈ σ_pT/p_T,
   and at |t| ≈ 1/B the transverse momentum is only **0.141 GeV**, so this is
   the assumption most likely to cost something. The finite-bin dilution of
   `cos2phi_fit_err` is small by comparison (1.0115 at 24 bins, 1.1107 at 8).
3. **ε_det flat in |t|, in φ_Δ, AND IN Q²**, and equal to arXiv:2511.05638's
   ⁷Li coherent-J/ψ number for all three mesons. Its |t| dependence is real
   and is not modelled here. **The Q² extrapolation is new and it is the
   largest unmeasured quantity in the chain**: that paper's sample was
   generated with 0.1 < Q² < 100 GeV², so **no point of its 17.75 % was ever
   evaluated below Q² = 0.1, where 85 % of the J/ψ rate sits.** Size: the
   whole photoproduction gain is ×6.75 in rate and ×2.60 in S, so a
   hypothetical ε_det that collapsed to zero below 0.1 would return the
   headline to 1.005 σ and 89 fb⁻¹/u. Two things bound it and neither is a
   measurement — the efficiency is a *recoil-nucleus* tagging acceptance
   ("tagging efficiency × acceptance", their Figs. 2 and 5) and so does not
   require the scattered electron, which at Q² < 0.1 is down the beam pipe;
   and the recoil |t| spectrum is unchanged to 0.4 % (§C2.3a), with ⟨W⟩ moving
   in the direction their Fig. 2 says *raises* the efficiency.

3a. **… and ε_det is NOT flat in BEAM ENERGY, which this list did not say
   until 2026-09-04, in the one direction the write-up never hedged.**
   0.1775 is a ⁷Li number at **18 × 117.9 GeV/u — that nucleus's own top
   energy** — and the chain applies it to a ⁶Li sample at **10 × 99.5**. The
   same paper measures the dependence (§V.B, on ³He: 32.23 % / 54.38 % /
   99.77 % at 183 / 100 / 41 GeV/u), and it is steep and one-directional.
   **Direction: UP. Size: ×1.122–1.158 for the 117.9 → 99.5 step**, from that
   scan's own two local power laws, d ln ε/d ln E = −0.866 and −0.681
   (`chang26_he3_energy_scan()`, §C2.3c(a)). The ⁷Li → ⁶Li substitution
   in the same clause is a **separate and undetermined** factor:
   ~~a further ×1.16–1.22 at fixed rigidity~~ is **retracted** (a non
   sequitur — fixing the rigidity fixes neither E/u nor Z), and read off the
   four `chang26_species_efficiency()` entries that share a beam energy it is
   **×0.99–1.33, straddling 1** (§C2.3c(b), fourth pass) — stated, not
   applied by default. **It is not the paper's own 1.687**: that is their 183 → 100
   ratio, a step 3.56× larger in ln E, and carrying it whole would give
   3.40 σ / 7.8 fb⁻¹/u by using the wrong lever arm. The ⟨W⟩ shift
   43.2 → 30.2 points the same way and is the *same* effect on another axis,
   not an independent one.

3b. **There is NO decay-lepton acceptance and NO reconstruction efficiency
   anywhere in the chain — this list did not carry it at all.** The chain is
   σ × BR(ℓ⁺ℓ⁻) × ε_recoil, and ε_recoil counts scattered **nuclei**: *"32.23 %
   of the scattered ³He nuclei occur within a safe distance from the beam"*,
   under a simulation that by its own p. 4 *"only accounts for the acceptance
   effect and does not incorporate the efficiencies of the detector.
   Additionally, we did not account for the efficiency and acceptance of the
   reconstructed distribution"*. **Direction: DOWN. Size: the geometric half
   is bounded here and is not the problem — A_geom = 0.994 at ⟨W⟩ = 30.2 and
   never below 0.89 anywhere the beams reach**, computed from this tree's own
   `Scenario::eta_max` = 3.5 and the J/ψ's lab rapidity (§C2.3c(c)), because
   at 10 × 99.5 the J/ψ is nearly at rest. **The per-lepton reconstruction
   efficiency is UNBOUNDED HERE**: no tracking or PID efficiency for a decay
   lepton exists anywhere in this repository. The 0.95/track used to draw the
   band's low end is the sibling `../PolarizedLithiumSim`'s own
   `HfsModel(eff_track=0.95)`, which that file labels a stand-in. This is what
   leaves the verdict band open below.

3c. **Items 3a and 3b PARTLY CANCEL** — ×1.122 against ×0.897 is 1.007 — which
   is why neither may be quoted alone. A reader given only 3a moves the
   verdict up; a reader given only 3b moves it down; together they move it by
   0.7 %.
4. **The coherent sample is identified as coherent** (breakup vetoed) at no
   further cost beyond ε_det.
5. **σ and B from eSTARlight**: one spherically symmetric Gaussian density, so
   no α+d clustering, no polarisation axis, no diffractive minimum, no
   saturation (`estarlight_li6.md` §6). The density choice alone moves S by
   **1.66×** between the two defensible radii (2.618 → 1.576 σ).
6. **a₂ from `a2_from_quadrupole`** — §C2.0. Not an amplitude. The map is
   anchored on a **photoproduction** calculation, which the new Q² window
   suits better than the old one, at an x_P an order of magnitude away from
   ours.
7. **Both lepton channels are reconstructed** with the same efficiency as
   e⁺e⁻ alone. `branching_all` = 0.11932; if μ⁺μ⁻ were unusable the headline
   falls to 1.85 σ and 3 σ moves to 26 fb⁻¹/u — still inside the band.
8. **The far-forward working point is the Yellow Report one**,
   `Optics::lumi_fraction` = 1, **not** LiPolGen's own de-squeezed ⁶Li tagging
   point. `estarlight_li6.md` §2e instructs to multiply by that fraction and
   this estimate does not; the reason is in §C2.2 and it is that the two are
   alternatives, not multipliers. **Size: 0.0781 at 10 × 100, i.e. ×0.279 on
   every S and ×12.8 on every 3 σ luminosity — the headline becomes 0.73 σ at
   168 fb⁻¹/u, and the whole §C2.3c band becomes 0.73–0.92 σ at
   106–167 fb⁻¹/u.** This is the largest single unmade choice in the item, and
   it is **the one correction that on its own restores the original NO** — at
   *both* ends of the band, which is what makes that statement survive the
   corrected efficiency chain.
9. **⟨P_zz²⟩ = 0.81, the background-immune difference**, not the optimal 0.90
   (§C2.2). The tables print both; the headline is the smaller.
10. **The programme luminosity carries no energy dependence.** 10 fb⁻¹/u is a
   scenario number with *no Li source at all* (`needs_survey.md` §3.8), so no
   penalty is applied for running at 10 × 99.5 rather than at a top energy —
   and none can be applied from anything in this tree. Direction and size both
   **unknown here**. It is listed because items 3a/8 make the beam energy and
   the optics explicit and this is the third thing the beam choice touches.


---

## C3. Open item O2 — the `li6.adr.fit` R₂ node is **real**, and it decides nothing

O2 as written asks whether the smoothed fit's R₂ node near 1.07–1.1 fm is a fit
artefact, and states that the answer "decides whether `FitRescaled` or
`OverlapRaw` is the honest default." **The first half is answerable and the
answer is "real"; the second half is a false premise** and is retracted below.

### C3.1 Where the node is, on each source

Measured on the sampler's own post-normalization, post-phase-fix tables
(the global phase flips both waves, so it cannot move a node):

| source | R₂ node [fm] | R₀ node [fm] |
|---|---|---|
| `FitRescaled` (default) | **1.06484697577** | 1.797 |
| `FitRaw` | 1.06484697577 | 1.797 |
| `OverlapRaw` (raw `li6.ad` r-block) | **1.11891335406** | 1.860 |

The raw block also has one **extra, noise-induced** crossing near 0.53 fm — its
R₂ there is +1.54e−3 ± 2.48e−3, consistent with zero — so "the first sign
change" is not well defined on the raw data. T22c takes the last crossing below
the outer positive lobe.

### C3.2 Does it survive the raw data's own error bars? Yes, at 3.3 σ

Averaging the raw `li6.ad` R₂ over the fit's inner negative lobe
(r < 1.065 fm), weighted by the file's own 1 σ MC column:

    weighted mean R₂ = −1.86947e−3 ± 5.62461e−4   →  −3.32 σ from zero
    restricted to 0.6–1.06 fm            = −2.16027e−3 ± 6.00601e−4  →  −3.60 σ

A pure fit artefact would sit at zero. It does not. Point by point, the raw R₂
in the lobe is (value ± σ, e−3): 0.45 → −5.38 ± 2.97, 0.65 → −3.23 ± 1.91,
0.75 → −4.71 ± 1.58, 0.85 → −3.20 ± 1.40, 1.05 → −1.64 ± 1.07 — five of the
seven points inside 0.4–1.06 fm are negative, three of them past 1.5 σ
individually.

Resampling the whole 0.55–1.75 fm window against the file's σ's (200 k draws),
the raw node position is **1.119 fm, 68 % CL [1.089, 1.153], 95 % CL
[0.983, 1.167]**. The fit's node at 1.0648 fm lies **1.6 σ inside** the raw
node's own error — a mild tension in *position*, none at all in *existence*.

**O2's physics question is closed: the node is a real feature of the AV18+UX
α–d overlap, not an artefact of the Forest et al. smoothing.**

### C3.3 …and it is irrelevant to every observable this module reports

The r⁴ weight of (G4) kills it, as the design guessed. Measured share of the
region r < 1.5 fm, on every source:

| source | share of q_int | share of the D-wave norm n₂ | share of the total norm |
|---|---|---|---|
| `FitRescaled` | **−0.092 %** | 0.500 % | 11.73 % |
| `OverlapRaw` | −0.085 % | 0.422 % | 11.24 % |
| `FitRaw` | −0.092 % | 0.500 % | 11.67 % |

Moving the node by the full 0.054 fm that separates the two tables cannot move
Q by more than a small fraction of 0.09 %. **The node does not decide the
default.** T22c pins both halves of this.

### C3.4 What actually separates the two sources

Not the node — the 2–9 fm shape, and it is a real, statistically significant
difference:

| | `FitRescaled` | `OverlapRaw` | gap |
|---|---|---|---|
| P_D(α–d) | 0.02011 | 0.02011 | — (the rescale matches norms by construction) |
| ⟨R²⟩_αd [fm²] | 16.9633 | 17.5482 | +3.4 % |
| q_int + q_dd [fm²] | −1.3204 | −1.4009 | +6.1 % |
| Q_charge [fm²] | −0.615448 | −0.669114 | **+8.7 %** |
| η | −0.0482161 | −0.0537547 | +11.5 % |
| `eps_b0_equivalent()` | −0.0506 | −0.0550 | +8.7 % |
| `a2_from_geometry(0.3, 1)` | 0.1976 | 0.2151 | +8.9 % |

Of the 0.5917 fm² ⟨R²⟩ gap, only **16 %** comes from r > 9.95 fm where the fit
table simply ends; **84 % is a shape difference inside the fit's own domain**.

And in the η window the smoothing is *not* a null operation. Propagating
`li6.ad`'s own MC errors through `asymptotic_ds_ratio()` (200 k draws):

    η(OverlapRaw) = −0.0537566 ± 0.0020800 (stat)
    η(FitRescaled) = −0.0482161            →  −2.66 σ of the raw's own error

The 6–8 fm window is **not** Monte Carlo noise, contrary to the wording that
was in `cluster_config.hpp` (corrected in this pass). Median (min) S/N of the
raw block's own columns:

| window [fm] | R₀ | R₂ |
|---|---|---|
| 6.0–8.0 (the η window) | **37.7 (24.9)** | **13.5 (8.9)** |
| 8.0–9.95 | 15.5 (9.0) | 6.8 (3.0) |
| 9.95–12.0 | 5.1 (2.8) | 2.8 (1.4) |

So the η window is well measured — it is the region *past* the fit's end that
is noise-dominated — and the fit and the raw genuinely disagree inside it at
−2.7 σ.

### C3.5 The decision: **`FitRescaled` stays the default — an author decision**

Recorded as an author decision, not a physics result, with the cost of each
choice measured:

**Why not `OverlapRaw`, given it is the direct measurement.**

1. **Neither is meaningfully closer to the one measured quantity in play.**
   Against GK's η = −0.025(12): `FitRescaled` is **−2.0 σ**, `OverlapRaw`
   **−2.4 σ**. The 2.7 σ by which they differ from *each other* is on the raw
   MC error, which is **5.6× smaller** (0.00208 against 0.011662) than the
   experimental error they are both being judged against. There is no measurement that prefers one.
2. **The gap is invisible under the systematic the module already carries.**
   The default choice moves Q by 8.7 %. `quadrupole_band_fm2()` already ships a
   factor-**7.5** band on that same number. An 8.7 % choice inside a 750 % band
   is not a physics decision.
3. **The raw block is genuinely rough at the level of its own error bars,
   everywhere.** Mean |second difference| in units of the point's own σ:

   | window [fm] | raw R₀ | raw R₂ | fit R₀ | fit R₂ |
   |---|---|---|---|---|
   | 2.0–6.0 | 0.86 | 0.64 | 0.17 | 0.06 |
   | 6.0–9.95 | 0.72 | 1.02 | 0.05 | 0.02 |
   | 9.95–15.0 | 0.74 | 0.88 | — | — |

   ≈1 is exactly what uncorrelated MC scatter gives. A **sampler** draws from
   these tables cell by cell; sampling `OverlapRaw` reproduces that scatter as
   structure in the density. The smoothed table is the better numerical object
   for the job, and `FitRescaled` already takes its *normalization* (hence
   P_D) from the raw block, so the raw data's one unambiguous content is kept.
4. **`OverlapRaw` extends to 19.95 fm**, where R₀ has S/N < 3, contributing
   0.63 % of the ⟨R²⟩ numerator from data that is noise.

**Cost of keeping `FitRescaled`, stated plainly.** A −2.7 σ (on MC error)
shape bias in the 6–8 fm Whittaker window relative to the raw data, which
propagates to an 8.7 % bias in Q_charge and 11.5 % in η — **in the direction
that makes the model's D-wave excess look smaller than the raw overlap says
it is**. Every §2.7/§8 number quoted for the default therefore sits at the
*low* end of the α–d source range, and the range −0.615…−0.730 fm² already
brackets it. `OverlapRaw` remains reachable as
`--alpha-d-source overlap-raw` / `AlphaDSource::OverlapRaw` and is pinned in
the §8 three-source table, so nothing is hidden.

**Cost of switching**, had it been taken: every `FitRescaled` row of design G
§8 and §2.7, T17, T22, T22b, T22c, the `USAGE.md` table, the
`PHYSICS_CHANNELS.md` row and the pytest anchors move. No
`validation/reference/*.json` depends on `cluster_config` (checked), so the
rtol 1e−12 gates would **not** have moved — the cost is entirely in pinned
documentation, not in the physics gates.

---

## C4. Open item O4 — ΔB is now defined, and `eps_b0` = −0.08 is a **deuteron** number

**Two findings, one of them load-bearing.** (1) The symbol ΔB is defined, once,
at the owning site, and the pre-existing label "relative slope modulation
ΔB₀/B" was off **by a sign**, not by a factor 2. (2) The shipped
`CoherentScenario::eps_b0` = −0.08 implies a ⁶Li charge quadrupole of
**−0.9345 fm²**, i.e. **11.42×** the measured −0.0818 fm² and 1.52× even the
α+d model's −0.615 — so every *generated* coherent tensor number in this
library is 11.4× the measured-quadrupole expectation. The default is not
changed (it is pinned bit for bit in `validation/reference/`); the consequence
is written down instead, in code as well as in prose.

### C4.1 The definition, once

`include/lipolgen/coherent.hpp` used ΔB in three docstrings and defined it
nowhere; `design_G_cluster_config.md` §2.8 listed three defensible readings
under which the docstring is wrong by −2, by −1, or not at all. The convention
is now fixed at the `eps_b0` declaration and reproduced nowhere else:

    |F_m(|t|, Φ)|² = exp(−|t| [B + ΔB_m cos 2(Φ − Φ_S)])

with Φ the momentum-transfer azimuth, Φ_S the spin azimuth and B = `slope_b`
the φ-averaged slope. For a Gaussian transverse profile that is exactly
ΔB_m = δ_m/2 with δ_m = ⟨x²⟩ − ⟨y²⟩ per nucleon in state m — design (G7) — and
first order against the anchor's 1 + 2a₂cos2Φ normalisation gives

    a₂(m) = −(ΔB_m/2)|t| = −(δ_m/4)|t| ,   δ₀ = −2δ_{±1}

which **is** `a2_m_state`, term by term. Hence

| relation | value |
|---|---|
| `eps_b0` ≡ δ_{±1}/B | the definition |
| `eps_b0` = +2 ΔB_{±1}/B | ΔB_{±1} = ½ `eps_b0` B = **−2.0 GeV⁻²** at the default |
| `eps_b0` = **−1** × ΔB₀/B | ΔB₀ = −`eps_b0` B = **+4.0 GeV⁻²** |

**So the old label was wrong by a sign and the code was right.** Reading (A) of
the design's table (the one the first draft asserted, factor −2) conflates two
conventions at once — it reads a₂ as the coefficient of cos 2Φ in the yield
rather than half of it — and reading (B) is the docstring's own literal, which
is right in magnitude and missing the minus sign. Reading (C) is the one above.

In code: `CoherentScenario::delta_b_m(m)` is the single home of ΔB_m and
`slope_at_azimuth(φ, m)` is the slope itself; **T10b** checks
a₂(m) ≡ −½ ΔB_m |t| to 1e−15 at every m and |t|, checks both relations above,
and checks that expanding exp(−|t| B(Φ)) really reproduces 2a₂ to 2 % at
|t| = 0.05 (the residual is the cubic term).

### C4.2 What `eps_b0` = −0.08 assumes about ⁶Li

`quadrupole_from_a2_slope` (new, `cluster_config.hpp`) is the exact algebraic
inverse of `a2_from_quadrupole` at m = ±1, so the scenario can be **asked**
rather than argued with. At the shipped `eps_b0` = −0.08, `slope_b` = 50:

| quantity | value |
|---|---|
| δ_{±1} = `eps_b0`·`slope_b` | **−4.0 GeV⁻²** |
| implied Q_matter(⁶Li) | −1.8690781872 fm² |
| implied **Q_charge(⁶Li)** | **−0.9345390936 fm²** |
| ÷ measured `LI6_QUADRUPOLE_FM2` = −0.0818 | **11.4246832958×** |
| ÷ the α+d model's −0.6154 | 1.5185× |
| a₂(±1, \|t\| = 0.3) | **+0.3000** |

That last row is the sanity check that settles it: the **deuteron's** own
digitized a₂(±1) at |t| = 0.3 is −0.28 (`mantysaari_a2_deuteron()`), and
Q_d = +0.2859 fm² is 3.5× *larger* in magnitude than Q(⁶Li). A ⁶Li scenario
that predicts a bigger azimuthal modulation than the deuteron is not a ⁶Li
scenario.

### C4.3 The band, re-derived for ⁶Li

`quadrupole_band_fm2()` through `a2_from_quadrupole` at the scenario's own
B = 50 GeV⁻² (nothing retyped; the values at B = 52.04, the *measured* point
radius `eps_b0_equivalent()` uses, are 4.1 % smaller in magnitude):

| Q_charge assumed | δ_{±1} [GeV⁻²] | **eps_b0 at B = 50** | eps_b0 at B = 52.04 | a₂(±1, 0.3) |
|---|---|---|---|---|
| measured −0.0818 | −0.350119 | **−0.0070023** | −0.0067278 | +0.026259 |
| GFMC AV18+IL7 −0.20(6) | −0.856037 | **−0.0171207** | −0.0164506 | +0.064203 |
|  … +1σ −0.14 | −0.599226 | −0.0119845 | −0.0115154 | +0.044942 |
|  … −1σ −0.26 | −1.112848 | −0.0222570 | −0.0213861 | +0.083464 |
| α+d model −0.61545 (`FitRescaled`) | −2.634232 | **−0.0526846** | −0.0506225 | +0.197567 |
| **shipped `eps_b0` = −0.08** | **−4.0** | **−0.08** | — | **+0.3000** |

**The honest ⁶Li band is −(0.0070 … 0.0527)** — it is §C1's factor-7.5
quadrupole budget and nothing else, no wider and no narrower — and the shipped
default sits **1.52× above the top of it**.

**The old documented band −(0.04 … 0.13) is a deuteron band**: it was set from
the deuteron's own eps_b0 = +0.105 at the deuteron's own
B_d = `gaussian_slope(1.967 fm)` = 33.1 GeV⁻², at the **opposite sign**, and
never rescaled. Said exactly, because "does not overlap" would be wrong: the
old band contains the α+d model row 0.0527 and **nothing else of ⁶Li** — the
GFMC 0.0171 and the measured 0.0070 are both below its floor of 0.04.

### C4.4 `eps_b0` and `slope_b` are **not** independent

Every observable uses the product δ_{±1} = `eps_b0`·`slope_b`
(`a2_deformation` is −(P_zz/4)·`eps_b0`·B·|t|), and δ is the physics: a
transverse second-moment difference fixed by the target's quadrupole, not by
how steep its form factor is. So a `slope_b` scan over its documented band
{40, 60} **at fixed `eps_b0`** silently scans the target's quadrupole too:

| `slope_b` | δ_{±1} [GeV⁻²] | implied Q_charge [fm²] | a₂(±1, 0.3) |
|---|---|---|---|
| 40 | −3.20 | −0.7476 | +0.2400 |
| 50 | −4.00 | −0.9345 | +0.3000 |
| 60 | −4.80 | −1.1214 | +0.3600 |

±20 % on a₂ for no physical reason. **Band δ, or band `eps_b0` and `slope_b`
together and anticorrelated; never `eps_b0` alone.** This is now stated at the
declaration.

### C4.5 The two consequences, priced

1. **The reach.** §C2's significance is linear in κ = a₂(±1)/|t|. The same
   J/ψ sample that gives **2.62 σ** at the measured quadrupole
   (κ = 0.08752978 GeV⁻², §C2.4's headline row) gives **29.9 σ** at the
   scenario's own κ = 1.0 GeV⁻². Had §C2 been run on `eps_b0` instead of on
   `a2_from_quadrupole`, it would have answered O5 the wrong way round by a
   factor 11.4 — an overstatement of the reach, where the Q²-window and
   lepton-channel errors §C2 shipped with were an understatement of it. The
   two are independent and both had to be fixed.
2. **`COHERENT_T_MAX_DEFAULT` = 0.2 is a consequence of the oversized
   `eps_b0`, not of the target.** The positivity edge |c₂| = 1 at P_zz = −2
   sits at |t| = (2·`amp` − 1)/(`eps_b0`·B):

   | `eps_b0` | |t| at \|c₂\| = 1 |
   |---|---|
   | −0.08 (shipped) | **0.245 GeV²** |
   | −0.0527 (α+d model) | 0.372 |
   | −0.0171 (GFMC) | 1.146 |
   | −0.0070 (measured Q) | **2.800** |

   At the measured quadrupole the azimuthal weight would stay a density out to
   |t| = 2.8 GeV², i.e. the |t| ceiling would be set by the digitisation range
   of the anchor and by nothing else.

### C4.6 The decision

**AUTHOR DECISION, recorded, not taken silently: the default stays −0.08 and
the band documentation is replaced.** The reasons, in order:

* `eps_b0` is pinned at rtol 1e−12 against `validation/reference/coherent.json`
  (`tests/test_coherent.cpp:57`, the `scenario_defaults` block) and this run may not move a reference gate.
  Changing it would move every generated coherent tensor number and every
  reference file that records one.
* The number that *should* replace it is not a single number: the ⁶Li band is a
  factor 7.5 wide (§C1) and its own 1σ on η spans a **sign change**, so
  −0.0070 would be as much a point-on-a-line as −0.08 is.
* What was actually missing was not a better default but the **statement of
  what the default assumes**, which now exists in code
  (`quadrupole_from_a2_slope`), in the declaration, and in **T10b**.

**Cost of that decision, stated:** every generated coherent tensor number —
`a2_deformation`, `cos2phi_coefficient`, the sampled azimuth, and any
pseudo-experiment built on them — is **11.4× the measured-quadrupole
expectation**, and `COHERENT_T_MAX_DEFAULT` = 0.2 inherits that. **Nothing
published from this channel may quote a single `eps_b0` row**: band it, and say
which quadrupole the row assumes.

*What would change the decision:* a dipole-model calculation of the ⁶Li
coherent amplitude (§C2's standing limitation), or a decision to re-pin the
reference files — the second is a one-line change plus a reference dump, and it
is the user's call because it moves published numbers.

---

## C5. The VMC inputs

Five threads, all offline. Four hold up and are answered; **C5.5's framing is
half wrong and the correction matters more than the original question**.

### C5.1 The α–d asymptotic D/S ratio: **no second dial** — a converter instead

**The question:** η = −0.0482 (model) against the measured
`LI6_ETA_DS_GK` = −0.025 ± 0.006 ± 0.010, with no dial that constrains it
directly. Does η deserve a first-class dial?

**No, and the reason is arithmetic, not taste.** §C1 measured that η is
**exactly linear** in `quadrupole_dial_s()` — the (G9) dial multiplies R₂ by s
and both waves by the same 1/√n(s), so the wave *ratio* carries s alone. An
`eta_target` option would therefore be a second name for the **same
one-parameter family**: two knobs onto one physics number, which
`docs/CONVENTIONS.md` forbids, and two knobs that could be set to
contradictory values in one options object with no way to decide which wins.

**What was actually missing was the conversion**, which existed only as a typed
number in a document (§C1.4's −0.18557 fm² at η = −0.025). It is now
`ClusterConfigSampler::quadrupole_for_eta(eta_target)`: the `quadrupole_target_fm2`
to request so that `asymptotic_ds_ratio()` lands on the target, computed from
the sampler's own **pre-dial** moments and its current dial
(s = dial_s × η_target/η_now), so it is exact on a dialled sampler too, and it
**throws** rather than clipping when the requested η is unreachable.

Re-derived from the tables, against §C1's typed values:

| η | `quadrupole_for_eta` [fm²] | §C1.3's typed value |
|---|---|---|
| −0.025 − σ_comb = −0.0366619038 | **−0.4004752268** | −0.4004752268 ✓ |
| −0.025 (GK central) | **−0.1842160147** | −0.1842160147 ✓ |
| −0.025 + σ_comb = −0.0133380962 | **+0.0297575059** | +0.0297575059 ✓ |

so the George–Knutson band — including the fact that its near edge is the
**wrong sign** — now re-runs from the data instead of rotting. `quadrupole_for_eta`
is not a default and changes nothing; **T23** pins it, the round trip through
the dial, and its exactness on an already-dialled sampler.

*Cost of the decision:* reaching a target η still costs a bisection inside the
dial (the converter returns a Q, the dial then bisects to it), and η remains
what §C1 measured it to be — a relabelling of the dial, not an independent
constraint on it. Anything derived from the central η alone is a point on a
line whose error bar does not pin its sign.

### C5.2 The Monte-Carlo errors are no longer thrown away

**The premise was exactly right.** `read_anl_momentum` / `read_anl_overlap`
parsed `AnlTable::err` and both factories (`vmc_from_momentum`,
`vmc_from_overlap_k`) dropped it on the floor, so `ClusterWaveSource::VmcAV18`
had no error band of any kind.

**What was added.** `VmcRadial` now carries `dpsi()` — the tables' own printed
1σ, **signed**, so it transforms with ψ under the global-phase flip the pair is
fixed with. For a momentum file ψ = s√ρ, hence dψ = s·dρ/(2√ρ) and 0 wherever
ρ ≤ 0. `shifted_by_sigma(n)` is ψ → ψ + n·dψ, `norm2_error()` reports the band
on ∫k²ψ²dk, and `li6_alpha_channel(..., vmc_mc_sigma)` /
`li7_alpha_channel(..., vmc_mc_sigma)` /
`PipelineConfig::cluster_vmc_mc_sigma` are the run-level knob. **0 is the
default and returns the cached tables bit for bit**; a non-zero band on a
family that carries no errors (Hulthén; the deuteron control's `fdeut.av18`,
which prints none) is **refused**, by the same rule the b₁ band scales are
refused under — a knob that did not run may not be recorded as if it had.

**It is not an independent-point error, and the report says so.** Every point
of an ANL table comes from one variational walk, so the point-to-point
correlation is unknown and is neither 0 nor 1. The knob takes the **fully
correlated** reading (the conservative envelope on any smooth functional);
`norm2_error()` returns both limits.

**The measured band** (`li6_ad1.momentum`, 1M samples):

| quantity | correlated 1σ | in quadrature |
|---|---|---|
| ∫k²ψ₀²dk (S wave) | **0.4705 %** | 0.1085 % |
| ∫k²ψ₂²dk (D wave) | **1.6176 %** | 0.3974 % |
| `li7_at3.momentum` P wave | 0.6129 % | 0.1343 % |

| n_σ (correlated) | P_D | Δ P_D | `tensor_dilution` | Δ |
|---|---|---|---|---|
| −2 | 0.018927265 | −2.210 % | 0.982960823 | +0.0392 % |
| −1 | 0.019139148 | −1.115 % | 0.982770077 | +0.0198 % |
| **0** | **0.019354933** | — | **0.982575817** | — |
| +1 | 0.019574536 | +1.135 % | 0.982378121 | −0.0201 % |
| +2 | 0.019797870 | +2.289 % | 0.982177065 | −0.0406 % |

Shapes move too: ⟨k⟩ per σ is +0.26 % (⁶Li S), +0.51 % (⁶Li D), +0.18 % (⁷Li P).

**The number that matters is the last column.** ±1σ of ANL statistics moves the
tagged tensor observable by **0.02 %**, against **6.6 %** for the Hulthén → VMC
wave-function change (§C5.5) and a **factor 7.5** for the quadrupole (§C1). The
Monte-Carlo error is now reportable, propagatable and **negligible**; recording
that it is negligible is the result, and it is the first time the band could be
computed at all. Pinned in **T24**.

### C5.3 N_αd 5 %, P_D 7 % — three numbers spanning three different amounts

**The premise needs a correction before it can be propagated.** The three
readings are

| reading | N_αd | N_D | P_D = N_D/N_αd |
|---|---|---|---|
| 2014 `li6_ad1.momentum` (**default**) | 0.819481 | 0.015861 | **0.0193549** |
| 2004 `li6.ad`, as printed | 0.855 (file prints `ndx` 0.856) | 0.017 | 0.0198830 |
| Wiringa et al., PRC **89** (2014) 024305 | 0.863 | 0.017 | 0.0196987 |
| **spread** | **5.311 %** | **7.181 %** | **2.729 %** |

so "N_αd spans 5 % and P_D 7 %" is quoting two of them about a third: the
**7.2 % belongs to the D-wave norm**, and P_D — the *ratio* every tagged
observable actually uses — spans only **2.7 %**.

**Propagated.** N_αd cancels exactly out of the tagged sector: `TaggedModel`
renormalises each wave to its own P_L, so only the ratio survives. What the
three readings do to the tagged observables:

| reading | P_D | `vector_dilution` | `tensor_dilution` |
|---|---|---|---|
| 2014 (default) | 0.0193549 | 0.970965994 | 0.982575817 |
| 2004 `li6.ad` | 0.0198830 | 0.970173789 | 0.982100390 |
| Wiringa 2014 | 0.0196987 | 0.970450277 | 0.982266320 |
| **spread** | 2.729 % | **0.0816 %** | **0.0484 %** |

**So the 5 % propagates to 0.05 % on the tagged tensor observable and to
nothing else outside b₁**, where it is exactly linear and therefore exactly
±5 % (`docs/open_items/run_2026-09-02/phase_D_numbers.md` § "The ±5 % N_αd
systematic (design §2.1) and the Q4 knob", already measured — the knob is
**Q4** and that file has no §Q3; this citation said Q3 until 2026-09-04). It is a real systematic in
one place and a rounding error in the other, and the honest statement carries
both. Pinned in **T26**.

### C5.4 The deuteron control: AV18 is available, and it is now opt-in

**The premise's blocker was false.** `pipeline.hpp` said the control "is always
Hulthén — there is no VMC d → p+n cluster table, the deuteron IS the cluster".
That is true of the `momenta/` files and **false of
`data/vmc/deuteron/fdeut.av18`**, whose u(k) and w(k) *are* the p–n relative S
and D waves and which the tree already reads (`read_fdeut_k`, used by the b₁
convolution). So `--cluster-wave vmc` was being **silently ignored** on that
channel.

`deuteron_channel(beta, p_d, source)` now takes `ClusterWaveSource`, and
`PipelineConfig::cluster_wave` reaches it. Measured:

| | Hulthén (default) | AV18 `fdeut` | change |
|---|---|---|---|
| P_D | 0.045 (`P_D_DEUTERON`, scenario) | **0.0575998919874** | **+28.00 %** |
| `vector_dilution` | 0.932494769 | 0.913594777 | **−2.027 %** |
| `tensor_dilution` | 0.959488074 | 0.948145618 | **−1.182 %** |

**Two validations came free.** (1) The k-block's own D fraction reproduces the
file's *r-space* header `dstate` = 0.057599 to **1.55e−5** relative — the
reader and the units convention checked by the file against itself. That
residual is the file's own r↔k quadrature spread and **not** print rounding:
0.05759989 would print as 0.057600 at six figures, not as 0.057599. (2) **The
relative S–D sign does not flip.** CDKS fix φ_L = i^L ψ_L with φ₂ = −W and
U, W ≥ 0 at low k, so the physical deuteron has ψ₂ = +W > 0 — exactly what the
positive-definite Hulthén forms already assume. This is the **opposite** of the
⁶Li α–d case, where the VMC overlap flips the sign below the S node. A first
draft of this section claimed the sign flips here too; it does not.

**AUTHOR DECISION: Hulthén stays the default, AV18 is opt-in.** Reasons: the
default is bit-for-bit pinned across `validation/reference/tagged.json`; the
control channel's job is to be the *Cosyn–Weiss tagged limit of this
machinery*, which is the analytic pair the ⁶Li Hulthén channel is built from,
so switching only the control would make the control and the channel it
controls two different wave-function families; and `fdeut.av18` prints no MC
errors, so the AV18 row carries no band of its own.

**AND THAT RATIONALE APPLIES TO THE OPT-IN PATH THIS SECTION SHIPPED, WHICH IS
WHAT §C5.5b HAD TO REPAIR.** "Switching only one thing would make it two
wave-function families" is precisely what `--cluster-wave vmc` then did one
level down: it made the **control's** deuteron AV18 while the ⁶Li α channel's
**embedded** deuteron stayed on the 0.045 scenario, so in one programme the
same flag meant "AV18 deuteron, P_D = 0.0576" on `--channel tagged-deuteron`
and "scenario deuteron, P_D = 0.045" on `--channel tagged-alpha`. It is fixed
in §C5.5b; the decision above (Hulthén is the *default*) is unchanged. **Cost, stated:** the
shipped control runs at a P_D that is 28 % below the AV18 deuteron's, worth
2.0 % on its vector dilution and 1.2 % on its tensor one — and, through
`P_D_DEUTERON`, the same 28 % sits inside `LI6_CLUSTER_POLARIZATION` (§C5.5).
Pinned in **T25**.

### C5.5 The inclusive-vs-tagged P_D split — the framing is half wrong, and the fix is worse than the defect

**First, the correction.** The item says "the same nucleus is described by two
different wave functions in **two channels of one run**". It is not one run: a
`Pipeline` builds **exactly one** channel (the `is_tagged` branch), a tagged run
never constructs the inclusive kernel and never reads `LI6().eff_pol_*`, and an
inclusive run never builds a `TaggedChannel`. So the *inclusive-vs-tagged*
drift is between **a tagged row and an inclusive row of one programme** — which
is worse, not better, because nothing in the code is in a position to notice
it.

> **CORRECTION, 2026-09-04 (§C5.5b).** The sentence this paragraph originally
> ended with — "so nothing *inside* one run is inconsistent" — was **false**,
> and it was shipped to five sites. A `--cluster-wave vmc` tagged-α run held
> **two deuteron wave-function families**: the ANL VMC AV18+UX α–d overlap for
> the relative motion, and the scenario Hulthén deuteron (P_D = 0.045) for the
> embedded deuteron, in both places a run reads it. Measured cost **+2.069 %**
> on every polarized tagged-α observable. §C5.5b has the numbers and the
> repair; what survives here — the 11.61 % between a tagged row and an
> inclusive row — is unaffected by it.

**Second, the size, measured.**

| | shipped | with the VMC α–d wave | drift |
|---|---|---|---|
| α–d vector dilution 1 − 1.5 P_D | 0.869950 (P_D = 0.0867) | 0.970966 (P_D = 0.0193549) | **+11.61 %** |
| `LI6_CLUSTER_POLARIZATION` | **0.811228375** | 0.905427287 | **+11.61 %** |
| rank-2 `LI6_B1_RANK2_TRANSFER` | **0.921947** | 0.982576 | **+6.58 %** |

**On the Hulthén default the CONVENTIONS.md claim is true and now measured**:
the tagged model's `vector_dilution()` = 0.869939 against the inclusive
constant's α–d factor 0.869950 — **1.22e−5**, a grid quadrature against a closed
form. Under `--cluster-wave vmc` the tagged side moves and the inclusive
constant does not.

**Third — and this is the finding — substitution is not the fix.** The one
*ab-initio* number for this very observable is the six-body VMC of Wiringa
et al., PRC **89** (2014) 024305 Table I: **0.848** (quoted in `beams.cpp`'s
`LI6()` since 2026-08-29, and now named `LI6_POLARIZATION_VMC_SIX_BODY`).

| reading of the cluster product | value | vs ab initio 0.848 |
|---|---|---|
| shipped (P_D = 0.0867, 0.045) | **0.811228** | **−4.34 %** |
| α–d from VMC (0.0193549, 0.045) | 0.905427 | **+6.77 %** |
| α–d from VMC + AV18 deuteron (0.0193549, 0.0576) | 0.887076 | +4.61 % |

**Making the inclusive constant consistent with the tagged VMC wave function
would move it AWAY from the only ab-initio anchor there is** — from 4.3 % below
to 6.8 % above. The product (1 − 1.5 P_D^{αd})(1 − 1.5 P_D^{d}) spans
**0.811 … 0.905** depending on which wave functions it is fed, the ab-initio
answer sits inside that span, and the **formula** carries the error, not the
choice of P_D. `P_D_LI6` = 0.0867 was never an α–d D-state probability — its own
docstring says it is chosen so that 1 − 1.5 P_D reproduces the 0.87 of
`b1_li6_from_deuteron` — so it is a *fitted effective depolarisation*, and
replacing a fitted effective parameter by a microscopic one inside a formula
the microscopic parameter does not obey is how a self-consistent answer becomes
a worse answer.

**AUTHOR DECISION: the default stays 0.811228, the drift is documented at
11.61 % and is not "closed".** What was added instead:

* `li6_cluster_polarization(p_d_alpha_d, p_d_deuteron)` (beams.hpp) — the ONE
  home of the product, so the alternatives can be written down without
  retyping (1 − 1.5 P). `LI6_CLUSTER_POLARIZATION` is bit for bit unchanged.
* `LI6_POLARIZATION_VMC_SIX_BODY` = 0.848 — the ab-initio anchor, named so the
  band can be checked instead of read out of a comment.
* `LI6_CLUSTER_POLARIZATION_VMC` (tagged.hpp) — the "consistent" reading, as a
  **diagnostic**, with its docstring saying in terms that it must never be
  substituted for the shipped one on the strength of being more self-consistent.
* The false half of the claim corrected at **every** site that carried it:
  `beams.hpp`, `tagged.hpp`, `docs/CONVENTIONS.md`, `docs/PHYSICS_CHANNELS.md`.

**Cost of the decision, stated:** a `--cluster-wave vmc` tagged row and an
inclusive row of the same programme describe ⁶Li with wave functions whose
vector dilutions differ by 11.61 % and whose rank-2 transfers differ by 6.58 %,
and nothing in the code refuses the combination (it cannot: they are different
runs). **No inclusive ⁶Li polarization may be quoted without the band
0.81 … 0.91.**

*The alternative, priced:* wiring an opt-in
`PipelineConfig` override that feeds `LI6_CLUSTER_POLARIZATION_VMC` into
`LI6().eff_pol_*` and the VMC `tensor_dilution` into `Li6B1`'s transfer is
about a day's work and would move every inclusive g₁ number by +11.6 % and
every inclusive b₁ number by +6.6 % — **away** from 0.848. It is not done, and
it should not be done until the cluster product is replaced by something that
reproduces 0.848 from a wave function rather than from a fit.

### C5.5b The claim §C5.5 shipped was false, and §C5.4 is what made it false

**§C5.5's repair introduced a new consistency claim** — *"a single run builds
exactly ONE channel, so nothing inside one run is inconsistent"* — at five
sites: `docs/CONVENTIONS.md`, `docs/PHYSICS_CHANNELS.md`, `docs/USAGE.md`,
`include/lipolgen/pipeline.hpp` (on the `cluster_wave` field itself, where a
user setting the flag reads it) and §C5.5 above. It was false, and **§C5.4's
own opt-in path is what made it false**.

#### C5.5b.1 The counterexample, measured

`src/core/tagged.cpp` set `c.dis_target = DEUTERON()` **unconditionally**,
outside the `VmcAV18` branch. So one `--channel tagged-alpha --cluster-wave
vmc` run used, simultaneously:

* the ANL VMC AV18+UX α–d overlap, P_D = 0.0193549, for the **relative
  motion**; and
* the scenario Hulthén deuteron, `P_D_DEUTERON` = 0.045 → eff_pol **0.9325**,
  for the embedded deuteron's **spin transfer to the struck nucleon**
  (`InclusiveKernel::g1a` → `PolSF::g1_nucleus`, on `channel.dis_target`).

§C5.4 measured that the AV18 deuteron belonging to that overlap has
P_D = 0.0575998919874, i.e. eff_pol **0.9136001620**. Reproduced here before
acting on it:

| | eff_pol of the embedded deuteron | source |
|---|---|---|
| what a `vmc` run used | **0.9325** | `P_D_DEUTERON` = 0.045, a *chosen* number |
| the AV18 deuteron of that overlap | **0.9136001620** | `deuteron_av18_p_d()` = 0.0575998919874 |
| ratio | **1.020687209533** | **+2.0687 % too high** |

**And it is exact, not approximate.** g₁A = Z·P_p·g₁p + N·P_n·g₁n is *linear*
in the effective polarization, so the whole polarized sector carried one common
factor. Measured on the kernel the pipeline actually builds
(`struck_cluster_kernel` → `InclusiveKernel(channel.dis_target)`), at four
(x, Q²) points and on both g₁ and A_∥:

| x | Q² | A_∥ (AV18) | A_∥ (as shipped) | ratio |
|---|---|---|---|---|
| 0.05 | 2 | 1.2022375054e−02 | 1.2271084446e−02 | 1.020687209533 |
| 0.10 | 5 | 3.8423691942e−02 | 3.9218570908e−02 | 1.020687209533 |
| 0.30 | 10 | 7.1440366345e−02 | 7.2918268173e−02 | 1.020687209533 |
| 0.50 | 20 | 2.7952367295e−01 | 2.8530623774e−01 | 1.020687209533 |

#### C5.5b.2 A SECOND site, not in the item as posed, found while checking it

`ClusterBreakup` — the **T1 tier, which is the default on tagged channels** —
carries a deuteron wave function of its own: it draws the struck nucleon's
momentum *and its spin projection* from `deuteron_channel(beta, p_d)`, at
`p_d` = `P_D_DEUTERON` = 0.045. `Pipeline` already forwarded `cluster_beta`
into it, with the comment *"the breakup's internal wave function … MUST be the
ones the rest of the run is made with, or the T1 tier would describe a
different nucleus from the T0 one"* — and did **not** forward the wave-function
family. So the same +2.0688 % sat in the RECORD as well as in the rate:

| | dilution the T1 spin draw implies | |
|---|---|---|
| Hulthén deuteron, P_D = 0.045 | 0.932494769 | |
| AV18 `fdeut` deuteron | 0.913594777 | **+2.0688 %** apart |

Note what this means for the record: the deuteron branch of `ClusterBreakup`
does **not** read `eff_pol_*` at all — the effective polarization *comes out of*
the CG sampling (`breakup.hpp`'s own `⟨2m₁⟩ = (1 − 1.5 P_D) m_S` identity), so
fixing only `dis_target` would have left the sampled spin labels and the rate
2.07 % apart **inside one event**.

#### C5.5b.3 The decision, and what it cost

**AUTHOR DECISION: the embedded deuteron follows the flag.** Of the three
options — follow, refuse the combination, or keep it and document it — *follow*
is the one §C5.4's own rationale requires (*"switching only the control would
make the control and the channel it controls two different wave-function
families"*), and it is the one that makes the flag mean **one wave function
everywhere it is read**. Refusing the combination would delete the α–d overlap,
which is the feature. Documenting a 2 % known-wrong number is not a repair.

What changed:

* `DEUTERON_AV18()` (tagged.hpp/tagged.cpp) — `DEUTERON()` with
  `eff_pol_{p,n}` = `vector_dilution_of(deuteron_av18_p_d())` = 0.913600 and
  **nothing else** different (same name, A, Z, spin), so the only thing that
  moves is the vector sector.
* `li6_alpha_channel` sets `dis_target` from `source`. `Hulthen` is
  `DEUTERON()`, bit for bit.
* `BreakupOptions::source` (breakup.hpp), forwarded by `Pipeline` next to
  `cluster_beta`; `ClusterBreakup` builds `deuteron_channel(beta, p_d, source)`.
* `vector_dilution_of(p_d)` (beams.hpp) — the ONE home of 1 − 1.5 P_D, which
  was retyped at four sites and would have become five.
  `ALPHA_D_VECTOR_POLARIZATION`, `DEUTERON_VECTOR_POLARIZATION` and
  `li6_cluster_polarization` are now written through it; all three are the same
  constexpr expression, so every value is bit for bit unchanged.

**Cost, stated.** Every polarized ⁶Li tagged-α observable made with
`--cluster-wave vmc` moves by **−2.027 %** (= 1/1.020687). Nothing on the
Hulthén default moves — no reference file, no `validation/reference/*.json`
gate, no assertion count. ⁷Li deliberately does **not** move: `TRITON()` is a
Faddeev-family per-nucleon slot and this tree has no AV18 A = 3 wave function
to switch it to, so on the α–t channel the flag selects the relative motion
alone — stated in the code so it is not a second silent version of the defect.
The d control's struck object is a free neutron on either setting.

**What one `--cluster-wave vmc` ⁶Li run now says the whole-nucleus
polarization is**, end to end, measured on 400 k generated events per setting
(⟨pol⟩/M of the struck nucleon, both M = ±1):

| run | sampled | closed form | which reading |
|---|---|---|---|
| `--cluster-wave hulthen` | 0.811253 ± 0.0011 | **0.811228** = `LI6_CLUSTER_POLARIZATION` | scenario α–d × scenario d |
| `--cluster-wave vmc`, **now** | 0.886169 ± 0.0011 | **0.887076** = `li6_cluster_polarization(VMC_P_D_LI6, deuteron_av18_p_d())` | AV18 α–d × AV18 d |
| `--cluster-wave vmc`, before | — | 0.905427 = `LI6_CLUSTER_POLARIZATION_VMC` | AV18 α–d × **scenario** d |

The `vmc` run is now **one Hamiltonian end to end**, and it lands on §C5.5's own
third table row. **§C5.5's conclusion is untouched by this**: 0.887076 is still
+4.61 % from the ab-initio `LI6_POLARIZATION_VMC_SIX_BODY` = 0.848 where the
shipped 0.811228 is −4.34 %, the cluster PRODUCT still carries the error, and
the inclusive default still stays at 0.811228.

#### C5.5b.4 The secondary defect: the flag was accepted where it is never read

`PipelineConfig::validate()` **accepted** `cluster_wave = VmcAV18` on
`Inclusive` and `CoherentLi6`. Verified by calling `validate()` on both: both
returned. That is the exact defect §C5.4 named — *"the flag was being SILENTLY
IGNORED on that channel, which is the same defect the b1 knobs are refused
for"* — and it is **the mechanism by which §C5.5's 11.61 % reaches a user**:
they set the flag wanting a "VMC ⁶Li" inclusive row, nothing refuses it, and
they get `LI6_CLUSTER_POLARIZATION` = 0.811228 with no warning. The same
`validate()` refuses `cluster_vmc_mc_sigma`, a **0.02 %** effect; the asymmetry
had no defence.

**It is now refused** on the two channels that never read it, with an error
naming the 11.61 % and the constant they would silently have got. What the
refusal does *not* claim: it does not make an inclusive VMC ⁶Li available —
there is none, and §C5.5 is the argument that substituting `VMC_P_D_LI6` into
`LI6_CLUSTER_POLARIZATION` would move it *away* from the only ab-initio anchor.
It converts a silent 11.61 % into an error message.

#### C5.5b.5 What is now true at the five claim sites

*"A single run builds exactly one channel, so nothing inside one run is
inconsistent"* is true of the code as it now stands, and every site says so
**with the correction attached** rather than as a standing claim: each records
that the sentence was false when written, by how much (2.069 %), where
(`dis_target` and the T1 breakup), and what closed it. The `--cluster-wave`
declaration in `pipeline.hpp` now lists *exactly where the flag reaches and
where it does not*, which is the list a user needs and the list whose absence
allowed the defect.

**Pinned:** `tests/test_tagged.cpp` **T27** (33 assertions: the two DIS
targets, the exact 1.020687209533 on g₁ and A_∥ at four (x, Q²), the breakup's
own model on both settings, the 0.887076 whole-nucleus reading, and the ⁷Li /
d-control non-moves), `tests/test_pipeline.cpp` (the `validate()` refusal on
all five channels), and the two pytest mirrors in `python/tests/test_module.py`.

*What would change this decision:* an off-shell or six-body embedded-deuteron
P_D. The free AV18 deuteron is not the deuteron inside ⁶Li — it is the same
approximation `P_D_DEUTERON` = 0.045 makes, taken at the right Hamiltonian
instead of at a chosen number. That is stated at `DEUTERON_AV18`.


---

## C6. Open item O3, offline — and the collaboration ask, reconciled

### C6.1 O3 — what is bounded offline, and what genuinely needs the configurations

O3 asks how much the incoherent/coherent split would move if the α core were
sampled from genuine correlated GFMC ⁴He configurations (as arXiv:2605.00454
uses) instead of the uncorrelated product of one-body densities
`ClusterConfigSampler` actually draws. **That question stays unanswerable
from this repository**: no configuration table for any nucleus is committed
here, and the split is a property of a Good–Walker amplitude's
configuration-to-configuration fluctuation under a dipole model, which
nothing in this tree computes.

What can be bounded with data already in the tree:
`data/vmc/bonus_other_clusters/he4.dd` is a genuinely correlated VMC overlap
of the SAME ⁴He wavefunction onto an α → d+d cluster channel — the identical
method (Forest *et al.*, PRC 54, 646 (1996)) and the identical 19-Apr-2004
run family as `li6.ad`/`li7.at`, just a different two-cluster decomposition
of the same nucleus. `validation/o3_alpha_correlation_bound.py` compares its
own relative-motion second moment to the analytic prediction of an
**uncorrelated** product of `he4.density` for the identical observable: the
separation D between the centroids of the α's two 2-nucleon halves.

**The null model, worked once.** `ClusterConfigSampler`'s alpha core already
IS this null model — four iid draws from a source density inflated by
λ² = 4/3 (`alpha_cm_inflate`'s own constant), recentred so the four sum to
zero, reproducing `he4.density`'s own ⟨s²⟩ = R₁² exactly after recentring
(`cluster_config.hpp` "the recentred single-nucleon density" note). For
D = (s₁+s₂)/2 − (s₃+s₄)/2, the sample mean cancels identically —
D = (v₁+v₂−v₃−v₄)/2 exactly, independent of recentring and of which of the
three 2-2 partitions of four iid draws is chosen — so

    <D²> = (1/4) Σᵢ <vᵢ²> = <v²> = (4/3) R₁² ,

a closed form with no free parameter, using only `he4.density`'s own
measured R₁.

**The measured comparison:**

| quantity | value |
|---|---|
| R₁ = √⟨r²⟩, `he4.density` (recomputed) | 1.44131 fm (file's own 1.4404) |
| `he4.dd` S+D norm (recomputed) | 0.97175 (file's own "ndx" 0.973) |
| `he4.dd` P_D (recomputed) | 0.02478 (file's own 0.024/0.973 = 0.0247) |
| r_rms(d–d separation), **correlated** (`he4.dd`) | **1.80948 fm** |
| r_rms(D), **uncorrelated** product of `he4.density` | **1.66429 fm** |
| ratio, variance ⟨r²⟩_corr / ⟨D²⟩_uncorr | **1.18210 (+18.2 %)** |
| ratio, rms | **1.08724 (+8.7 %)** |

**Reading.** Genuine 4-body correlation widens this ONE position-space
moment, for this ONE clustering channel, by a several-to-twenty-percent
effect — not a factor of several. That is small next to §C1's factor-7.5
quadrupole gap and comparable to the wave-function-choice systematic already
measured on the tagged tensor observable (6.6 %, §C5.2 note). **This is not a
measurement of how the incoherent/coherent split itself would move** — it is
a different observable (α → d+d, not α → 4 free nucleons; a position-space
moment, not a Good–Walker amplitude's fluctuation) standing in for the
general SIZE of correlation effects in this nucleus, on the only genuinely
correlated data this repository holds. What genuinely needs the actual GFMC
configurations: any statement about the split itself.

### C6.2 The collaboration ask, reconciled

Two committed copies of the Mäntysaari-group collaboration request
disagreed: `docs/OPEN_ITEMS_SOLUTIONS.md` §11.2 and
`docs/open_items/run_2026-09-02/design_G_cluster_config.md` §10, by one
number ("3 %" vs "4 %" on the α+d model's ⁶Li point-radius agreement — both
correct, against different reference radii: 3.0 % against the measured
`LI6_R2_POINT_FM2` = 2.4655 fm, 3.9 % against the *ab initio* VMC
`LI6_R_POINT_VMC_FM` = 2.4433 fm the design's own T6 gate quotes). More than
that number, both drafts predate §C2's answer to O5.

The reconciled draft — **the single canonical copy** —
`docs/open_items/run_2026-09-03/mantysaari_collaboration_draft.md` folds in
§C1 (the factor-7.5 decomposition), §C2 (O5) and §C6.1 above. Both earlier
copies now point to it rather than restating the ask.

**What was done to those two copies, exactly.** In each of them the ask
blockquote and both editor's notes were **DELETED** — in
`OPEN_ITEMS_SOLUTIONS.md` §11.2 and in `design_G_cluster_config.md` §10 — and
replaced by a pointer to the canonical file plus a one-line statement of what
each note had found. They were **not** annotated in place, and an earlier
version of this section, of the two change lists below and of
`run_2026-09-03/PLAN.md`'s C6 bullet said they were; that claim is withdrawn.
The deletion is what "reconcile into ONE canonical draft" means: an ask
standing in three places in three versions is three asks, and annotating would
have left the "3 %" / "4 %" disagreement sitting in the tree next to the file
whose whole purpose is to resolve it. The removed text is in git at 31ed181,
is quoted nowhere else, and carried no measured number that is not in the
canonical draft.

**The rule this run follows for a dated record**, which is what makes that
consistent with the opposite treatment everywhere else in the run:

* A **finding, conclusion or premise** is never dropped. It stays where it
  stood, struck through or standing, with a dated correction beside it; where a
  struck-through bullet's body is elided to a stub, the closure that replaces it
  restates the finding in full. §C5.5 above keeps its
  own retracted *"nothing inside one run is inconsistent"* and carries the
  correction as a block under it; `run_2026-09-03/PLAN.md`'s **C1 bullet keeps
  the premise §C1 measured false** (the ±0.012 band's near edge is **+0.030
  fm², the wrong sign**, not −0.05) with the retraction annotated below it,
  because a brief is a record of what was *asked*, not of what is true.
* A **duplicated artefact** — one text committed in more than one copy, the
  copies disagreeing — is reconciled into one canonical copy and the superseded
  copies are **replaced by a pointer**, because keeping them keeps the
  disagreement. What each copy said, and what its editor's notes found, is
  carried in one line at the pointer.

The ask is the second case. C1's premise and C5.5's retraction are the first.
Nothing in this run deletes a result.

The original draft's request (c) asked the group to run their existing
polarized-J/ψ machinery on this repository's ⁶Li configuration tables.

> *(Revised 2026-09-04.)* **This paragraph previously said that request (c)
> was "a request to spend person-months confirming a result already in this
> tree: the channel is blind (0.749 σ at one EIC year, 160 fb⁻¹/u for 3 σ)".
> That is withdrawn.** §C2 now gets a **band** for the same channel —
> 2.63 σ at its low edge, 2.84–3.29 σ at its top, 3 σ at 8.3–13.0 fb⁻¹/u —
> once the photoproduction region and the
> second lepton channel are included, so request (c) is a request for a
> calculation that a plausible measurement would need — which is what it was
> drafted as. The draft has been rewritten a second time accordingly.

Requests (a) (the GFMC ⁴He configurations — the only route to closing
O3 above) and (b) (the code generalization to any A) are unaffected, since
neither is J/ψ-specific. Request (c) now states the deciding number, the
kinematic window it belongs to, the two-factor quadrupole budget and the two
things that are **not** established (the efficiency below Q² = 0.1 and the
far-forward working point) plainly rather than omitting them.

**Nothing has been sent.** The canonical draft says so at its own top and
leaves that decision to the author.

---

## C7 The anchor gate could not fail, and 53 anchors had drifted behind it

**The gate's rule was the defect.** `validation/check_physics_channels_links.py`
accepted a reference `` `path:line` `name` `` if `name` appeared *anywhere*
within **±2 lines** of `line`. Every way an anchor can go wrong by one or two
lines — sliding onto the doc comment above its declaration, onto the blank line
before it, onto the closing `};` of the struct above, or onto the NEXT
declaration in the file — is inside that window, and a symbol's own doc comment
almost always names it. So the gate reported **0 broken** on a document in which
it could not have reported anything else. That is three phases in a row.

### C7.1 What had drifted (measured, not asserted)

Criterion: the reference's value at `31ed181` pointed at a line that DECLARES or
DEFINES the name in that file's `31ed181` content, and its value in the working
tree does not, in the working tree's content. **53 anchors**, all of them
pointing *earlier* than the declaration:

> **What 53, 113 and 166 were measured against, and why they cannot be
> re-derived (added 2026-09-04, C8).** The document those three counts describe
> was the *pre-fix* working-tree `PHYSICS_CHANNELS.md` — an UNCOMMITTED
> intermediate that this pass's own `--fix` run then overwrote. It no longer
> exists in the tree or in git, so re-running the gate today reproduces none of
> them. They are kept because they are the finding; they are not reproducible
> figures.
>
> **CORRECTION (2026-09-04, C9): the replacement figures were not reproducible
> either, and are withdrawn.** This box previously offered `1071 | 111` and
> `1071 | 373` "under today's strict gate". Following its own recipe gives
> neither: no run of the gate prints 111 or 373, and 1071 is not that gate's
> count of the document — it is the count printed by the **`31ed181` gate**,
> whose consuming name-group swallowed one citation. The two numbers came from
> two different definitions (111 excludes the allow-list-pin failures, 373
> includes them but excludes the unnamed one) and neither is a figure the gate
> prints. §C9.1 re-measures with **one** definition and **one** recipe.

| file | offset | count | what the anchor now lands on |
|---|---|---|---|
| `include/lipolgen/cluster_config.hpp` | **+1** | **14** | 12 doc comments, one `};`, and `VMC_LI6_DENSITY` on **`VMC_HE4_DENSITY`'s declaration** |
| `include/lipolgen/coherent.hpp` | **+2** | **19** | 9 doc comments / blanks / `};`, and **10 on another symbol's declaration** |
| `src/core/tagged.cpp` | +1 | 3 | two blank lines, one doc comment |
| `src/core/tagged.cpp` | +2 | 1 | a `// ----- TaggedModel` banner |
| `include/lipolgen/breakup.hpp` | +2 | 1 | the doc comment of `triton_sf` |
| everything else | +3 … +57 | 15 | comments, banners, call sites |

The `+1` / `+2` split across three files is the signature of a **uniform line
shift applied to the document** that under-counted inserted header lines — not
of `--fix`, which moves one reference at a time. **13** of the 53 land on a
different symbol's declaration line, which is the worst shape: it reads as
correct. `include/lipolgen/coherent.hpp:228` `x_coh` pointed at
`double f0 = 0.04;`, `:230` `slope_b` at `double x_coh = 0.01;`, `:421`
`x_pom_max` at `double m_x_min = …;`, `:508` `set_optics` at
`xpom_model()`, `:509` `set_weighted_azimuth` at `t_max()`.

*(An adversarial review of the same document reported **29**, twelve of them in
`cluster_config.hpp`, and put the `+2` group in `tagged.hpp`. Re-measured here:
it is **53**, **14** in `cluster_config.hpp`, and the `+2` group is
`coherent.hpp` (19) plus `src/core/tagged.cpp` (1) — `tagged.hpp`'s four are off
by +27, +42, +55 and +57. The review's count was low and its file attribution
for the `+2` group was `tagged.hpp` where the tree says `src/core/tagged.cpp`;
its description of the mechanism was right.)*

### C7.2 The gate, rebuilt

**Strict is the default.** A reference must satisfy both, on the referenced line
itself, with no window:

* **S1** the name occurs there as **code** — not inside a `//`, `///`, `/* */`
  or `#` comment, and not inside a string that merely mentions it. A comment
  that names a symbol is documentation *about* it, not where it lives.
* **S2** if the file declares or defines the name **anywhere**, the reference
  must be on one of those lines. A call is not a declaration.

S2 bites only in the declaring file, and that is the carve-out that replaces the
±2 window rather than widening it: the document cites call sites on purpose
(`src/core/xsec.cpp:88` `g1_nucleus` is the one place a kernel calls it), and
those live in translation units that do not declare the symbol, so they pass
with no exemption. Recognised declaration forms are listed in the script's
docstring; they include the two things the old window was really covering —
multi-line declarations (the name is on the first line, which is where the
anchor goes) and macro wrappers (`LIPOLGEN_CHANNEL(LI6_ALPHA_TAG, …)` counts,
as does a name that exists only as a string: a CLI flag, or a binding key in
first-argument or subscript position).

**Six exemptions, each with its reason**, printed on every run so an exemption
cannot hide: `xsec.hpp:220` `g1_rank3` (no such symbol exists — the row cites
the comment that says why there is no rank-3 slot); `pipeline.cpp:567`
`coherent` (the English word inside `validate()`'s refusal message, which is the
text the row quotes); and four use sites cited on purpose *inside* their
defining file — `pipeline.cpp:1085` `optics_lumi_factor`, `rc.cpp:1003`
`is_tagged_channel`, `breakup.cpp:332` `proton_fraction`,
`pythia_bridge.cpp:679` `dis_parton_fraction`. An exemption is reached only
after the declaration rule has already rejected the reference, so the *other*
references to the same name stay strictly checked.

> **CORRECTION (2026-09-04, C8): the sentence that stood here — "an exempted
> reference must still carry its name somewhere on the exact line, so a shift
> is still caught" — was FALSE, and it was false in two ways.** The test was
> `mentions()`, which searches the raw line including comments and strings, so
> an allow-listed anchor could drift anywhere in its file that happens to name
> the symbol — a doc comment above the declaration being the likeliest landing
> spot, which is the exact shape C7 existed to kill. And the key was
> `(path, name)`, not the reference: the *second* claim above, that the other
> references to the same name stay checked, held only for names cited in a
> DIFFERENT file. Cite the same name twice in the file that owns the
> exemption and both citations were exempted by the one entry.
>
> Fixed in C8. An exemption is now keyed on **(path, line, name)** — one
> citation, never a name — and carries a **PIN**: a distinctive substring of
> the line it names (`"cfg_.lumi_pb * optics_lumi_factor()"`,
> `"is_tagged_channel(channel_)"`, …). The pin must be on that exact line or
> the gate breaks, so a shift moves the code out from under the pin and is
> caught; a second citation of the same name in the same file gets no
> exemption and is told so by name.

`--fix` re-anchors onto the **nearest declaration** (refusing an exact tie
between two, which is how the two `finite_q_delta` members were separated by
hand) instead of onto the nearest textual match. `--loose` keeps the historical
±2 rule for editing a file whose line numbers are in flux; it is not the gate.

### C7.3 What the strict gate found, and what was done

**166 violations** on first run (measured against the pre-fix uncommitted
document described in C7.1's box — not reproducible today): the 53 above plus
**113 pre-existing**, which the ±2 window had hidden since before this run — 16 CLI flags in `cli.py`
anchored one line above their own `add_argument`, ~30 in `pipeline.hpp` on the
doc comment above the field, the `cluster.cpp` reader functions on their banner
comments, and so on. **160 were re-anchored** onto the declaration (across 27
files, 57 lines of the document) and **6 allow-listed** with the reasons above.
Ten needed judgement rather than nearest-match and were done by hand: the three
`src/core/tagged.cpp` `dis_target` references now point at the three assignments
they describe (`:265` α-tag, `:295` ⁷Li, `:319` d-control) instead of all three
at one comment; `b1_nuclear.hpp:611` `finite_q_delta` goes to `:689`
(`DeuteronConvolutionB1::Options`, default **true**) and not to the equidistant
`:533` (`Li6ConvolutionOptions`, default **false**) which is cited beside it;
and `tagged.hpp` `n_of_kc` goes to the three-argument overload the row's own
sentence names.

`python/tests/test_doc_link_gate.py` — **16 tests** — pins the document against
the strict gate and the rule against fixtures: a comment naming the symbol, a
blank line and the neighbour's declaration are each rejected; the same anchor is
*accepted* under `--loose`, which is the hole itself, written down; `--fix`
walks past a nearer comment to the declaration while `--fix --loose` stops at
the comment; and every `ALLOW` key must still be reached by the real run, so a
dead exemption fails the suite instead of covering the next drift.

**Tally after C7: 1092 references checked, 0 broken, 6 allow-listed, strict.**
No source file, no test expectation and no reference JSON was touched by C7;
the C++ suite is unchanged at 388 cases / 17 217 290 assertions.
(C8 raises the checked count to **1094 point citations + 97 ranges**; see §C8.)

### C7.4 Minor: `shifted_by_sigma` had no pytest naming it

`VmcRadial::shifted_by_sigma` is bound in `python/bindings.cpp` and was covered
only through the `li6_alpha_channel(vmc_mc_sigma=…)` path. Two tests in
`python/tests/test_module.py` now name it: that a Python-built table (no error
column) is returned unchanged by every shift — which is what makes the band's
zero row bit for bit — and, on the real ANL table, that psi → psi + n·dpsi
elementwise **in the same direction at every k**, that `dpsi` survives the shift
so the band can be re-applied, that the provenance records it, and that the
`vmc_mc_sigma` knob IS this function rather than a second path.

---

## C8 The hardened gate was blind to two other citation shapes, and to itself

C7 made `path:line` `name` strict and left three shapes in the same document
checked by nothing or by almost nothing. Measured 2026-09-04 on the post-C7
document.

### C8.1 What was blind

| shape | how many | what checked it before C8 |
|---|---|---|
| `` `path:first-last` `` a RANGE | **93 citations** | **nothing.** The point regex requires a backtick straight after the digits, so a range never matched it — before OR after C7 |
| `` `:first-last` `` a bare continuation range | **5** | nothing; not even the path was resolved |
| `` `path:line` `` with NO name | **17** | file exists, line in range. Nothing else |
| a citation ADJACENT to another | **10** (1 point, 9 ranges) | **nothing.** The point regex's name group `` [^`\n]{0,12}`(...)` `` CONSUMED the next citation as this one's "name", so the swallowed citation was never matched in its own right and the swallowing citation was checked against a name that was really a path |
| an ALLOW-listed citation | 6 | see C7.2's correction box: much less than the docstring claimed |

### C8.2 The 20 stale ranges, and how they were found

Criterion, the same shape as C7.1's: the citation's `first-last` is unchanged
from `31ed181`, and the block of file content that lived at those numbers at
`31ed181` no longer lives there. **20 of the 89 distinct blocks**, offsets +1 to
+319:

| file | offset | blocks |
|---|---|---|
| `include/lipolgen/coherent.hpp` | **+319** | 2 |
| `include/lipolgen/tagged.hpp` | **+100** | 2 (one of which was already mis-anchored before `31ed181` — see C8.4) |
| `docs/USAGE.md` | **+140**, +88 | 4 |
| `docs/CONVENTIONS.md` | **+48**, +35 | 6 |
| `src/core/coherent.cpp` | +78 | 1 |
| `include/lipolgen/pipeline.hpp` | +60 | 1 |
| `include/lipolgen/beams.hpp` | +37 | 1 |
| `python/lipolgen/__init__.py` | +3 | 1 |
| `tests/test_tagged.cpp` | +1 | 1 |
| `include/lipolgen/tagged.hpp:53-59` | 0 | 1 — the block was REWRITTEN in place, still at 53-59; **no change needed** |

Nineteen needed their numbers moved; the twentieth did not. The task that
commissioned this pass predicted **six**. The criterion above is its own
criterion applied to ranges and it gives **20**; the six are presumably the
subset whose *source file* is C++ or Python (10 by that reading) or the subset
with a visible structural tell (4 — the ones with a blank edge). The number is
20, and every one is listed above by file and offset.

**Only 4 of the 20 had a blank edge**, so a rule that checks bounds and edges —
"the minimum" — would have found four of them. That is why C8 fingerprints
instead.

### C8.3 R4: the block is fingerprinted

`validation/physics_channels_ranges.json`, 89 entries, one per distinct cited
block: `sha256` of the block text plus its first and last line in clear so the
sidecar's own diff is readable. Written and rewritten only by
`check_physics_channels_links.py --record-ranges`. On a mismatch the gate hashes
every window of the same length in the file, so it reports which of two things
happened:

* **MOVED** — the recorded block is found elsewhere: the message names the new
  range, and `--fix` rewrites the citation and carries the fingerprint across
  (the block is byte-identical, only its address changed);
* **EDITED in place** — found nowhere: re-read the row, confirm the citation
  still says something true, then `--record-ranges`.

**The running cost is the point.** Editing a block that `PHYSICS_CHANNELS.md`
cites breaks this gate until a human confirms the row. That is a treadmill, and
it is the only mechanism that would have caught 16 of the 20 above.

### C8.4 What was corrected in the document — 32 citations, 38 occurrences

* **19 stale ranges re-pointed** (C8.2).
* **11 sloppy range edges tightened**: 10 ranges whose first or last line was
  blank were trimmed inward, and `physics_literature.md:28-29` — a one-line
  paragraph plus the blank after it — became the point citation `:28`. These
  are pre-existing, not phase-C drift.
* **1 range that was mis-anchored before `31ed181`**: the row citing [CW20]/
  [CW26] "for the spectator-on-shell/struck-off-shell impulse approximation"
  pointed at `tagged.hpp:255-257`, which at `31ed181` was `double e_lab…;`,
  `};` and a blank — the tail of `SpectatorLab`. It now points at
  `tagged.hpp:364-375`, the CONVENTION paragraph that states exactly that
  ("The spectator is put ON SHELL … so the struck cluster is the OFF-SHELL
  one. This is the standard impulse approximation of spectator tagging").
* **1 range wrong since before `31ed181` in the other direction**:
  `docs/CONVENTIONS.md:119-121` (the hadronic-X balance row) held the X-must-be-
  TIMELIKE rule at `6f67d84`, had drifted to `:189-191` by `31ed181` and to
  `:237-239` now. Re-pointed to `:237-239`.
* **2 unnamed citations that landed on blank lines**:
  `include/lipolgen/breakup.hpp:153` → `:160` (the `Nucleon,  ///< … the d /
  3He control channels` enumerator the row calls "the ³He control") — a **phase
  C regression**, the block moved +7; and `tests/test_rc.cpp:588` → `:495`
  (`SUBCASE("the first C0 zero lies in [2.9, 3.3] fm^-1")`, which is the refit
  window the row quotes) — **blank at `31ed181` too**, so pre-existing.

### C8.5 The other three fixes to the rule

1. **The name group no longer consumes.** `REF` and `RANGE` both match the name
   by lookahead, and the name must be **adjacent** — one space. All 1075 named
   point citations are written that way; the twelve-character window bought
   nothing and cost ten citations that were never checked plus one false
   attachment (`physics_literature.md:28`, whose "name" was the next sentence's
   `tensor_harmonics_gamma`).
2. **A base class is not a declaration.** `class ToyF2 : public UnpolSF {`
   counted as declaring `UnpolSF`. `TensorSF` had five declaration lines in
   `sf.hpp` — `:359` and four derived classes' headers `:421`, `:432`, `:444`,
   `:494`. No citation sat on one, but `--fix` could have moved one there. The
   base-clause guard leaves `TensorSF` with `:359` and `:361` (the class and its
   virtual destructor).

   > **CORRECTION (2026-09-04, C9): the sentence that stood here — "and drops
   > the count of named citations accepting more than one line from 156 to
   > 120, and of those surviving a uniform +2 drift from 36 to 34" — is
   > WITHDRAWN. The guard's effect on both counts is 0.** Toggling `CLASS_HEAD`
   > changes neither figure on the working-tree document (120 / 34 either way)
   > nor on the `31ed181` one (119 / 38 either way against working-tree sources,
> 120 / 42 either way against `31ed181` sources — the re-review of 2026-09-04
> found the table had not said which). Its measured reach over
   > the cited set is **6 citations** — `sf.hpp:77` `UnpolSF`, `:161` `PolSF`,
   > `:359` `TensorSF`, `triton_sf.hpp:226` `TritonSpectralFunction`,
   > `fsi.hpp:214` `FsiWeight`, `rc.hpp:356` `Spin1ElasticFF` — and every one
   > of the six **already accepted more than one line and already survived a
   > +2 drift with the guard ON**, because each is a `class X {` whose
   > `virtual ~X() = default;` sits two lines below it. The guard is real and
   > is kept: it closes a `--fix` hazard (`--fix` could have re-anchored
   > `TensorSF` onto a derived class's header). It just does not move those
   > two counts, and 156 reproduces from no state measured here. See §C9.2.
3. **ALLOW is per citation and pinned** (C7.2's correction box). Measured: the
   six exemptions covered **33 lines** between them under the old
   `(path, name)` + `mentions()` rule — every line of their file that names the
   symbol, comments and strings included. `src/core/pipeline.cpp` `coherent`
   alone covered **16** (`:57`, `:269`, `:409`, `:567`, `:972`, `:973`, `:985`,
   `:988`, …, most of them other refusal messages); `src/core/breakup.cpp`
   `proton_fraction` covered **4** — the definition `:178`, the DEUTERON draw
   `:258`, and the two triton draws `:285` and `:332` — so the row's citation
   could have moved from the triton branch to the deuteron branch, a different
   physics statement, and been reported as "allowed" with the triton reason
   printed beside it. Under the new rule each exemption covers exactly ONE
   citation; moving that one now says so by name ("`proton_fraction` is
   allow-listed in this file at another line only; an exemption covers one
   citation, not a name"). The pin is a *content* check on the pinned line, not
   a uniqueness test — `is_tagged_channel(channel_)` occurs on three lines of
   `rc.cpp` and `proton_fraction(in.x, in.q2, 1, 3)` on two of `breakup.cpp`;
   uniqueness comes from the line key, drift detection from the pin.

   Five mutations were run against the real document to check the gate can
   fail: inserting three lines at the top of `coherent.hpp` (**41 broken**, two
   of them "the cited block MOVED to :473-474 / :516-522"); a +1 typo on one
   range (**2 broken** — the fingerprint is orphaned and the new range is
   unrecorded); an unnamed citation moved onto a blank line (**1**); one
   fingerprint deleted from the sidecar (**1**); and `src/core/rc.cpp` shifted
   by one line, which takes the code out from under an ALLOW pin (**32**, the
   pin among them).

   > *Re-run 2026-09-04 against the gate as C9 leaves it: **41** (2 of them
   > MOVED), **2**, **1**, **1** and **35**. Four reproduce exactly. The fifth
   > was 32 and is now 35 because `rc.cpp` carries three S3-pinned use sites
   > (`:953` `a_transfer_frac`, `:983` `lam_e`, `:994` `CoherentLi6`) that the
   > same one-line shift also moves; one ALLOW pin failure, as before.*

### C8.6 Residuals — what the gate still cannot see, measured

1. An **unnamed** citation is pinned to "non-blank and not bare punctuation".
   Of the 19, **14 survive a uniform +1 drift, 12 a uniform −1, 10 both**.
   Naming the symbol is the only real fix.
2. A **named** citation is pinned to the set of lines that declare its name:
   **120 of 1075** accept more than one line, **34** still land on a
   declaration after a uniform +2 drift — nearly all `class X {` on N with
   `virtual ~X() = default;` on N+2. Invisible by construction.

   > **UNDER-STATED, and corrected 2026-09-04 (C9).** That "120 of 1075" was
   > over the wrong denominator: it counts only the **1027** citations that
   > are pinned to declaration lines at all. **42 of the 1075 are not** — they
   > pass by the use-site carve-out, which pinned them to EVERY line of the
   > file carrying the name as code. **29 of those accepted more than one
   > line**, the widest **14** (`python/lipolgen/__init__.py:249` `isotope`),
   > then 11 (`src/core/pipeline.cpp:531` `B1Model`), 9, 6, 6, 6. Over all
   > 1075 the honest figure was **149**, not 120. §C9.3 closes that class
   > (rule S3, a fingerprint of the cited line): the 42 use-site citations now
   > accept exactly one line each, so over all 1075 it is **120**, which is
   > the declaration-pinned residual and nothing else.
3. A **range**, a pinned use site or an external line with a valid fingerprint
   that a maintainer re-records without re-reading the row is blessed. The gate
   can force the question; it cannot answer it.
4. Files outside the path set (`include/ src/ python/ data/ validation/ docs/
   tests/ examples/` + `README.md`) are unchecked. `pyproject.toml` and
   `CMakeLists.txt` are cited by line range in the packaging row and are not in
   the set; that row already says so. `README.md` was added in C8, which is why
   the point-citation count is 1094 and not 1093.

   > **INCOMPLETE, and corrected 2026-09-04 (C9).** The path set is not the
   > only thing that decides what is checked: **5 point citations (naming 6
   > lines) and 1 range name PYTHIA's own sources** — `ProcessLevel.cc:645`,
   > `Pythia.cc:662`, `BeamRemnants.cc:662,935`, `BeamSetup.cc:645`,
   > `BeamSetup.cc:869-875`, `BeamParticle.cc:178` — which resolve outside the
   > repository, matched no rule, and were **silently unchecked**, not
   > declared unchecked. Rule D (§C9.4) reads them from the dependency tree
   > `env.sh` sets up and skips them **by name** when it is absent.
5. R3 (a range's name must be declared inside it) is implemented and currently
   fires on **0 of 97**: no range in the document supplies an adjacent name.
   Attaching the prose name instead was measured — wrong **9 times in 11** by
   backward proximity, **3 in 4** by forward proximity — and is not done.

### C8.7 Tally

`python/tests/test_doc_link_gate.py` goes from **16 tests to 37**: a moved
block, an edited block, a dead fingerprint, a bare continuation range, both
blank range edges, out-of-bounds and reversed ranges, R3, an unnamed citation on
a blank line and on `};`, a base class, an adjacent citation read as a name, the
ALLOW pin, and an exemption refusing to cover a second citation of the same
name in the same file.

**1094 point citations + 97 range citations (89 distinct blocks) checked,
0 broken, 6 allow-listed, strict.** No source file, no test expectation and no
reference JSON was touched by C8.

*(C9 adds the 6 external citations to the printed line and the S3 use-site
pins to the sidecar: the gate now prints `1094 references checked (strict), 97
ranges, 6 external, 0 broken, 6 allow-listed`, and
`python/tests/test_doc_link_gate.py` goes from 37 tests to 48.)*

---

## C9 The gate is sound; its RECORD was not (2026-09-04)

C7 and C8 hardened the mechanism and it holds: a uniform one-line insertion at
the top of all **77** cited files is caught on **1075 of 1075** named point
citations and **97 of 97** ranges; `--fix` cannot land on a base-class mention;
five of the six exemptions were per-citation. What did not hold is the
arithmetic written down about it. **Four recorded measurements reproduce from
nothing** — the counts **111**, **373** and **1071**, and the **156 → 120 /
36 → 34** attributed to the base-clause guard — one residual was stated over
the wrong denominator, one exemption's escape was in the carve-out beside it
rather than in the exemption, and one citation class was unchecked without
saying so. This section re-measures each and closes the two mechanism holes the
re-measurement exposed. Nothing here changes a physics number, and no C++
source, test expectation or reference JSON is touched.

### C9.1 The "reproducible replacement" counts were not reproducible

**Withdrawn: `1071 | 111` and `1071 | 373`** (C7.1's box). Following that box's
own recipe reproduces neither figure, and 1071 is not the reference count of
the gate the box says it used. Diagnosed:

* **1071** is what the **`31ed181` gate** prints. Today's gate counts **1072**
  citations in the same document — C8's lookahead fix un-swallowed one citation
  that the old consuming name-group had eaten. The table mixed the old gate's
  count with the new gate's verdict.
* **111** and **373** come from two *different* definitions. Neither is a
  number any run prints: 111 is the named-point failures **excluding** the
  allow-list-pin failures, 373 is a count **including** them but **excluding**
  the unnamed one.

**ONE definition, used everywhere below.** *The `broken` count that
`validation/check_physics_channels_links.py` prints* — the whole of it, ranges
and points and all — from the gate **as it stands at the end of this pass**,
run with **no sidecar** (`validation/physics_channels_ranges.json` absent, as
it is at `31ed181`) and with `$LIPOLGEN_DEPS`, `$LIPOLGEN_DEPS_PREFIX` and
`$LIPOLGEN_PYTHIA_SRC` all **unset**, so the external citations of rule D are
skipped rather than counted and the figure does not depend on which
dependencies the reader has installed. The gate's rules move, so
a figure like this is only meaningful with the gate version named; that is why
it is named.

**ONE recipe. It runs from the repository root and writes nothing into the
working tree:**

```sh
T=$(mktemp -d); U=$(mktemp -d); Z=$(mktemp -d)
git archive 31ed181 | tar -x -C $Z && (cd $Z && python3 validation/check_physics_channels_links.py | sed -n 1p)   # row 0
git archive 31ed181 | tar -x -C $T                       # row 1: 31ed181 sources
cp validation/check_physics_channels_links.py $T/validation/
tar -c --exclude=./build --exclude=./.git . | tar -x -C $U   # row 2: worktree sources
git show 31ed181:docs/PHYSICS_CHANNELS.md > $U/docs/PHYSICS_CHANNELS.md
rm -f $U/validation/physics_channels_ranges.json         # the sidecar is not at 31ed181
for d in $T $U; do (cd $d && env -u LIPOLGEN_DEPS -u LIPOLGEN_DEPS_PREFIX \
    -u LIPOLGEN_PYTHIA_SRC python3 validation/check_physics_channels_links.py \
    | sed -n 1p); done
```

| # | gate | doc | sources | what it prints |
|---|---|---|---|---|
| 0 | `31ed181` | `31ed181` | `31ed181` | `1071 references checked, 0 broken` |
| 1 | today's | `31ed181` | `31ed181` | `1072 references checked (strict), 98 ranges, 6 external, 239 broken, 3 allow-listed, 6 skipped` |
| 2 | today's | `31ed181` | working tree | `1072 references checked (strict), 98 ranges, 6 external, 483 broken, 3 allow-listed, 6 skipped` |

Row 0 is the finding itself, and it is the only row that needs nothing but git:
the gate that shipped at `31ed181` **could not fail**. Row 1 is the debt the
±2 window was hiding; row 2 is that debt plus everything phase C's line motion
added.

**What the two totals are made of**, from the same output — save it and run:

```sh
grep -cE '^   [^ ]+:[0-9]+  '        # point-citation failures
grep -c  'names it as code on'       #   ... of which S3 use-site pin demands
grep -cE '^   [^ ]+:[0-9]+-[0-9]+  ' # range failures
grep -c  'block is not fingerprinted' ; grep -c 'line of the range is blank'
```


| | row 1 | row 2 |
|---|---|---|
| point citations | **141** (27 of them S3 pin demands a sidecar-less run cannot satisfy) | **385** (10 S3) |
| ranges | **98** — 80 unfingerprinted, 18 with a blank edge | **98** — 82 unfingerprinted, 16 blank edge |
| external | 6 skipped, not counted | 6 skipped, not counted |
| **total printed** | **239** | **483** |

The range column is the same 98 in both rows — every range of a document that
predates the sidecar is unfingerprinted — so the difference **483 − 239 = 244**
is exactly the difference in point failures, **385 − 141**.

### C9.2 The base-clause guard's drop was 0, not 156 → 120

See the correction box in §C8.5. Measured by toggling `CLASS_HEAD` between a
real pattern and one that never matches, and recomputing every citation's
acceptance set: **6 citations change, and the two counts do not.**

| | guard ON | guard OFF |
|---|---|---|
| accept more than one line (working-tree doc) | 120 | 120 |
| survive a +2 drift | 34 | 34 |
| accept more than one line (`31ed181` doc, working-tree sources — row 2's tree) | 119 | 119 |
| survive a +2 drift (`31ed181` doc, working-tree sources — row 2's tree) | 38 | 38 |
| accept more than one line (`31ed181` doc, `31ed181` sources — row 1's tree) | 120 | 120 |
| survive a +2 drift (`31ed181` doc, `31ed181` sources — row 1's tree) | 42 | 42 |

The 6 citations it does reach are `sf.hpp:77` `UnpolSF` (accepts `:77, :79`
with the guard, `+ :113, :201` without), `sf.hpp:161` `PolSF` (`+ :181, :216`),
`sf.hpp:359` `TensorSF` (`+ :421, :432, :444, :494`), `triton_sf.hpp:226`
`TritonSpectralFunction` (`+ :291`), `fsi.hpp:214` `FsiWeight` (`+ :294`) and
`rc.hpp:356` `Spin1ElasticFF` (`+ :592, :631`). Every one is a `class X {`
whose `virtual ~X() = default;` is two lines below, so every one **already**
accepted more than one line and **already** survived +2 with the guard on. The
guard is kept for what it actually does — it stops `--fix` re-anchoring a base
name onto a derived class's header — and 156 reproduces from no state measured
here.

### C9.3 The residual class, over the right denominator — and rule S3

The 1075 named point citations split three ways, not one:

| class | how many | pinned to | accept > 1 line |
|---|---|---|---|
| declaration-pinned | **1027** | the lines that declare the name | **120** |
| use site (file declares it nowhere) | **42** | *every* line naming it as code | **29** |
| allow-listed | 6 | one line, by pin | 0 |

So over all 1075 the honest figure was **149**, and the widest citations in the
document were all in the class the residual did not mention: **14** accepting
lines for `python/lipolgen/__init__.py:249` `isotope`, **11** for
`src/core/pipeline.cpp:531` `B1Model`, then 9 (`__init__.py:424` `b1_model`),
6, 6, 6 (`sampler.cpp:41` `q2_edges`, `pipeline.cpp:957` `InclusiveSampler`,
`triton_sf.cpp:90` `HBARC_GEV_FM`).

**The decision.** *Land on the FIRST use* was rejected: it is wrong on the
document as written — `isotope`'s first use is `__init__.py:10`, and the row
cites `:249` because that is the site it is talking about — so the rule would
demand 42 re-pointings, most of them onto lines that say something else.
*Within N lines of a named anchor* was rejected for the same reason and because
N is a free parameter with nothing to set it from. What distinguishes use site
249 from use site 251 is **content and nothing else**, so the only rule that
can catch the drift is a content pin — which the gate already has, for ranges
(R4) and for exemptions.

**S3, implemented.** A use-site citation whose name appears as code on **more
than one line** of the file must carry a fingerprint of the cited LINE in the
same sidecar; one that appears on exactly one line needs no entry, because the
rule already pins it. Cost: **27 new sidecar entries** (29 citations, two pairs
of which cite one line); 13 lone use sites stay free. The gate reports MOVED
(and `--fix` re-points) or EDITED, exactly as for a block. After it, all 42
use-site citations accept exactly one line, so **120 of 1075** is now the whole
residual and it is entirely the declaration-pinned class. The widest remaining
are 6 (`fsi.hpp:238` `sigma_xn_mb`), 5 (`rc.cpp:695` `fc`), 5
(`rc.hpp:758` `n_eta`), then 4s.

### C9.4 The one exemption that was not per-citation — and where the leak was

The sixth is **`src/core/pipeline.cpp:567` `coherent`**. Its `ALLOW` entry is
per-citation and pinned like the other five; the leak was **around** it, in the
carve-out. The other five names are *declared* in the file that holds their
exemption, so S2 governs every other citation of them there and says so
("an exemption covers one citation, not a name"). `pipeline.cpp` declares no
symbol `coherent`, so the use-site carve-out accepted **`:996` and `:1009`** —
both `cfg_.coherent`, a member access — with **nothing printed at all**. The
exemption's effective reach over that file was **3 lines, of which 2 were
invisible**.

S3 closes it without a new mechanism: those two lines are an ambiguous use
site, so a citation of either is now reported as needing a fingerprint. Checked
directly (`test_an_exempted_name_cannot_escape_through_the_carve_out`, and by
running the gate on a one-line document citing each of `:567`, `:996`,
`:1009`): `:567` allowed, `:996` and `:1009` broken. Each of the six exemptions
now covers exactly one citation, and nothing else in its file passes unpinned.

### C9.5 The PYTHIA upstream citations were silently unchecked — rule D

Five point citations naming six lines, and one range, name PYTHIA's own
sources: `ProcessLevel.cc:645` and `Pythia.cc:662` (the `doVarEcm` gate and the
abort that make `Beams:allowMomentumSpread` useless here), `BeamRemnants.cc:662`
and `:935` (`wPosRem = eCM`, the reason the final state is forced to the
initialisation beam total), `BeamSetup.cc:645` (the on-shell PDG-mass forcing),
`BeamSetup.cc:869-875` (990 counts as a hadron) and `BeamParticle.cc:178`
(meson-like beam handling). They match no in-repo path, so **the gate never saw
them** — it did not report them as unchecked either.

Rule **D** now resolves them against the dependency tree `env.sh` sets up.
*Two names exist for that prefix and both are honoured*: `env.sh` exports
`$LIPOLGEN_DEPS`, while `CMakeLists.txt`, `README.md:89` and `docs/USAGE.md:19`
call the same directory `$LIPOLGEN_DEPS_PREFIX` (it is a CMake cache variable
that CMakeLists.txt also reads from the environment). Both name
`<...>/deps/install`; the unpacked sources sit beside it at
`<...>/deps/src/pythia8*/src`, and `$LIPOLGEN_PYTHIA_SRC` overrides. Each citation
is checked for bounds, a non-blank edge, and the same sha256 fingerprint a
range gets, so a PYTHIA upgrade that moves the cited line says so instead of
leaving the row quietly wrong (verified by injecting one line at the top of a
copy of `BeamSetup.cc`: `the text recorded against pythia8317… is no longer
there -- it MOVED to :646`). **All seven anchors were re-read against
`pythia8317` and all seven are correct.** With no PYTHIA tree on the machine
the gate asserts nothing, prints each citation as `skipped` **by name**, leaves
the exit status alone, and `--record-ranges` carries the recorded external
entries across instead of dropping them.

### C9.6 The drift experiment, end to end, and the tally

Not a model of a drift — the drift itself: insert (or delete) *k* lines at the
top of every one of the **77** cited files in a mirror tree, run the real gate,
and count the citations it does **not** report. A "**+k** drift" is *k* lines
lost above the citation (the document's number now points *k* lines later in
the old content); "−k" is *k* lines inserted.

| drift | named point citations surviving (of 1075) | before S3 |
|---|---|---|
| −1 (one line inserted) | **0** | 0 |
| −2 | **0** | 0 |
| +1 (one line deleted) | **2** | 4 |
| +2 | **34** | 36 |

Every survivor is now declaration-pinned and is RESIDUALS 2 exactly: a
`class X {` cited on N whose `virtual ~X() = default;` sits on N+2. The
use-site survivors S3 removed were `python/lipolgen/cli.py:604` and
`python/lipolgen/__init__.py:424` at +1, and `__init__.py:249` and `:424` at
+2. Ranges: all **89 distinct cited blocks** — hence all **97** range
citations — are reported at every one of the four drifts.

C8.5's five mutations were re-run against the gate as this pass leaves it:
**41 / 2 / 1 / 1 / 35** broken. Four reproduce exactly; the fifth (a one-line
shift of `src/core/rc.cpp`) was 32 and is 35 because that file now carries
three S3-pinned use sites the shift also moves.

Sidecar: **89 → 123 entries** — 89 blocks, 27 pinned use-site lines, 7 upstream
lines. `python/tests/test_doc_link_gate.py`: **37 → 48 tests** (a lone use site
needing no pin, an ambiguous one needing one, a drifted one caught and fixed, an
edited one, the exempted-name escape, and six for rule D — recorded, changed,
comma list and range, either name for the deps prefix, skipped by name, and
`--record-ranges` not dropping what it cannot check).

**Measured after this pass:** **390 doctest cases / 17 217 416 assertions /
1 skipped / 0 failed**, **345 pytest passed**, and the strict docs gate at
**1094 references checked, 97 ranges, 6 external, 0 broken, 6 allow-listed**
(`, 6 skipped` when the dependency tree is absent, exit status 0 either way).
No source file, no test expectation, no reference JSON and no physics number
was touched: the C++ suite is unmoved at 390 / 17 217 416, and pytest is
334 → 345 by the eleven new gate tests alone.

---

## C2 FOURTH PASS (2026-09-04) — the band's TOP was an unquantified factor presented as quantified

*Everything under this heading is what the fourth adversarial pass over §C2
moved. The verdict (MARGINAL, as a band), the band's LOW edge (2.63 σ, 3 σ at
13.0 fb⁻¹/u) and "3 σ inside {1, 10, 100} at both ends" all survived it
unchanged. The band's TOP did not.*

### A1. The non sequitur, retracted at every site that carried it

*"Every entry is at the same rigidity, so its A-ordering is a species lever on
its own."* Sites: `validation/o5_a2_reach.py` (the `species_scaling` docstring
and the §3b(b) report text), `include/lipolgen/coherent.hpp`
(`chang26_species_efficiency()`'s comment and `COHERENT_JPSI_EFF_IR8_LI7`'s),
**`python/bindings.cpp`'s docstring for the same table** (a seventh site the
first sweep missed, found by grepping for the phrase rather than for the
number), `phase_C_numbers.md` §C2.3c(b) and §C2.8 item 3a,
`OPEN_ITEMS_SOLUTIONS.md` §11.3b, `STATUS.md`, and the **name of a
`tests/test_coherent.cpp` TEST_CASE** ("… so it is a lever" → "… is fixed
RIGIDITY, NOT fixed E/u").

**Why it is a non sequitur.** Fixing R = A·E/Z eliminates rigidity and nothing
else. At fixed R the per-nucleon energy is **E/u = R·Z/A** and the total beam
momentum is **p_z = Z·R**, so both still vary down the list:

| | across `chang26_species_efficiency()` |
|---|---|
| E/u | **118 GeV/u (⁷Li) → 183 (³He)** — a ×1.5508 spread |
| Z | **1 → 8** |
| p_z = A·E | **274 → 2192 GeV** |

**And the confound is the size of the effect.** The chain's *own other leg*
measures d ln ε/d ln E_ion = **−0.8656 and −0.6807** (`he3_energy_slopes`), so
the list's ×1.5508 spread in E/u is worth **×1.3481–1.4620** in ε — as large
as the entire ×1.16–1.22 that was being read off the list as "species". The
list is a **joint (A, Z, E/u) lever**; nothing in it decomposes that, and the
chain asserted the decomposition without testing it. Pinned in T10d
(`eu_lo`/`eu_hi`/`z`/`p_z`/`eu_worth`) and in
`test_the_species_list_is_fixed_rigidity_not_fixed_energy`.

### A2. The one validation was mis-specified — the wrong Z — and its directional instruction is deleted

The form is ε = exp(−B(A)·p_T,cut²), B = `gaussian_slope(1.2 A^{1/3})`,
p_T,cut = **0.19577 GeV** inverted from the ⁷Li row. arXiv:2511.05638's own
criterion is *"within a safe distance **from the beam**"* — a cut on the
**angle** — so **p_T,cut = θ·p_z = θ·Z·R and it carries the charge.** ³He and
⁴He are **Z = 2**; ⁷Li is **Z = 3**. The shipped test applied ⁷Li's Z = 3 cut
unchanged to the Z = 2 pair, and its own comment said *"same Z"* without
noticing it was the **wrong** Z.

| reading of p_T,cut | ³He/⁴He, measured **1.0955** | ²D | ⁴He | ⁹Be | ³He |
|---|---|---|---|---|---|
| **∝ Z** — (2/3)² on the pair | **1.0967**, i.e. **0.1 %** | 1.953 | 2.003 | 0.214 | 2.005 |
| **Z-independent** | 1.2309 — *the shipped one*, "overstates by 12 %" | **1.003** | **1.034** | **1.047** | **1.161** |

(The last four columns are pred/meas on the list's **absolute** values.) So:
reading (1) reproduces the ³He/⁴He *ratio* to 0.1 % but wrecks the absolute
values; reading (2) reproduces the absolute values to **0.3–4.7 %** for ²D, ⁴He
and ⁹Be and misses **³He by 16 %** — and ³He is precisely the one nucleus at
**E/u = 183** rather than the family's 137, for which the chain's own energy
scan predicts a **17.9–22.2 % shortfall** (0.7783–0.8211). **Reading (2) is
not a description of the whole list**, and this record does not claim it is:
the form degrades with A and misses ¹²C by 32 % (1.3217) and ¹⁶O by a factor
3.1 (3.1320). It is a *local* interpolation over the light end, and A = 6 is
inside the range where it works. **Three readings of one table, and nothing in the table picks
between them.** The third pass took the third — the shipped one, with the
wrong Z — because it is the only one that yields the directional instruction
***"Read the low end"***, which was printed at three sites
(`o5_a2_reach.py`'s docstring and report, `tests/test_coherent.cpp`). **That
instruction is deleted**, and T10d now carries `pred_z` = 1.0967 and
`pred_wrong_z` = 1.2309 side by side.

### A3. The four same-energy entries the chain never used, and the top

⁶Li's own fixed-rigidity energy is Z/A × 275 = **137.5 GeV/u**. **Four** of the
seven entries — ²D (A = 2), ⁴He (4), ¹²C (12), ¹⁶O (16) — sit at **137 GeV/u**,
so they bracket A = 6 with **no energy step at all**, and ⁴He → ¹²C is the
adjacent pair. The same closed forms on that set:

| construction | ε(⁶Li, 137 GeV/u) | factor on 0.1775 | S at the band's top | 3 σ at |
|---|---|---|---|---|
| `tag_acceptance`, p_T,cut from ¹²C | 0.17629 | **×0.9932** | **2.842 σ** | 11.1 fb⁻¹/u |
| log-linear in R_G | 0.17823 | ×1.0041 | 2.858 σ | 11.0 fb⁻¹/u |
| log-linear in A | 0.20059 | ×1.1302 | 3.032 σ | 9.8 fb⁻¹/u |
| `tag_acceptance`, p_T,cut from ⁴He | 0.20124 | ×1.1338 | 3.037 σ | 9.8 fb⁻¹/u |
| linear in A | 0.23655 | **×1.3327** | **3.292 σ** | 8.3 fb⁻¹/u |

**→ ×0.99–1.33: the species leg STRADDLES 1.** Against the ×1.16–1.22 of the
confounded reading, whose two anchors are 19 GeV/u apart in E/u — worth
×1.107–1.138 on the chain's own slopes, i.e. most of what it returned.

> **"The band's TOP just crosses 3 σ" is NOT ESTABLISHED, and is withdrawn.**
> The honest form is a span: **2.84 … 3.29 σ**, two forms below 3 σ and three
> above. It was printed as a single **3.15 σ / 9.1 fb⁻¹/u** in
> `OPEN_ITEMS_SOLUTIONS.md`, `PHYSICS_CHANNELS.md`, `PLAN.md`, `STATUS.md`,
> `phase_C_numbers.md`, `physics_literature.md`, `design_G_cluster_config.md`,
> `mantysaari_collaboration_draft.md`, `python/lipolgen/configs.py` and the
> two pins; every one now carries the span.

**What is unaffected and stays.** The band's LOW edge (**2.63 σ**, 3 σ at
**13.0 fb⁻¹/u**) contains no species leg at all, so it does not move. MARGINAL
stands. And **3 σ is inside {1, 10, 100} at both ends on every form: 8.3 …
13.0 fb⁻¹/u.** The de-squeezed-optics rung moves with the top —
**0.73–0.92 σ, 3 σ at 106–167 fb⁻¹/u** — and is still outside the band at both
ends, so `Optics::lumi_fraction` is still the single correction that on its own
restores a NO.

### B. Three sites still quoted the point instead of the band

The run's own ruling is *"QUOTE IT AS A BAND … neither may be quoted alone"*
(§C2, `OPEN_ITEMS_SOLUTIONS.md` §11.3b, `configs.py`, `PHYSICS_CHANNELS.md`,
`PLAN.md`, `physics_literature.md`). Three sites violated it, all in lines the
third pass itself added:

1. **`include/lipolgen/cluster_config.hpp`, the `a2_from_geometry` docstring**
   — the code-level home of the number, and one of the five sites the third
   pass *named* as fixed. Three defects: (a) it quoted **"2.62 sigma, 3 sigma
   needs 13 fb^-1/u — MARGINAL"**, the forbidden point alone, citing only
   §11.3a; (b) it said **"Two further things are UNESTABLISHED"** where every
   other site says three or four, and the one it omitted was **the
   decay-lepton reconstruction efficiency** — precisely the factor that is
   unbounded below and sets the band's open bottom; (c) it gave the
   de-squeezed optics as a bare **"the headline is 0.73 sigma"** where that
   too is a band. All three fixed, and the count is now **four**.
2. **`docs/OPEN_ITEMS_SOLUTIONS.md` §C6's live paragraph** — *"marginally
   measurable — 2.62 σ at one EIC year, 3 σ at 13 fb⁻¹/u, inside the band"*.
   Now the band.
3. **A whole-tree sweep** for `2.62`, `13 fb`, `13.1`, `3.15` and `9.1 fb` near
   σ. Every remaining occurrence is now either explicitly labelled as the
   band's **uncorrected** middle (with the band beside it) or is inside a
   dated historical record of a superseded pass.

### What the pins do now

`_self_check` and `python/tests/test_o5_reach.py` no longer carry a scalar
top. The keys `significance_hi`, `lumi_for_3sigma_lo` and `band_factor_hi` are
**removed from the returned dict** — a `KeyError` is the cheapest guard
against the point coming back — and both files assert their absence. The pin
is

```python
assert res["significance_hi_min"] < 3.0 < res["significance_hi_max"]
```

which **fails if "the top crosses 3 σ" is ever re-asserted as established**,
because that means the top's low reading has been pushed above 3. Beside it:
the species leg must straddle 1 (`species_lo < 1.0 < species_hi`); the
confounded reading must still clear 1 at both ends, so the retraction stays
visible; the Z-corrected ³He/⁴He prediction must agree within 1 %; the
Z-independent absolute values must fit ²D/⁴He/⁹Be better than the Z-scaled ones
do **and must still break down at ¹²C/¹⁶O**, so "it fits the absolute values"
cannot be widened to the list.

**Seven mutations were run against those pins, and all seven fail:**

| mutation | caught by |
|---|---|
| species leg restored to the confounded ×1.16–1.22 | `species_lo < 1.0 < species_hi` |
| the top's **low** reading pushed above 3 σ ("the top crosses 3 σ" re-asserted) | the straddle assertion, by name |
| the top's **high** reading pulled below 3 σ ("the top cannot reach 3 σ") | the same assertion, from the other side |
| a scalar `significance_hi` re-added to the returned dict | the "no scalar top" assertion |
| the **wrong Z** put back into the ³He/⁴He test | `abs(he3_he4_ratio_z − 1) < 0.01` |
| the band's low edge moved off 2.63 σ | `2.0 < significance_lo < 3.0` |
| 3 σ pushed outside {1, 10, 100} at the top | the per-key band check |

### Files changed by this pass

* `include/lipolgen/coherent.hpp` — the retraction and the E/u–Z–p_z
  arithmetic in `chang26_species_efficiency()`'s comment; the ×1.16–1.22 in
  `COHERENT_JPSI_EFF_IR8_LI7`'s comment replaced by ×0.99–1.33 "straddles 1".
* `python/bindings.cpp` — the same retraction in the binding docstring for
  `chang26_species_efficiency`.
* `include/lipolgen/cluster_config.hpp` — `a2_from_geometry`'s docstring: the
  band instead of the point, four unestablished things instead of two, the
  optics as a band.
* `tests/test_coherent.cpp` — the TEST_CASE renamed; the E/u, Z, p_z and
  same-energy-count checks; the same-energy interpolation; the Z-corrected
  ³He/⁴He prediction beside the wrong-Z one; the Z-independent absolute-value
  readings; *"Read the low end"* deleted.
* `validation/o5_a2_reach.py` — `species_list_is_not_a_species_lever()`,
  `_interp_forms()`, `he3_he4_readings()`, `species_scaling_same_energy()`;
  `species_scaling()` demoted to "the confounded reading"; §3b(b) rewritten;
  §7 quotes the top as a span; the scalar top keys removed; `_self_check`
  re-keyed.
* `python/tests/test_o5_reach.py` — the species test renamed and rewritten,
  a new three-readings test, the two-leg and verdict pins re-keyed.
* `python/lipolgen/configs.py`, `docs/OPEN_ITEMS_SOLUTIONS.md` (new §11.3c),
  `docs/PHYSICS_CHANNELS.md`, `docs/open_items/physics_literature.md`,
  `docs/open_items/run_2026-09-03/{PLAN,STATUS,phase_C_numbers,mantysaari_collaboration_draft}.md`,
  `docs/open_items/run_2026-09-02/design_G_cluster_config.md` — the band, at
  every site that carried the point or the 3.15 σ top.

**The shipped default is bit for bit.** `COHERENT_JPSI_EFF_IR8_LI7` still
holds the paper's 0.1775, no correction is folded into any constant,
`validation/reference/` is untouched and the rtol 1e−12 gates do not move.

**Measured after this pass:** **390 doctest cases / 17 217 416 assertions /
1 skipped / 0 failed**, **334 pytest passed**, and the strict docs gate at
**1094 references checked, 97 ranges, 0 broken, 6 allow-listed** *(the gate
printed no `external` field yet; C9 adds it — the line is now `1094 references
checked (strict), 97 ranges, 6 external, 0 broken, 6 allow-listed`)*. Before
it:
390 / 17 217 385 / 333 and the same gate line — so the pass moved no C++ test
case, added **31 assertions** (all inside the rewritten T10d) and **one
pytest**, and changed no citation count.

## What was changed in the tree by this pass

* `include/lipolgen/cluster_config.hpp` — new `LI6_ETA_DS_GK` /
  `_STAT` / `_SYST` (single home for the George & Knutson measurement); the
  `LI6_QUADRUPOLE_GFMC_FM2` attribution nit (it lives in *this* header, not
  `rc.hpp`); the two-vs-three-factor budget note; the exact-linearity note on
  `asymptotic_ds_ratio()`; and the corrected description of the raw block above
  9 fm ("run out of signal", with the measured S/N, replacing "is Monte Carlo
  noise", which is false for R₀ out to ~11 fm).
* `tests/test_cluster_config.cpp` — **T22b** (35 assertions) and **T22c**
  (14 assertions); **T22**'s two `-0.025` literals and its MESSAGE now read
  through `LI6_ETA_DS_GK` (+49 assertions in total, 377 → 379 cases).
* `python/bindings.cpp` — the three new `LI6_ETA_DS_GK*` constants, plus
  `LI6_QUADRUPOLE_FM2` (already used by the cluster tests as a bare `-0.0818`
  literal in five places; it now reads through the constant, beside the
  existing `LI6_R2_POINT_FM2` binding).
* `python/tests/test_cluster_config.py` —
  `test_eta_rides_the_dial_and_the_two_factor_q_budget`,
  `test_the_gk_eta_band_spans_a_sign_change_in_q`.
* Prose, at **every** site that carried the old one-line version:
  `design_G_cluster_config.md` §2.7 (the new η → dial → Q derivation), its
  up-front framing paragraph and §10 (O1 and O2 both marked CLOSED);
  `OPEN_ITEMS_SOLUTIONS.md` §11.2 (the quadrupole paragraph, the "what remains
  open" O1/O2 entries, and the section's closing rule);
  `docs/USAGE.md` (the dial paragraph, the η table row and the caveat block);
  `docs/PHYSICS_CHANNELS.md` (the tensor-signal row and the `OverlapRaw`
  description); `python/lipolgen/configs.py`'s printed-on-every-run caveat.
  The two collaboration-letter drafts were **annotated, not rewritten** by
  *this* pass. The later §C6 pass then **deleted** both ask blockquotes and
  their editor's notes and replaced them with a pointer to the canonical
  draft — §C6.2 has the rule and the reason.

### …and by the §C2 pass (2026-09-04)

* `include/lipolgen/coherent.hpp` / `src/core/coherent.cpp` — new
  `EstarlightLi6Row` + `estarlight_li6_coherent()`, on the same pattern as
  `mantysaari_a2_deuteron()`: the **single code home** of the eSTARlight ⁶Li
  cross sections, their Q² > 0.7 restriction, their measured-radius variants,
  the fitted slopes and the branching fractions, which until now existed only
  as prose in three documents. Inert — nothing in the generator reads it.
* `python/bindings.cpp` — `EstarlightLi6Row` and `estarlight_li6_coherent()`.
* `validation/o5_a2_reach.py` — **new**: the whole O5 estimate, importable and
  self-checking, reading all three of its inputs from the bindings.
* `python/tests/test_o5_reach.py` — **new**, 8 tests: the table, the signal,
  the |t|-slope suppression, the flip plan, the agreement of the Fisher form
  with `estimators::cos2phi_fit_err`, and both halves of the verdict.
* `tests/test_coherent.cpp` — **T10a** (41 assertions), the same arithmetic on
  the C++ side (379 → 380 cases).
* Prose, at **every** site that carried the old O5 wording or the a₂ = +0.026
  claim without its reach: `OPEN_ITEMS_SOLUTIONS.md` row 11, §11.1 (the "to
  4 %" correction) and §11.2 (the new §11.3, the O5 entry, and the editor's
  note on the drafted letter); `design_G_cluster_config.md` §10's O5 bullet;
  `docs/PHYSICS_CHANNELS.md`'s eSTARlight row and tensor-signal row;
  `include/lipolgen/cluster_config.hpp`'s `a2_from_geometry` docstring;
  `estarlight_li6.md` §2a **annotated, not rewritten** (it is a dated record).

### …and by the §C2 SECOND pass (2026-09-04) — the verdict reversed

**The first §C2 pass answered NO and the NO was wrong.** What was actually
done, beyond the prose:

* **A new eSTARlight run**, not an estimate — same build
  (`939b11a24499398392d959db81c7502aeec91046`), same beams, same seed
  5574531, same everything, with `MIN_GAMMA_Q2` lowered. The Q² > 0.1 rows
  reproduce `estarlight_li6.md` §2a/§2b to every printed digit (σ, B, ⟨|t|⟩,
  ⟨W⟩, ⟨Q²⟩), which is what makes the scan comparable. Written up as
  **`estarlight_li6.md` §2f**.
* `include/lipolgen/coherent.hpp` / `src/core/coherent.cpp` — new
  `EstarlightLi6Q2Row` + `estarlight_li6_q2_floors()` (the scan, 9 rows);
  `EstarlightLi6Row::b_rmeas` and `::branching_all` (two numbers that had
  been **retyped into the validation script**, one of them with no code home
  at all); and `COHERENT_JPSI_EFF_IR8_LI7` = 0.1775, likewise. Inert.
* `python/bindings.cpp` — all four of those.
* `validation/o5_a2_reach.py` — rewritten. The docstring's "three in-tree
  inputs and NOTHING ELSE" was **false** and is replaced by the five that
  actually enter; the Q² windows and the second lepton channel are priced;
  `PZZ2_DIFFERENCE` = 0.81 is the new headline ⟨P_zz²⟩ and `PZZ2_OPTIMAL`
  = 0.90 is printed beside it; `Optics::lumi_fraction` is stated with its
  size; and `_self_check`'s unfailable
  `assert lumi_for_3sigma_jpsi > 100.0` is replaced by pins in **both**
  directions.
* `tests/test_coherent.cpp` — **T10c** (new case: the scan, its reproduction
  of §2a/§2b, the rate multipliers and the slope invariance) and T10a
  extended with `b_rmeas`, `branching_all`, the efficiency constant, the
  0.81 difference estimator and the new headline. 385 → 386 cases.
* `python/tests/test_o5_reach.py` — 8 → 12 tests; the two "the verdict is NO"
  tests are replaced by four that pin the restricted row, the headline, the
  pessimism ladder and the perfect-detector ceiling.
* Prose, at **every** site that carried the NO or quoted the number bare:
  `OPEN_ITEMS_SOLUTIONS.md` row 11, §11.3 and the new **§11.3a** (the record
  of what moved); `PHYSICS_CHANNELS.md`; `physics_literature.md` §11 and its
  cos 2φ warning; `PLAN.md`'s C2 bullet; `STATUS.md`; `configs.py`;
  `a2_from_geometry`; `design_G_cluster_config.md` §10; and the
  **Mäntysaari draft, rewritten a second time** — request (c) restored from
  "we are not asking that anymore" to a request for the polarized J/ψ
  calculation *in photoproduction*, with the two unestablished factors
  written inside the ask. `estarlight_li6.md` §2e is **annotated, not
  rewritten**.
* Tally after the second pass: **386 doctest cases / 17 217 245 assertions /
  1 skipped / 0 failed** (385 / 17 217 169 before it), **287 pytest** (283),
  **1088 doc references / 0 broken**, **0 reference JSONs moved**. No shipped
  default changed: `estarlight_li6_coherent()`'s six existing columns are
  bit-identical, both new tables are inert, and `CoherentScenario`'s values
  are untouched.

### …and by the §C4 + §C5 pass (2026-09-04)

* `include/lipolgen/coherent.hpp` / `src/core/coherent.cpp` — **ΔB defined**,
  once, on the `eps_b0` declaration, plus `delta_b_m(m)` (its single code home)
  and `slope_at_azimuth(phi, m)`; the `eps_b0` and `a2_m_state` docstrings
  repaired (the old "ΔB₀/B" was off by a **sign**); the ⁶Li band, the implied
  quadrupole, the `eps_b0`–`slope_b` degeneracy and the author decision written
  at the declaration.
* `include/lipolgen/cluster_config.hpp` / `src/core/cluster_config.cpp` —
  `quadrupole_from_a2_slope` (the exact inverse of `a2_from_quadrupole`, so a
  scenario can be **asked** what it assumes) and
  `ClusterConfigSampler::quadrupole_for_eta` (C5.1's converter, plus the stored
  pre-dial moments it needs).
* `include/lipolgen/cluster.hpp` / `src/core/cluster.cpp` — `VmcRadial` carries
  `dpsi()`, `has_errors()`, `shifted_by_sigma()` and `norm2_error()`; both
  momentum/overlap factories now propagate the files' own 1σ columns instead of
  dropping them (C5.2).
* `include/lipolgen/tagged.hpp` / `src/core/tagged.cpp` — `vmc_mc_sigma` on the
  two lithium channels (refused on Hulthén); `deuteron_channel(..., source)`
  with the AV18 `fdeut` waves and `deuteron_av18_p_d()` (C5.4);
  `LI6_CLUSTER_POLARIZATION_VMC` as a **diagnostic** (C5.5).
* `include/lipolgen/beams.hpp` — `li6_cluster_polarization()` as the one home
  of the cluster product (`LI6_CLUSTER_POLARIZATION` bit for bit unchanged) and
  `LI6_POLARIZATION_VMC_SIX_BODY` = 0.848, the ab-initio anchor.
* `include/lipolgen/pipeline.hpp` / `src/core/pipeline.cpp` —
  `cluster_vmc_mc_sigma` with its `validate()` refusal; `cluster_wave` now
  reaches the deuteron control instead of being silently ignored there.
* `python/bindings.cpp` — every one of the above.
* `tests/test_coherent.cpp` **T10b** (46 assertions),
  `tests/test_cluster_config.cpp` **T23** (18),
  `tests/test_tagged.cpp` **T24** / **T25** / **T26** (159 across the three);
  `python/tests/test_module.py` — five mirrors (380 → 385 doctest cases,
  278 → 283 pytest).
* Prose, at **every** site that carried the old version: `docs/CONVENTIONS.md`
  (the ⁶Li effective-polarization entry and the N_αd entry),
  `docs/PHYSICS_CHANNELS.md`, `docs/USAGE.md`,
  `docs/OPEN_ITEMS_SOLUTIONS.md`, and `design_G_cluster_config.md` §2.8 / §10's
  O4 bullet **annotated, not rewritten** (it is a dated record). Stated exactly,
  since §C6.2's rule is about precisely this: §2.8 keeps its full original text
  with the 2026-09-04 annotation inline, while §10's **O4** and **O5** bullets
  are struck through with their bodies elided to a stub — the closure that
  follows each restates the finding in full, so no finding is lost, but the
  original wording of those two bullets is in git at 31ed181, not in the file.

**No shipped default changed** by either pass, and
`validation/reference/*.json` is untouched (verified: 0 files).

### …and by the §C6 pass (2026-09-04)

* **New file** `docs/open_items/run_2026-09-03/mantysaari_collaboration_draft.md`
  — the single reconciled collaboration draft (C6). Not sent.
* **New file** `validation/o3_alpha_correlation_bound.py` — the O3 offline
  bound (§C6.1), self-checking, no pytest gate (an offline analysis on
  `bonus_other_clusters/he4.dd`, not a shipped-physics number).
* `docs/OPEN_ITEMS_SOLUTIONS.md` §11.2 — the inline ask paragraph and its two
  editor's notes **deleted** and replaced by a pointer to the canonical draft
  (§C6.2's rule; nothing else in §11.2 was removed); new §11.5
  (O3 offline bound + C6 summary); top summary-table row 11 updated; the O3
  bullet in §11.2's "what remains open" list updated.
* `docs/open_items/run_2026-09-02/design_G_cluster_config.md` §10 — the
  inline "ask, in one paragraph" and its two editor's notes **deleted** and
  replaced by the same pointer; the O3 bullet updated with the offline bound.
  Everything in the design that is a *finding* is **annotated, not rewritten**:
  the O1–O5 bullets are struck through with their closures beside them and
  §2.7 / §2.8 keep their original text under dated annotations. The duplicated
  ask is the one thing removed, under §C6.2's rule.
* No code, no test, no reference JSON touched by this pass; 0 files under
  `validation/reference/` changed (verified).

### …and by the §C5.5b pass (2026-09-04) — the false consistency claim, repaired

* `include/lipolgen/beams.hpp` — `vector_dilution_of(p_d)`, the ONE home of
  1 − 1.5 P_D (it was retyped at four sites and about to become five).
  `ALPHA_D_VECTOR_POLARIZATION`, `DEUTERON_VECTOR_POLARIZATION` and
  `li6_cluster_polarization` are written through it; the same constexpr
  expression, so **every value is bit for bit unchanged**.
* `include/lipolgen/tagged.hpp` / `src/core/tagged.cpp` — `DEUTERON_AV18()`,
  the embedded deuteron at `deuteron_av18_p_d()` = 0.057600 (eff_pol 0.913600
  against the scenario 0.9325), and `li6_alpha_channel` sets `dis_target` from
  `source`. The ⁷Li and d-control `dis_target` lines carry, in code, the
  statement of why they do **not** move.
* `include/lipolgen/breakup.hpp` / `src/core/breakup.cpp` —
  `BreakupOptions::source`; `ClusterBreakup` builds
  `deuteron_channel(beta, p_d, source)`, so the T1 struck-nucleon spin draw is
  the run's deuteron. The `⟨2m₁⟩ = (1 − 1.5 P_D) m_S` identity in the header is
  now stated as the consistency gate it is.
* `src/core/pipeline.cpp` — `bo.source = cfg_.cluster_wave` next to
  `bo.beta = cfg_.cluster_beta`; and `validate()` **refuses** `cluster_wave` on
  `Inclusive` and `CoherentLi6`, where it is never read (§C5.5b.4).
* `include/lipolgen/pipeline.hpp` — the `cluster_wave` declaration rewritten as
  an explicit *reaches / does not reach* list, which is the list whose absence
  allowed the defect.
* `python/bindings.cpp` — `vector_dilution_of`, `deuteron_av18()`,
  `BreakupOptions.source`.
* `tests/test_tagged.cpp` **T27** (33 assertions) and the stale
  "there is no VmcAV18 overload of `deuteron_channel` at all" comment corrected;
  `tests/test_pipeline.cpp` — the `validate()` refusal across all five channels
  (12); `python/tests/test_module.py` — two mirrors.
* Prose, at all five sites that carried the false claim: `docs/CONVENTIONS.md`,
  `docs/PHYSICS_CHANNELS.md` (§3 dilution row **and** §5's channel bullets),
  `docs/USAGE.md`, `include/lipolgen/pipeline.hpp`, and §C5.5 above (a dated
  correction block, not a rewrite). §C5.4's decision paragraph records that its
  own rationale is what forbade the path it shipped.
* **MINOR, also fixed:** `docs/CONVENTIONS.md` and §C5.3 cited
  "`phase_D_numbers.md` §Q3" for the b₁ ±5 % linearity. That file's section is
  "The ±5 % N_αd systematic (design §2.1) and the Q4 knob" — the knob is **Q4**
  and there is no §Q3. Both citations now name the file's real path and heading.

**The shipped default is bit for bit.** `validation/reference/` is untouched
(0 files), the rtol 1e−12 gates do not move, and the C++ assertion count is
unchanged on the pre-existing cases. What moved is the **opt-in**
`--cluster-wave vmc` tagged-α path, by −2.027 %, deliberately.

### …and by the §C2 THIRD pass (2026-09-04) — the efficiency chain, and a claim that was false at two of three sites

The MARGINAL verdict survives; the chain it rested on had **two defects
pulling opposite ways, neither of them in §C2.8's "every assumption, in one
list"**, and one printed self-contradiction.

* **`include/lipolgen/coherent.hpp` / `src/core/coherent.cpp`** —
  `COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV` = 117.9 and
  `COHERENT_JPSI_EFF_IR8_LI7_W_MEAN_GEV` = 43.2 record *the beam and the ⟨W⟩
  the 0.1775 was measured at*, beside the constant, because it is not flat in
  either and open item O5 applies it at 10 × 99.5 / ⟨W⟩ = 30.2. The doc
  comment, which said only "none at 10 × 99.5", now carries the **direction
  and the size** the same paper supplies, and states in as many words that the
  number is a recoil-nucleus acceptance carrying **no decay-lepton acceptance
  and no reconstruction efficiency**.
* **`chang26_he3_energy_scan()`** — arXiv:2511.05638 §V.B's three e+³He points
  (18 × 183 → 32.23 %, 10 × 100 → 54.38 %, 5 × 41 → 99.77 %), verified against
  the PDF. The only measurement of the beam-energy dependence this tree has.
* **`chang26_species_efficiency()`** — that paper's p. 4 species list with
  Fig. 2's beam energies, which exists to make **one fact checkable**: every
  entry is at `A/Z × E` = 274.0–275.3 GeV/e, i.e. the same rigidity, so the
  A-ordering is a species lever separate from the energy one.
* **`EstarlightLi6Q2Row::w_mean_gev`** — ⟨W⟩ per (meson, Q² floor), which both
  transfer arguments run on and which lived only in `estarlight_li6.md` prose.
* **`validation/o5_a2_reach.py`** — `beam_energy_scaling`, `species_scaling`,
  `jpsi_lab_rapidity`, `lepton_pair_acceptance`, `jpsi_w_mean` and
  `eps_det_chain`; a new report §3b printing the whole corrected chain; the
  ladder in §5 with **the optics row inside it** and a second ladder on the
  corrected ε_det; and a verdict rewritten as a band. Nothing is retyped:
  the beams come from `default_configs`, the efficiency lever from the two new
  tables, the acceptance edge from `Scenario::eta_max`, ⟨W⟩ from the row.
* **The self-check is re-keyed.** It used to read
  `assert 2.0 < res["jpsi"]["significance"] < 3.0` on a point estimate. It now
  pins the **band**, both **directions** (beam factor > 1; lepton factor ≤ 1),
  the **cancellation**, the fact that the beam factor is **not** the paper's
  1.687, and that the optics still restores a NO at **both** ends. Nine
  mutations of the claim were checked to make it fail, including pasting 1.687
  in, dropping the lepton factor, and flipping either direction.
* **`tests/test_coherent.cpp` T10d** (two cases) and five new pytests pin the
  same things from both languages.
* **The false claim, at three sites.** *"Every single conservatism on its own
  leaves 3 σ inside the band"* is not true — the de-squeezed optics row is a
  counterexample, and `o5_a2_reach.py` printed the claim **three lines below**
  its own bullet saying so. §C2.4a had the exception buried in the same
  sentence; the script had none; `mantysaari_collaboration_draft.md` had none
  **and omitted the optics row from its ladder**. All three now state the
  exception as an exception, and the draft's ladder carries the row.

**What moved and what did not.** The verdict is still MARGINAL and 3 σ is
still inside the {1, 10, 100} fb⁻¹/u band — the two omissions cancel to 0.7 %
— but it is now quoted as a band, **2.63 σ at the low edge and (after the
fourth pass) 2.84 … 3.29 σ at the top, 3 σ at 8.3 … 13.0 fb⁻¹/u**, open
below, with ε_det named as the dominant unquantified factor and
`Optics::lumi_fraction` as the single correction that on its own restores a NO.
**The shipped default is bit for bit**: `COHERENT_JPSI_EFF_IR8_LI7` still holds
the paper's 0.1775 (the corrections are printed, never folded into the
constant), `validation/reference/` is untouched, and the rtol 1e−12 gates do
not move.
