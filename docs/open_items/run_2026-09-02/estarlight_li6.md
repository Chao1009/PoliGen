# eSTARlight unpolarized coherent baseline for e + ⁶Li — open item 11, step (1)

Run date 2026-09-02.  Scope: the *unpolarized* coherent vector-meson rate for
e + ⁶Li asked for in `docs/OPEN_ITEMS_SOLUTIONS.md` row 11 / §11 and
`docs/open_items/physics_literature.md` §11 ("Recommended path (1)"), i.e.
replace the hand-tuned `CoherentScenario` numbers with citable simulated ones.
⁷Li is run alongside as the cross-check against arXiv:2511.05638.

Nothing outside this directory was modified and nothing was committed.

---

## 1. What was run

**Generator.** eSTARlight, <https://github.com/eic/estarlight>, commit
`939b11a24499398392d959db81c7502aeec91046` ("Changes to satisfy newer Mac C++
compiler and CMake Version 4.0", 2025-05-02), version string `trunk`.
This is the same code arXiv:2511.05638 cites (its footnote 1, p. 2, and refs. [19, 20]).

**Build.** Plain out-of-source CMake, no options, no ROOT and no PYTHIA
needed for the default ASCII output:

```
git clone https://github.com/eic/estarlight
mkdir build && cd build && cmake ../estarlight && make -j8
```

cmake 3.22.1 / g++ 11.4.0 (Ubuntu 22.04) is enough — the upstream
`cmake_minimum_required` is 3.6.  The build defaults to `CMAKE_BUILD_TYPE=Debug`
and `-Wall -Wextra -Werror`; it compiled clean.  ROOT *is* on this machine but
is not used: `OUTPUT_FORMAT = 0` writes the plain-text `slight.out` event
record, which already carries everything needed (`t:`, `GAMMA: E_γ Q²`,
`TARGET:` = scattered nucleus four-vector, `SOURCE:` = scattered electron).
The single executable produced is `e_starlight`.

**⁶Li needs no code change — it runs today.**  There is no per-nucleus table
to extend.  `nucleus::init()` (`src/nucleus.cpp`) has explicit `case`s only for
Z = 82, 79, 29, ±1; everything else falls to `default:`, which sets
`_Radius = 1.2·A^{1/3}` and, for `_Z < 7`, `_rho0 = _A`.  `nucleus::rws()` and
`nucleus::formFactor()` then take the light-nucleus branch:

```cpp
// src/nucleus.cpp, formFactor(), Z < 7
const double R_G = nuclearRadius();
return exp(-R_G*R_G*t/(6.*hbarc*hbarc));            // Gaussian form factor
```

So Z = 3 lands in the Gaussian branch `physics_literature.md` predicted it
would.  The only visible symptom is a warning, which is cosmetic:

```
Warning: density not defined for projectile with Z = 3. using defaults.
```

Setting `TARGET_BEAM_Z = 3`, `TARGET_BEAM_A = 6` is all that is required.

**One optional patch, used for a sensitivity variant.**  The `1.2·A^{1/3}`
default gives R_G = 2.1805 fm for ⁶Li and 2.2955 fm for ⁷Li, both well below
the measured radii.  Since the Z < 7 branch reads `_Radius` as the Gaussian
**rms** radius (ρ(r) ∝ exp(−3r²/2R_G²) has ⟨r²⟩ = R_G²), the measured rms
charge radii can be substituted directly.  Exact diff (also saved as
`estarlight/nucleus_li_radii.patch`):

```diff
--- a/src/nucleus.cpp
+++ b/src/nucleus.cpp
@@ -111,6 +111,20 @@ void nucleus::init()
 		  }
 		}
 		break;
+	case 3:
+		{
+		  // Lithium: explicit rms charge radii, Angeli & Marinova,
+		  // At. Data Nucl. Data Tables 99 (2013) 69:
+		  //   R_rms(6Li) = 2.589(39) fm, R_rms(7Li) = 2.444(42) fm.
+		  // The Z < 7 branch below reads _Radius as the Gaussian rms
+		  // radius R_G of rho(r) ~ exp(-3r^2/(2 R_G^2)), i.e. exactly
+		  // the rms radius, so these go in unconverted.
+		  if(_A == 6)      _Radius = 2.589;
+		  else if(_A == 7) _Radius = 2.444;
+		  else             _Radius = 1.2*pow(_A, 1. / 3.);
+		  _rho0 = _A;
+		}
+		break;
 	default:
 		printWarn << "density not defined for projectile with Z = " << _Z << ". using defaults." << endl;
```

The **default (unpatched) build is what the main results table below uses**, so
that the ⁷Li numbers are directly comparable to arXiv:2511.05638, which ran the
code unmodified.  The patched build supplies the `R_rms` variant rows.

**Beams.**  LiPolGen configuration 1 = `ELECTRON_ENERGIES[1] = 10.0` GeV ×
99.5 GeV/u (`src/core/beams.cpp:119`; `docs/USAGE.md:153` "10 GeV e × 99.5 GeV/u").
eSTARlight parameterizes beams by Lorentz γ and internally uses
`E_ion = γ·A·m_p` (`src/inputParameters.cpp:197`), so γ is fixed by the
per-nucleon *momentum*:

| beam | γ used | note |
|---|---|---|
| e⁻ 10 GeV | `ELECTRON_BEAM_GAMMA = 19569.5` | 10 / 0.000511 |
| e⁻ 18 GeV | `ELECTRON_BEAM_GAMMA = 35225.1` | 18 / 0.000511 |
| ion 99.5 GeV/u | `TARGET_BEAM_GAMMA = 106.0507` | m_p·√(γ²−1) = 99.5 |
| ion 117.9 GeV/u | `TARGET_BEAM_GAMMA = 125.6605` | LiPolGen ⁷Li top (paper writes "18 × 118") |

**Q² range.**  `MIN_GAMMA_Q2 = 0.1`, `MAX_GAMMA_Q2 = 100.` — this is
arXiv:2511.05638's range verbatim: *"Each sample, consisting of 10 million
events with a kinematic range of 0.1 < Q² < 100 GeV², was passed through the EIC
afterburner"* (p. 3, §IV, right column).  A second ⁶Li set was run with
`MIN_GAMMA_Q2 = 0.7` to match LiPolGen's own inclusive analysis window (§5).

**Other settings** (identical in every run): `PROD_MODE = 12` (e+A),
`BREAKUP_MODE = 5`, `QUANTUM_GLAUBER = 1`, `SELECT_IMPULSE_VM = 0`,
`W_GP_MIN/MAX = -1` (no W cut), `CUT_PT = 0`, `OUTPUT_FORMAT = 0`,
`RND_SEED = 5574531`.  `PROD_PID` = `443011` (J/ψ → e⁺e⁻), `333` (φ → K⁺K⁻),
`113` (ρ → π⁺π⁻).  200 000 events per main run (§2a), 100 000 per
measured-radius variant (§2b), and only 2 000 per `MIN_GAMMA_Q2 = 0.7` run
(§5) — those are used for σ alone, which does not depend on `N_EVENTS`.

Representative input files (9 of the 19 runs; the others differ only in the
fields listed in §7) are in `estarlight/*.in` next to this note.

**How the cross section is read.**  eSTARlight prints two numbers:

```
 Total cross section: 1.773 nanob.                                  <- sigma(e+A -> e'+A'+VM)
 The cross section of the generated sample is 105.842 picobarn.     <- x branching ratio
```

`Total cross section` is the **production** cross section in the stated Q²
range and is computed from the luminosity table, so it does not depend on
`N_EVENTS`; the second line is that number times the decay branching fraction of
the chosen `PROD_PID` (1.773 nb × 0.0597 = 105.8 pb for J/ψ → e⁺e⁻).  The table
below quotes the first.

---

## 2. Results

`B` is the fit slope of dN/d|t| ∝ exp(−B|t|), maximum likelihood on a truncated
exponential over 0 < |t| < 0.10 GeV² (the |t| < 0.05 and < 0.20 fits agree to
within 5 % / 1 % respectively — from the `ana.py` output over `slight.out`,
not from eSTARlight's own log).  `⟨W⟩` is
√(m_p² + 2m_p E_γ − Q²) averaged over generated events, the same W the paper
defines in its Eq. (2).

### 2a. Default eSTARlight radii (R_G = 1.2·A^{1/3}: ⁶Li 2.1805 fm, ⁷Li 2.2955 fm), 0.1 < Q² < 100 GeV²

| run | σ_coh | B [GeV⁻²] | ⟨\|t\|⟩ [GeV²] | ⟨W⟩ [GeV] | median W | ⟨Q²⟩ [GeV²] | N |
|---|---|---|---|---|---|---|---|
| ⁶Li J/ψ, 10 × 99.5 | **1.773 nb** | 38.9 | 0.0254 | 32.2 | 30.0 | 0.906 | 2×10⁵ |
| ⁶Li φ,   10 × 99.5 | **30.16 nb** | 39.2 | 0.0253 | 20.5 | 15.6 | 0.292 | 2×10⁵ |
| ⁶Li ρ,   10 × 99.5 | **506.4 nb** | 38.6 | 0.0256 | 17.8 | 12.2 | 0.240 | 2×10⁵ |
| ⁷Li J/ψ, 10 × 99.5 | 2.156 nb | 43.3 | 0.0229 | 32.4 | 30.2 | 0.907 | 2×10⁵ |
| ⁷Li φ,   10 × 99.5 | 36.49 nb | 43.5 | 0.0229 | 20.7 | 15.7 | 0.292 | 2×10⁵ |
| ⁷Li ρ,   10 × 99.5 | 608.5 nb | 43.1 | 0.0231 | 18.0 | 12.4 | 0.241 | 2×10⁵ |
| ⁷Li J/ψ, 18 × 117.9 | 3.193 nb | 43.9 | 0.0227 | 43.2 | 39.3 | 0.913 | 2×10⁵ |
| ⁷Li φ,   18 × 117.9 | 44.07 nb | 43.8 | 0.0227 | 27.5 | 19.8 | 0.294 | 2×10⁵ |
| ⁷Li ρ,   18 × 117.9 | 707.8 nb | 43.4 | 0.0230 | 24.1 | 15.5 | 0.242 | 2×10⁵ |

The fitted slopes reproduce the analytic Gaussian value B = R_G²/(3ħc²) to 4 %
(⁶Li 40.70 predicted vs 38.6–39.2 fitted; ⁷Li 45.11 vs 43.1–43.9); the small
deficit is the |t_min| ≈ (E_pom/γ)² floor and the Q² dependence folded in by
the electroproduction flux, not a fit artefact.  **The slope is set by the
density alone and is identical for all three vector mesons to 1.5 %** — there
is no VM-dependent slope in eSTARlight for a coherent target.

### 2b. Measured rms charge radii (patched build), 0.1 < Q² < 100 GeV²

| run | σ_coh | B [GeV⁻²] | ⟨\|t\|⟩ | ⟨W⟩ | N |
|---|---|---|---|---|---|
| ⁶Li J/ψ, 10 × 99.5, R = 2.589 fm | 1.255 nb | 55.0 | 0.0181 | 32.3 | 10⁵ |
| ⁶Li φ,   10 × 99.5, R = 2.589 fm | 22.08 nb | 54.8 | 0.0182 | 20.6 | 10⁵ |
| ⁶Li ρ,   10 × 99.5, R = 2.589 fm | 379.5 nb | 54.1 | 0.0184 | 17.9 | 10⁵ |
| ⁷Li J/ψ, 18 × 117.9, R = 2.444 fm | 2.821 nb | 49.9 | 0.0199 | 43.2 | 10⁵ |
| ⁷Li φ,   18 × 117.9, R = 2.444 fm | 39.44 nb | 49.6 | 0.0201 | 27.5 | 10⁵ |
| ⁷Li ρ,   18 × 117.9, R = 2.444 fm | 639.4 nb | 49.0 | 0.0203 | 24.1 | 10⁵ |

Enlarging ⁶Li from 2.18 to 2.59 fm raises B by 41 % and **lowers** σ by 29 % —
the larger nucleus squeezes the same coherent amplitude into a narrower |t|
range.  This is the dominant model uncertainty on the ⁶Li rate, ±30 %, and it
is entirely a statement about the assumed density, not about the dynamics.

### 2c. \|t\| distribution — no diffractive minimum

**There is no first minimum, at any A, at any energy.**  A Gaussian density has
a Gaussian form factor with no zeros, so dσ/d|t| is a pure exponential all the
way out.  The generated distributions are single exponentials that simply run
out of statistics: the largest generated |t| is 0.27–0.29 GeV² (⁶Li:
0.270 J/ψ, 0.292 φ, 0.292 ρ) and 0.25–0.26 GeV² (⁷Li: 0.246–0.264) at
2×10⁵ events.  Anything that looks like a
diffractive minimum in a Woods–Saxon nucleus is *absent by construction* here,
which is the single most important caveat for using these numbers as an imaging
baseline (see §6, and `physics_literature.md`'s note on arXiv:2608.23445,
"enhancing the secondary diffractive lobes for light nuclei").

### 2d. A-scaling sanity check

σ(⁷Li)/σ(⁶Li) at 10 × 99.5, same code path, same energy:
J/ψ 1.216, φ 1.210, ρ 1.202, against the coherent expectation
(7/6)^{4/3} = 1.228.  All three within 2 % of A^{4/3} — the generator is
behaving as a coherent A-scaling calculation should, which is the cheapest
available validation of the Z = 3 configuration.

### 2e. Rates

Per **fb⁻¹** of e+⁶Li at 10 × 99.5 GeV/u, 0.1 < Q² < 100 GeV² (produced VMs,
before any branching fraction, acceptance or efficiency):

| channel | σ | events / fb⁻¹ | events / fb⁻¹ in the quoted decay |
|---|---|---|---|
| coherent ρ⁰ | 506 nb | 5.1 × 10⁸ | 5.1 × 10⁸ (ρ → ππ, ≈100 %) |
| coherent φ | 30.2 nb | 3.0 × 10⁷ | 1.5 × 10⁷ (φ → K⁺K⁻, 49 %) |
| coherent J/ψ | 1.77 nb | 1.8 × 10⁶ | 1.1 × 10⁵ (J/ψ → e⁺e⁻, 5.97 %) |

LiPolGen carries **no default integrated luminosity** — `PipelineConfig::lumi_pb`
is 0.0 (`include/lipolgen/pipeline.hpp:323`) and every example takes `--lumi`
from the command line — so these are quoted per fb⁻¹.  Multiply by
`Optics::lumi_fraction` for the far-forward working point actually used
(1 at the Yellow Report envelopes, 0.1467 at the ⁶Li 5 × 41 tagging point,
`docs/USAGE.md`).  With the measured-radius density the ⁶Li rates drop 25–29 %.

---

## 3. ⁷Li cross-check against arXiv:2511.05638

**The paper quotes no absolute cross sections.**  It is an acceptance study:
Chang et al., PRD 113 (2026) 032018, report detection efficiencies, Q²–x
correlations, |t| distributions and their Fourier transforms, but no σ appears
anywhere in the text or the figures (checked by full-text search for
"cross section", "nb", "µb").  So the cross-check has to be on configuration
and shape, not on a number-to-number comparison of σ.  What can be checked:

1. **Configuration is identical.**  Same generator, same commit-era code, same
   `0.1 < Q² < 100 GeV²` (p. 3, §IV right column), same VM set (J/ψ, φ, ρ),
   same ⁷Li beam (their Fig. 2 legend, p. 4: "e ⁷Li 18 × 118"; LiPolGen's
   rigidity-derived top energy is 117.9 GeV/u, a 0.1 % difference).  Our ⁷Li
   runs use the **unpatched** code, so they are the same density the paper used.
2. **|t| slope and reach.**  Their Fig. 6 (p. 6) plots the generated
   coherent-J/ψ |t| distribution for e⁷Li at 18 × 118 (open black squares)
   and e³He at 18 × 183 (open red circles) as *counts* per bin on a log axis
   running from 10² to ≈ 4×10⁶, over 0 < |t| < 0.5 GeV², from 10⁷ events.
   Read off the figure (digitized by eye, ±10 %): the ⁷Li generated curve
   starts at ≈ 2×10⁶ per bin at |t| = 0 and hits the 10² floor at
   |t| ≈ 0.22–0.23 GeV², i.e. 4.3 decades over 0.22 GeV², **B ≈ 44–45 GeV⁻²**,
   against our fitted 43.9.  With 10⁷ events and B = 44 the first-bin content
   implies a 0.005 GeV² bin width (2×10⁶ / (10⁷ · 44) = 0.0045).  The ³He
   curve falls ≈ 2.5 decades over 0.2 GeV² (B ≈ 28; eSTARlight's default
   1.2·3^{1/3} = 1.73 fm gives 25.7 analytic) and reaches the floor at
   ≈ 0.40–0.45 GeV².  Both curves are straight on the log axis apart from the
   first bins, as a Gaussian form factor requires.  Consistent — and this is
   the closest thing to a quantitative check the paper allows.
3. **Fourier transform width.**  Their Fig. 7 (p. 7, right panel) is F(b)
   for ⁷Li from |t| < 0.4 GeV², Eq. (3): F(b) = (1/2π)∫dΔ Δ J₀(Δb)
   √(dσ/d|t|), Δ = √(−t), plotted as F(b)/∫F(b)db.  For dσ/d|t| ∝ e^{−B|t|}
   the integrand is e^{−BΔ²/2} and F(b) ∝ e^{−b²/(2B)}, i.e. 1/e at
   b = √(2B)·ħc = **1.85 fm** with B = 43.9 GeV⁻² = 1.71 fm² (the same
   Gaussian as eSTARlight's projected density, T(b) ∝ e^{−3b²/2R_G²}).  The
   paper's generated (red squares) ⁷Li curve peaks at ≈ 0.38 fm⁻¹ at b = 0
   (an untruncated Gaussian normalized to unit integral would peak at
   1/√(2πB) = 0.31), is at ≈ 0.1 by |b| ≈ 2 fm, ≈ 0 by |b| ≈ 3 fm, and shows
   the small negative lobes at |b| ≈ 3–4 fm expected from the |t| < 0.4
   truncation.  Consistent at the level a by-eye reading allows.
4. **No shape structure.**  Neither their figure nor ours shows a diffractive
   minimum, as required by the shared Gaussian form factor.  Consistent, and
   for the same reason in both cases.
5. **Not reproducible here:** their global detection efficiencies
   (⁷Li 17.75 % for coherent J/ψ at top energy, p. 4) and the tagging-efficiency
   curves.  Those come from the IR-8 far-forward geometry through the EIC
   afterburner, which is downstream of eSTARlight and outside this exercise.

**Verdict:** the ⁷Li setup here is configuration-identical to the paper and
agrees with every shape statement the paper makes.  No absolute-σ validation is
possible against that paper because it publishes none; if an absolute check is
wanted, the reference to use is eSTARlight's own published ep/eA benchmarks
(Lomnitz & Klein, PRC 99 (2019) 015203), not 2511.05638.

---

## 4. Implications for `CoherentScenario::slope_b`

**This is the clean, directly transferable result.**

`gaussian_slope(r_rms)` in `src/core/coherent.cpp:53` is

```cpp
const double r = r_rms_fm / GEV_PER_FM_INV;   //  r_rms / ħc
return r * r / 3.0;                           //  B = R_rms²/(3 ħ²c²)
```

which is **algebraically the same object** as eSTARlight's light-nucleus form
factor, `exp(−R_G² t / 6ħ²c²)` squared → `exp(−R_G² t / 3ħ²c²)`.  The two codes
already share a convention; only the radius differs.  So:

| R_rms (⁶Li) | source | B analytic | B fitted from events | inside the 40–60 band? |
|---|---|---|---|---|
| 2.1805 fm | eSTARlight default 1.2·A^{1/3} | 40.70 | 38.6 – 39.2 | at/just below the lower edge |
| 2.417 fm | ← LiPolGen's own `slope_b = 50.0` | 50.00 | — | centre by construction |
| 2.589 fm | Angeli–Marinova rms **charge** radius | 57.38 | 54.1 – 55.0 | inside, near the upper edge |

Read the other way, the scenario band `slope_b ∈ {40, 60}` corresponds to
R_rms ∈ **{2.162, 2.647} fm**, and the default 50.0 to **2.417 fm**.

**Conclusion: `slope_b = 50 ± 10` survives contact with eSTARlight unchanged,
and is now citable rather than hand-tuned.**  The band's two ends are almost
exactly the two defensible ⁶Li densities — eSTARlight's own `1.2·A^{1/3}` at
the bottom and the measured charge radius at the top — and the default 2.417 fm
sits close to the ⁶Li point-*matter* rms radius (≈2.45 fm), which is the right
proxy for a gluon density.  If anything the band could be tightened to
45–58 GeV⁻², but there is no reason to: the ±30 % rate swing in §2b shows the
density assumption is the leading systematic and the wide band is honest.

One genuinely new fact: **the slope is VM-independent** in this model (1.5 %
spread across ρ/φ/J/ψ), so `slope_b` does not need a channel index.

## 5. Implications for `CoherentScenario::f0`

`f0` is *"the coherent fraction of the DIS rate at x → 0"*, with
f_coh(x) = f0/(1 + (x/x_coh)²), and the LiPolGen coherent channel carries rate
only above `COHERENT_MX_MIN_DEFAULT = 1.2 GeV`.  **eSTARlight cannot measure
that quantity**, and it is important to say so plainly rather than quote a
number that looks like it does:

* eSTARlight generates **exclusive** VM production only — M_X is the vector
  meson mass, exactly.  It has **no diffractive continuum**.
* Two of the three channels are *below* LiPolGen's M_X floor by construction:
  ρ(775) and φ(1019) both sit under 1.2 GeV.  `coherent.hpp` already says this
  — the ρ/ω/φ window "is the EXCLUSIVE VECTOR-MESON channel — a genuinely
  different process with its own t-slope, spin-density matrix and decay, not a
  low-mass limit of this one … It is Phase 2."
* Only J/ψ (3.097 GeV) is inside LiPolGen's coherent window, and it is one
  narrow resonance out of a continuum.

With that stated, the ratio that *can* be formed is still informative.  Running
⁶Li at `MIN_GAMMA_Q2 = 0.7` to match LiPolGen's own inclusive analysis window,
against `generate_inclusive --isotope 6Li --config 1`, which reports
σ_incl = 5.92 × 10⁵ pb = **592 nb** (Q² ≥ 0.70, y ∈ [0.004, 0.985], W² ≥ 8,
|η| ≤ 3.8, E′ ≥ 0.30):

| channel (⁶Li, 10 × 99.5, Q² > 0.7) | σ | σ / σ_incl |
|---|---|---|
| coherent ρ⁰ | 15.32 nb | 2.59 × 10⁻² |
| coherent φ | 1.897 nb | 3.20 × 10⁻³ |
| coherent J/ψ | 0.605 nb | 1.02 × 10⁻³ |
| **sum, all three** | **17.82 nb** | **3.01 × 10⁻²** |
| **sum, M_X ≥ 1.2 GeV (J/ψ only)** | **0.605 nb** | **1.02 × 10⁻³** |

So:

* **f0 = 0.04 is not refuted and looks reasonable at the order-of-magnitude
  level.**  The *total* exclusive coherent VM fraction of the DIS rate is
  3.0 × 10⁻², squarely inside the {0.02, 0.08} band — but that number is
  90 % ρ⁰ and therefore counts events LiPolGen deliberately excludes.
* **The part eSTARlight can see inside LiPolGen's own window is
  1.0 × 10⁻³, i.e. 40× below f0 = 0.04.**  This is a *lower bound*, not a
  measurement: exclusive J/ψ is a small fraction of coherent diffraction at
  M_X ≥ 1.2 GeV, and the continuum that dominates it is simply not in the
  generator.
* Caveat on the ratio itself: the two windows are not identical.  eSTARlight
  applies no y, W², η or E′ cut, so the denominator is the more restricted of
  the two and the ratio is if anything an overestimate.  Also, f0 is the
  x → 0 *limit* of f_coh; the table is x-integrated over the window, so it maps
  onto f0 only through the assumed 1/(1+(x/x_coh)²) shape.

**Recommendation:** leave `f0 = 0.04` and its {0.02, 0.08} band as a scenario,
and record in `coherent.hpp` that eSTARlight brackets it from below at
1.0 × 10⁻³ (exclusive J/ψ) and from above at 3.0 × 10⁻² (all exclusive VMs) but
does not determine it.  Determining f0 needs a coherent *diffractive-DIS*
calculation (a coherent-A analogue of the H1/ZEUS diffractive PDFs), which is
neither eSTARlight nor Sartre; it is the same gap `physics_literature.md` §11
identifies.  `slope_b`, by contrast, is now settled (§4).

---

## 6. Caveats

1. **No clustering.**  eSTARlight has one spherically symmetric density per
   nucleus.  The α + d structure that carries the entire physics case for ⁶Li is
   absent, so nothing here bears on the α-tagged channel, on S/D interference,
   or on the near-null quadrupole moment Q(⁶Li) = −0.0806 fm².
2. **Gaussian density → no diffractive structure.**  §2c.  There are no
   secondary lobes and no first minimum, so these |t| distributions must not be
   used as an imaging baseline where lobe positions matter
   (cf. arXiv:2608.23445).  They are a *rate* and *slope* baseline only.
3. **Unpolarized.**  There is no polarization axis anywhere in eSTARlight, so
   nothing here constrains `amp`, `eps_b0`, `a2_deformation`, or any tensor
   observable.  Item 11's polarized half is untouched; arXiv:2408.13213 remains
   the only polarized coherent calculation in existence.
4. **No saturation.**  The paper says so explicitly of its own runs (p. 3:
   "no saturation effect is included").  Same here.
5. **Radius is the leading rate systematic**, ±30 % (§2b), and eSTARlight's
   own default (2.18 fm) is *not* the physically preferred value.
6. **⁶Li vs ⁷Li at the same energy is not the physical comparison.**  Rows 1–6
   of §2a share 99.5 GeV/u so that the A-scaling check in §2d is clean; the
   ⁷Li rows at 18 × 117.9 are the ones to compare with the paper.  LiPolGen's
   own ⁷Li top energy is 117.9 GeV/u, not ⁶Li's 137.5.
7. **eSTARlight's per-nucleon mass is m_p**, not the AME2020 nuclear mass per
   nucleon LiPolGen boosts with (`nuclear_mass(3, 6)`/6 = 5.60152/6 =
   0.93359 GeV for ⁶Li, `src/core/spectator.cpp:32`).  The γ values in §1
   were chosen to make the per-nucleon *momentum* match at 99.5 GeV/u
   exactly, which is the quantity every cross section depends on; the
   per-nucleon mass differs by 0.50 % (0.93827 vs 0.93359) and at fixed
   momentum the per-nucleon energy by < 10⁻⁵, both below every other
   uncertainty here.
8. **Statistics.**  2 × 10⁵ events per main run vs the paper's 10⁷.  σ is
   independent of `N_EVENTS` (it comes from the luminosity table), so the
   cross sections are exact to the 4 digits printed; only the |t| tail reach
   and ⟨W⟩ carry the MC error, and ⟨W⟩ is stable to < 0.5 % between the
   2 × 10⁵ and 10⁵ samples.

---

## 7. Reproduction

```bash
S=/tmp/estarlight            # any scratch dir
mkdir -p $S && cd $S
git clone https://github.com/eic/estarlight
cd estarlight && git checkout 939b11a24499398392d959db81c7502aeec91046 && cd ..
mkdir build && cd build && cmake ../estarlight && make -j8      # -> ./e_starlight

# optional: measured-radius variant
cd $S && cp -r estarlight estarlight_li
patch -p1 -d estarlight_li < <LiPolGen>/docs/open_items/run_2026-09-02/estarlight/nucleus_li_radii.patch
mkdir build_li && cd build_li && cmake ../estarlight_li && make -j8

# one run: 6Li J/psi at LiPolGen config 1
mkdir -p $S/run && cd $S/run
cp <LiPolGen>/docs/open_items/run_2026-09-02/estarlight/li6_jpsi_10x99.5.in slight.in
$S/build/e_starlight | tee run.log
grep 'Total cross section' run.log
```

The other eight input files in that directory differ only in
`TARGET_BEAM_A`, `TARGET_BEAM_GAMMA`, `ELECTRON_BEAM_GAMMA`, `PROD_PID` and
(for `q7_li6_jpsi.in`) `MIN_GAMMA_Q2`.

Slope, ⟨W⟩ and the |t| histogram come from `slight.out`, whose per-event lines
are `t: <t>` (signed, so |t| = −t) and `GAMMA: <E_gamma_target_frame> <Q2>`;
W = √(m_p² + 2 m_p E_γ − Q²).  The MLE slope on a truncated exponential over
0 < |t| < t_max is a two-line bisection on
n/B − n·t_max e^{−B t_max}/(1 − e^{−B t_max}) − Σ|t| = 0.

The inclusive denominator of §5:

```bash
<LiPolGen>/build/generate_inclusive --isotope 6Li --config 1 --events 200 --out /tmp/x.txt
# -> "azz0  0.33333  0  592031  0"   i.e. 5.92e5 pb
```

---

## 8. Verification log (2026-09-02, adversarial re-check)

Every σ in §2 and §5 was re-read from `runs/*/run.log` (`Total cross
section:` lines) and every B, ⟨|t|⟩, ⟨W⟩, median W, ⟨Q²⟩ and N was
regenerated with `ana.py` over the corresponding `slight.out`; all match the
tables as printed.  The unpatched build (`build/`, compiled 19:54:57 from a
clean checkout of `939b11a`) produced §2a and §5; the patched build
(`build_li/`, 19:59:25, `git status` shows only `src/nucleus.cpp` modified)
produced §2b — confirmed from the `printCompilerInfo()` line of each
`run.log`.  `generate_inclusive --isotope 6Li --config 1` reproduces
592 031 pb.  LiPolGen configuration 1 is 10 GeV × 99.5 GeV/u for both ⁶Li
and ⁷Li, and ⁷Li configuration 2 is 18 × 117.9 GeV/u (`generate_inclusive`
banner).  Paper attributions re-checked against the local PDF (v2, 25 Feb
2026; journal ref. PRD 113, 032018 confirmed on the arXiv abstract page).

Corrected in this pass: footnote-1 page (p. 2, not p. 3); event count of
the `MIN_GAMMA_Q2 = 0.7` runs (2 000, not 100 000); Fig. 7 page (p. 7, not
p. 6); the F(b) width (e^{−b²/(2B)}, 1/e at 1.85 fm — the earlier e^{−b²/B},
1.32 fm dropped the factor 2 from √(dσ/dt)); the description of Fig. 6's
axes and the |t| reach of its curves (counts 10²–4×10⁶ in ≈ 0.005 GeV² bins,
⁷Li drawn to ≈ 0.23 GeV², ³He to ≈ 0.45), replaced by a digitized-slope
comparison; the ⁶Li per-nucleon mass (0.9336 GeV, not 0.9316) and the size
and nature of the m_p mismatch (0.50 % in mass, not "0.7 % in energy"); the
maximum generated |t| per run; and the statement of where the B fit-range
comparison comes from (`ana.py`, not the eSTARlight log).
