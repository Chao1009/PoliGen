All three prototypes are built and produce numbers. Here is the report.

---

# 1. Coherent diffractive final state for e+⁶Li → e′ + X + ⁶Li(g.s.)

## What I found in PYTHIA 8.317

Two facts settle it, both verified in source and by running code:

**(i) `Diffraction:doHard` cannot be driven externally and is off for DIS.** `HardDiffraction::isDiffractive` (`src/HardDiffraction.cc:144-226`) invents the Pomeron *itself* a posteriori by comparing flux⊗PomPDF to the inclusive PDF; the Pomeron beams `beamPomA/B` are constructed only inside `BeamSetup.cc:563-569` and swapped in by `PartonLevel::setupHardDiff` — there is no setting pointing them at user input. Worse, hard diffraction is only initialised for `gammaMode != 4` (`PartonLevel.cc:269`) and each side is gated on `isHadron()` **or a resolved photon** (`PartonLevel.cc:401-414`). A virtual photon is neither, so γ*+p hard diffraction simply never fires. The one ep example, `examples/main342.cc`, is photoproduction (`Photon:Q2max = 1.`), not DIS. **Option (a) as stated is dead.**

**(ii) But `Beams:idA = 990` IS a legal user beam, and it works through LHAup.** `BeamSetup.cc:869-875` counts 990 as a hadron explicitly; `BeamParticle.cc:178` gives it meson-like beam handling (so the remnant is a single antiquark); `getPDFPtr(990)` returns a real Pomeron PDF. I initialised PYTHIA with `Beams:idA=990, idB=11, frameType=5` and it printed *"We collide Pomeron with e- at a CM energy of 63.246 GeV"*. All `PDF:PomSet` values I tried (1,2,3,6,11) init fine.

That opens a third option nobody listed, and it is the right one.

## The recommendation: γ*–Pomeron DIS through the existing bridge

The identification is **exact**, not an analogy:

| DIS bridge (today) | coherent bridge |
|---|---|
| target `p_N`, id 2212/2112 | `P_IP = P_ion − P_recoil`, id **990** |
| `W² = (q+p_N)²` | `M_X² = (q+P_IP)²` — the pipeline already computes it (`CoherentEvent::m_x2`) |
| `ξ = x_Bj` | `ζ = β = Q²/(M_X²+Q²)` — **exact**, because PYTHIA's Pomeron has `m0 = 0`, so `P_IP² = 0` |
| beam remnant = proton remnant (diquark) | beam remnant = Pomeron remnant (antiquark) → automatically colour- and charge-neutral |

So the whole surrogate + frame-map + conservation argument of `docs/PYTHIA_BRIDGE.md` §1–5 carries over verbatim. `PythiaBridge` grows a **third PYTHIA instance** (`idA = 990`) beside its proton and neutron ones, `hadronize()` grows one target branch, and `Pipeline::make_coherent` writes one new `Role::Pomeron` particle. Nothing else changes.

## Options table

| option | verdict | why |
|---|---|---|
| **(a) `Diffraction:doHard` / `PomFlux` / `hardDiffSide`** | **no** | not externally drivable (`HardDiffraction.cc:144`), and structurally off for virtual photons (`PartonLevel.cc:269, 401-414`). Resolved-photon photoproduction only. |
| **(a′) Pomeron BEAM through LHAup** ← *new* | **YES** | `BeamSetup.cc:870` accepts 990; verified init + 0-veto running. Reuses ~all existing bridge code, gets DGLAP ISR/FSR, real β dependence, gluon-dominated Pomeron, neutral remnant for free. |
| **(b) hand-filled q q̄ string + `forceHadronLevel`** | **fallback** | works perfectly (numbers below) but the flavour/topology is hand-imposed: no β dependence, no q q̄ g from backward evolution, no Pomeron PDF. Keep as the low-M_X / no-PDF safety net. |
| **(c) exclusive VM (ρ, φ, J/ψ)** | **yes, separately** | a genuinely different channel, not a competitor. It is the only thing that can fill `M_X < 1.2 GeV`, where both (a′) and (b) veto. Phase 2. |

## Numbers the prototypes produce

**`pom_dis.cc`** (option a′) — Q² = 5 GeV², 2000 ev/point, default PomSet:

| M_X | β | ⟨n⟩ | ⟨n_ch⟩ | ⟨Q_tot⟩ | M_had/M_X | veto | gluon frac. of IP |
|---|---|---|---|---|---|---|---|
| 1.5 | 0.690 | 5.04 | 2.24 | **0.0000** | 1.035 | 0 | 0.545 |
| 3.0 | 0.357 | 6.75 | 3.18 | **0.0000** | 1.008 | 0 | 0.755 |
| 5.0 | 0.167 | 9.70 | 4.54 | **0.0000** | 1.000 | 0 | 0.866 |
| 12 | 0.034 | 16.4 | 7.76 | **0.0000** | 1.000 | 0 | 0.920 |
| 40 | 0.0031 | 28.8 | 13.7 | **0.0000** | 1.000 | 0 | 0.934 |

Charge is exactly zero on every event; the hadronic invariant mass reproduces M_X. Throughput **≈ 40 000 ev/s** (20 k events in 0.53 s), matching the existing DIS bridge's 33 k/s. Veto rate is 0 for M_X ≥ 1.4, 6 % at 1.2, 29 % at 1.0, 75 % at 0.8 → **raise `COHERENT_MX_MIN_DEFAULT` from 1.0 to 1.2** (or route below it to (c)).

Two implementation traps I hit and fixed, both silent:
- **colour-tag orientation must follow the initiator's sign** — an incoming antiquark needs `(0,101)`, not `(101,0)`. Getting it wrong is a **50 % silent veto rate** (measured exactly 49–52 %).
- at `Q² = 1, β < 0.1` the LO Pomeron grid has **literally no quarks** (gluon fraction = 1.000, flavour weights all zero). The bridge's existing `q2_pdf_min` floor is mandatory here, plus an e_q² fallback.

**`mx_string.cc`** (option b) — γ = 60 boost, 2000 ev/point, then boosted onto P_X:

| M_X | ⟨n_ch⟩ qq̄ | ⟨n_ch⟩ qq̄+FSR | ⟨n_ch⟩ gg | max \|ΔP\|/E_X | max \|ΔQ\| |
|---|---|---|---|---|---|
| 1.5 | 1.85 | 1.85 | 1.80 | 8e-16 | 0 |
| 3 | 3.19 | 3.16 | 3.06 | 7e-16 | 0 |
| 8 | 6.10 | 6.52 | 7.15 | 1e-15 | 0 |
| 20 | 8.24 | 10.28 | 12.05 | 1e-15 | 0 |
| 40 | 10.05 | 14.25 | 15.73 | 1e-15 | 0 |

Recipe (from `examples/main234.cc:70-76`): `ProcessLevel:all = off`, `Check:event = off`, `event.reset()`, `append(id,23,101,0,…)` / `append(-id,23,0,101,…)`, optional `forceTimeShower(1,2,M_X)`, then `pythia.next()`. Momentum closes at **1e-15 relative**, charge exactly 0, 130 k ev/s. Below M_X ≈ 1.1 the ss̄ strings fail (17–19 % ≈ the 1/6 strange weight) — a flavour-dependent mass threshold is needed. Note `forceHadronLevel()` itself is documented as *not* handling resonance decays; `pythia.next()` with `ProcessLevel:all=off` is the better entry point (it routes to `forceHadronLevel` at `Pythia.cc:1057` but also fixes up `event[0]` and info).

## Effort

- **(a′) Pomeron-beam bridge: ~2–3 days.** Third PYTHIA instance + one `Role::Pomeron` in `make_coherent` + one target branch + flavour sampler on `getPDFPtr(990)` with the two traps above + tests. The conservation tests already exist; they just need the coherent channel enabled.
- **(b) as fallback: ~1 day**, mostly the flavour-threshold table.
- **(c) exclusive VM: ~1 week**, a genuinely new channel (VM t-slopes, spin-density-matrix, decay).

---

# 2. Triton remnant realism

## The headline result: the CS parameterization hands you the two-body/three-body split

BeAGLE's `DT_KFERMI` (`/home/cpeng/Projects/polli/BeAGLE/source/src/dpmjet3.0-5F-new.f:17216-17403`) implements

```
n(k) = Σ_{i=0,1,2}  A_i exp(−B_i k²) / (1 + C_i k²)²      [k in fm⁻¹, n in fm³]
```

| A | A₀ | B₀ | C₀ | A₁ | B₁ | C₁ | A₂ | B₂ | C₂ |
|---|---|---|---|---|---|---|---|---|---|
| 2 (d) | 157.4 | 1.24 | 18.3 | 0.234 | 1.27 | 0 | 0.00623 | 0.220 | 0 |
| 3 (³He **and** ³H) | 31.7 | 1.32 | 5.98 | 0.00266 | 0.365 | 0 | — | — | — |
| 4 (α) | 4.33 | 1.54 | 0.419 | 5.49 | 4.90 | 0 | — | — | — |

Cited in-code to *"Claudio Ciofi & S. Simula, PRC 53 (1996) 1689"* (line 17237, 17104). **There is no separate ³H set** — the branch is `IF (ANUCLEUS .EQ. 3)`, Z is never passed, so ³He and ³H are identical.

BeAGLE renormalizes numerically and so never notices what the A_i mean. I computed it (`nk_triton.py` §1):

```
∫₀^∞ n(k) k² dk  =  1.0031 (A=2)    0.6525 (A=3)    0.7996 (A=4)
```

The CS convention is `∫ n k² dk = 1` **without the 4π**, and the A=3 / A=4 sets are only the **n₀(k)** piece — the part where the residual (A−1) system is left in its **ground state**. That is exactly what the source comment says (*"These are n0k parametrization not including n1k"*), and it makes the deficit a physics number:

> **A=3: 65.3 % two-body breakup (bound remnant), 34.7 % three-body continuum.**
> A=4: 80.0 % / 20.0 %. A=2: 100.3 % — the deuteron has no excitable remnant, as it must be.

0.653 lands right on the ≈ 2/3 the ³He(e,e′p)d literature quotes for the p–d spectroscopic factor. **You do not need to go outside the parameterization you were already transcribing to get the split.** (I could not pin an exact published figure by search — Ciofi degli Atti & Kaptari PRC 71 024005 (2005) and the E89-044 papers do not put the number in an abstract — so treat 0.653 as the CS-internal value and 0.60–0.70 as the band.)

**Two BeAGLE consequences worth recording:** because `DT_KFERMI` renormalizes n₀ to unity, BeAGLE samples the A=3/A=4 Fermi momentum from the *ground-state-remnant distribution alone* and drops the 35 % / 20 % correlated continuum — its high-k tail is too soft by construction. Separately, `DT_DFERMI` (A>4) has an `.OR.`-for-`.AND.` chain making the C/O/Ca/Fe branches dead code, and uninitialised `SAVE`d `X0`/`CDF` that corrupt every call after the first.

## Tail comparison — the requested table

`P(k > k_cut)`, k in GeV (CS converted with ħc = 0.19733):

| distribution | k>0.1 | k>0.2 | k>0.3 | k>0.45 |
|---|---|---|---|---|
| **CS A=3 (³He/³H)** | 0.4266 | 0.06276 | 0.009621 | 0.002365 |
| CS A=2 (deuteron) | 0.3124 | 0.08288 | 0.03185 | 0.01405 |
| CS A=4 (α) | 0.6863 | 0.1665 | 0.01734 | 0.00012 |
| **Hulthén t\*→n+d** (LiPolGen) | 0.5633 | 0.1680 | 0.05246 | 0.01181 |
| **Hulthén t\*→p+nn** (LiPolGen) | 0.6246 | 0.2035 | 0.06571 | 0.01508 |
| Hulthén nn split (κ=10.4 MeV) | 0.06196 | 0.01255 | 0.003516 | 0.00075 |

**Ratio Hulthén(n+d) / CS(A=3): 1.32, 2.68, 5.45, 5.0.** ⟨k⟩ = 133 MeV (Hulthén n+d), 145 MeV (p+nn) vs **102 MeV** (CS A=3).

So the current Hulthén branch is **too hard by a factor 2.7 at 0.2 GeV and 5.5 at 0.3 GeV**. The cause is structural: κ = √(2μS) is 88.5 / 103.0 MeV for the two triton channels, driven by tiny separation energies (6.26 / 8.48 MeV), and the β = 0.30 GeV short-range term then supplies a power-law tail the real n₀ does not have. Note the CS n₀ tail is *itself* too soft relative to the true n(k) (it is missing n₁), so the true discrepancy at high k is smaller than the table shows — but the direction and the ⟨k⟩ mismatch stand.

## Recommended `TritonSpectralFunction`

Signatures are in `interfaces_sketch.hpp` §2. The design:

- **Model S_N(k, E)**, not a decay chain. Three channels, and the branching is **k-dependent**, `p_two_body(k) = n₀(k)/(n₀(k)+n₁(k))`, integrating to 0.653:
  - struck **n** → remnant `d` (bound, E=0) with weight n₀(k)/n(k) ≈ 0.65
  - struck **n** → remnant `(pn)` continuum with weight n₁(k)/n(k) ≈ 0.35 ← **this channel does not exist today at all**
  - struck **p** → remnant `(nn)` continuum, always (no bound nn)
- Keep the species draw as `Z·F2p : N·F2n` — that part is right.
- Keep the impulse-approximation rule (partners on shell, struck nucleon absorbs the difference) unchanged: conservation is independent of all of this.
- **Two implementations behind one interface**, matching the library's `Backend` convention: `CiofiSimulaTriton` (analytic, the coefficients above) and a table-driven Faddeev/AV18 `S(k,E)` when one is available.
- **`sample()` must consume a fixed number of uniforms whatever the branch**, so the RNG stream does not depend on the channel — the rule `CoherentXpomModel::draw` already follows.
- `BreakupOptions` grows one slot `std::shared_ptr<const TritonSpectralFunction> triton_sf`; null keeps today's behaviour bit-for-bit.

**Effort: ~2 days** for `CiofiSimulaTriton` + the third channel + tests; the n₁ shape is the one input you have to choose (the CS paper's n₁ coefficients are not in the BeAGLE source — only n₀ was transcribed).

---

# 3. FSI on the spectator α/d

## The formula

Cosyn–Weiss transplanted to X–α, in the form you asked for. The Glauber statement `Ψ_FSI(b) = ψ(b)[1 − Γ(b)]` is a convolution in momentum space:

```
Ψ_FSI(α_s, k_T) = ψ(α_s, k_T) − ∫ d²k′_T/(2π)²  Γ̃(k_T − k′_T)  ψ(α_s, k′_T)

Γ̃(q) = (σ_Xα (1 − iε) / 2) exp(−B_α q²/2)          [GeV⁻²],  Γ̃(0) = σ_Xα(1−iε)/2
Γ(b) = (σ_Xα (1 − iε) / (4π B_α)) exp(−b²/2B_α)
w_FSI = |Ψ_FSI|² / |ψ|²
```

Two structural points:
- The rescattering is **eikonal**: it transfers transverse momentum only, so `α_s` (equivalently `k_z`) is a spectator of the FSI too. The convolution is 2-D at fixed `k_z`. This is why the ratio → 1 at θ_k → 0, π and is worst at 90°.
- The longitudinal integral is done by residue (pole part). The principal-value piece is O(Δ/k) with **Δ ≈ k²/(2m_α)** = 0.3 / 1.3 / 5.4 / 12.1 MeV at k = 0.05 / 0.1 / 0.2 / 0.3 GeV, against κ(α–d) = **60.7 MeV**. So Δ ≪ κ and the eikonal limit is good — which is also why the ratio comes out symmetric under θ → π−θ. Independent support: Ciofi degli Atti & Kaptari (nucl-th/0407024) find GEA and plain Glauber differ by **< 3–4 %** for p_m < 0.6 GeV/c.

**σ_Xα is not 4σ_XN.** `fsi_alpha.py` builds it by Glauber from the X–N profile, `Γ_α(b) = 1 − [1 − (Γ_N ⊗ T_α)(b)]⁴`, with the α point-nucleon Gaussian `a² = (r_ch²(α) − r_ch²(p))/3 = 0.700 fm²`:

| σ_XN [mb] | σ_Xα [mb] | 4σ_XN | shadowing | B_α [GeV⁻²] |
|---|---|---|---|---|
| 10 | 38.1 | 40 | 0.953 | 24.8 |
| 20 | 72.5 | 80 | 0.907 | 25.6 |
| 30 | 103.5 | 120 | 0.862 | 26.4 |
| **40** | **131.0** | 160 | **0.819** | **27.2** |
| 50 | 155.5 | 200 | 0.778 | 28.0 |

## FSI/IA ratio for the ⁶Li α tag (W ≈ 10 GeV)

σ_XN = 40 mb, ε = −0.5, B_XN = 6 GeV⁻²:

| k [GeV] | θ=0 | 30 | 60 | 90 | 120 | 150 | 180 |
|---|---|---|---|---|---|---|---|
| 0.02 | 0.787 | 0.786 | 0.782 | 0.780 | 0.782 | 0.786 | 0.787 |
| 0.05 | 0.732 | 0.721 | 0.694 | 0.678 | 0.694 | 0.721 | 0.732 |
| 0.10 | 0.607 | 0.572 | 0.469 | 0.381 | 0.469 | 0.572 | 0.607 |
| 0.15 | 0.491 | 0.435 | 0.260 | 0.119 | 0.260 | 0.435 | 0.491 |
| 0.20 | 0.398 | 0.328 | 0.130 | 0.178 | 0.130 | 0.328 | 0.398 |
| 0.30 | 0.270 | 0.194 | 0.078 | **1.618** | 0.078 | 0.194 | 0.270 |
| 0.40 | 0.193 | 0.128 | 0.106 | **2.595** | 0.106 | 0.128 | 0.193 |

and on the variables the kernel actually uses:

| k_z \ k_T | 0.00 | 0.02 | 0.05 | 0.08 | 0.12 | 0.20 | 0.30 |
|---|---|---|---|---|---|---|---|
| 0.00 | 0.800 | 0.780 | 0.678 | 0.510 | 0.257 | 0.178 | 1.618 |
| 0.10 | 0.607 | 0.597 | 0.546 | 0.459 | 0.310 | 0.085 | 0.344 |
| 0.30 | 0.270 | 0.268 | 0.259 | 0.241 | 0.209 | 0.132 | 0.065 |

**Where pole dominance holds.** Not where the question's framing expects. Because the kernel is a function of `k_T` alone, the FSI-safe corner is **small k_T**, not "backward spectator" — the ratio is identical at θ and π−θ. Concretely, at σ_XN = 40 mb the distortion is ≤ 20 % only for `k ≲ 0.02 GeV`, ~27 % at 0.05, ~40–60 % at 0.10, and the amplitude *flips sign* around k_T ≈ 0.25 GeV where the rescattering feed-in from the large low-k strength overwhelms the steeply falling Hulthén tail (ratio > 1). At σ_XN = 20 mb everything is roughly halved (0.877 / 0.844 / 0.764 at k = 0.02 / 0.05 / 0.10, θ=0).

Two caveats I want on the record: the ~20 % suppression that survives at k → 0 is nuclear transparency (σ_Xα ≈ 131 mb makes the α nearly black), and the **σ_XN(W) formation-length ramp is the weakest input in the whole model** — at EIC energies the DIS debris has a formation length of many fm, so the effective σ is well below the free-hadron 40 mb and the 20 mb row is arguably the realistic one. This must be banded, not fixed.

**The useful operational conclusion:** the Roman-Pot α tag selects small p_T, i.e. small k_T, which is exactly the FSI-suppressed corner. The α-tag acceptance and the FSI-safe region coincide.

## The C++ hook

`interfaces_sketch.hpp` §3. **FSI must enter as a weight, never as a momentum shift** — moving the spectator would break both the "never recoil-correct the light spectator" rule and the exact conservation the whole record rests on.

```cpp
class FsiWeight {
  virtual double weight(const FsiKinematics& kin) const = 0;             // FSI/IA
  virtual double weight_normalised(const FsiKinematics& kin) const = 0;  // shape only
  virtual double sigma_eff_mb(double w) const = 0;                       // what was used
};
struct FsiKinematics { double k, cos_theta_k, w, q2, x; int spectator_z, spectator_a; };
```

`TaggedSampler` gains `set_fsi(std::shared_ptr<const FsiWeight>)` and a new `TaggedEvent::weight` field (default 1); `Pipeline` multiplies it into `Event::weight` exactly as the coherent channel already does with its azimuthal weight. Null pointer = today's PWIA, bit-for-bit. Evaluate on a `(k, cos θ_k)` grid at construction and read by bilinear interpolation per event — the same shape as `acceptance_weights`, so the throughput target is untouched. Expose `weight_normalised` separately so a run can take the shape distortion without changing the total rate.

**Effort: ~2–3 days** for `GlauberFsiWeight` + the hook + tests, given the Python is already the reference implementation. The σ_XN(W) ramp is the open physics input.

---

# Prototype files

All in `/tmp/claude-1000/-home-cpeng-Projects-polli/ad581d75-f90b-4f90-9f7b-196cb9e9bda5/scratchpad/design/`:

| file | what it produces |
|---|---|
| `pom_beam_probe.cc` | proves `Beams:idA=990` + `frameType=5` initialises; dumps the Pomeron PDF (xg = 0.746 at β=0.05, Q²=10) |
| `pom_dis.cc` | the **recommended** coherent T2: ⟨n_ch⟩, ⟨Q⟩=0, M_had/M_X=1, veto rate, gluon fraction vs M_X at Q² = 1 / 5 / 30; 40 k ev/s |
| `mx_string.cc` | the fallback: qq̄ / qq̄+FSR / gg string hadronization, ⟨n⟩ vs M_X, ΔP = 1e-15, ΔQ = 0, 130 k ev/s |
| `nk_triton.py` | CS coefficients transcribed + the 0.6525 / 0.7996 normalizations (= the 2-body/3-body split) + the P(k>cut) table |
| `fsi_alpha.py` | σ_Xα Glauber table, FSI/IA vs (k, θ_k) and (k_z, k_T), the dropped Δ vs κ |
| `interfaces_sketch.hpp` | proposed `CoherentBridgeOptions`, `TritonSpectralFunction`, `FsiWeight` signatures |

Build the C++ ones with `source LiPolGen/env.sh && g++ -O2 -std=c++17 X.cc -o X $(pythia8-config --cxxflags --ldflags)`. Nothing outside the scratchpad was touched.

**Sources consulted for §2/§3 literature:** [Ciofi degli Atti & Kaptari, GEA analysis of ²H/³He(e,e′p)](https://arxiv.org/abs/nucl-th/0407024) · [JLab E89-044, exclusive ³He(e,e′p) below the QE peak](https://arxiv.org/pdf/nucl-ex/0310021) · [microscopic spectral function with 2N/3N SRC vs ab initio for A=3](https://ar5iv.labs.arxiv.org/html/1701.08211)