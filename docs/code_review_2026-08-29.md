# Adversarial code/physics review — 2026-08-29 (v0.1 pre-release)

Independent review of the C++ core against the Python originals. Findings C1–C6/P1–P8 were dispatched for fixing the same day; see git log.

## Verdict

The port is careful and unusually well documented; the doctest suite (198 cases, 6.1 M assertions) passes, throughput is 1.5–3.5 M ev/s (target was 10⁵), and the cross-language reference gate (`validation/reference/*.json` at rtol 1e-12) is real, not a self-comparison. The defects below are mostly **integration-layer and stale-default** problems, not kernel algebra. I found nothing wrong with the m-ordering, φ' = φ − φ_S, `TENSOR_LL_SIGN`, the θ_S/φ_S rotation, the λ_e sign, the J = 3/2 branch, the inverse-CDF indices, or the PYTHIA Lorentz map.

---

# CONFIRMED

### C1 — Coherent channel emits a **spacelike** hadronic X in 100 % of events; the negative mass² is silently clipped
`src/core/pipeline.cpp:764` (`add_hadronic_x` call), `src/core/pipeline.cpp:620` (`x.mass = std::sqrt(std::max(p_x.m2(), 0.0))`), `include/lipolgen/coherent.hpp:197` (`x_pom_ = 0.0`).

`PipelineConfig` has **no `x_pom` knob** and `CoherentSampler::set_x_pom` is never called by `Pipeline`, so x_P ≡ 0 for every coherent event. Then P_recoil differs from P_ion only by p_T = √|t|, so
X = k + P_ion − k′ − P_recoil ≈ q, and M_X² ≈ −Q².

Measured (6Li, 10 × 99.5 GeV/u, 200 events, defaults):

```
[coherent] 200 events, 200 with M_X^2 < 0, worst M_X^2 = -35.67 GeV^2 (Q2 there 35.4)
           event0: Q2=6.805  M_X^2=-6.502  stored mass=0  t=0.00479  x_pom=0
```

Failure scenario: `examples/generate_coherent.cpp --out ev.hepmc` writes a status-1 `pdg 92` particle with `set_generated_mass(0)` (`src/hepmc/hepmc_writer.cpp:130`) whose four-vector satisfies E² − p² = −6.5 GeV². Any consumer that recomputes the mass gets a tachyon; a diffractive-mass analysis gets nonsense. `tests/test_pipeline.cpp:108` passes because it only checks Σp and Σq, and `tests/test_pipeline.cpp:530` actively *asserts* `ev.kin.x_pom == 0.0`, freezing the bug in as expected behaviour.

Fix: sample x_P from the diffractive mass. Given a target M_X (or a β = x/x_P draw), x_P = (M_X² + Q²)/(W² + Q²); set `p_z,recoil = (1 − x_P) A p_u` accordingly. Minimum viable guard: add `x_pom` to `PipelineConfig`, wire `set_x_pom`, and assert `p_x.m2() >= 0` in `add_hadronic_x` instead of clipping.

---

### C2 — `InclusiveKernel` default `target_mass` is `false`; the Python default has been `True` since 2026-08-29
`include/lipolgen/xsec.hpp:123` vs `evgen/polligen/xsec.py:149`.

Both paths exist and both are pinned against the reference JSON (`tests/test_reference.cpp:52`), so this is invisible to the suite — but `InclusiveKernel(ion)` in C++ and `InclusiveKernel(ion)` in Python are different physics. `default_inclusive_kernel` (`src/core/pipeline.cpp:280`) also leaves it off, so the shipped generator runs the massless A∥. Measured ratio A∥(exact)/A∥(massless), 6Li at 10 × 99.5:

```
  x      Q2      ratio
0.005    2.0    1.00005
0.050    2.0    1.00438
0.200    2.0    1.07046
0.200   10.0    1.01411
0.500    2.0    1.44045
0.500   10.0    1.08823
```

Well beyond the "0.1–1 % over the published x range" the Python quotes, because the *generator* window (W² ≥ 8, Q² ≥ 0.7) reaches further into high-x/low-Q² than the analysis window. Fix: flip the default to `true` and mirror the Python's g2 handling (next item), or state the divergence in `CONVENTIONS.md`.

### C2b — C++ refuses a configuration Python explicitly permits
`src/core/xsec.cpp:63-67` throws for `target_mass && g2_mode == kZero`. `xsec.py`'s docstring (lines 126-129) states this is *"legitimate and is no longer refused"* — it is exactly the `g2_scale = 0` twist-3 variation `evgen/scripts/target_mass_bound.py` needs. The C++ also has **no `g2_scale`** at all, so the 1.5 × g2^WW half of that systematic is unreachable. Fix: allow the combination (Python zeros g2 in `tables()` and proceeds), and add `g2_scale` to `Options`.

---

### C3 — `TaggedModel`'s three `mutable` caches have no mutex; safety rests entirely on the Pipeline warm-up
`include/lipolgen/tagged.hpp:157-159` (`amp2_`, `n_`, `cdf_`), lazily filled in `src/core/tagged.cpp:174, 186, 313`.

`InclusiveSampler` guards its equivalent cache with `cache_mutex_` (`sampler.cpp:262`); `TaggedModel` does not. `Pipeline`'s constructor (`src/core/pipeline.cpp:408-427`) warms every reachable `(M, m_S)` key, and I verified the coverage is exact: `rates()` skips `p_ms[b] <= 0` and the warm loop warms every `p_ms[b] > 0`, so the threaded pipeline path is safe *today*.

It is not safe for anyone using the module directly: `TaggedSampler::sample_category` / `sample_one` on a fresh `TaggedModel` from two threads inserts concurrently into `std::map` → UB. Nothing in the header says the model must be warmed first, and the comment that documents this is buried in `pipeline.cpp`, not in `tagged.hpp`. One added wave, one changed `j_ion`, or one direct API user reintroduces the race. Fix: add a `std::mutex` (the cost is one uncontended lock per event on a warmed map) or build all tables eagerly in the constructor.

---

### C4 — The T2 hadronizer hook is invoked on **every** channel, including coherent
`src/core/pipeline.cpp:791`: `if (cfg_.hadronizer) cfg_.hadronizer(out, rng);` with no channel guard.

A coherent event carries no `Role::StruckNucleon` and no `Role::StruckCluster`, so `PythiaBridge` falls into its inclusive branch (`src/pythia/pythia_bridge.cpp:365-381`), invents a nucleon at P_ion/A, and appends a γ*N final state on top of the intact recoil. The whole-nucleus balance then fails by ≈ P_ion(1 − 1/A) — a ~500 GeV residual, not a rounding error. No test binds a hadronizer to the coherent channel. Fix: refuse (or skip) the hook for `PipelineChannel::CoherentLi6` in `PipelineConfig::validate()`.

---

### C5 — `unpolarized_emc_ratio` / `tmt_valence_scale` are frozen at the pre-2026-08-29 Python baseline, presented as fact
`src/core/sf.cpp:237` (12-point table), `sf.cpp:245` (`cbt_valence_scale()` returns a hard-coded `1.0`), `sf.cpp:247-259`.

Python's `unpolarized_emc_ratio` default mode is now `'epps21'`, not `'table'`; `valence_scale`'s default baseline is now `'epps21'`, giving **0.5106** for CBT and **0.2027** for TMT, against the C++'s hard-coded 1.0 and 0.397009 (`polarized.py:484-513` documents both sets). Measured C++ output:

```
unpolarized_emc_ratio(0.5) = 0.936667      (python mode='table'; python default 'epps21' differs)
tmt_valence_scale() = 0.397009   cbt_valence_scale() = 1.000000
tmt_polarized_emc_ratio(0.5) = 0.940021
```

That is a factor ~2 on the transferred polarized-EMC depletion — the headline observable of `money_polemc.py`. Freezing at the table baseline is defensible (the epps21 mode needs LHAPDF nuclear grids and the core must not link LHAPDF), but it is neither documented as a divergence nor selectable, and `cbt_valence_scale()` returning a literal `1.0` is exactly the "physics number in two places" `CONVENTIONS.md` forbids. Fix: add an `EmcBaseline` enum with the legacy value as an explicit choice, and compute `cbt_valence_scale` rather than asserting it.

---

### C6 — `Optics::lumi_fraction` never reaches the event count
`src/core/pipeline.cpp:235` sets it; nothing in `Pipeline` reads it (verified by grep — only `examples/generate_tagged.cpp:243` applies it, in a printed column).

Running at fixed `lumi_pb` with `OpticsChoice::Tagging` gives the de-squeezed *acceptance* without the de-squeeze's luminosity penalty (`lumi_fraction` = 0.067–0.147 in `kTaggingPerConfig`). A user comparing tagging vs Yellow-Report optics through the `Pipeline` API over-states the tagged yield by 7–15×. `tests/test_pipeline.cpp:302` only checks the *value* of the field, never that it is applied. Fix: either fold it into `lumi_[k]` or document loudly in `pipeline.hpp` that it is the caller's job.

---

### C7 — The PYTHIA frame map is correct (verified), including the cos 2φ-relevant structure
Not a defect — recorded because you asked. I generated the same event twice with an identical RNG stream, changing only the e′ azimuth (0.3 → 2.4 rad), and measured each hadron's azimuth about q̂ referred to the lepton plane:

```
case 0: phi_e'=0.300  det(e1,e2,e3)=+1.000000000000000  lambda=0.999999999999988
case 1: phi_e'=2.400  det(e1,e2,e3)=+1.000000000000000  lambda=0.999999999999988
n_had=18  same pdg order=1  worst |d(phi_h - lepton plane)| = 4.44e-16 rad
```

Both triads are right-handed (det = +1), so R = Σ e′ᵢ ⊗ eᵢ is a **proper** rotation — the azimuthal sign is not flipped, and the hadronic azimuthal structure about γ* relative to the lepton plane is transported exactly. I also re-derived the surrogate algebra (`pythia_bridge.cpp:400-406`), the ζ solve, and `dis_parton_fraction`; all three are right, and the rationalised ξ form is genuinely the stable branch. **No tensor-weight double counting**: `Impl::run` never touches `ev.weight` or `ev.xsec_pb`; the polarized physics stays entirely in the core sampler.

**Neutron handling is correct** — I checked PYTHIA 8.317 itself: `PDF::resetValenceContent` gives 2112 `nu=1, nd=2 → beamType = -1`, and `PDF::xf` then returns `xu` for `id == 1` (`src/PartonDistributions.cc`), so the u↔d swap the doc claims really happens, and keeping two instances is the right call for the remnant diquark.

---

# PLAUSIBLE

### P1 — The struck nucleon is drawn **flat** Z : N, not cross-section weighted
`src/pythia/pythia_bridge.cpp:378` (`rng.uniform() * ion.A < ion.Z`) and `:360` for the cluster case. The inclusive rate is Z·F2p + N·F2n; with F2n/F2p ≈ 0.63 at x = 0.5 the true split for 6Li is 61:39, not 50:50. The *cross section* is right (the kernel builds F2A correctly); only the T2 final-state flavour/charge composition is biased, growing with x. `docs/PYTHIA_BRIDGE.md` §9 says "the nucleus enters only through (a) which nucleon is struck" but never says the draw is unweighted. Fix: weight by `f2p : f2n` at the event's (x, Q²), or document.

### P2 — Charm is offered at the wrong ζ
`pythia_bridge.cpp:428, 441`: the flavour weight uses `xf(id, zeta_light, q2pdf)` for every flavour, but the ζ actually handed to PYTHIA is `zeta_of(m_q)`, which for charm is (Q² + m_c²)/(Q²) ≈ 2.1× larger at Q² = 2 GeV². Since xf falls steeply with x, charm is over-offered at low Q²; when the true ζ_c ≥ 1 the retry loop just `continue`s (`:475`) and burns a retry without removing charm from the pool. Fix: evaluate each flavour's weight at its own ζ.

### P3 — Two different "implicit struck nucleon" masses in one library
`src/core/generator.cpp:34` uses `M_NUCLEON` (0.9383, free) while `src/pythia/pythia_bridge.cpp:369` uses `ion.mass_per_nucleon()` (0.9338 for 6Li), and `docs/PYTHIA_BRIDGE.md` §6 argues explicitly for the bound mass. Only the bridge branch fires when no `Role::StruckNucleon` exists, so the Pipeline inclusive path is self-consistent — but the two definitions coexist and disagree by 4.5 MeV. This is the exact class of defect `plans/08 C1` was written about.

### P4 — No positivity guard on the coherent azimuthal density
`src/core/coherent.cpp:89-91`. c₂ = −(P_zz/2)·ε_B0·B·|t| + amp·P_zz is unbounded in |t|, and `coherent_t_max` defaults to 0.5 while the model's own docstring says the linear-in-|t| form is valid only for |t| ≲ 0.2. Measured:

```
defaults, Pzz=1.0                    c2(|t|=0.5) =  1.0100   min weight = -0.0100
band edge eps=-0.13 B=60, Pzz=1      c2(|t|=0.5) =  1.9600   min weight = -0.9600
default (eps=-0.08,B=50,Pzz=1): c2 crosses 1 at |t| = 0.4950
```

In weight mode this yields **negative event weights**; in `weighted_azimuth` mode the rejection silently samples max(W, 0). In practice the exponential |t| makes it a ~2 × 10⁻⁷ event at the band edge (200 000 sampled events gave zero), but `InclusiveSampler::build_state` (`sampler.cpp:247-260`) *throws* on exactly this condition — the two channels are inconsistent. Fix: mirror the inclusive guard, or clamp `coherent_t_max` to the model's validity.

### P5 — `pot_levers("5x41").r34` is transcribed as −1.0; Python has 4.56
`src/core/spectator.cpp:223` vs `farforward.py:365`. `r34` is unused today (`over_rigid_route` destructures it away in both languages), so no number moves — but it is a wrong data value sitting in a table, and `separation_at_pots` is the obvious next consumer.

### P6 — Per-event heap allocations in the tagged hot loop
`src/core/tagged.cpp:554` (`m_values(...)` builds a `std::vector` per event) and `:569` (four `std::vector<double>` constructed inside `sample_one`, each resized to 1). Five allocations per tagged event; tagged throughput is 1.51 M ev/s against inclusive's 3.47 M ev/s. Well above the 10⁵ target, so cosmetic — fix by hoisting `m_values` into the `TaggedModel` and using a small fixed buffer in `KinematicsSource::sample`.

### P7 — Minor `np_interp` fidelity gap
`src/core/numerics.cpp:116` falls back on `!isfinite`; NumPy's `compiled_interp` tests `isnan` only. A ±inf interpolant would take the fallback branch here and not in NumPy. Unreachable with the current tables.

### P8 — `amplitudes()` gates on `state.j` but `tensor_moments()` uses `ion_.spin`
`src/core/xsec.cpp:167` vs `:150-155`. For the deuteron control channel (S_c = 1, DIS target = free neutron, spin 1/2) the rank-2 branch is entered and then returns (0, 0) — correct by accident. This is a faithful port of the same structure in `xsec.py:332/309`, so it is not a regression, but a channel with S_c ≠ target spin and a non-zero b1 would produce a silently wrong geometry.

---

# What the tests do not cover

1. **Nothing asserts `HadronicX.p.m2() >= 0`** on any channel — the reason C1 survived. `tests/test_pipeline.cpp:530` asserts `x_pom == 0.0`, i.e. it pins the bug.
2. **No hadronizer test on the coherent channel** (C4).
3. **No test that the C++ kernel default equals the Python kernel default** — `tests/test_reference.cpp:52` builds `target_mass` explicitly per reference block, so C2 is structurally invisible.
4. **No thread-safety test outside the Pipeline.** `tests/test_pipeline.cpp:637` passes only because the constructor warms the maps; a direct `TaggedSampler::sample_category` from two threads is untested UB. No TSAN run in the suite.
5. **No test that `Optics::lumi_fraction` is applied**; `test_pipeline.cpp:302` checks only its value.
6. **No test pinning the EMC baseline against Python's current default** (C5).
7. **No positivity test for the coherent azimuth**, though the inclusive one has `tests/test_sampler.cpp:236`.
8. Tautological by construction (fine as sampler self-consistency, but they prove no physics): `test_tagged.cpp:587` "sample_kc reproduces its own density", `test_pipeline.cpp:185` and `:341` (generated sample vs the same `TaggedModel` tables), `test_sampler.cpp:92`, `test_xsec.cpp:394` ("r_func = r_sigma_lt is bit-for-bit r_func = none"). `test_pipeline.cpp:594` is a refactor guard, not a physics check. The genuine external anchors are `test_reference.cpp`, `test_lhapdf.cpp` and `test_pythia.cpp:515`.

---

# Solid

- **`spin.cpp`** — a line-for-line port including the Racah k-ranges, the `exp(β(m − m_max))` maxent stabilisation, and NumPy's clip semantics. I found no divergence.
- **`numerics.cpp`** — NumPy's block-128 pairwise summation and `np.interp`'s exact-node short-circuit reproduced correctly; this is what makes the rtol-1e-12 gate meaningful rather than aspirational.
- **Sign/convention discipline.** `TENSOR_LL_SIGN` appears in exactly two places and they are mutually consistent (Azz thirds → +(2/3)b1/F1); φ' = φ − φ_S in both `dsigma` and `draw_event`; R_z(φ_S)R_y(θ_S) maps ẑ → n̂ correctly; m = +J…−J everywhere; the J = 3/2 rank-2 branch uses one Q_NN geometry for both spins and is test-pinned.
- **`density_min`** — I re-derived it; the vertex/endpoint case split is exact, and it is a genuine improvement on the 1 + |A| + |B| envelope.
- **`struck_cluster`** — the spectator is on shell to machine precision (E²−p² = m_spec² identically under the longitudinal boost), P_X = P_ion − p_spec closes the whole-nucleus balance exactly, and charge closes as −1 + Z_spec + Z_partner.
- **The PYTHIA frame map** — verified proper, azimuth-preserving, λ = 1 to 1.2 × 10⁻¹⁴ (C7). The surrogate algebra, the ζ solve, `dis_parton_fraction`, and the two-instance neutron design are all correct; the c/b mass fix and `SpaceShower:QEDshowerByQ = off` are real, non-obvious catches.
- **Luminosity bookkeeping** — the share genuinely reaches counts once and cross sections never (`sampler.cpp:154`, `pipeline.cpp:536`), and `test_pipeline.cpp:693` proves it.
- **Determinism** — `Pipeline::event(i)` is a pure function of the index; the threading construction is correct and the counter-based RNG has no stream reuse anywhere I could find (`kCountStreamEvent = ~0ULL` is properly reserved out of the event range).