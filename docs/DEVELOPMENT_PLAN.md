# LiPolGen — development plan (2026-08-29)

**LiPolGen** = doubly polarized e + ⁶Li / e + ⁷Li DIS event generator for the
EIC, C++ core with a Python interface, built as a **library on top of stock
PYTHIA 8.317** (no fork) and designed as a drop-in, faster replacement for the
numpy generator `PolarizedLithiumSim/evgen/polligen`.

## 0. Feasibility verdict (from the three surveys in `docs/surveys/`)

| candidate | verdict | reason |
|---|---|---|
| fork **BeAGLE** | **no** | unlicensed F77, closed-source FLUKA link, A ≥ 5 is a spherical Woods–Saxon ball (no α+d / α+t), silent Fermi-motion bug for 5 < A ≤ 56 (`.OR.` for `.AND.`), coherent recoil impossible by construction, zero spin support. Mine it for parameters only (Ciofi–Simula n(k), LF spectator prescription, "never recoil-correct the light spectator" rule, output columns). |
| fork **PYTHIA 8 core** | **no** | GPL allows it but MCnet guidelines discourage forks; nothing we need requires a core patch once the hard process is injected through `LHAup`. Core has no spin, no nucleus API for e+A (any nucleus id diverts into Angantyr, which has no lepton physics), and `Beams:allowMomentumSpread` runs the hard process at the *initialisation* √s (`ProcessLevel.cc:645`), which is wrong for ±20 % Fermi-smeared per-nucleon momenta. |
| **new C++ library + stock PYTHIA 8 as hadronizer** | **yes** | the physics that does not exist anywhere (polarized-nucleus vertex, spin ⊗ cluster spectator correlation, coherent channel) is already *specified and validated* in `polligen` (305 tests, rtol 1e-12 identities, external anchors). PYTHIA supplies showers/hadronization through `LHAup` (`Beams:frameType = 5`) with exact per-event four-momenta; spectators are appended post-`next()` (precedent: `Pythia8Plugins/PythiaCascade.h`). |

So: **doable**, and the dominant risk is physics input (VMC overlaps, spin-3/2
rank-2 basis, coherent amplitude for A > 2), not engineering. The generator is
designed so those inputs are swappable tables behind interfaces.

## 1. Architecture

```
 lipolgen (C++17, libLiPolGen.so)                                python: import lipolgen
 ┌────────────────────────────────────────────────────────────┐
 │ core/    no external deps                                  │
 │   spin        ρ(m,m′) J=1, 3/2 · Wigner-d · CG · moments   │
 │   sf          F1,F2,R,g1,g2(WW),b1,b2,Δ backends            │◄── LHAPDF6 (CT18, NNPDFpol, EPPS21) optional
 │   xsec        master formula dσ/dx dQ² dφ (HJM ⊕ CW)        │
 │   beams       γ-matched EIC configs, physical masses        │
 │   sampler     grid inverse-CDF (x,Q²) · φ accept/reject     │
 │   cluster     wave functions S/P/D (Hulthén) · VMC tables   │
 │   tagged      |A_{m_S}(M;k,k̂)|² · struck populations · boost │
 │   coherent    f_coh(x)·e^{-B|t|}·cos2φ(deformation, glue)   │
 │   bookkeeping SpinCategory/RunPlan · per-(run,bunch) RNG    │
 │   event       Event record T0/T1: e′, spectators, struck X  │
 │   io          HepMC3 Asciiv3 writer (+ spin attributes)     │◄── HepMC3
 ├────────────────────────────────────────────────────────────┤
 │ pythia/  T2 tier                                           │
 │   LhaupBridge  per-event e + q(x P_N) hard process → LHAup  │◄── libpythia8
 │   Hadronizer   showers/hadronization, spectator append      │
 ├────────────────────────────────────────────────────────────┤
 │ python/  pybind11 module + numpy exporters                 │
 │   HFSSample .npz · inclusive/tagged dict schemas (polligen) │
 └────────────────────────────────────────────────────────────┘
```

Rules:
- Core has **no PYTHIA/HepMC dependency**; everything heavy is C++, Python is
  orchestration + numpy views (zero-copy).
- All conventions are **named constants with a single definition**
  (`TENSOR_LL_SIGN = +1`, `LI6_EFF_POL`, rank-2 geometry `Q_NN`), mirrored
  from `polli_fastsim/asymmetries.py`. Flipping is a deliberate act with a test.
- Every physics input (b₁ curves, Δ ansatz, EMC ratio, cluster radial
  functions) is a `Backend` interface with a *toy* implementation and a
  *table* implementation, so VMC overlaps / new theory land as data files.
- Reproducibility: counter-based RNG stream per (run, bunch, event) — same
  seed gives the same event on any thread count.

## 2. Physics scope (what is implemented, by tier)

| tier | content | source of truth |
|---|---|---|
| T0 | e′ (x,Q²,y,φ), spin labels (λ_e, M, m_S), spectator 4-vector, struck cluster as pseudo-particle | `polligen/{spin,xsec,sample,tagged,coherent,bookkeeping}.py` |
| T1 | struck-cluster internal nucleon + partner spectators (d* → p+n? no: α+p from d*; α+d/nn from t*), physical AME masses | `plans/05` 5.D, `spectator.py` masses |
| T2 | PYTHIA 8 hadronization of γ*–nucleon with the nucleon's *actual* Fermi-smeared momentum; spectators appended | `tools/pythia8/gen_dis_hfs.py` settings incl. the two silent-cut fixes |

Explicitly out of scope for v0.1 (documented interfaces left): FSI, tensor-sector
radiative corrections, FLUKA-grade evaporation, ⁷Li coherent scenario, polarized
lepton QED radiation (PYTHIA dipole-recoil limitation).

## 3. Work breakdown and agent assignment

**Status 2026-08-29 (end of day 1):** P0–P7 done and merged (215 doctest cases / 9.0 M
assertions, 114 pytest cases); adversarial review done (`docs/code_review_2026-08-29.md`),
all 12 findings fixed with tests. T1 tier (struck cluster → nucleon + partner spectator)
in progress; see §6 for what remains open.

Model policy: **Fable 5** = planning, architecture, review gates, merges.
**Opus 5** = the physics-bearing C++ modules. **Sonnet 5** = tooling, bindings,
scripts, docs, reference dumps.

| # | task | owner | depends on | gate |
|---|---|---|---|---|
| P0 | plan, repo skeleton, CMake, deps (HepMC3 3.3.0, LHAPDF 6.5.5, PYTHIA 8.317 → `../deps/install`) | Fable | – | builds |
| P1 | `validation/dump_polligen_reference.py`: JSON tables of ρ moments, kernel amplitudes (w_avg, a₁, a₂), asymmetries, tagged densities, boost outputs at fixed points, from the Python code | Sonnet | – | script runs in the existing env |
| P2 | core: `spin`, `sf` (toy + LHAPDF), `xsec`, `beams` + doctest suite vs P1 tables at rtol 1e-12 | Opus | P0,P1 | identities pass |
| P3 | core: `sampler` (inclusive, weighted mode), `bookkeeping`, RNG, `event` record | Opus | P2 | estimator closure vs analytic errors (15 %) |
| P4 | core: `cluster`, `tagged`, `coherent` + Cosyn–Weiss deuteron gate, ⁷Li P₂ = −T/5 polarimeter, ⁷Li P_p = 0.866 | Opus | P2 | CW TABLE II (+1/−2) |
| P5 | `io`: HepMC3 Asciiv3 writer with spin attributes (status-4 beams, 10-digit ion codes, omitted origin vertex) + reader smoke test | Sonnet | P3 | pyhepmc reads it back; 4-momentum/charge conserved |
| P6 | `pythia`: LHAup bridge (struck nucleon at Fermi-smeared momentum, SPINUP set), hadronizer, spectator append, HFS truth identities Σ(E−p_z) = 2E_e y + m²/(E_N+p_N), Σp_T = p_T,e | Opus | P3, PYTHIA built | identities to 1e-3 / 1e-10 |
| P7 | pybind11 module, numpy exporters (`HFSSample` .npz, inclusive/tagged dicts), CLI `lipolgen-run` + YAML/JSON config | Sonnet | P3–P6 | polligen scripts consume the output |
| P8 | review, closure vs `polligen` end-to-end (FOM maps), performance target ≥ 10⁵ T0 ev/s single core, docs | Fable | all | closure report |

## 4. Validation matrix (must hold before v0.1)

1. ρ moments exact (atol 1e-12), J = 1, 3/2, arbitrary axis; max-entropy fills
   reproduce (P_z, P_zz) = (8/13, 4/13) and (0.7, 0.4).
2. Master formula sector identities at rtol 1e-12 against P1 tables:
   (w₊−w₋)/(2+w₊+w₋) = P_e A∥; thirds → A_zz; a₂ → A_cos2φ; F₂ cancels in a₂;
   magic angle kills tensor rate shift; unpolarized limit.
3. Estimator closure: pulls unbiased, spreads within 15 % of
   err_A∥ = 1/(P_e P_z √N), err_Azz = √(2/N)/P_zz, err_cos2φ = √(2/N)/P_zz;
   relative-luminosity biases −(2/3)δ/P_zz and δ/(2P_eP_z) reproduced and removed.
4. Tagged: Cosyn–Weiss deuteron limit (P₂ factorization < 1e-5, peak k = 0.31 GeV,
   A_T∥ extremes +1/−2); ⁶Li S-wave = inclusive deuteron; ⁷Li ⟨P₂⟩ = −T/5;
   P_p = 0.866 and P_n ≈ −0.037 (the neutron half is an *open* gate, report it).
5. Spectator boost matches `spectator.py` quantile-by-quantile at P_D = 0;
   rigidity R(⁶Li α) = 0.99813, R(⁷Li α) = 0.85571.
6. Coherent: ⟨t⟩ = 1/B, acceptance = exp(−B c²), m-state relation a₀ = −2a₁.
7. HFS truth identities on T2 output; HepMC3 round trip; PYTHIA cross section
   matches `gen_dis_hfs.py` at identical settings (same 8.3 physics).
8. Deterministic: same (seed, run, bunch) → identical events at any thread count.

## 5. Conventions carried over verbatim

- Frame: head-on, ion along +z at p_u GeV/u, electron along −z; crossing angle
  25 mrad applied only in the lab transform (`reco.py`).
- Beam energies γ-matched: ⁶Li 40.8 / 99.5 / 137.5 GeV/u, ⁷Li 40.8 / 99.5 / 117.9.
- Populations ordered m = +J … −J. `TENSOR_LL_SIGN = +1` (plans/08 D1 pending).
- Generator window looser than analysis: Q² ≥ 0.7, y ∈ [0.004, 0.985], W² ≥ 8.
- PYTHIA: `WeakBosonExchange:ff2ff(t:gmZ)`, `SpaceShower:dipoleRecoil=on`,
  `PDF:lepton=off`, `TimeShower:QEDshowerByL=off`, `PhaseSpace:pTHatMinDiverge=0.5`,
  `PhaseSpace:mHatMin=0.5` (both silent-cut fixes); with LHAup these become
  shower settings only, the hard phase space is ours.
- Ion spin in HepMC3: attributes `spin_J`, `spin_M`, `spin_axis_theta/phi`,
  `P_e`, `lam_e`, `P_z`, `P_zz`, `struck_cluster_m` on the GenEvent (proposed
  convention, plans/04 #17).

## 6. Open items after day 1

- **T1 breakup realism**: deuteron internal wave function is Hulthén S+D; the triton
  remnant (d vs nn) is crude and flagged (plans/05 risk table).
- **Coherent T2**: no coherent-diffractive final-state model; the hadronizer is refused
  on the coherent channel by default (`hadronize_coherent` opt-in reproduces v0).
- **Physics inputs still external** (unchanged from plans/04): VMC α+d / α+t overlaps,
  spin-3/2 rank-2 basis, b₁ for A > 2, coherent amplitude for polarized A > 2,
  tensor-sector radiative corrections, polarized nuclear PDFs. All are `Backend`s.
- **Conventions awaiting the author**: `TENSOR_LL_SIGN = +1` (plans/08 D1),
  ⁶Li effective polarization 1/3 vs 0.81 (plans/04 #6), `EmcBaseline` default = Epps21
  (mirrors the Python default of 2026-08-29).
- **HepMC3 → abconv → npsim smoke test** not yet run (needs the eic-shell container).
- Packaging: no `pyproject.toml` yet (PYTHONPATH route via `env.sh`).
