# LiPolGen — development plan (2026-08-30)

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
  (`TENSOR_LL_SIGN = -1`, `LI6_CLUSTER_POLARIZATION`, rank-2 geometry `Q_NN`), mirrored
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
lepton QED radiation (PYTHIA dipole-recoil limitation). *(That is the v0.1 scope
as of 2026-08-30; FSI and the RC weight family have landed since, both **opt-in**
— `--fsi glauber-cluster|glauber-nucleon` and `--rc tensor-band` — and neither
changes a default run. See §6.)*

## 3. Work breakdown and agent assignment

**Status 2026-08-30:** P0–P7 done and merged; the **T1 tier is done and default**
(252 doctest cases / 15.2 M assertions, 136 pytest cases), and the library is
synced with `PolarizedLithiumSim` runs 15 and 16 (§6). Adversarial review done
(`docs/code_review_2026-08-29.md`), all 12 findings fixed with tests. See §6 for what
remains open.

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
- Populations ordered m = +J … −J. `TENSOR_LL_SIGN = −1` — the literature
  convention (Cosyn Eq. 27 / HERMES), decided 2026-08-29 (plans/08 D1).
- Generator window looser than analysis: Q² ≥ 0.7, y ∈ [0.004, 0.985], W² ≥ 8.
- PYTHIA: `WeakBosonExchange:ff2ff(t:gmZ)`, `SpaceShower:dipoleRecoil=on`,
  `PDF:lepton=off`, `TimeShower:QEDshowerByL=off`, `PhaseSpace:pTHatMinDiverge=0.5`,
  `PhaseSpace:mHatMin=0.5` (both silent-cut fixes); with LHAup these become
  shower settings only, the hard phase space is ours.
- Ion spin in HepMC3: attributes `spin_J`, `spin_M`, `spin_axis_theta/phi`,
  `P_e`, `lam_e`, `P_z`, `P_zz`, `struck_cluster_m` on the GenEvent (proposed
  convention, plans/04 #17).

## 6. Open items after day 1

*Updated 2026-09-03: six more deliverables landed — see
`docs/OPEN_ITEMS_SOLUTIONS.md` for the full write-ups. Closed since: the
physics-channels reference document (`docs/PHYSICS_CHANNELS.md`, gated by
`validation/check_physics_channels_links.py`); the tensor-sector
radiative-correction band (`include/lipolgen/rc.hpp`, `src/core/rc.cpp`,
`--rc tensor-band`); the spin-3/2 finite-γ theory note
(`docs/theory/SPIN32_FINITE_GAMMA.md`); packaging (`pyproject.toml` with
scikit-build-core, `pip install -e .`); the eSTARlight unpolarized coherent
⁶Li baseline (`docs/open_items/run_2026-09-02/estarlight_li6.md`) and the
polarized ⁶Li α+d configuration sampler (`include/lipolgen/cluster_config.hpp`,
console script `lipolgen-configs`). Landed, and its A = 2 gate CLOSED on
2026-09-03: the four-term α–d convolution backend for b₁(⁶Li)
(`include/lipolgen/b1_nuclear.hpp`, `--b1-model li6-convolution`) ships opt-in
and band-mandatory — see the "Physics inputs still external" bullet below and
`docs/open_items/run_2026-09-03/phase_A_numbers.md`.*

*Updated 2026-09-02: every item below was investigated — see
`docs/OPEN_ITEMS_SOLUTIONS.md` for the ranked solutions and the
implemented-with-numbers status of rows 5–7. Closed since: VMC cluster wave
functions (implemented), the ePIC chain gate (passed), coherent T2
(implemented — the γ*–Pomeron tier, `docs/PYTHIA_BRIDGE.md` §12), the
Ciofi–Simula triton spectral function (implemented, `--triton-sf`), the
Glauber spectator-FSI weight (implemented, `--fsi`).*

**Run-15/16 sync, 2026-08-30 (the Python moved after the port was made;
`PolarizedLithiumSim` `bb636b5..568ff40`).** Everything below is ported and
gated at rtol 1e-12 against regenerated `validation/reference/*.json`
(252 doctest cases / 15.2 M assertions, 136 pytest cases):

- **D1 CLOSED — `TENSOR_LL_SIGN = −1`**, the literature convention.  The guard
  test now asserts the literature relation itself (`A_zz(1 + εR) = −(2/3)b₁/F₁`)
  rather than the constant, so flipping back fails it.  `A_zz^tag(k)` does NOT
  move: it is a ratio of cluster-wave populations and carries no b₁.
- **plans/04 #6 CLOSED — ⁶Li effective polarization = the cluster picture.**
  `LI6_CLUSTER_POLARIZATION = 0.81123` whole-nucleus, built in `beams.hpp` from
  the SAME `P_D_LI6` / `P_D_DEUTERON` the tagged sector uses (they moved out of
  `tagged.hpp` into `beams.hpp`, mirroring `beams.py`); `LI6_NAIVE_ONE_THIRD`
  keeps the retired Cloet value reachable.  ⁶Li g₁ is multiplied by 0.81123 and
  the deuteron slot became exact (0.9325, was a rounded 0.93).
- **D2 LANDED, off by default — the exact finite-γ tensor sector.**
  `InclusiveKernel::Options::tensor_gamma`, `b3_func`, `b4_func`,
  `theta_q_cos_sin`, `cosyn_tensor_sfs`, `cosyn_unpolarized_sfs`,
  `tensor_harmonics_gamma`; `tests/test_tensor_gamma.cpp` is the full port of
  `evgen/tests/test_tensor_gamma.py`, including both finite-γ rows of Cosyn's
  Table 1 (1e-10), the γ → 0 identity for any b₂ with b₃/b₄ cancelling (1e-12),
  and the leakage coefficient Δ_fake/(γ² b₁) = 0.14–0.16 that retires the
  "γ² b₁ × 1.15" bound.  `xsec.json` gained a `kernel_tensor_gamma` block.
- **The finite-γ kinematics moved to `asymmetries.hpp`** (one implementation for
  both halves of the library), gaining `a_parallel_exact`, the
  `a_parallel(..., g2)` overload and `depolarization_effective` — the divisor
  that inverts A_par with no O(γ²) bias, which `fom.project_observables` now
  uses on the Python side.
- **`EMC_VALENCE_DEPLETION_EPPS21` = 0.031052** (was 0.029789): the free-nucleon
  denominator is CT18ANLO, EPPS21's own proton baseline, not CT18NLO.  Both
  transferred camps move (CBT 0.5322, TMT 0.2113).
- Already in place before the sync and re-verified against the new Python: the
  `target_mass = true` default, the per-configuration tagging optics
  (`OpticsChoice::Tagging` vs `TaggingLegacyLevers`), and R34 = 4.56 m at 5×41.

- **T1 tier: DONE** (`breakup.hpp` / `src/core/breakup.cpp`, default on the three
  tagged channels through `PipelineConfig::tier`). The struck cluster is resolved into
  a struck nucleon + on-shell partner spectator(s), so
  `k + P_ion = k' + p_spec + Σ p_partner + hadrons` holds exactly and the T2 bridge
  needs no caller-side hook (`docs/T2_CHAIN.md` §1a: worst residual 1.4e-13 relative,
  charge 0, no-surrogate tail 0/300 on each tagged channel).
- **T1 breakup realism** (what stays open inside it): the deuteron internal wave
  function is Hulthén S+D at P_D = 0.045 with the full m_S angular correlation, which
  is as good as the cluster model gets without VMC overlaps; the **triton remnant
  default is the crude, flagged sequential model** (plans/05 risk table) — a
  sequential two-body decay at the AME2020 separation energies with an isotropic
  S-wave relative direction and the unbound nn split at its virtual-state pole, no
  Faddeev/AV18 three-body correlation and no tensor structure — and since
  2026-09-01 the opt-in `--triton-sf ciofi-simula` (`PipelineConfig::triton_sf`,
  `triton_sf.hpp`) replaces it with the three-channel Ciofi–Simula spectral
  function: k-dependent n₀/(n₀+n₁) branching (S₀ = 0.6525) and the
  struck-n → (p n) continuum channel the sequential model has no room for. **No FSI moves any fragment's four-vector,
  on any channel** — since 2026-09-01 the tagged channels carry the Glauber
  spectator FSI as an optional per-event **weight** (`PipelineConfig::fsi`, off by
  default; `fsi.hpp`, `docs/USAGE.md` §3): it distorts the tagged spectator's
  spectrum statistically but never a vector, and the breakup fragments' own
  rescattering stays open.
- **Coherent T2 — closed 2026-08-30**: the coherent channel hadronizes by default
  through the γ*–Pomeron tier (`docs/PYTHIA_BRIDGE.md` §12): a third PYTHIA instance
  on a Pomeron beam (`Beams:idA = 990`) runs the same surrogate with `W² → M_X²` and
  `ζ = β` exactly, so the whole record conserves (worst 4.8e-14 relative, charge
  exact, veto 0 at M_X ≥ 1.4 GeV). The old refusal and its `hadronize_coherent`
  opt-in are gone; the knob is `PythiaBridgeOptions::coherent_t2 = {Pomeron, Off}`.
  Still open in the coherent sector: no exclusive-VM channel below M_X = 1.2 GeV,
  and `PDF:PomSet` (default 6) is the tier's largest systematic.
- **Physics inputs still external** (updated 2026-09-03 from plans/04): tensor-sector
  radiative corrections and b₁ for A > 2 have moved from wholly external to
  opt-in `Backend`s (`--rc tensor-band`; `--b1-model li6-convolution`) — the RC
  band is unrestricted, and the b₁(⁶Li) convolution passed its A = 2 magnitude
  gate on 2026-09-03 **for one unpolarized nucleon input** (G3b peak ratio
  0.843 with MSTW2008 LO at CDKS Eq. (21) — `--b1-unpol mstw`; **0.440 on the
  shipped `toy` default, outside the [0.5, 2] window**, so quote numbers made
  with the selector set to `mstw`; 1.000338 with CD-Bonn as well, which is a
  residual below the error of digitizing a published figure and not three-digit
  agreement with CDKS; `docs/open_items/run_2026-09-03/phase_A_numbers.md`),
  which lifted its publication ban — the mandatory ±100 % band stays, and the
  gate is A = 2 and says nothing about the α–d step. The coherent amplitude
  for polarized A > 2 has a first rung — the eSTARlight unpolarized ⁶Li
  baseline plus the polarized α+d configuration sampler — but not yet a full
  polarized coherent amplitude. The spin-3/2 rank-2 basis has its theory note
  (`docs/theory/SPIN32_FINITE_GAMMA.md`); the rank-3 sector it describes stays
  off by default, so no code behaviour changed. Still wholly external:
  polarized nuclear PDFs.
- **Conventions, all decided**: `TENSOR_LL_SIGN = −1` (plans/08 D1, closed),
  ⁶Li effective polarization = the cluster picture's 0.81123 (plans/04 #6,
  closed; the 0.81–0.85 band whose top is the Wiringa VMC 0.848 is the
  remaining uncertainty, not the convention), `EmcBaseline` default = Epps21
  on the CT18ANLO denominator (mirrors the Python default).
- **NOT ported, deliberately** (no C++ counterpart exists): the `reco.py` /
  `recopseudo.py` reconstruction chain and its published seven-bin coherent |t|
  window, `fom.project_observables`' D_eff extraction (the analytic divisor IS
  ported, as `depolarization_effective`; the projection that uses it is
  fastsim-only), `polarized.unpolarized_emc_ratio`'s grid modes and
  `structure.NuclearF2Ratio` (they need LHAPDF at call time — the C++ carries
  the one number they produce, `EMC_VALENCE_DEPLETION_EPPS21`), the
  `ToyG1.g2_nucleus` per-grid g2 cache (a Python performance device; the C++
  computes the quadrature per point), `miller_b1_q2_scale` /
  `toy_b1(q2_evolve=True)` (off by default in the Python and used only by
  `money_b1.py`), and the figure/report scripts of run 14's addendum.
- **HepMC3 → abconv → npsim smoke test: PASSED 2026-09-02** (10/10 events
  through `npsim` directly and via `abconv`; see `docs/OPEN_ITEMS_SOLUTIONS.md`
  §2 and the README's ePIC-chain paragraph).
- **Packaging: DONE 2026-09-02** — `pyproject.toml` / scikit-build-core,
  `pip install -e .` (`docs/OPEN_ITEMS_SOLUTIONS.md` §12–13; the PYTHONPATH
  route via `env.sh` still works and is unchanged).
