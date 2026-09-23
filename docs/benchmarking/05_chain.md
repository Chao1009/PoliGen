<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 05 — CHAIN: MC-internal closures, the ePIC software chain, and other EIC
# generators' validation practice

Companion to `00_in_tree_checks.md`. That file inventories what this tree
already checks (its Table F is the chain row set, F-1 … F-8). **This file does
not re-list those rows**; it says, for each of them and for the closures the
tree does *not* have, what EXTERNAL resource exists, whether it is obtainable,
what it would actually validate, and — crucially — **what it would not**.

Everything below marked *verified* was checked on this machine on **2026-09-06**
by fetching the page/repository, listing the remote store, or running the tool.
Nothing was downloaded into either repository; the only files written are under
`docs/benchmarking/`.

---

## 0. The frame

Three separate things get called "validation" in a Monte Carlo project, and the
chain domain contains all three. Keeping them apart is the whole point of this
file:

| level | what it can prove | what it can never prove |
|---|---|---|
| **MC-internal closure** | the code computes what its own equations say; estimators are unbiased; the record is conserved, reproducible and well formed | that the equations are right |
| **software-interface check** | a downstream tool accepts the file and produces the object it is supposed to | anything physics-side |
| **generator/data comparison** | the modelled final state resembles a measured one | anything about the tensor-polarized sector, for which **no data and no second generator exist** |

The honest headline for this domain: **the chain level is the ONE place where
LiPolGen can be checked as thoroughly as any mature MC project**, because none
of it depends on tensor polarization. Everything in §2–§4 below is achievable.
It is also the level that proves the least about the physics the generator
exists to produce. A reader who takes "the chain gates all pass" as evidence
that the ⁶Li tensor observables are right has misread this file.

---

## 1. MC-internal closures

### 1.1 What exists in `tests/` today

Compact map only — the details are in `00_in_tree_checks.md` F-3 … F-7 and in
the TEST_CASE names themselves.

| closure | present? | where |
|---|---|---|
| per-event 4-momentum + charge conservation, every channel, T0/T1/T2 | **yes** | `tests/test_pipeline.cpp`, `tests/test_breakup.cpp`, `tests/test_t2.cpp` (9 cases), `tests/test_rc_pipeline.cpp` |
| seed / thread / order determinism, incl. byte-identical HepMC3 | **yes** | `tests/test_rng.cpp` (5), `tests/test_pipeline.cpp` "threaded generation is bit-identical", `tests/test_t2.cpp` determinism case, `tests/test_cluster_config.cpp` T13/T14 |
| estimator closure vs the analytic error formulas | **yes** | `tests/test_sampler.cpp`, 6 `closure:` cases — pulls unbiased, spread within 15 % of 1/(P_e P_z √N) and √(2/N)/P_zz, rel-luminosity biases reproduced *and removed* |
| weighted vs unweighted agreement | **partial** | `tests/test_sampler.cpp` "Mode-W weights reproduce the polarized rate ratio"; `tests/test_rc.cpp` T17 (weighted RC slots); `tests/test_fsi.cpp` "weight_normalised … averages to 1". **There is no single test that generates the same physics point unweighted and weighted and compares the two σ estimates with their MC errors** |
| spin-density-matrix trace / positivity | **yes** | `src/core/spin.cpp` `trace()`; `tests/test_spin.cpp` "multipole operators are orthonormal under the trace product", "unphysical populations are rejected"; `tests/test_bookkeeping.cpp` "pzz_true … through the operator trace"; cross-section positivity in `tests/test_xsec.cpp` "density_min is the exact minimum over phi" / "production scenarios keep a healthy positivity margin", coherent-sector positivity in `test_pipeline.cpp` / `coherent.cpp` |
| HFS identities Σ(E − p_z), Σp_T | **yes** | `python/tests/test_hfs.py` (through the `.npz` export), `docs/PYTHIA_BRIDGE.md` §8, `docs/T2_CHAIN.md` §1a — exact form holds at 1.5e-13 GeV, the closed-form collinear `hfs_sigma_empz_truth` is approximate at the 1.4 % level at T1 and that is documented as the formula's assumption failing |
| Close–Kumano ∫b₁ dx = 0 | **reported, not enforced** | `close_kumano_integral` in `sf.hpp` / `b1_nuclear.hpp`; `00_in_tree_checks.md` E-4. The digitized tables stop before x → 2, so the integral cannot close on the data that exists |
| ∫Δ_T f dy = 0 (the convolution's own rank-2 sum rule) | **yes, enforced** | `tests/test_b1_nuclear.cpp` "gate layer 0: int delta_T f dy = 0 (G0e)" — this is the one sum rule the tree *asserts* |
| **Burkhardt–Cottingham ∫g₂ dx = 0** | **absent** | see §1.2 |
| **Bjorken sum rule** | **absent** | see §1.2 |
| unitarity of the FSI operator (∫dΓ S[FSI+FSI²] = 0, Strikman–Weiss) | **absent** | `fsi.hpp` checks the σ→0 limit and ⟨w⟩ = 1, not the unitarity integral |

### 1.2 The two sum rules that are missing, and what they are worth

**Burkhardt–Cottingham, ∫₀¹ g₂(x) dx = 0.** LiPolGen builds g₂ by the
Wandzura–Wilczek relation (`numerics.hpp`, `g2_ww`, pinned in
`tests/test_sf.cpp` against its analytic power-law form). **∫g₂^WW dx = 0 is an
algebraic identity of the WW form**, so an in-tree BC test would be a cheap and
worthwhile *self-check of the quadrature* — it would catch a wrong integration
range or a lost Jacobian — and nothing more. It would **not** test the physics:
the twist-3 piece that BC is actually sensitive to is set to zero by
construction here. Comparing to measured g₂ (SLAC E155x, JLab Hall A / RSS /
SANE) would test the WW *approximation*, which is a structure-function-model
question, not a chain question, and belongs in the theory survey rather than
here.

**Bjorken.** ∫(g₁^p − g₁^n)dx = g_A/6g_V × (1 + α_s corrections). In LiPolGen
g₁^p and g₁^n come from an external polarized PDF set (NNPDFpol1.1 through
`LhapdfG1`, or `ToyG1`). The sum rule is then a property of the **PDF set**, not
of this generator, and NNPDFpol1.1 already imposes/reports it. Adding the test
would check that `LhapdfG1`'s flavour handling and the finite-flavour convention
(`tests/test_lhapdf.cpp` "NNPDFpol11 grid carries c/b, PartonG1 drops them")
have not silently broken the isovector combination — again a plumbing check, of
real but limited value. **Neither sum rule touches the tensor sector at all.**

### 1.3 What a standard MC validation suite would add — and the templates for it

| addition | template, verified | value here |
|---|---|---|
| a **generation / analysis / plotting / timing** validation package that reruns a fixed set of physics points across versions and diffs the output | The HepMC3 paper's own validation tooling (Buckley *et al.*, *Comput. Phys. Commun.* **260** (2021) 107310, [arXiv:1912.08005](https://arxiv.org/abs/1912.08005)) — a generation tool writing a common ROOT n-tuple for Pythia 8 / Herwig++ / Pythia 6+Tauola, an analysis tool, a plotting tool and a **timing tool** to compare per-event cost | LiPolGen already prints events/s per channel in `generate_full` and has throughput TEST_CASEs, but there is **no cross-version regression of physics distributions** — only the rtol-1e-12 `validation/reference/*.json` pins, which are self-referential (see `00`'s Table A) |
| a **shipped example/regression corpus** run in CI | PYTHIA 8.317 in `deps/src/pythia8317/examples`: **141 `mainNNN.cc` programs + `runmains` + `Makefile`** — *verified present on this box*. This is PYTHIA's de-facto test suite and the convention LiPolGen's `examples/` directory imitates | LiPolGen has 5 example binaries (`generate_full`, `generate_inclusive`, `generate_tagged`, `generate_coherent`, `hadronize_example`, `write_hepmc_example`) and no `runmains` equivalent |
| an **I/O round-trip + format-corner test corpus** | HepMC3 3.3.0's own `test/` directory: **86 files** — `testIO1…IO26`, `testAttributes`, `testBoost`, `testDelete`, `testMultipleCopies`, `testThreads`, plus 20+ reference input files in Asciiv3, HepMC2, HEPEVT, LHEF, ROOT and protobuf. *Verified present in `deps/src/HepMC3-3.3.0/test/`* | `tests/test_hepmc.cpp` covers 6 cases and `test_t2.cpp` two round trips. The HepMC3 corpus is the template for the corners LiPolGen does not test (multiple copies, thread safety of the writer, attribute survival across formats) |
| a **spec conformance validator** for the event record | **NuHepMC** — S. Gardiner, J. Isaacson, L. Pickering, *SciPost Phys. Codebases* **57** (2025), [arXiv:2310.13211](https://arxiv.org/abs/2310.13211). Conventions are named `<Component>.<Category>.<Index>` over G (run metadata) / E (event) / V (vertex) / P (particle), split into **R**equirements, **C**onventions and **S**uggestions, with `NuHepMC/ReferenceImplementation` and **`NuHepMC/Validator`**, a tool that checks a HepMC3 file against the standard. *Verified from the spec repository* | **This is the exact template for `docs/HEPMC3_CONVENTION.md`.** LiPolGen's ion-spin attribute schema is a proposed convention with no machine-checkable conformance test; NuHepMC shows what "propose it to the ePIC MC group" should look like as a deliverable |
| **MCnet-style independent replay of published measurements** | Rivet — see §3, which is the largest single finding of this survey | none today |

**MCPLOTS (mcplots.cern.ch) — could NOT be verified**: the site's TLS
certificate has expired as of 2026-09-06 and the page would not load. Do not
cite it as an available resource on the strength of this file.

---

## 2. The ePIC chain

### 2.1 What is actually installed on this box (verified 2026-09-06)

| item | path | verified contents |
|---|---|---|
| `singularity` | `/usr/bin/singularity` | present |
| **eic_xl-nightly.sif** | `~/Projects/eic-2026/local/lib/` (4.5 GB, 2026-08-28) | **Rivet 4.1.2**, `rivet-mkhtml`, **abconv 0.1.3**, `npsim`, `eicrecon`, `HepMC3-config` → **3.3.0**, YODA 2.1.2, **DD4hep 1.37**, **Geant4 11.4.1**, FastJet 3.5.0, LHAPDF 6.5.5, `sanitize_hepmc3`, `runBeamShapeHepMC` |
| **jug_xl-nightly.sif** | `~/Projects/eic/local/lib/` (3.2 GB, 2024-09-13) | the legacy container; `xrdfs`/`xrdcp` work here, and the sibling records that its `pyHepMC3.rootIO.ReaderRootTree` reads the official `.hepmc3.tree.root` files that the new container **segfaults** on |
| official ePIC generator store | `root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/` | **reachable, listed live 2026-09-06** — see §2.4 |

So `00_in_tree_checks.md`'s row "F-1, F-2 — **not runnable from this tree**" is
precise about the repository's dependency tree but should be read with this
addendum: **the containers exist on this machine and the chain is runnable
here today.** What is not runnable is the gate *from `pytest` / `lipolgen_tests`*.

### 2.2 What the 2026-09-02 gate verified — and the leg it did not

The record (`docs/T2_CHAIN.md` §3 item 1, `docs/OPEN_ITEMS_SOLUTIONS.md` §2,
the `DEVELOPMENT_PLAN.md` validation matrix) is: **10/10 events through
`npsim --compactFile epic_craterlake_10x100.xml` directly, and 10/10 via
`abconv -p ip6_hiacc_100x10` first**, with EDM4hep `MCParticles` carrying the
full role chain, and one cosmetic fix found (`generated_mass` for a massless
electron).

What that establishes: the HepMC3 Asciiv3 file is a legal ePIC-chain input; the
10-digit ion codes survive; status-4 beams are accepted; the afterburner will
process it under a manual preset.

**What it does not establish, and should be stated in the same breath:**

1. **The EICrecon leg was never run.** `PolarizedLithiumSim/plans/05` step 5.D
   asks for "100-event abconv → npsim → **EICrecon** smoke per plans/03 step
   2.1.4". The 2026-09-02 pass stopped at `npsim`/EDM4hep. `eicrecon` is in the
   container (`/opt/local/bin/eicrecon`). **This is the single cheapest open
   chain gate in the project** — see CH-14 below.
2. **10 events is a smoke test, not a gate.** It cannot see a 1-in-10³ malformed
   vertex, a rare role that Geant4 refuses, or a channel-dependent failure. Only
   the inclusive ⁶Li channel was exercised; the tagged (α / partner spectators),
   T1 triton-remnant and **coherent (Pomeron, PDG 990, status 3, negative
   `generated_mass` = −√|t|)** records went through no chain at all.
3. **`abconv -p 1` cannot decode a ⁶Li ion** (per-nucleon energy comes out 0) —
   a known, unfixed limitation of the downstream tool, worked around by naming
   the preset by hand. See §2.3 for how deep that goes.
4. Nothing physics-side. It is a file-format and Geant4-primary-generator check.

### 2.3 `abconv` has no light-ion preset at all (verified 2026-09-06)

`abconv --help` in the container lists, verbatim:

```
-p,--preset flag, values [0,1,2] set config and auto determine energy from source file:
  0: IP6, High divergence, auto read energy [default],
  1: IP6, High acceptance, auto read energy
  2: IP6, eAu, auto read energy
  3: IP8, High divergence …  4: IP8, High acceptance …  5: IP8, eAu …
The other options sets energy settings manually, not checking the source file:
  ip6_hidiv_41x5, ip6_hidiv_100x5, ip6_hidiv_100x10, ip6_hidiv_275x10, ip6_hidiv_275x18
  ip6_hiacc_41x5, ip6_hiacc_100x5, ip6_hiacc_100x10, ip6_hiacc_275x10, ip6_hiacc_275x18
  ip6_eau_41x5,   ip6_eau_110x5,   ip6_eau_110x10,   ip6_eau_110x18
  ip8_… (the same three families)
```

The manual presets are **ep** (`hidiv`, `hiacc`) or **eAu**. `ip6_hiacc_100x10`,
the one the 2026-09-02 gate used for a ⁶Li beam at 99.5 GeV/u, is a **proton**
preset: its crossing angle, divergence, crab-cavity kick and vertex smearing are
the ep high-acceptance numbers. Using it for ⁶Li is a *placeholder that lets the
file through*, not a beam-effects model for a light ion — which is exactly the
subject of `PolarizedLithiumSim/plans/10_beam_divergence_light_ions.md`. Any
far-forward acceptance number taken through this chain inherits a proton beam
envelope. **Say so wherever such a number is quoted.**

### 2.4 The official ePIC generator store — listed live, 2026-09-06

`xrdfs dtn-eic.jlab.org ls /volatile/eic/EPIC/EVGEN/` returns
`BACKGROUNDS  CI  DDIS  DIS  Djangoh  EW_BSM  EXCLUSIVE  SIDIS  SINGLE  Test`.
The four sub-trees that matter here, all verified by listing:

| tree | contents (verified) | what it is for |
|---|---|---|
| `DIS/pythia8.316-1.0/NC/noRad/ep/` | `5x41 5x100 5x130 9x130 9x275 **10x100** 10x250 10x275 18x275`, each split `q2_1to10 … q2_1000toINF` | **the ep 10×100 point is LiPolGen's own per-nucleon energy** (10 × 99.5 GeV/u). The community reference final state at exactly the point the T2 bridge runs |
| `DIS/DJANGOH4.6.21-1.0/NC/` | **`Rad` and `noRad`**, `ep/18x275`, `q2_{1to10,10to100,100to1000,1000to10000}` | a **paired radiative / non-radiative** sample from the community RC code. The Rad/noRad ratio is the **inelastic** QED correction factor. The production ran with the elastic radiative tail OFF (IEL2 = IEL31..33 = 0, `eic/InclusiveDjangohSamples` logs; `../open_items/run_2026-09-23/phase_B3_chain_rc.md` §1.2) |
| `DIS/BeAGLE1.03.02-3.1/` | `eAg eAu eH2 eHe3`; `eH2/{en,ep}/9x130/q2_1to1000/…_ab_run0NN.hepmc3.tree.root` | the **deuteron** spectator control, already afterburned. The sibling ran 20 000 events of this on 2026-08-26 |
| `DDIS/rapgap3.310-{1.0,2.0}/{noRad,ewRad}/ep/{10x100,9x130,9x275}` | `rapgap3.310-1.0_ddis_ep_noRad_10x100_ab.hepmc3.tree.root` etc. | **RAPGAP** diffractive DIS on a proton — the closest external analogue of the coherent γ*–Pomeron T2 tier's hadronic final state, and it uses the same H1 2006 DPDF family as `PDF:PomSet = 6` |
| `EXCLUSIVE/DIFFRACTIVE_RHO/Sartre/sartre_bnonsat_Au_rho_*.hepmc` | coherent ρ on **Au** | coherent nuclear diffraction, wrong nucleus, unpolarized |

There is **no lithium anywhere in the store**, and nothing polarized. That is
the expected result and it is worth writing down: the ePIC campaign has no
light-ion DIS sample against which a ⁶Li/⁷Li generator could be normalized.

### 2.5 The ePIC benchmark repositories — what they are NOT

All three were fetched. None of them is a generator benchmark suite:

- **`github.com/eic/detector_benchmarks`** — "a maintained set of performance
  plots for individual detector subsystems", driven by single-particle
  simulations. It benchmarks the *detector*, from guns, not generators.
- **`github.com/eic/physics_benchmarks`** — evaluates "the ePIC detector as a
  whole with respect to a given physical process", working on reconstructed
  `EDM4eic` output (example: `Exclusive-Diffraction-Tagging/diffractive_vm`).
  It benchmarks *reconstruction*, and it is the right home for a far-forward
  tagging benchmark, but it does not validate a generator's physics.
- **`github.com/eic/epic-analysis`** — an *analysis framework* for SIDIS/DIS
  (kinematic reconstruction methods, binning, Delphes/ePIC/ATHENA/ECCE inputs).
  Explicitly "not a benchmark suite with reference plots".

**Conclusion, stated plainly: there is no ePIC "generator validation benchmark"
to plug into.** An ePIC milestone benchmark for LiPolGen would have to be built,
and the two realistic shapes for it are (a) a `physics_benchmarks`-style
reconstructed far-forward tagging benchmark, and (b) a Rivet analysis — which is
where the community is heading (§3.4).

### 2.6 `sanitize_hepmc3`: the EIC's own format checker, run on a LiPolGen file

The eic_xl container ships `/opt/local/bin/sanitize_hepmc3`, a Python filter
that reads Asciiv3 from stdin and validates the header, the
`START_/END_EVENT_LISTING` framing, the per-event vertex and particle **counts
declared in the `E` line against the `V`/`P` lines actually present**, and the
displaced-vertex `@` annotation — skipping any event whose counts do not match.

**Run 2026-09-06** on 300 inclusive ⁶Li events from
`./build/generate_full --channel inclusive --events 300` (seed 20260713):
**exit 0, no event skipped, one warning:**

```
WARNING: Ignoring lines after END_EVENT_LISTING was reached
```

That warning is **not a LiPolGen defect**: HepMC3's own writer emits the footer
as `"HepMC::Asciiv3-END_EVENT_LISTING\n\n"` —
`deps/src/HepMC3-3.3.0/src/WriterAscii.cc:353` — so the trailing blank line the
sanitizer objects to is present in *every* file HepMC3 writes, and
`sanitize_hepmc3` warns on all of them. Recorded here so nobody chases it.

This is a **free, runnable, external format gate** that the tree does not
currently invoke. See CH-15.

### 2.7 The status-code trap that the chain nearly walks into

Verified from upstream: **AIDASoft/DD4hep issue #918** — *"the hepmc3 reader
does not check for the particles status codes and passes all particles to
Geant4"*, in contrast to the HepMC2 reader which did have a status checker; the
reported symptom was a Geant4 fatal decay exception on a virtual particle with a
wrong mass. Closed by **PR #920** (merged 2023-02-23), which also handles excited
ions and prints a warning for excitation levels Geant4 cannot represent.

Why this matters for **this** generator specifically: LiPolGen deliberately
writes documentation-only rows at status 3 — `VirtualPhoton` (mass written
negative), `HadronicX` (PDG 92, forced to 3), `StruckNucleon` / `StruckCluster`
(off shell), and in the coherent channel the **`Pomeron`, PDG 990, status 3,
`generated_mass` = −√|t|**. A reader that ignores status codes hands Geant4 a
PDG-990 object with a negative mass. The 2026-09-02 pass (DD4hep 1.37, i.e. long
after PR #920) did **not** exercise the coherent channel at all. Running the
coherent record through `npsim` is therefore not a formality — it is the one
place where the tree's own convention and a documented downstream failure mode
meet. See CH-13.

---

## 3. Rivet — the largest finding in this survey

### 3.1 Rivet 4.1.2 is installed here, with 46 H1/ZEUS DIS analyses

*Verified 2026-09-06*: `rivet --version` → **4.1.2** inside
`eic_xl-nightly.sif`; `rivet --list-analyses` returns **23 H1** and **23 ZEUS**
analyses. Rivet reads HepMC3 Asciiv3 directly. The reference data are HEPData
records, shipped with Rivet — no download needed. The DIS hadronic-final-state
ones most relevant to LiPolGen's T2 tier:

| analysis | measurement | HEPData / journal | beams |
|---|---|---|---|
| `H1_1994_I372350` | energy flow (lab and hCM) + energy–energy correlations, Q² ∈ [10, 100] GeV²; historically the "target" analysis for shower models | `hepdata.net/record/ins372350`; *Z. Phys.* **C63** (1994) 377 | 820 × 26.7 GeV |
| `H1_1996_I422230` | **charged-particle multiplicity distributions and their moments** in the current fragmentation region of the hCM, vs W and Q² | `hepdata.net/record/ins422230` | 27.5 × 820 GeV |
| `H1_1996_I424463` | transverse-momentum spectra of charged particles in DIS | INSPIRE 424463; hep-ex/9610006, *Nucl. Phys.* **B485** (1997) 3 | HERA |
| `H1_2006_I699835` | event shapes — thrust, jet broadening, jet mass, C-parameter, 14 < Q < 200 GeV | `hepdata.net/record/ins699835`; hep-ex/0512014 | 27.6 × 820/920 |
| `H1_2013_I1217865` | charged-particle densities in η*, p_T* in bins of x and Q², 5 < Q² < 100 GeV² | `hepdata.net/record/ins1217865`; arXiv:1302.1321 | **"Beam energies: ANY"** |

`H1_2013_I1217865` is singled out because its Rivet metadata declares **ANY**
beam energies — at fixed (x, Q²) the hadronic W² = Q²(1−x)/x is fixed
independently of the beam energies, so this analysis is in principle applicable
at EIC energies inside the overlapping (x, Q²) bins. It is the natural first
target.

`H1_1996_I422230` is the direct external replacement for the tree's
**in-tree-only** multiplicity check: `tests/test_pythia.cpp` "the charged
multiplicity agrees with a stock WeakBosonExchange run in the same (x, Q²)
window" compares the bridge to a *locally built stock PYTHIA*, i.e. two
configurations of the same library. `H1_1996_I422230` compares it to **data**.

### 3.2 …and Rivet cannot read a ⁶Li beam. Measured, not inferred.

This is the load-bearing result of this section and it was obtained by running
the tools, so it is stated with the numbers.

**Setup.** 300 events, `./build/generate_full --channel inclusive --events 300
--out …/li6_incl.hepmc` (seed 20260713, e 10 GeV × ⁶Li 99.5 GeV/u, T2 on).
Median x = 6.4e-3, median Q² = 1.78 GeV². A three-line Rivet plugin (`XProbe`)
declaring `DISKinematics(DISLepton())` was compiled with `rivet-build` in the
container and run on the file.

**Result 1 — Rivet reads the file.** `rivet` parses it, counts 300 events, and
picks up the cross section: `Cross-section = 5.920306e+05 pb`. The Asciiv3
writer and the `GenCrossSection` attribute are fine.

**Result 2 — `DISKinematics` returns nothing.** With the real beam code
`1000030060`, every event gives

```
x = -1   Q2 = -1   y = -1   s = -1   W2 = -1
beamHadron().pid() = 10000   beamHadron().E() = 0   beamLepton().E() = 0
```

i.e. the projection's "unset" sentinels. **`--ignore-beams` changes nothing** —
it was tried both ways. `scatteredLepton().E()` *is* correct (9.717 GeV on event
1), so `DISLepton` finds the scattered electron unaided; it is the **beam pair**
that Rivet fails to resolve when one beam is a nucleus. Note that Rivet does
understand nuclear PDG codes — `/opt/local/include/Rivet/Tools/ParticleIdUtils.hh`
has `isNucleus(pid)` "implementing the 2006 Monte Carlo nuclear code scheme" —
so this is a beam-hadron identification gap in the DIS projections, not a PID
parsing gap.

**Consequence:** running `rivet -a H1_… li6.hepmc` **silently produces empty
histograms**. It does not error. The only outward sign in the run above was four
`WARN … invalid scale factor = nan` lines at finalize, and a `.yoda` in which
exactly one object (`/RAW/_EVTCOUNT`) had entries out of 114. **Anyone who wires
Rivet into this project without checking `numEntries()` will publish a plot of
nothing.**

**Result 3 — with the beam relabelled a proton, the projection runs, and then
the *nuclear* x is what you get.** Substituting PDG `1000030060 → 2212` in the
file (four-vector untouched) makes `DISKinematics` produce numbers:
`beamHadron().E() = 597.026`, `s = 23911.9 GeV²`. An independent recomputation
of the exact invariants straight from the same file gives
`s = (k+P)² = 23911.9` and, event by event,

| event | `dis_Q2` (attr) | exact x_A | `dis_x` (attr) | `dis_x` / x_A | exact y_A | `dis_y` (attr) |
|---|---|---|---|---|---|---|
| 1 | 5.50938999588257 | 0.0054839425 | 0.0329041417948201 | **6.000000** | 0.042069442 | 0.0420697446 |
| 2 | 1.0109232228603 | 0.0027072646 | 0.0162438871956878 | **6.000000** | 0.015636646 | 0.0156367012 |
| 3 | 1.22294657663993 | 7.9828277e-05 | 0.000478980150752076 | **6.000000** | 0.64151505 | 0.6415151207 |

So **LiPolGen's record is internally exact**: its `dis_x` is the per-nucleon x
and equals A × x_A to seven digits, and its `dis_y` equals the exact invariant.
A Rivet DIS analysis, which builds x against the whole beam four-vector, would
bin these events at **x/6**. Any H1/ZEUS comparison must therefore be done with
a per-nucleon beam in the file, or with a LiPolGen-specific Rivet projection
that reads `Role::StruckNucleon` — it cannot be done by pointing `rivet` at a
production ⁶Li file.

**One unresolved discrepancy, recorded rather than explained.** In the
proton-relabelled run Rivet's own `Q2` matched `dis_Q2` to all printed digits,
but Rivet's `x` and `y` differed from the exact invariants above by a
*non-constant* factor — 1.031 on event 1, 0.923 / 1.084 on event 2, 0.998 /
1.002 on event 3 (Rivet's x·y product matched the exact x·y in each case). The
cause was **not determined**, and it is a prerequisite for trusting any number
that comes out of a Rivet run on this chain. Reproduce with: the `XProbe` recipe
above, printing `k.x()`, `k.y()`, `k.s()` beside the file's own `dis_*`
attributes.

### 3.3 What Rivet would and would not validate

**Would**: the T2 tier's hadronic final state — charged multiplicity, its
moments, η*/p_T* densities, energy flow, event shapes — against real HERA data,
replacing an in-tree comparison of PYTHIA against PYTHIA with a comparison
against measurement, with the analysis code written and validated by other
people (`00`'s F-row on the multiplicity check is *upgraded*, not replaced —
keep both).

**Would not**: anything nuclear or polarized. The reference data are **ep**.
Between LiPolGen's ⁶Li record and an H1 ep measurement sit Fermi motion, the
absence of any nuclear FSI/formation-time modelling for the DIS debris
(`T2_CHAIN.md` §3 item 2 lists this as open), the EMC effect, and the fact that
the bridge hands PYTHIA one off-shell nucleon rather than a free proton. A
disagreement would not localize to any of them. **Rivet on this chain is a
regression net for the hadronization tier, not a physics validation of the
nucleus.**

### 3.4 The community is moving this way

*Verified*: the **MC4EIC 2025** workshop (9–11 July 2025, JLab Indico event 930)
ran a **Rivet hackathon** to "provide computer codes within the Rivet framework"
for **ep and ed** data, "which will then be used to enhance comparisons and
constraints on event generators", alongside a session on "comparisons and tuning
to relevant ep and ed data". The **ed** half is the interesting one for this
project: the deuteron is the only nucleus where LiPolGen's physics has external
anchors at all. The precedent for how HERA measurements get into Rivet is
[arXiv:2112.12598](https://arxiv.org/abs/2112.12598) — 19 HERA measurements coded
during the DESY 2021 summer-student programme, validated against the legacy
HZTool implementations.

Rivet's own citations, for the record: **v4** — C. Bierlich *et al.*,
*SciPost Phys. Codebases* **36** (2024), [arXiv:2404.15984](https://arxiv.org/abs/2404.15984);
**v3** — *SciPost Phys.* **8**, 026 (2020),
[arXiv:1912.05451](https://arxiv.org/abs/1912.05451), whose abstract names
"heavy-ion and **ep** physics" among the v3 additions.

---

## 4. How the other EIC generators document their validation

The template question: what does a generator in this space actually publish as
its evidence? Answered from the primary sources, with the honest caveat attached
to each.

| generator | what its own documentation claims as validation | verified from | the caveat |
|---|---|---|---|
| **PYTHIA 8.3** | The manual (*"A comprehensive guide to the physics and usage of PYTHIA 8.3"*, Bierlich *et al.*, *SciPost Phys. Codebases* 8-r8.3 (2022), [arXiv:2203.11601](https://arxiv.org/abs/2203.11601)) states the guiding principle as reproducing measured collision properties, and the release ships **141 example programs**. The *DIS-specific* validation is in the follow-up literature and is explicitly against HERA: multi-jet merging in DIS with Vincia compared to **H1** data ([arXiv:2410.20950](https://arxiv.org/abs/2410.20950)), and NLO multiplicative matching for NC DIS used "to describe reduced cross-sections measured at the HERA collider" ([arXiv:2605.00502](https://arxiv.org/abs/2605.00502)), both by Helenius *et al.* | abstracts fetched | These validate PYTHIA's **shower and matching** on a proton. LiPolGen supplies the hard process itself through `LHAup` (`src/pythia/lhaup_dis.hpp`, Strategy 3, unweighted) and uses PYTHIA only for shower + hadronization, so it inherits the HERA pedigree *for that part only* — and only if the bridge's surrogate reproduces a standard configuration, which is what the in-tree stock comparison, and §3's Rivet plan, are for |
| **BeAGLE** | W. Chang *et al.*, *"BeAGLE: Benchmark eA Generator for LEptoproduction in high energy lepton-nucleus collisions"*, **Phys. Rev. D 106 (2022) 012007**, [arXiv:2204.11998](https://arxiv.org/abs/2204.11998) — the abstract's own framing is "model and data comparisons in particle production in both ep and eA collisions" with "tuning of the parameters in BeAGLE based on available experimental data" | abstract fetched | **The abstract does not name E665, HERMES or NMC.** They are widely believed to be among the comparison sets, but this survey did not verify which datasets appear in the figures, and it should not be asserted without opening the paper. Independently: `PolarizedLithiumSim/tools/beagle/README.md` records that **A > 4 uses the C-12 Fermi-momentum parameterization with Woods–Saxon geometry and no α+d/α+t clustering**, so BeAGLE cannot validate the cluster wave function even in principle |
| **eSTARlight** | The generator is `github.com/eic/estarlight`; the physics paper is M. Lomnitz & S. Klein, *"Exclusive vector meson production at an electron-ion collider"*, **Phys. Rev. C 99 (2019) 015203**, [arXiv:1803.06420](https://arxiv.org/abs/1803.06420); the parent is STARlight, Klein *et al.*, **Comput. Phys. Commun. 212 (2017) 258**, [arXiv:1607.03838](https://arxiv.org/abs/1607.03838) | both abstract pages and the repository fetched | **Neither abstract, nor the repository page, advertises a validation against HERA data, and the repository shows no test suite or reference output.** The in-tree D-row uses eSTARlight as a *baseline number*, and this survey found no published validation pedigree to attach to it. If one is needed, it must be dug out of the paper bodies — do not assert it from the abstracts |
| **DJANGOH** | HERACLES + DJANGO6; "the emphasis is put on the inclusion of QED radiative corrections (single photon emission from the lepton or the quark line, self energy correction, complete set of one-loop weak corrections)", plus radiative elastic background; interfaces to LEPTO, ARIADNE, PYTHIA/JETSET, SOPHIA. The COMPASS note (COMPASS 2018-1, N. Pierre) devotes whole sections to **"Self-consistency of DJANGOH"** and **"Comparison between DJANGOH and RADGEN"** — i.e. its documented validation is *code-vs-code against another RC generator*, plus inclusive and semi-inclusive RC factor tables | the COMPASS note PDF was fetched and read | The DJANGOH homepage (`wwwthep.physik.uni-mainz.de/~hspiesb/djangoh/`) **refused connection** on 2026-09-06 and `desy.de/~hspiesb/mcp.html` returned 403; the version/manual details here come from the COMPASS note, not from the author's site. **DJANGOH's RC are the unpolarized/vector sector**; it says nothing about the tensor sector, which is where LiPolGen's RC band lives |
| **RAPGAP** | H. Jung, *"Hard diffractive scattering in high-energy ep collisions and the Monte Carlo generator RAPGAP"*, **Comput. Phys. Commun. 86 (1995) 147**; v3.3 released Dec 2021; the ePIC campaign ships `rapgap3.310` DDIS samples | reference confirmed by search; the ePIC samples listed live | Proton target, unpolarized. It is the right external analogue for the γ*–Pomeron tier's *final state*, not for the ⁶Li coherent amplitude |
| **RADGEN** | I. Akushevich, H. Böttcher, D. Ryckbosch, *"RADGEN 1.0: Monte Carlo Generator for Radiative Events in DIS on Polarized and Unpolarized Targets"*, [hep-ph/9906408](https://arxiv.org/abs/hep-ph/9906408) — built on **POLRAD 2.0**, "analytical and numerical tests are performed and discussed in detail"; the community tool used by HERMES/COMPASS | search + ADS record | **Vector (g₁/A_∥) polarization only.** It is the natural external partner for the tree's POLRAD gate and would give an *event-level* RC cross-check where the tree has only compiled numbers — but the tensor sector stays exactly as uncovered as `00`'s gap list says |

**The template, distilled.** None of these generators publishes a formal
conformance suite. What they publish is: (i) a manual with the physics written
out; (ii) a set of comparisons to measured distributions, usually through Rivet
or a bespoke script; (iii) code-vs-code comparison where no data exists (DJANGOH
vs RADGEN is the closest analogue to LiPolGen's situation); (iv) shipped
examples. LiPolGen already has (i) at unusual depth and (iv); (iii) exists only
against its own Python predecessor, which is not an independent code; (ii) does
not exist and §3 is the cheapest route to it.

---

## 5. Format and convention checks

### 5.1 Ion PDG codes — settled, and the tree is right

The PDG *Monte Carlo Particle Numbering Scheme* (Review of Particle Physics,
current revision at `pdg.lbl.gov/2024/reviews/rpp2024-rev-monte-carlo-numbering.pdf`)
gives nuclear codes as **±10LZZZAAAI**, with A the total baryon number, Z the
charge, L the number of Λ's and I the isomer level (0 = ground state); the
deuteron is `1000010020`. LiPolGen's `1000030060` (⁶Li), `1000030070` (⁷Li),
`1000020040` (α), `1000010020` (d) all conform. Rivet's `isNucleus()` and
DD4hep's ion handling both implement the same scheme. **No open item here.**

### 5.2 Status codes — the tree conforms; the risk is downstream

HepMC3's status-code table (0 = null, **1 = undecayed physical particle**,
2 = decayed, **3 = documentation line**, **4 = incoming beam**, 5–10 reserved,
11–200 generator-dependent, 201+ simulation-dependent) is what
`docs/HEPMC3_CONVENTION.md` maps `lipolgen::Status` onto, including forcing
`Role::HadronicX` to 3 so a status-1 sum is never double-counted. That is
correct and conservative. The exposure is §2.7: a *reader* that ignores status
codes. DD4hep 1.37 is well past the fix, but the coherent record has never been
put through it.

### 5.3 Spin in HepMC3 — no convention exists, and that is upstream's position too

*Verified from the HepMC3 author's own EIC-audience talk* (A. Verbytskyi, HepMC3
library, BNL Indico event 8345): among "obsolete and removed features (vs.
HepMC2)", **"Flow and Polarization classes are removed and the information
should be represented by generic attributes."** HepMC2 carried polarization as
(θ, φ) class members on `GenParticle`; HepMC3 deliberately dropped them, judging
them "rarely filled and meaningful".

So LiPolGen's approach — named `GenEvent` attributes `spin_J`, `spin_M`, `P_z`,
`P_zz`, `spin_axis_theta`, `spin_axis_phi` — is **the mechanism HepMC3 itself
prescribes**. What does not exist is a *named, standard* attribute class, the
way `GenCrossSection`, `GenPdfInfo` and `GenHeavyIon` are standard attributes.
Two consequences worth stating in `HEPMC3_CONVENTION.md`:

1. The proposal to the ePIC MC group (plans/05 step 5.D, plans/04 #17) has a
   concrete, verified precedent to point at: **NuHepMC** (§1.3) standardized
   neutrino-generator semantics entirely inside HepMC3's attribute system, with
   an R/C/S tiering and a `Validator`. That is the shape the ion-spin proposal
   should take, and it is a much stronger ask than "here are six attribute
   names".
2. The `abconv`/`npsim` chain **does not carry event attributes into EDM4hep**
   (already recorded in `T2_CHAIN.md` §2 and `00`'s F-4). So the spin label is
   readable only from the HepMC3 file itself. Any analysis that needs it must
   join back on run/event index. **Nothing downstream of `npsim` can be checked
   against the spin labels**, and no amount of chain validation changes that.

---

## 6. The resource table

Priority is "what buys the most confidence per unit of work", not importance of
the physics.

| # | resource | kind | obtainable | in tree | effort | priority |
|---|---|---|---|---|---|---|
| CH-1 | Rivet 4.1.2 in `eic_xl-nightly.sif` (46 H1/ZEUS DIS analyses) | chain | in-tree *(container on this box)* | no | days | **high** |
| CH-2 | `H1_1996_I422230` charged multiplicities in DIS | data | HEPData (`ins422230`, ships with Rivet) | no | hours once CH-1 works | **high** |
| CH-3 | `H1_2013_I1217865` charged-particle densities, beams ANY | data | HEPData (`ins1217865`) | no | hours once CH-1 works | **high** |
| CH-4 | `H1_1994_I372350` energy flow + EEC | data | HEPData (`ins372350`) | no | hours | medium |
| CH-5 | `H1_2006_I699835` DIS event shapes | data | HEPData (`ins699835`) | no | hours | medium |
| CH-6 | ePIC `pythia8.316-1.0` ep **10×100** EVGEN samples | generator | xrootd, verified live | no | days | **high** |
| CH-7 | ePIC `DJANGOH4.6.21-1.0` **Rad / noRad** ep 18×275 | generator | xrootd, verified live | no | days | **high** |
| CH-8 | ePIC `BeAGLE1.03.02-3.1` eH2 deuteron spectator samples | generator | xrootd, verified live; used once in the sibling | no (sibling only) | days | **high** |
| CH-9 | ePIC `DDIS/rapgap3.310` ep diffractive samples | generator | xrootd, verified live | no | days | medium |
| CH-10 | Sartre coherent ρ on Au | generator | xrootd, verified live | no | days | low |
| CH-11 | local BeAGLE build for A = 6, 7 | generator | **license-blocked** (FLUKA) | no | weeks + license | low |
| CH-12 | `abconv` 0.1.3 — no light-ion preset | chain | in-tree *(container)* | yes (F-2) | hours | medium |
| CH-13 | `npsim` on the **coherent** and **tagged** records | chain | in-tree *(container)* | partly (F-1) | hours | **high** |
| CH-14 | **`eicrecon` leg** — the third leg of plans/05 5.D | chain | in-tree *(container)* | **no** | hours–days | **high** |
| CH-15 | `sanitize_hepmc3` as a CI-able format gate | chain | in-tree *(container)*, run 2026-09-06 | no | hours | medium |
| CH-16 | HepMC3 3.3.0 `test/` corpus (86 files) as a template | chain | in-tree *(`deps/src`)* | no | days | medium |
| CH-17 | PYTHIA 8.317 `examples/` + `runmains` as a template | generator | in-tree *(`deps/src`)*, 141 mains | no | days | low |
| CH-18 | **NuHepMC** spec + `Validator` as the model for the spin convention | chain | public repositories | no | days (proposal), weeks (validator) | **high** |
| CH-19 | PDG MC numbering scheme (10LZZZAAAI) | theory | public | yes (conforms) | none | low |
| CH-20 | HepMC3 status codes + attribute mechanism | chain | public | yes (conforms) | none | low |
| CH-21 | DD4hep #918 / #920 status-checker history | chain | public | no | hours | medium |
| CH-22 | RADGEN 1.0 (POLRAD-based RC event generator) | generator | code on request; paper public | no | weeks | medium |
| CH-23 | `eic-smear` (generator-text ↔ ROOT ↔ HepMC3) | chain | public repository | no | days | low |
| CH-24 | MC4EIC 2025 Rivet hackathon (ep **and ed** analyses) | chain | public agenda; analyses' status unverified | no | — | medium |
| CH-25 | Burkhardt–Cottingham ∫g₂ dx = 0 on `g2_ww` | theory | analytic | **no** | hours | medium |
| CH-26 | `epic-analysis` / `detector_benchmarks` / `physics_benchmarks` | chain | public | no | — | low — **none of them is a generator benchmark** |

---

## 7. What none of this validates

Repeating the frame, because the volume of green in §6 is misleading:

1. **The tensor sector.** Not one resource above is sensitive to b₁, b₂, A_zz,
   the cos 2φ amplitude, P_zz, or the ion spin state. Every HERA analysis, every
   ePIC EVGEN sample, every conformance check is unpolarized (or, at best,
   vector-polarized in the case of RADGEN). If **all 26 rows passed**, the
   statement earned would be: *"LiPolGen writes a well-formed, conserving,
   reproducible event record that the ePIC chain accepts, and its hadronic final
   state is consistent with HERA ep measurements and with other EIC generators
   at the same per-nucleon kinematics."* Nothing about lithium and nothing about
   tensor polarization.
2. **The nucleus.** The store has no lithium and no light ion of any kind; the
   only nuclear samples are ³He, Ag and Au (BeAGLE) and Au (Sartre), and
   BeAGLE's A > 4 treatment has no α+d clustering. The α–d cluster wave
   function, the T1 breakup, the ⁷Li triton remnant and the Glauber FSI weight
   have **no external counterpart at any level of the chain**.
3. **Beam effects for a light ion.** `abconv` has proton and Au presets only
   (§2.3); a ⁶Li run borrows a proton envelope.
4. **The spin label after `npsim`.** Event attributes do not reach EDM4hep. The
   chain is validatable; the *label the chain exists to carry* is not.
5. **The residual Rivet kinematics discrepancy of §3.2** is open. Until it is
   understood, a Rivet number from this chain is not quotable.

## 8. One-line summary

The chain domain is the one place where LiPolGen can be validated to the
standard of a mature MC project, and the tooling to do it — Rivet 4.1.2, the
full ePIC container stack, and a live xrootd store of official PYTHIA 8, DJANGOH
Rad/noRad, BeAGLE and RAPGAP samples — is **already installed on this machine**;
the two highest-value gaps are the never-run `eicrecon` leg and the fact that
**Rivet silently produces empty histograms on a ⁶Li beam**, both measured here
on 2026-09-06; and none of it says anything whatsoever about tensor-polarized
lithium.
