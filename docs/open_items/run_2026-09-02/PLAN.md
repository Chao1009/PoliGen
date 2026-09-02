# Scheduled run 2026-09-02 — remaining open items + physics-channels document

Brief for the deferred (post-usage-reset) run. Authored 2026-09-02 17:05 CDT by the
supervising session; executed phase by phase by the cron kick-off. Everything a
fresh context needs is here or in the files it names. Progress lives in
`STATUS.md` next to this file — read it first, execute the first phase not marked
DONE, mark it, commit, continue.

## Ground rules (apply to every phase)

- Repo: `/home/cpeng/Projects/polli/LiPolGen` (own git repo). Env: `source env.sh`
  (deps in `../deps/install`). Build: `cmake --build build -j` ; C++ tests:
  `build/tests/lipolgen_tests` (doctest; 275 cases green at start); Python:
  `pytest python/tests` (145 green at start). Both suites must stay green; the
  rtol 1e-12 gates against `validation/reference/*.json` must not move.
- Model policy (user instruction, memory `model-assignment-policy`): the supervising
  session (Fable) plans, reviews and verifies; physics-bearing implementation →
  `model: "opus"`; docs, scripts, bindings, tables, reference dumps → `model: "sonnet"`.
  Verification / adversarial-review stages: omit `model` (inherit Fable).
- Ultracode is on: every substantive phase runs as a Workflow (parallel readers,
  Opus implementers with worktree isolation only when they would collide,
  adversarial verify with ≥3 refuters or 3 distinct lenses, a completeness critic).
- Every phase ends with: tests green, docs updated, ONE git commit
  (`Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` + the
  `Claude-Session:` line), and `STATUS.md` updated. Do NOT push; the user pushes.
- Conventions: `docs/CONVENTIONS.md` (`TENSOR_LL_SIGN = −1`, ⁶Li cluster
  polarization 0.81123, `EmcBaseline::Epps21` on CT18ANLO). New physics is opt-in
  behind a `Backend`/option with the current default bit-for-bit unchanged, and
  gets a doc section in `docs/USAGE.md` plus a row in `docs/OPEN_ITEMS_SOLUTIONS.md`.
- No number defined twice: constants go in `constants.hpp` or the owning header;
  Python bindings (`python/bindings.cpp`) expose every new option/constant and a
  pytest gates it.
- Literature already surveyed with equation numbers: `docs/open_items/physics_literature.md`,
  `docs/open_items/code_designs.md`, `docs/open_items/engineering.md`,
  `docs/open_items/vmc_*.md`, prototypes in `docs/open_items/prototypes/`,
  upstream questions in `../PolarizedLithiumSim/plans/04_open_questions.md`
  (#6, #9, #10, #14, #16, #18) and `plans/08_simulation_chain_completion.md`.
- If a usage limit interrupts a phase: the next cron fire resumes from `STATUS.md`
  (partial work is on disk; agents must write files early, not only at the end).

## Phase A — `docs/PHYSICS_CHANNELS.md` (explicit user request)

Deliverable: one document listing every physics channel the generator handles,
the approach used for each, literature references, and the functions/classes
that implement it (as `file:line` links, one per row). Table of contents shape:

1. Beam and target: e⁻ (λ_e) on polarized ⁶Li (J=1: vector + tensor) / ⁷Li
   (J=3/2: vector + rank-2 (+ rank-3 moment)), spin-density matrices, ⟨P₂⟩ = −T/5
   polarimeter — `spin.hpp`, `beams.hpp`.
2. Inclusive DIS kernel (T0): unpolarized F₁/F₂ + R, vector g₁/g₂ (A_∥, A_⊥,
   finite γ), tensor b₁…b₄ (HJM master formula, Cosyn finite-γ tensor sector,
   `tensor_gamma`), effective polarizations, EMC/EPPS21 baseline, target mass —
   `xsec.hpp`, `asymmetries.hpp`, `sf.hpp`, `lhapdf_sf.hpp`.
3. Structure-function backends: toy, LHAPDF (CT18NLO / CT18ANLO, NNPDFpol11_100,
   EPPS21nlo_CT18Anlo_Li6), toy b₁ (Miller-type), where each number enters.
4. Sampling and bookkeeping: inclusive sampler (weighted / unweighted), RNG
   discipline, estimator closure — `sampler.hpp`, `bookkeeping.hpp`, `rng.hpp`.
5. Tagged (incoherent) channels: ⁶Li → α + d* (struck d), ⁶Li → d + α*, ⁷Li → α + t*
   (+ whichever else `tagged.hpp` enumerates): cluster wave functions (Hulthén S+D,
   VMC AV18 α+d / α+t tables), spin ⊗ cluster correlation, A_zz^tag(k), spectator
   kinematics and light-front prescription — `cluster.hpp`, `tagged.hpp`,
   `spectator.hpp`.
6. T1 breakup: struck cluster → struck nucleon + partner spectators (deuteron
   Hulthén S+D with m_S correlation; triton sequential default; Ciofi–Simula
   spectral-function option) — `breakup.hpp`, `triton_sf.hpp`.
7. Spectator FSI weight (Glauber cluster / nucleon variants, σ band 20–40 mb) — `fsi.hpp`.
8. Coherent channel: elastic-like scattering on the whole ⁶Li (T0 scenario), the
   γ*–Pomeron T2 hadronization tier (`Beams:idA = 990`, ζ = β), M_X floor — `coherent.hpp`,
   `pythia_bridge.hpp` (§12 of `docs/PYTHIA_BRIDGE.md`).
9. T2 hadronization: LHAup bridge (`Beams:frameType = 5`), surrogate construction,
   mass repair, silent-cut fixes, spectator append, HFS identities —
   `pythia_bridge.hpp`, `src/pythia/lhaup_dis.hpp`.
10. Output: HepMC3 (status codes, 10-digit ion ids, spin attributes, weights),
    numpy/npz exporters, CLI — `hepmc_writer.hpp`, `python/lipolgen/export.py`, `cli.py`.
11. What is NOT handled (explicit): exclusive VM below M_X = 1.2 GeV, RC (until
    Phase C), b₁ for A > 2 (until Phase D), reconstruction, polarized nuclear PDFs.

For each row: physics approach in 2–5 sentences with formulas where they define
the implementation (not a derivation); references with arXiv/journal ids (pull
from `docs/CONVENTIONS.md`, `docs/open_items/physics_literature.md`, header
comments, `docs/PYTHIA_BRIDGE.md`; known anchors: Hoodbhoy–Jaffe–Manohar NPB 312
(1989) 571; Jaffe–Manohar NPB 321 (1989) 343; Cosyn et al. arXiv:2410.12764;
Cosyn–Weiss arXiv:2603.23699; HERMES hep-ex/0506018; Wiringa et al. PRC 89 (2014)
024305; Piarulli et al. 2023; Ciofi degli Atti–Simula PRC 53 (1996) 1689;
Akushevich et al. POLRAD 2.0 CPC 104 (1997) 201; Cloët–Bentz–Thomas; EPPS21
(Eskola et al. EPJC 82 (2022) 413); CT18 (Hou et al. PRD 103 (2021) 014013);
NNPDFpol1.1 (Nocera et al. NPB 887 (2014) 276); PYTHIA 8.3 (Bierlich et al.
SciPost Phys. Codebases 8 (2022)); HepMC3 (Buckley et al. CPC 260 (2021) 107310);
LHAPDF6 (Buckley et al. EPJC 75 (2015) 132); H1 2006 diffractive PDFs fit B
(EPJC 48 (2006) 715) for `PDF:PomSet`); and the implementing symbols with
`file:line` (verified to exist at that line at commit time).

Workflow shape: parallel readers per subsystem (Sonnet, one per header group) →
structured `{channel, approach, formulas, refs, symbols[{file,line,name}]}` →
per-row adversarial verify (Fable, 3 lenses: symbol-exists-at-line via grep,
physics-claim-matches-code, reference-is-real-and-correctly-attributed) →
Opus writer assembles the document → completeness critic (what channel/option in
`PipelineConfig`, the CLI `--help`, or `docs/USAGE.md` is missing from the doc?)
→ fix → final Fable read. Gate: every `file:line` resolves (`sed -n` check
script), every CLI option appears in some row, README links the doc.

## Phase B — VMC follow-through (item 1 leftovers) + cosmetic writer fix (item 2)

- Re-run the published tag fractions and tagged A_zz^tag(k) curves with
  `--cluster-wave vmc` (script `validation/vmc_tag_fractions.py` exists — extend or
  rerun it; both YR high-acceptance and tagging optics; 5×41, 10×100, 18×275 if
  configured in `beams.hpp`), write the table into
  `docs/OPEN_ITEMS_SOLUTIONS.md` §1 and `docs/USAGE.md`, and retire the Hulthén β
  band language in the docs in favour of "VMC vs Hulthén default" (the β option
  stays in code). Sonnet for scripts/tables, Fable checks numbers against the ones
  already quoted (⁶Li α tag 0.0249 → 0.0348 at 10×99.5; ⁷Li 0.973 → 0.998).
- Item 2 cosmetic: the writer already calls `set_generated_mass` (`src/hepmc/hepmc_writer.cpp:136`);
  verify the electron row gets m_e = 0.51099895 MeV and the test in
  `tests/test_hepmc.cpp` pins it; fix if not.

## Phase C — Tensor-sector radiative-correction band (item 9)

Opt-in `rc.hpp` / `src/core/rc.cpp` (`RcOptions{off, band}`), never a momentum
shift: (i) spin-blind ISR/Born ratio used as a per-event weight or, if the
polligen chain already carries an RC weight, the same interface; (ii) POLRAD 2.0
tensor elastic tail (Eq. A.4 as transcribed in `docs/open_items/physics_literature.md`)
with ⁶Li charge/quadrupole form factors (VMC, Wiringa–Schiavilla PRC 58 (1998)
2698 — tabulate the parametrisation in `data/`), and (iii) the band: 1.5 % on
A_zz at x ≳ 0.05 (E12-13-011 precedent), 10–30 % at x ≲ 0.01. Output: an extra
weight column / HepMC weight name (`rc_tensor_lo`, `rc_tensor_hi`), CLI
`--rc {off,tensor-band}`, `PipelineConfig::rc`, binding + pytest, doctests for
the tail's known limits (→ 0 as Q² → ∞ faster than the inelastic part; symmetry
in the polarization sign; band monotone in x). Opus implements; Fable
adversarial review of formulas against the transcription; doc section in
`docs/USAGE.md` and a §9 update in `OPEN_ITEMS_SOLUTIONS.md` with measured numbers.

## Phase D — b₁ for ⁶Li: three-term α–d convolution (item 10)

`b1_nuclear.hpp` (`B1Li6Convolution`): embedded-deuteron term (b₁ᵈ ⊗ f_{d/⁶Li}(z)
from the VMC α–d relative momentum distribution, light-front z), the α–d D-wave
term (P_D(⁶Li) with F₁ᵈ), and the Clebsch–Gordan depolarization; kernel per
Cosyn–Dong–Kumano–Sargsian (PRD 95 (2017) 074036) as transcribed in
`physics_literature.md`. Validate on A = 2 first: reproduce the deuteron b₁ shape
against the Cosyn figures/tables already tabulated in the docs (state which
numbers were available); then ⁶Li with a mandatory 100 % band (quadrupole
puzzle). Wire as an opt-in `Backend` for the b₁ slot of the inclusive kernel
(`--b1-model {toy,li6-convolution}`), keep `toy_b1` default. Opus implements
in two agents (kernel/convolution; wiring+tests) with Fable design first and
review after. Measured numbers into `OPEN_ITEMS_SOLUTIONS.md` §10.

## Phase E — Spin-3/2 finite-γ theory note (item 8)

`docs/theory/SPIN32_FINITE_GAMMA.md`: the Jaffe–Manohar J = 3/2 basis, the four
leading-twist functions of arXiv:2209.12161 Eqs. 19a–d (rename their "g₂" as
`g1_rank3`), the time-reversal argument that odd-rank multipoles are
unobservable with an unpolarized beam in inclusive DIS (so rank ≤ 2 is exact for
the generator's observables — state the theorem and proof sketch), the
finite-γ decomposition for J = 3/2 by analogy with the Cosyn J = 1 tensor sector
(`tensor_gamma`), and the list of what would change in `spin.hpp` /
`xsec.hpp` if the rank-3 sector were ever switched on (the rank-3 moment already
exists: `spin.hpp:93`). Opus writes, Fable checks every equation against the code
conventions (`c_m = 3m² − 2`, ⟨P₂⟩ = −T/5). Also add `g1_rank3` naming in code
comments where the rank-3 slot is mentioned. No default behaviour change.

## Phase F — Packaging (item 12)

`pyproject.toml` with scikit-build-core (prototype recipe in
`docs/open_items/engineering.md`): `pip install -e . --config-settings=cmake.define.LIPOLGEN_DEPS_PREFIX=…`
must build and `python -c "import lipolgen"` + `pytest python/tests` pass from the
installed module; document in README (build section) with the GPL-3 note and
the `auditwheel` caveat (non-relocatable RPATH). Sonnet.

## Phase G — Coherent ⁶Li amplitude on-ramp (item 11), best effort

(i) An α+d configuration sampler module (`cluster_config.hpp`): α core from the
VMC density, p–n pair from AV18 u/w oriented by the polarization axis, α–d
separation from the VMC momentum distribution; unit tests (normalisation, ⟨r²⟩
against the VMC value, polarization-axis alignment). (ii) eSTARlight ⁶Li
unpolarized coherent J/ψ, φ, ρ rate estimate only if the code builds offline
within ~1 h wall clock (skip and record otherwise). (iii) A doc section in
`docs/OPEN_ITEMS_SOLUTIONS.md` §11 describing the graft into
`subnucleondiffraction` and the collaboration ask. Opus for (i), Sonnet for (iii).

## Phase H — Close-out

Full suites green; `README.md` status paragraph, `docs/OPEN_ITEMS_SOLUTIONS.md`
table, `docs/DEVELOPMENT_PLAN.md` §6 updated with dated status lines; final
adversarial code review of the whole diff since commit `2e404b8` (Fable, 3
lenses); memory file `lipolgen-project.md` updated; commit. Then delete the cron
job named in `STATUS.md` and write the user-facing summary (what was done,
numbers, what stayed open, that nothing was pushed).
