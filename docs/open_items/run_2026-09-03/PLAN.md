# Run 2026-09-03 — the residual open items

Brief for the run that follows the 2026-09-02 close-out (`../run_2026-09-02/`,
commits `8c41899..a94fd6e`, all pushed). Authored by the supervising session
after a 14-agent inventory of the whole tree; the inventory found **44** items
still open across six domains, and two of its feasibility probes **overturned a
documented blocker**. Progress lives in `STATUS.md` next to this file — read it
first, execute the first phase not marked DONE, mark it, commit, continue.

## Baseline at the start of this run

`a94fd6e`, tree clean, both suites re-run and green on this machine:
**363 doctest cases / 17 203 863 assertions / 1 skipped** and
**211 pytest passed** (35.75 s). Every number this run publishes is measured
against that baseline.

## Ground rules (apply to every phase)

- Repo `/home/cpeng/Projects/polli/LiPolGen`. Env `source env.sh`
  (deps in `../deps/install`). Build `cmake --build build -j`; C++ tests
  `build/lipolgen_tests`; Python `pytest python/tests`. Both suites must stay
  green and the rtol 1e-12 gates against `validation/reference/*.json` must not
  move unless a phase says so and re-pins them deliberately.
- Model policy (memory `model-assignment-policy`): the supervising session plans,
  reviews and verifies; physics-bearing implementation → `model: "opus"`; docs,
  scripts, bindings, tables, reference dumps → `model: "sonnet"`; adversarial
  verification stages inherit the session model.
- Ultracode is on: every substantive phase runs as a Workflow — parallel
  readers, implementers under worktree isolation only when they would collide,
  adversarial verification with ≥3 refuters or 3 distinct lenses, and a
  completeness critic.
- New physics is **opt-in** behind an option or `Backend` with the current
  default bit-for-bit unchanged, gets a `docs/USAGE.md` section, a row in
  `docs/OPEN_ITEMS_SOLUTIONS.md`, and a row in `docs/PHYSICS_CHANNELS.md`.
- No number defined twice: constants live in `constants.hpp` or the owning
  header. Python bindings expose every new option and a pytest gates it.
- `validation/check_physics_channels_links.py --fix` runs at the end of every
  phase that moves a line number, and its output is reviewed, not trusted.
  A phase that edits a block `docs/PHYSICS_CHANNELS.md` cites by RANGE, a
  PINNED USE-SITE line (rule S3) or a PYTHIA upstream line (rule D) also
  re-reads the citing row and re-records the fingerprint (`--record-ranges`);
  the gate refuses to bless an edited block or line on its own. Run the gate
  with `env.sh` sourced, or the six external citations are skipped rather than
  checked — it says which, by name.
- Every phase ends with: suites green, docs updated, ONE git commit, and
  `STATUS.md` updated. **Do NOT push; the user pushes.**
- Author decisions are **collected with their evidence and surfaced to the
  user in one batch**, not taken silently and not used to block a phase.
  Nothing outward-facing (no email, no PR, no upload) is sent by this run.

## What the inventory changed about the picture

Three findings reorder the work relative to `OPEN_ITEMS_SOLUTIONS.md`:

1. **The b₁ gate's blocker is not real.** Every document says MSTW2008 LO "is
   not installed" (`phase_D_gate.md:66,313,370`, `b1_nuclear.hpp:17`,
   `USAGE.md:287`, `cli.py:388`). That is true only of LHAPDF's set store. The
   MSTW2008 LO central grid is **on disk** as
   `deps/install/share/Pythia8/pdfdata/mstw2008lo.00.dat`, and PYTHIA ships
   `MSTWpdf` to read it. A probe compiled against it measured MSTW/CT18NLO on
   F₂ᵖ = 1.17 / 1.42 / 1.56 at x = 0.5 / 0.7 / 0.8 — and the G3b peak sits at
   x = 0.736, so the gate's measured ratio 0.7193 should move to ≈ 0.95–1.05.
   **G3b plausibly passes**, which is what lifts the ⁶Li b₁ publication ban.
   (Landmine, measured: `sizeof(Pythia8::MSTWpdf)` is 7 988 336 bytes — a
   stack-local instance segfaults; it must be heap-allocated.)
2. **A published RC number's sign is not stable.** `F_point`'s shape parameters
   are unfitted starting values (`rc.hpp:307-324`, "TO BE REFIT"). The in-tree
   VMC density `data/vmc/density/li6.density` — already committed, already read
   by the cluster-config tests — has a Fourier transform 2.6× the shipped
   `F_point` at q = 2 fm⁻¹. Propagated through `polrad_sigma_el_t`, the
   published `(1/6)σ^el_T/σ^el_U = +1.5595e−04` at x = 0.10 becomes
   **negative**, so the §9 headline "the tensor fraction changes SIGN with x"
   rests on an unfitted guess. This is a correctness problem in a number that
   is already written down, and it outranks the items it was filed behind.
3. **CD-Bonn needs no download of data.** It is a published closed-form 11-term
   parameterisation (Machleidt, PRC 63 (2001) 024001, Tables XVII/XVIII;
   e-print nucl-th/0006014), transcribable into a header — not a table to
   obtain.

## Phase A — lift the ⁶Li b₁ publication ban (open item 10)

The single highest-leverage cluster: it is the only gate in the tree that
forbids publishing a number. `OPEN_ITEMS_SOLUTIONS.md:556-573` lists seven
conditions; A1–A7 are those seven.

- **A1** `include/lipolgen/mstw_sf.hpp` + `src/pythia/mstw_sf.cpp`:
  `class MstwSF : public UnpolSF` over a **heap-allocated**
  `Pythia8::MSTWpdf(2212, iFit=3, pdfdata)`. Charge weights and the
  `f2_from_weights` construction copied verbatim from
  `src/lhapdf/lhapdf_sf.cpp:40-70` so MSTW and CT18 are compared like with
  like. `LIPOLGEN_PYTHIA8_PDFDATA` added beside `LIPOLGEN_PYTHIA8_XMLDOC`
  (`CMakeLists.txt:135`). Bound beside `LhapdfSF` (`bindings.cpp:1155`); the
  `UnpolSF` binding has no trampoline, so a Python-only prototype is impossible
  and the C++ class is mandatory. Doctest + pytest.
- **A2** Two MSTW rows in `checklist_item4()`
  (`validation/b1_li6_table.py:158-170`) and the gate rerun. **G3a will abort
  before G3b is reached**: `tests/test_b1_nuclear.cpp:475` is a `REQUIRE(z.size()
  == 2)` inside a counting window `[0.02, 1.0]`, and the docs already record
  that a realistic PDF drops the low-x zero below the scan floor. The design's
  own tolerance for that clause is Δx = ±0.08 about 0.0656, i.e. it admits a
  zero anywhere in [0, 0.146] — widening the counting window to the tolerance
  the design already grants is defensible; redefining the clause is not. The
  choice is written down in the design document before the test moves.
- **A3** Miller's normalisation (Q5(b)). `refs/1311.4561.pdf` and
  `refs/1702.05337.pdf` are on disk and extract cleanly. The premise of the
  question is already known to be **false** — Miller's Fig. 5 caption carries
  no per-nucleon/per-deuteron label and the phrase never appears in the paper —
  so the answer comes from his Eqs. (1)/(5)/(6)/(24) and Table I instead, and
  the retraction of the Fig. 5 claim ships in the same edit. The arbiter for
  both camps is the HERMES b₁ paper (hep-ex/0506018), which is **not** in
  `refs/`; fetch it.
- **A4** Q5(a) — whether `CdksB1` stops halving the per-nucleon column.
  A3 shows Q5(a) and Q5(b) are **not independent**: Miller Eq. (1) has no 1/A,
  CDKS Eq. (10) has one, and both curves sit on the same HERMES points. Resolve
  the pair together, or state precisely why it cannot be resolved.
  `B1_PER_DEUTERON_TO_PER_NUCLEON` (`constants.hpp:86`) is shared by both camps,
  so neither may be changed through it alone.
- **A5** CD-Bonn u(p), w(p) transcribed into a header from Machleidt's own
  tables, replacing the AV18 D-state rescaling proxy in the gate; the proxy
  stays reachable and the difference is measured, not asserted.
- **A6** `Li6ConvolutionOptions`'s r-default (`r1998` vs `r_sigma_lt`).
- **A7** Re-open G3b. If it passes: lift the ban — **fourteen** files carry it
  (`README.md`, `docs/{USAGE,PHYSICS_CHANNELS,OPEN_ITEMS_SOLUTIONS,CONVENTIONS}.md`,
  the five `run_2026-09-02/*.md`, `b1_nuclear.hpp`, `cli.py`,
  `test_b1_nuclear.cpp`, `b1_li6_table.py`), and `PHYSICS_CHANNELS.md` carries
  it twice. If it does not: restate the residual budget with the new number.

## Phase B — the radiative tail, where it is wrong rather than merely incomplete

- **B1** `F_point`'s shape (finding 2 above). Ship a **band**, not a refit: a
  `C0Shape { Ho, VmcFt }` edge on `HoSpin1FFOptions` whose VMC branch is the j₀
  transform of `li6.density`, r-rescaled so ⟨r²⟩_point = 6.0788 fm² is exact and
  T11's closed forms still pass on both edges. Republish every σ^el_T row and
  the "changes SIGN with x" headline as sign-indeterminate at x = 0.10 if that
  is what the band says. Two policy points must be argued in the header, not
  asserted: `design_G_cluster_config.md:588` declares `li6.density` a validation
  target rather than an input (T6 becomes partly circular — needs a written
  override), and the same file is being called too wrong for the quadrupole and
  right enough for the monopole. Split off the one-line contradiction between
  T11's q₀ ∈ [2.9, 3.3] fm⁻¹ and `PHYSICS_CHANNELS.md:370`'s |t| ≈ 0.31 GeV²
  (⇒ q = 2.82 fm⁻¹, outside the gate).
- **B2** Promote the s-/p-peaks. Six functions (`ll_radiator`,
  `rosenbluth_spin1`, `rosenbluth_nucleon`, `dsigma_el_dq2`, `ll_peaks_spin1`,
  `ll_peaks_qe`, `tests/test_rc.cpp:659-733`) move into `src/core/rc.cpp`; T8(c)
  then tests the shipped path instead of a test-local construction, which is a
  green, independently commit-able first step before any table moves. Then widen
  `TailTriple`/`tail_sigma_at`/`build_tail_tables` and add the tail-model knob
  so the t-peak-only number stays reproducible. **There is no tensor s/p peak**:
  `ll_peaks_spin1` returns the unpolarised Rosenbluth only, so the s+p term goes
  into a column deliberately kept out of the tensor numerator, and that is
  stated rather than allowed to silently dilute the tensor fraction.
- **B3** The polarised quasi-elastic tail — unpriced on a tail that is 99.9 %
  quasi-elastic at x = 0.30. A bounded `qe_tensor_scale` on
  `RcModel::tail_ratio_at` (`rc.cpp:1014-1016`), where scale = 1 is exactly the
  bound "QE tensor fraction ≤ elastic tensor fraction". Cap the scope at a
  bound: the full step is a POLRAD Eq. (44) tensor analogue for an A = 6 spin-1
  nucleus, which does not exist in the literature.
- **B4** The tagged channels' quasi-elastic half — `rc_tail ≡ 1` there is only
  half a kinematic fact.
- **B5** The `Im^el_6` transcription error still printed in the design document.
- **B6** The A = 2 → A = 6 transfer band on δ (`RcModel::delta`, `rc.cpp:756`),
  added in quadrature and named as Q8. The Gakh–Shekhovtsova *shape* is
  deliberately **not** adopted: its δ is unbounded and changes sign near the
  Born zero (−0.537 at x = 0.169, +1.026 at x = 0.214 for Q² = 4, with
  δ ≡ (Δσ_RC − Δσ_Born)/Δσ_Born as fixed in `phase_B_numbers.md` §B6.3), it is not
  monotone in x (so it falsifies T3 by construction), its Q² support is
  non-rectangular, and its Born is not this generator's Born. The reasons are
  recorded so the option is not re-litigated.

## Phase C — coherent ⁶Li and the VMC inputs (open items 11 and 1)

- **C1** O1, the 7.5× quadrupole budget, closed **quantitatively**. The
  mechanism is real and undocumented: `quadrupole_target_fm2` bisects a dial that
  mutates `ad_f2_` in place (`cluster_config.cpp:398-435`), so η moves with it —
  at the measured η = −0.025 the dial is 0.5200 and Q = −0.18557 fm². The budget
  is **two** factors, not three: 3.317 × 2.269 = 7.524 exactly. The
  1/√S_αd = 1.17 is *already inside* the −0.615 baseline
  (`cluster_config.cpp:369-374`) and must be written as a ceiling on what a
  coherent missing-component model could add, never as a third factor. The
  George–Knutson error ±0.012 on η spans Q from −0.05 to past −0.40, so the
  3.3 leg is really 1.5×–12×: pin the numbers as a regression anchor, do not
  dress them as a physics uncertainty.
- **C2** O5 — whether ⁶Li's a₂ survives EIC statistics and separates from the
  linearly-polarised-photon cos 2φ background. The document says this is
  answerable today; it is the question that decides whether the collaboration is
  worth proposing, so it comes before the ask is finalised.
  *(Answered 2026-09-04, `phase_C_numbers.md` §C2 / `OPEN_ITEMS_SOLUTIONS.md`
  §11.3: **MARGINAL — it survives, barely, and as a BAND** — S = 2.63 σ at
  the band's low edge and 2.84 … 3.29 σ at its top, 3 σ at 8.3 … 13.0 fb⁻¹/u,
  inside the {1, 10, 100} band at both ends, over the whole Q² range with both
  lepton channels. Whether the top crosses 3 σ is **not established**
  (§11.3c): the ⁷Li → ⁶Li efficiency substitution that sets it is undetermined
  in direction, ×0.99–1.33. a₂ ∝ |t| on an e^{−B|t|} sample is only a 0.25 %
  modulation, but that is 2.6 σ on 4.2e5 events. The separation is **free**
  and survives Q² → 0 — the P_zz flip is a 1.50× gain, not a cost. **Read the
  number with its limitation: `a2_from_quadrupole` is a CLOSED FORM, not a
  Good–Walker dipole-model amplitude** (§C2.0) — no amplitude, no saturation,
  matter quadrupole as a proxy for the gluon one; a dipole run could move it
  by ×1.15 either way across the 3 σ line. And what is still not established
  is not statistics: no detection efficiency exists below Q² = 0.1, where
  85 % of the rate is; **no decay-lepton reconstruction efficiency exists in
  this tree at all**; and the far-forward working point is unchosen — the last
  is the **single correction that on its own restores a NO**. So C6 finalises
  the ask **on photoproduction**, with all three gaps stated.
  The first pass of this bullet read "it does not survive, in J/ψ — 0.75 σ,
  160 fb⁻¹/u … C6 should not finalise the ask as drafted"; that was computed
  from one lepton channel in one Q² window and is withdrawn,
  `OPEN_ITEMS_SOLUTIONS.md` §11.3a.
  **The second pass then shipped its MARGINAL verdict on an efficiency chain
  with two defects of comparable size and opposite sign, in neither case in
  its own "complete" assumption list** — ε_det taken at ⁷Li's TOP energy
  18 × 117.9 and applied to a ⁶Li sample at 10 × 99.5 (×1.12–1.16 UP, by that
  paper's own ³He energy scan), and no decay-lepton acceptance or
  reconstruction efficiency in the chain at all (DOWN, geometry 0.99 and
  bounded, reconstruction efficiency unbounded here). **They cancel to 0.7 %**,
  which is why the point barely moved and why the verdict is now a band;
  `OPEN_ITEMS_SOLUTIONS.md` §11.3b, `phase_C_numbers.md` §C2.3c. That pass
  also repaired a claim — "every single conservatism on its own leaves 3 σ
  inside the band" — that was **false at two of the three sites printing it**,
  the de-squeezed optics row being the counterexample.)*
  *(And the C1 bullet above carries a premise §C1 measured false: the ±0.012
  band's near edge is **+0.030 fm², the wrong sign**, not −0.05, so the leg is
  "3.32× (1 σ: 1.54× … sign change)", not 1.5×–12×. Left as written — this is a
  record of the brief, not of a result.)*
- **C3** O2 — the `li6.adr.fit` R₂ node near 1.07–1.1 fm, real or fit artefact;
  it decides the honest default between `FitRescaled` and `OverlapRaw`.
- **C4** O4 — `coherent.hpp` uses ΔB without defining it; define it and
  re-decide `eps_b0` and its band for ⁶Li.
  *(DONE 2026-09-04, `phase_C_numbers.md` §C4. ΔB is defined once at the
  declaration — |F_m|² = exp(−|t|[B + ΔB_m cos2(Φ−Φ_S)]), ΔB_m = δ_m/2 — and the
  old label was off **by a sign**, not a factor 2. The band was re-decided:
  eps_b0 = −0.08 implies Q_charge(⁶Li) = **−0.9345 fm², 11.42× the measured**,
  the ⁶Li band is −(0.0070…0.0527), and the default is **kept** with the cost
  written down because it is a reference gate.)*
- **C5** The VMC inputs: the α–d asymptotic D/S ratio (η = −0.048 model against
  −0.025 measured) with a dial that actually constrains it; the parsed-then-
  discarded Monte-Carlo errors; the 5 % N_αd and 7 % P_D spreads; the Hulthén-
  pinned deuteron control channel; and the inclusive-vs-tagged P_D split under
  `--cluster-wave vmc`.
  *(DONE 2026-09-04, `phase_C_numbers.md` §C5. **C5.1** no η dial — η is exactly
  linear in the quadrupole dial, so a converter (`quadrupole_for_eta`) was added
  instead. **C5.2** the MC errors are carried and are **negligible** (0.02 % on
  the tagged tensor observable). **C5.3** propagated — and this bullet's own
  framing needs the correction: N_αd spans 5.311 %, the D-wave **norm** 7.181 %,
  and P_D (the ratio) only **2.729 %**. **C5.4** the deuteron control reaches
  AV18, opt-in; the flag had been silently ignored there. **C5.5** the split is
  real (+11.61 % / +6.58 %) but "two channels of one run" is **false** — one run
  builds one channel — and substituting the VMC P_D moves the inclusive number
  6.8 % **above** the ab-initio 0.848, so it is documented, not closed.
  **C5.5b (2026-09-04)** that correction was itself half wrong: "nothing inside
  one run is inconsistent" was **false** — a `--cluster-wave vmc` α-tag run held
  the VMC α–d overlap and the 0.045 scenario deuteron together, worth
  **+2.069 %** on every polarized tagged-α observable, in `dis_target` and in
  the T1 breakup. Both now follow the flag (`DEUTERON_AV18`,
  `BreakupOptions::source`): opt-in path −2.027 %, default bit for bit, `vmc`
  whole-nucleus reading **0.887076**. `validate()` now also **refuses**
  `cluster_wave` on `inclusive`/`coherent`, where it is never read.)*
- **C6** Reconcile the two disagreeing committed copies of the Mäntysaari ask
  into one, informed by C2. **Drafted, not sent** — sending is the user's.
  *(DONE 2026-09-04, `phase_C_numbers.md` §C6.* The two copies disagreed by
  one number (3 % vs 4 % on the α+d model's point-radius match — both correct,
  against different reference radii). Reconciled into ONE canonical draft,
  `docs/open_items/run_2026-09-03/mantysaari_collaboration_draft.md`; both
  earlier copies now point to it instead of restating the ask: in each of them
  the ask blockquote and its two editor's notes were **DELETED**, not annotated
  — a duplicated artefact is reconciled and the superseded copies are replaced,
  or the disagreement stays in the tree beside the file that resolves it. The
  rest of each document is untouched, and the rule (with why it differs from
  how the retracted C1 premise annotated above is handled) is stated in
  `phase_C_numbers.md` §C6.2. Informed by
  C2/O5: the letter's original request (c) asked for exactly the channel O5
  prices, and O5's second pass makes that channel **marginally measurable**
  (2.63 σ at the band's low edge and 2.84–3.29 σ at its top, 3 σ at
  8.3–13.0 fb⁻¹/u), so (c) stays a request to
  perform the calculation and moves to **photoproduction**, with the four
  unestablished factors written inside the ask (the Q² < 0.1 efficiency, the
  ⁷Li → ⁶Li species transfer — undetermined in direction, and what makes the
  band's top a span — the missing decay-lepton factor, and the working
  point) — requests (a) and (b) are
  unaffected, being channel-independent. *(An intermediate version rewrote (c)
  into "is this worth pursuing at all", on the first pass's 0.75 σ; withdrawn,
  `OPEN_ITEMS_SOLUTIONS.md` §11.3a.)* **O3 bookkeeping**
  (§11's O3, blocked on the same collaboration): still unanswerable in the
  form asked (no GFMC configuration table exists in this tree, and the split
  is a property of a Good–Walker amplitude, not a one-body density), but a
  genuinely correlated input already committed (`he4.dd`, an α → d+d overlap)
  bounds the SIZE of 4-body correlation on a comparable position-space moment
  at **+18.2 % (variance) / +8.7 % (rms)** against the uncorrelated-product
  null, `validation/o3_alpha_correlation_bound.py` — a several-to-twenty-
  percent effect, not a factor of several, though not a measurement of the
  split itself. **Still not sent** — the canonical draft says so at its own
  top. No code, no test, no reference JSON changed.)*

## Phase D — what the repository-wide sweep found that the documents missed

- **D1** A ⁷Li inclusive run's entire tensor and cos 2φ sector is **identically
  zero** — the rank-2 structure-function input is never set. Decide and
  implement an input, or make the zero explicit and loud at the run surface.
- **D2** Toy F₂/g₁/R are the shipped defaults everywhere and the tagged channels
  have no injection point for a real backend.
- **D3** The per-nucleon Glauber FSI variant is algebraically identical to the
  cluster one — the two "variants" are one.
- **D4** `PDF:PomSet` is called the coherent T2 tier's largest systematic and has
  never been scanned.
- **D5** The coherent channel's hard |t| ≤ 0.2 GeV² ceiling and c₂ positivity.

## Phase E — release hygiene

- **E1** Widen the documentation link gate beyond one document and wire it into
  a suite instead of leaving it hand-run.
- **E2** Release metadata: no copyright holder, authors or citation exist
  anywhere although GPL-3 is chosen. **Needs the author's name** — collected in
  the author-decision batch.
- **E3** `T2_CHAIN.md`'s stale "unrun" ePIC gate; validation-matrix row 7.
- **E4** CI and the portable wheel, written but not pushed; both are the user's
  call because both are outward-facing on their repository.

## Phase F — close-out

Whole-diff adversarial review, `PHYSICS_CHANNELS.md` refreshed and its link gate
green, `OPEN_ITEMS_SOLUTIONS.md` rewritten against what this run actually
measured, `STATUS.md` final, both suites green, and the author-decision batch
put to the user with its evidence.
