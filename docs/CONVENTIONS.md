# Coding and physics conventions

- Namespace `lipolgen`; headers in `include/lipolgen/<module>.hpp`; sources in
  `src/core/`, `src/lhapdf/`, `src/hepmc/`, `src/pythia/`; tests `tests/test_<module>.cpp`
  (doctest, one `TEST_CASE` per identity, reference JSON in `validation/reference/`).
- Units: GeV, GeV², fm only where stated; angles rad; azimuth φ ∈ [0, 2π).
- Frame: head-on. Ion +z, electron −z. Lab (25 mrad crossing) only in an explicit transform.
- Spin: populations ordered m = +J … −J. Quantization axis n̂(θ_S, φ_S) in the head-on frame.
  `TENSOR_LL_SIGN = -1.0` — one constant, defined once (`constants.hpp`), test-guarded.
  **On a TAGGED channel that axis is the ION fill's and reaches the event through
  the spectator rotation alone**: the struck cluster's |S_c m_S⟩ is evaluated with
  its own axis along the BEAM (`InclusiveKinematicsSource::plan_for` builds the
  pure category at θ_S = φ_S = 0 whatever the fill says, which is what `polligen`
  does). That is exact while the struck-cluster DIS kernel is spin-blind, so
  `Pipeline`'s constructor **refuses** a tilted fill together with a term that is
  not — P_e ≠ 0 on any tagged channel, or `--inclusive-b1` on `tagged-6Li-alpha`,
  the one channel that reads it. The two shipped tilted plans (`transverse-tensor`,
  `tensor-flip`) carry P_e = 0 by construction and are untouched, bit for bit.
- Structure-function inputs are `Backend` interfaces: toy implementation always available,
  table/LHAPDF implementations optional. No physics number is hard-coded in two places.
- Randomness: counter-based stream keyed by (seed, run, bunch, event); never a global RNG.
- Every double compared to `polligen` is compared at rtol 1e-12 unless the reference itself is MC.
- No exceptions for control flow in the event loop; errors in setup throw `std::runtime_error`.

## A knob that did not run may not be recorded as if it had

**The rule.** If a selector, scale or switch cannot affect *this* run's output,
the run may not write its value into the npz `meta` or print it in the banner
as though it had. It is either **refused** or **labelled**, never silent and
never bare.

**The mechanism, and there is exactly one.**
`Pipeline::knob_provenance(KnobRunContext)` (`include/lipolgen/pipeline.hpp`)
returns one row per user-settable knob — `{name, flag, value, status ∈ {read,
not-read, refused}, reason, label, at_default}` — and the three surfaces read
that table and nothing else: `meta["knob_provenance"]` plus the individual keys
that go through `KnobProvenance::meta_value()` — `pol_sf`, `rc_scope`,
`coherent_t2`, `pom_set`, `pom_rescale`, `optics`, `pot_config`, and (since
2026-09-05) `b1_model` and `b1_unpol`, whose ⁷Li value is the same
`rank2_none_label()` the row carries and which used to be computed a second
time in `python/bindings.cpp` — the CLI's `KNOB PROVENANCE` banner block
(`cli.knob_provenance_lines`), and `python/tests/test_knob_provenance.py`,
which rebuilds the (spec × knob) matrix — 12 (isotope, channel, plan) specs
× 81 knob variants = 627 cells (re-measured 2026-09-15; 76 / 592 was Phase A's,
left behind by Phase B), not the full (channel × plan) product; see
`USAGE.md` §7c — and asserts the table against
the **output hash**. **A knob added without a row fails that test.** Do not add
a per-knob reach sentence to the banner or a bare per-knob key to the `meta`:
add the row.

**Refuse or label — the criterion, once.** REFUSE when the value would name a
*variation of a piece that did not run* (a scale, a band edge or a shape on a
term the run computes as identically 1): such a value claims a systematic was
PRICED, and no label makes a priced systematic un-priced — `b1_band_scale` on
Miller, `rc_qe_tensor_scale` with no quasi-elastic tail, `rc_sp_tensor_scale` with no s-/p-peaks (2026-09-06), every rc tail sub-knob on a tagged channel. **A REFUSAL CAN ALSO HAVE THE OPPOSITE CAUSE, and the `reason` must say which:** `rc_sp_tensor_scale` is refused on `--rc-tail-model polrad-full` too (2026-09-06) — not because the s/p tensor term did not run but because it **DID**, since POLRAD Eq. (18) carries it and the stand-in would double-count a term that ran. Same status, opposite sentence.
LABEL when it names a *backend, an axis
or a member of a family that a channel-, plan- or set-scan sets uniformly
across runs*: refusing one cell of such a scan costs more than it buys —
`--pol-sf` on coherent and under the unpolarised-beam plans, `--pzz-mode`
under the three tensor plans, `--pom-set` off the coherent channel. (`--pzz`
under `helicity-flip` was the second example until 2026-09-06, when
`--pzz-mode typed` made that flag reachable there and the label moved onto the
mode.) A rule that depends on
the RUN PLAN can only be labelled: `validate()` has no plan (`rc_scope` under a
θ_S = 0 fill).

Why it needs a mechanism at all: the rule was enforced knob by knob and then
broke five times in one run, each time on an axis the previous fix had not
looked at — the channel, the run plan, the T2 tier, the rc sub-knobs, and
knobs recorded nowhere at all
(`docs/open_items/run_2026-09-03/phase_D_numbers.md` §D6).

## Physics defaults that are a CHOICE, and where the single copy lives

- **Tensor-sector sign.** `TENSOR_LL_SIGN = -1.0` since 2026-08-29 (author
  decision, plans/08 D1): the LITERATURE convention, Cosyn et al. EPJ A 61
  (2025) 83 Eq. (27) and HERMES,

      A_zz(θ_S = 0) (1 + ε(y) R) = −(2/3) b₁/F₁   exactly, at every y,

  with P_zz = n₊ + n₋ − 2n₀, i.e. b₁ > 0 means the m = 0 state has the LARGER
  cross section.  It was `+1` — the repository's own transcription of
  Hoodbhoy–Jaffe–Manohar — until that date, and setting the constant back is
  the whole of the change: nothing else in the library knows the sign.  The
  guard test (`tests/test_xsec.cpp`, "Cosyn Eq. 27: A_zz(theta_S=0)(1 + eps R)
  = -(2/3) b1/F1") is written against the LITERATURE relation itself, with no
  reference to the constant, so flipping it back fails it; the neighbouring
  "the program sign IS the literature sign, deliberately" additionally pins
  the constant's VALUE (`CHECK(TENSOR_LL_SIGN == -1.0)`), so the flip is a
  visible act.  This bullet named the second case as the constant-free one
  until 2026-09-16, which is the one case that does reference the constant.  What flips with it: A_zz at
  fixed b₁, the by-product κ of the spin-state ratio, and any b-sector
  subtraction built on κ — including the O(γ²) tensor leakage into cos 2φ.
  What does NOT: |A_zz|, the whole Δ (cos 2φ) sector, and `A_zz^tag(k)`
  (`azz_tensor_curve*`), which is a ratio of cluster-wave populations and
  carries no b₁ at all.
- **⁶Li effective polarization.** The CLUSTER PICTURE since 2026-08-29
  (plans/04 #6, closed).  `LI6_CLUSTER_POLARIZATION = (1 − 1.5 P_D_LI6)
  (1 − 1.5 P_D_DEUTERON) = 0.81123` whole-nucleus, and `LI6().eff_pol_p =
  eff_pol_n = LI6_CLUSTER_POLARIZATION/3`, so `Z·P_p = N·P_n = 0.81123`.  The
  two D-state probabilities live in **`beams.hpp`** — one source of truth, as
  in `polli_fastsim.beams` — and `tagged.hpp` uses those names rather than
  keeping copies, so **on the Hulthén default** the INCLUSIVE effective
  polarization and the TAGGED S/D interference of `li6_alpha_channel` are the
  same wave function seen in two experiments — measured 2026-09-04 and
  re-measured 2026-09-06 on the fixed S–D interference sign, they agree
  to **7.84e−6** (0.869950 against the tagged model's **0.8699431789**, a closed
  form against a grid quadrature; 6.82e−6 absolute).  *(The pre-fix pair was
  0.869939 and 1.22e−5: the fix moves the 96-cell quadrature residual by
  +4.3e−6 and nothing else here — the dilution is angle-integrated.)*  **Under `--cluster-wave vmc` they do not**: the
  tagged α–d wave becomes the ANL VMC overlap at `VMC_P_D_LI6` = 0.019355 and
  the inclusive constant does not follow, so the two differ by **+11.61 %** in
  the vector sector and **+6.58 %** in the rank-2 one (`LI6_B1_RANK2_TRANSFER`
  = 0.921947 against the VMC channel's `tensor_dilution` = 0.982576).  It is
  **not closed by substitution and that is the finding**: feeding
  `VMC_P_D_LI6` into the product gives `LI6_CLUSTER_POLARIZATION_VMC` =
  0.905427, **6.8 % above** the one ab-initio number for this observable
  (`LI6_POLARIZATION_VMC_SIX_BODY` = 0.848, Wiringa PRC 89 (2014) 024305
  Table I), where the shipped 0.811228 is 4.3 % **below** it; adding the AV18
  deuteron's own P_D = 0.057600 gives 0.887076, still 4.6 % above.  The
  product 1 − 1.5 P_D × 1 − 1.5 P_D therefore spans **0.811 … 0.905** on the
  wave functions in this tree, the ab-initio answer sits inside that, and the
  FORMULA carries the error, not the choice of P_D.  Author decision
  2026-09-04: the default stays 0.811228, the drift is documented at 11.61 %
  rather than "fixed", and no inclusive ⁶Li polarization is quoted without the
  band 0.81–0.91 (`docs/open_items/run_2026-09-03/phase_C_numbers.md` §C5.5,
  pinned in `tests/test_tagged.cpp` **T26**).  The drift is between a tagged
  row and an inclusive row of one programme, never inside one run — **and that
  sentence was false when it was first written (§C5.5b, 2026-09-04).**  A
  `--cluster-wave vmc` tagged-α run did hold two deuteron wave functions:
  the ANL VMC AV18+UX α–d overlap for the RELATIVE motion, and the scenario
  Hulthén deuteron (P_D = 0.045) for the embedded deuteron in both places a
  run reads it — `TaggedChannel::dis_target`, i.e. the struck cluster's g₁,
  and `BreakupOptions`, i.e. the T1 struck-nucleon MOMENTUM **and** spin draw
  (`ClusterBreakup` samples k from that same deuteron).  Every polarized
  tagged-α observable was **+2.069 %** high against the AV18 deuteron
  (P_D = 0.057600) belonging to that overlap.  Since 2026-09-04 both follow
  the flag (`DEUTERON_AV18()`, `BreakupOptions::source`), the opt-in path moved
  by −2.027 %, the Hulthén default is bit for bit, and a `vmc` run's
  whole-nucleus reading is one Hamiltonian's: 0.887076.  `validate()` now also
  **refuses** `cluster_wave` on `Inclusive` and `CoherentLi6`, where it is
  never read — accepting it unread was how the 11.61 % above reached a user
  who thought they had asked for a VMC ⁶Li (**T27**, and the pipeline test).
  The deuteron slot carries `DEUTERON_VECTOR_POLARIZATION =
  1 − 1.5 P_D_DEUTERON = 0.9325` verbatim, so per-nucleon
  g₁(⁶Li)/g₁(d) = (1 − 1.5 P_D_LI6)/3 = 0.290 exactly (the deuteron's own D
  state cancels between the two isoscalar ions).  The retired Cloet
  convention stays reachable as `LI6_NAIVE_ONE_THIRD` — a whole-nucleus 1.0,
  1.233× this one and above the **0.81–0.91** band of this tree's wave
  functions, INSIDE which the ab-initio Wiringa VMC 0.848 sits (0.811 … 0.905;
  this bullet said "the 0.81–0.85 band whose top is the Wiringa VMC 0.848"
  until 2026-09-05, which contradicted the same file's own C5.5 bullet).
- **Exact finite-γ tensor sector.** `InclusiveKernel::Options::tensor_gamma`
  defaults to FALSE, matching `xsec.py`.  True replaces the massless
  Hoodbhoy–Jaffe–Manohar b-sector with the Cosyn Eqs. (9)/(10)/(14)/(16)/(17)/
  (24) kernel (`theta_q_cos_sin`, `cosyn_tensor_sfs`,
  `cosyn_unpolarized_sfs`, `InclusiveKernel::tensor_harmonics_gamma`), which
  carries the unmeasured `b3_func`/`b4_func` slots and leaks the rate sector
  into cos 2φ at O(γ²).  The two paths agree IDENTICALLY at γ = 0 for any b₂,
  with b₃ and b₄ cancelling, which is what makes the switch reversible; it is
  off by default because that leakage is carried as a SYSTEMATIC of the Δ
  extraction rather than as a correction the extraction subtracts, and
  because switching it on would silently move every published tensor number.
  Its size is Δ_fake/(γ² b₁) = 0.14–0.16 — the twist-4 Eq. (17e) term almost
  alone, the leading-twist T_LL and twist-3 T_LT channels standing 3 : −3 : 1
  and cancelling — which RETIRES the old bound "γ² b₁ × 1.15".
  `tests/test_tensor_gamma.cpp` is the port of
  `evgen/tests/test_tensor_gamma.py`, anchored on the two finite-γ rows of the
  paper's own Table 1 at 1e-10.

- **Struck nucleon mass.** The implicit struck nucleon of the per-nucleon
  subsystem is on shell at the FREE nucleon mass `M_NUCLEON` (0.9383), never at
  `Ion::mass_per_nucleon()` (0.9338 for 6Li).  Every per-nucleon label on the
  record — `w2_from_xq2`, `Kinematics::nu` — is built on `M_NUCLEON`, so a
  target at the ion's mass per nucleon would describe a different object from
  its own kinematics.  The binding energy is not lost: it is carried by the
  (A−1) remnant, which the per-nucleon balance deliberately does not write.
  Both places that build that nucleon —
  `InclusiveGenerator::target_nucleon` and `PythiaBridge`'s implicit-target
  branch — use `M_NUCLEON` (P3).
- **Struck nucleon species.** Drawn `Z F2p(x, Q²) : N F2n(x, Q²)` at the
  event's own kinematics, because that is what the inclusive rate
  `Z F2p + N F2n` is made of — never flat `Z : N`, which is its x-independent
  limit and is wrong wherever `F2n/F2p ≠ 1` (0.5 against a true 0.615 for 6Li
  at x = 0.5).  `InclusiveGenerator::proton_fraction` is the one rule;
  `NucleonChoice::ByStructureFunctions` applies it in the bridge (P1).
- **Target mass.** `InclusiveKernel::Options::target_mass` defaults to TRUE,
  matching `xsec.py`.  Identities written against the massless
  `A_par = D(y) g1/F1` must construct the massless kernel explicitly.  The
  finite-γ kinematics themselves (`gamma_squared`, `epsilon_gamma`,
  `depolarization_d_gamma`, `eta_gamma`, `a_parallel_exact`,
  `depolarization_effective`) live in **`asymmetries.hpp`** — ONE
  implementation for both halves of the library, mirroring their move into
  `polli_fastsim.asymmetries`; `xsec.hpp` re-exports only the alias
  `depolarization_gamma`, exactly as `polligen.xsec` does.  `M_NUCLEON`
  (0.9383, the FREE nucleon mass) stays in `constants.hpp`, which is this
  library's single-definition rule for a convention constant.
  `target_mass = true` with `G2Mode::kZero` is legitimate — it is the
  `g2_scale = 0` twist-3 variation, and `g2_scale` is the knob that spans it.
- **Polarized-EMC baseline.** The CBT and TMT curves are quoted on nuclei that
  are not ours, so each is transferred by a valence scale
  `⟨1 − R_unpol,baseline⟩ / ⟨1 − R_unpol,table⟩` over
  `POLEMC_VALENCE_WINDOW`.  That baseline is a CHOICE and is named:
  `EmcBaseline::Epps21` (the default, matching
  `polli_fastsim.polarized.POLEMC_BASELINE`) or `EmcBaseline::LegacyTable`
  (CBT on itself — the pre-2026-08-29 constants 1 and 0.397009).  Both scales
  come out of one code path; the EPPS21 depletion is the single stored number
  `EMC_VALENCE_DEPLETION_EPPS21` = 0.031052077003862335, because computing it
  needs LHAPDF.  Its free-nucleon denominator is **CT18ANLO**, EPPS21's own
  proton baseline, so the fit cancels and the ratio is the nuclear
  modification alone; against CT18NLO it was 0.02979, 4.2 % shallower.  On
  that baseline the two transferred camps are CBT 0.5322 and TMT 0.2113.
- **α–d spectroscopic factor N_αd.** ONE home: `VMC_N_ALPHA_D_LI6` =
  0.80362 + 0.015861 = **0.819481** in **`tagged.hpp`**, the ANL momentum
  file's own two printed norms, and `VMC_P_D_LI6` is now *expressed through it*
  (`0.015861 / VMC_N_ALPHA_D_LI6`) rather than repeating them — the file's
  S- and D-wave norms are defined once.  It is a CHOICE and it carries a stated
  **±5 %** systematic: three tabulations span that (0.819481 from the 2014
  `li6_ad1.momentum`, 0.856 from the 2004 `li6.ad`, 0.863 from Wiringa et al.,
  PRC 89 (2014) 024305), entries 1 and 3 being the same year and Hamiltonian
  family, so the spread is not a version difference anyone can name.
  **THREE NUMBERS SPAN THREE DIFFERENT AMOUNTS — measured 2026-09-04 (§C5.3),
  and "N_αd 5 %, P_D 7 %" was quoting two of them about a third.**  N_αd spans
  **5.311 %** (0.819481 → 0.863), the D-wave NORM spans **7.181 %** (0.015861 →
  0.017), and P_D — the *ratio* every tagged observable actually uses — spans
  only **2.729 %** (0.0193549 → 0.0198830).  And the ±5 % PROPAGATES TO ALMOST
  NOTHING outside b₁: `TaggedModel` renormalises each wave to its own P_L, so
  N_αd cancels out of the tagged sector exactly and only the ratio survives,
  worth **0.048 %** on `tensor_dilution` and **0.082 %** on `vector_dilution`
  across the three readings.  In b₁ it is exactly linear and therefore exactly
  ±5 % (`docs/open_items/run_2026-09-02/phase_D_numbers.md` § "The ±5 % N_αd
  systematic (design §2.1) and the Q4 knob" — the knob is **Q4**, and that
  file has no §Q3; the citation said Q3 until 2026-09-04).  Pinned in
  `tests/test_tagged.cpp` **T26**.  It is a
  knob (`Li6ConvolutionOptions::norm_target`), and the default reading is the
  conservative one: the 18 % of ⁶Li that is not α+d gets b₁ = 0.
  `VMC_S_ALPHA_D_LI6` = 0.81971 is the file's *total* block and is a THIRD,
  different value — do not "unify" them; P_D must be S₂/N with the same N.
- **What each published b₁ curve is PER — the two camps stopped sharing one
  constant on 2026-09-03.**  Every consumer here pairs b₁ with a per-NUCLEON
  F₁, but the two digitized deuteron curves do not arrive normalised the same
  way.  `B1_MILLER_TABLE_TO_PER_NUCLEON` = **0.5** (and it *is*
  `B1_PER_DEUTERON_TO_PER_NUCLEON`, not a second copy of the number) because
  Miller's Eqs. (1)/(5)/(6)/(20) are per deuteron with every ½ in the chain
  spoken for; `B1_CDKS_TABLE_TO_PER_NUCLEON` = **1** because CDKS Eq. (10)
  carries an explicit 1/A and the text under their Eq. (16) says "b₁ is defined
  by the one per nucleon".  The arbiter is neither paper but HERMES, whose
  published b₁ᵈ — the data both camps plot against — is per nucleon by their
  Eq. (5).  The CDKS half is CERTAIN and `CdksB1` doubled; the Miller half is
  LIKELY, not certain, because his own Table I and Fig. 5 compare an unhalved
  curve to HERMES at face value, so the paper is self-inconsistent by exactly
  this factor.  Keeping the 0.5 is the status quo and is an author decision
  (`docs/OPEN_ITEMS_SOLUTIONS.md` §10; argument in
  `docs/open_items/run_2026-09-03/phase_A_miller_normalisation.md`).  A
  consequence worth stating: `close_kumano_integral(false)` and
  `close_kumano_integral(true)` integrate the RAW columns and are therefore on
  DIFFERENT scales — never compare them to each other.
- **The b₁ convolution's δ-function — and the two objects default the OTHER WAY
  from each other, on purpose.**  `DeuteronConvolutionB1::Options::
  finite_q_delta` defaults to **TRUE** (since 2026-09-03), i.e. CDKS Eq. (21)'s
  κ = |q⃗|/ν = √(1 + γ²); `Li6ConvolutionOptions::finite_q_delta` defaults to
  **FALSE**, i.e. Eq. (17)'s κ = 1 as printed.  **Both are CHOICES, not
  derivations**, and the asymmetry is the point: Eq. (21) is CDKS's exact
  definition (their Eq. (18), of which Eq. (17) is the "≃"), and the A = 2
  object exists to reproduce *their* figure, where κ is worth a factor
  **1.57–1.69** at x ≥ 0.5 and **1.97** on ∫b₁ dx.  It is also load-bearing for
  the gate's verdict: on 2026-09-03 the G3b peak ratio came out **0.843** at
  Eq. (21) against **0.520** at Eq. (17), both with CDKS's own MSTW2008 LO, so
  the pass is comfortable at Eq. (21) and clears the window by 4 % at Eq. (17).
  In ⁶Li the same
  switch is worth only −1 % to +7 % (the term κ multiplies is the small orbital
  one), so there the κ = 1 form is kept: it is what makes f(y) a function of y
  alone, one cached table reused at every x, against a 0.24 s density rebuild
  per x.  Measured both ways in
  `docs/open_items/run_2026-09-02/phase_D_gate.md` checklist item 0.
  **R differs between the two as well, and that was decided on 2026-09-03**:
  the gate defaults to `r1998` (CDKS's SLAC world fit) and
  `Li6ConvolutionOptions` to `r_sigma_lt`.  The deciding argument is that the
  observable is a RATIO — the tensor weight is K/D_φ, and D_φ's F₁ comes from
  `InclusiveKernel`'s own `UnpolSF`, whose R is `r_sigma_lt` because nothing
  sets it — so a different R in the numerator would not cancel.  Measured cost
  of the choice on the shipped ⁶Li x·b₁: −5.5 % / +24.6 % / +5.4 % / +2.7 % /
  −1.9 % at x = 0.05 / 0.10 / 0.20 / 0.30 / 0.50, all of it in the two orbital
  terms (term (1) is bit-identical, because it carries a b₁ᵈ table with no R in
  it); on the gate the same swap is +2.7 % on G3b.  Author decision, recorded
  in `docs/OPEN_ITEMS_SOLUTIONS.md` §10 and stated in `docs/USAGE.md` §2a.
  **Since 2026-09-06 the "one R in both" option EXISTS as an opt-in and the
  defaults still do not agree by design**: `--r-source {unset, sigma-lt,
  r1998}` (`PipelineConfig::r_source`) threads ONE R object into
  `Li6ConvolutionOptions::r_func` and `InclusiveKernel::Options::r_func`
  together; `unset`, the default, does not enter that branch, so nothing
  moved.  It is NOT a free repair: at y = 0.5, Q² = 2.5 sharing `r1998` moves
  K/D_φ by −3.8357 % / +26.3988 % / +3.3518 % at x = 0.05 / 0.10 / 0.30
  against the numerator-only −5.5148 % / +24.6099 % / +2.6629 % — larger at
  x = 0.10 — and the shared shift is y-dependent where the numerator-only one
  is not.  Priced in
  `docs/open_items/run_2026-09-06/phase_A_numbers.md` §A1.
  The same file records the quadrature choices: the inner k integral is
  **Simpson**, not `numerics.hpp::trapezoid`, because b₁ at small x is a
  three-decade cancellation of the exact ∫δ_T f dy = 0 and the trapezoid's
  endpoint bias survives it (22 % at x = 0.05 and n_k = 2001) — so `n_k` must be
  ODD and an even value is bumped at construction rather than silently falling
  back to the trapezoid.  Every y-grid integral stays on `trapezoid`, on a
  **2400 / 2400 / 3200** three-segment grid (the design's 600 / 2400 / 800 was
  2.3 % out on the struck-deuteron orbital term at x = 0.1, because that
  constituent's z-width scales with M_α/M_d = 1.987).
- **Coherent |t| range.** `COHERENT_T_MAX_DEFAULT = 0.2 GeV²`, and it STAYS
  0.2.  Since 2026-09-04 the reason is the **anchor range**:
  `mantysaari_a2_deuteron()` is digitized over |t| ≤ 0.30 in four rows and is
  linear in |t| only as |t| → 0, so 0.2 is inside the input and 0.5 is not.
  That reason is a property of the input table and does not move with
  `CoherentScenario::eps_b0`.  The **second** reason is contingent and must be
  quoted as such: the cos 2φ coefficient is linear and unbounded in |t| and
  crosses −1 at |t| = 0.245 for P_zz = −2 — a number now DERIVED by
  `CoherentScenario::t_positivity_edge` rather than repeated, and one that
  becomes 2.80 GeV² at the measured ⁶Li quadrupole — three figures, because
  that is arithmetic on a ROUNDED `eps_b0` = −0.0070 (the derived −0.0070024
  gives 2.7990); never write 2.8000 — i.e. 9.3× outside the anchor.  `CoherentSampler` still checks `CoherentScenario::positivity_margin`
  at SETUP and throws, mirroring `InclusiveKernel::positivity_margin`: that is
  the GUARD on the truncated weight the sampler uses, not the derivation of
  the constant.
- **Coherent x_P.** Per-nucleon pomeron fraction,
  `x_P = (M_X² + Q²)/(W² + Q²)` with the per-nucleon W², drawn log-uniform on
  `[x_P(M_X,min), 0.1]`; the nucleus loses `x_P/A` of its own light-cone
  momentum.  `β = x/x_P` and `M_X²` are on the record.  X must be TIMELIKE on
  every channel — a hard check in `Pipeline::add_hadronic_x` and in
  `InclusiveGenerator`, never a clip.
- **Cluster radial forms.** `ClusterWaveSource::Hulthen` is the DEFAULT
  everywhere; its radial forms are bit-compatible with `polligen` and with
  every number published **since** the 2026-09-06 S–D sign fix (which is every
  ⁷Li number ever published).  The spin-1 tagged numbers published BEFORE that
  date carry the inverted S–D sign and are not reproduced, by design —
  A_zz^tag(0.20 GeV) went +0.845 → −1.207, and `tagged.json`'s `li6_alpha` and
  `deuteron` blocks were re-pinned from the fixed C++ from 2026-09-06 until `polligen` took the same phase (PolarizedLithiumSim `1066555`);
  since 2026-09-23 every block is dumped from `polligen` again (`validation/README.md` §`tagged.json`; `tests/test_tagged.cpp`, and eleven bullets below).  This
  bullet said "bit-compatible with every published number" until 2026-09-16.
  The analytic
  two-parameter forms and their `beta` band (0.20–0.40, default 0.30) are in
  `cluster.hpp`.  `ClusterWaveSource::VmcAV18` replaces them — on all three
  tagged channels: the two lithium α tags and, since 2026-09-04, the deuteron
  control and the ⁶Li α tag's EMBEDDED deuteron (one deuteron per run) — with
  the ANL VMC tables:

  | wave | magnitude |ψ_L(k)| | sign |
  |---|---|---|
  | ⁶Li α+d, L = 0 and L = 2 | `momenta/li6_ad1.momentum` S/D blocks, ψ_L = √ρ_L | `overlap_old/li6.ad`, k-space amplitudes |
  | ⁷Li α+t, L = 1 | `momenta/li7_at3.momentum` (3/2⁻ ground state) | none needed — a single wave has no interference |

  Two files per wave because they carry different things.  A momentum density
  is |ψ_L|² and has no phase, but the **S–D relative sign is observable**: the
  interference term of `n_M(k, k̂)` goes as ψ₀ψ₂ and, at P_D ≈ 0.019, it
  DOMINATES the tensor asymmetry (negating the D table flips A_zz^tag by 0.97
  at k = 0.20 GeV).  The signed `overlap_old` amplitudes supply
  it.  **The two numbers swapped roles on 2026-09-06.**  This sentence read
  "flips A_zz^tag from +0.45 to −0.52"; the shipped answer is now **−0.5191**
  and negating the D table gives **+0.4518** — measured at the k = 0.1979 GeV
  cell, acceptance-weighted at the YR high-acceptance optics.  The claim
  itself did not move; what moved is which side of it the code is on, because
  the tagged amplitude was summing ψ₂ where it needed φ₂ = i²ψ₂ = −ψ₂
  (`src/core/tagged.cpp` `build_amp2`;
  `docs/benchmarking/07_cw_sign_investigation.md`;
  `docs/open_items/run_2026-09-06/phase_CW_numbers.md`).  `vmc_from_momentum` does not copy sign(A_L(k)) point by point — it
  reduces the reference to its zero CROSSINGS below 3 fm⁻¹ (⁶Li: one S node at
  0.678 fm⁻¹ = 0.134 GeV, one D node at 2.25 fm⁻¹ = 0.444 GeV, both confirmed
  bin-for-bin by minima of the momentum file's own ρ_L) and anchors the phase
  at the reference's largest |A|, where its Monte Carlo sign is beyond doubt.
  Above ~3 fm⁻¹ both columns are at the noise floor and carry ~1e-4 of the
  norm, so a noise-driven flip there would be all cost and no signal.  The
  global phase is then fixed to ψ₀(k → 0) > 0, which is unobservable but makes
  sign(ψ₂) read as the relative sign.
  On the VMC path `beta` and `p_d` are IGNORED: the shape is the table's and
  P_D is a property of the wave function (`VMC_P_D_LI6` = 0.01935, the file's
  own 0.015861/(0.80362+0.015861), against the 0.0867 SCENARIO placeholder).
  Tables are ZERO outside 5 fm⁻¹ = 0.9866 GeV, so the sampler cannot draw a
  spectator past it.  The deuteron control channel is Hulthén **by default and
  follows `--cluster-wave vmc`** (C5.4, 2026-09-04): `fdeut.av18`'s u(k), w(k)
  ARE the p–n relative waves, so the sentence that used to stand here — "the
  deuteron control channel is always Hulthen … no d → p+n two-cluster table
  exists" — was false, and the flag had been *silently ignored* there.
  `deuteron_channel(beta, p_d, source)` now takes the flag and selects
  P_D = 0.057600.  **The Cosyn–Weiss tensor gate runs on the AV18 control
  since a7b3d18** (2026-09-06; `tests/test_tagged.cpp`, registration-skipped
  without `data/vmc/deuteron/fdeut.av18`), because CW quote TABLE II *for*
  AV18 — the Hulthén pair's f₂/f₀ never reaches √2 anywhere on the grid (it
  tops out at 1.2866) and cannot carry the landmarks
  (`docs/benchmarking/07_cw_sign_investigation.md` §8); the Hulthén pair keeps
  only the two statements that do not depend on the wave function. The three
  cell-centre rows are pinned at −1.937 / +0.999 / +0.967 (atol 1e−4).
  This bullet said "pinned on the Hulthén DEFAULT, which is bit for bit what
  it was" until 2026-09-16 — false since a7b3d18, and the pins had moved from
  +0.617 / −0.318 / −1.656 with the S–D sign fix.
  `AUTHOR_DECISIONS.md` §B9(c) and `run_2026-09-06/STATUS.md` both recorded
  the change; this was the site that did not.
  Reconciliation, provenance and every number:
  `docs/open_items/vmc_reconciliation.md`, `data/vmc/README.md`.
- **Triton spectral function.** `TritonSfChoice::Hulthen` is the DEFAULT
  everywhere and keeps the sequential two-body triton decay of `breakup.hpp`
  bit-compatible with every published ⁷Li number.  `TritonSfChoice::CiofiSimula`
  (CLI `--triton-sf ciofi-simula`) replaces it, on the ⁷Li α tag only, with
  the Ciofi degli Atti–Simula spectral function (`triton_sf.hpp`).
  Provenance, all transcribed and none refit: n₀(k) is CS PRC 53 (1996) 1689
  Eq. (74) with Table A.1 (A = 2, 3, 4 carried as data — the same numbers
  BeAGLE's `DT_KFERMI` holds, checked column for column against the paper);
  n₁(k) is Eq. (76) for A = 3 and Eq. (75)/Table A.3 for A = 4; k in fm⁻¹ is
  converted at ħc = `HBARC_GEV_FM` = 0.19733 GeV·fm — a `constants.hpp`
  single definition since 2026-09-01, with `coherent.hpp`'s `GEV_PER_FM_INV`
  now an alias of it.  S₀ = ∫dk k²n₀ = 0.6525 is the ³He(e,e′p)d
  spectroscopic factor (the literature's ~2/3), S₀ + S₁ closes on 1 to
  3.4 × 10⁻⁴ untuned, and the 2-body/3-body branching is the RATIO
  n₀/(n₀ + n₁) at the event's own k, never an assumed constant.  Table A.1's
  A = 3 column is the ³He **proton** distribution; under the isospin mirror
  ³He ↔ ³H it is the triton's struck-**neutron** distribution and the bound
  remnant of the mirror is the same deuteron (BeAGLE branches on A alone with
  no mirror argument).  The continuum pair is split at its own virtual-state
  pole: nn at `KAPPA_NN_VIRTUAL` = ħc/|a_nn| = 10.44 MeV (a_nn = −18.9 fm),
  pn at `KAPPA_PN_SINGLET` = ħc/|a_pn| = 8.31 MeV (a_pn = −23.74 fm, the
  ¹S₀ channel — the ³S₁ pn channel IS the bound deuteron and is the other
  branch); both constants live in `triton_sf.hpp`, one definition each.
  **The BeAGLE n₀-only caveat**: `DT_KFERMI` renormalizes its A = 3 and A = 4
  tables — the n₀ piece ALONE — to unity, silently dropping the 34.7 %
  (A = 3) / 20.0 % (A = 4) correlated continuum: its high-k tails there are
  too soft by construction and it has no three-body breakup channel to put
  that strength in.  LiPolGen does not renormalize n₀; it adds n₁ and lets
  the deficit BE the branching.  `sample()` consumes a fixed 6 uniforms on
  every branch (the breakup branch a fixed 9), so the stream position never
  depends on which channel came out.  The `Pipeline` constructs the
  `CiofiSimulaTriton` itself at the run's own `cluster_beta`, so the
  continuum pair's q-shape and every other radial form share one β.
  Full derivation and every measured number: the `triton_sf.hpp` header and
  `tests/test_triton_sf.cpp`.
- **Spectator FSI.** `PipelineConfig::fsi = PipelineFsi::Off` is the DEFAULT:
  every tagged number is the plane-wave impulse approximation bit for bit
  unless a run turns the Glauber weight on.  Three rules when it is on
  (`fsi.hpp`, `tests/test_fsi.cpp`):
  (1) **FSI is a WEIGHT, never a shift** — the rescattering multiplies
  `Event::weight` (the FSI/IA density ratio at the drawn (k, cos θ_k)) and
  moves NO four-vector; the spectator is the measurement, and the tagged
  record of an FSI-on run is bit-identical to the FSI-off run.
  (2) **The weight is SPIN INDEPENDENT by construction** — numerator and
  denominator are m-summed, so it is the unpolarized-shape distortion, which
  is what the literature supports (Cosyn–Weiss VI C leaves the spin
  dependence open).  Quote its effect as an **unpolarized-shape systematic of
  the tagged spectrum, never as a correction to A_zz**.
  (3) **σ_XN is a BAND, 20–40 mb, never one number** — 40 mb (the default,
  `--fsi-sigma-mb`) is the free hadron, 20 mb the formation-length end; run
  both.  The survival probability (∫w dΓ/∫dΓ ≈ 0.52 at 40 mb on the ⁶Li α
  tag) is logged by the run summary and exposed as
  `Pipeline::fsi_weight()->survival()`, and since 2026-09-04 the whole FSI
  block is in the npz `meta` (`fsi`, `fsi_sigma_mb`, `fsi_sigma_cluster_mb`,
  `fsi_sigma_cluster_el_mb`, `fsi_survival`, …) — without it the two ends of
  the mandatory band were indistinguishable files, which is the same failure
  the `b1_*` block exists to prevent.  Variant (a) `GlauberCluster`
  shadows σ_Xα over the α's own profile (131.0 mb at σ_XN = 40, not
  4 × 40); variant (b) `GlauberNucleon` is the unshadowed A σ_XN = 160 mb
  single-scattering bracket — a POINTWISE low-k bracket, not an integrated
  one (its survival lands above (a)'s; `fsi.hpp` header).  **The two are NOT
  one**, which the open-items inventory claimed until 2026-09-04: measured on
  ONE event stream reweighted both ways — `--channel tagged-6Li-alpha
  --events 20000 --seed 1234`, `tensor-thirds` at P_z = 0.7 / P_zz = 0.6 /
  P_e = 0.7, at the default σ_XN = 40 mb, with `k`, `cos_theta_k`, `phi_k`,
  `x` and `q2` bit-identical across the three runs — they differ on
  **99.50 %** of events by more than **1 %** (|w_nucleon/w_cluster − 1| > 0.01;
  Σw 10419.07 vs 11632.10, i.e. Σw/Σw_off **0.520954** vs **0.581605**) and by
  up to a factor **68.52** per event.  Σw/Σw_off is the SAMPLE mean weight and
  is not the model's own grid-integrated `GlauberFsiWeight::survival()`, which
  is **0.520239** vs **0.582899** — the number `meta["fsi_survival"]` carries
  and `tests/test_fsi.cpp` pins.  Quote that window with the two numbers; the
  whole claim is about one stream and one σ_XN.  The weight table
  is built on a (k_z, k_T) grid — a deliberate deviation from the design
  note's literal (k, cos θ_k): the eikonal kernel transfers k_T only, so the
  table is smooth and even in k_z there — and read per event by bilinear
  interpolation, consuming no randomness, so stream discipline and
  event-index determinism are untouched.
- **The WEIGHT-FAMILY rule — which weight a new effect goes on.** There are
  now two kinds of per-event weight in the library and they are NOT
  interchangeable:
  (1) A **correction to the model** multiplies `Event::weight`. FSI is the
  only one: the plane-wave impulse approximation is *wrong* and the Glauber
  ratio *fixes* it, so the nominal sample must carry it.
  (2) A **systematic variation** or a **background** goes in its own named
  weight family, never on `Event::weight`. The RC band
  (`rc_tensor_lo`/`rc_tensor_hi`) is a variation and `rc_tail` is a
  background; both leave the Born sample alone and an analysis multiplies one
  in on purpose. They travel in `Event::rc_weights`, appear under their own
  HepMC3 names (APPENDED after `nominal` and every `spin_weight_k`, so no
  existing index moves — `docs/HEPMC3_CONVENTION.md`) and their own npz
  columns, and those columns exist **only when the family is on**, on BOTH
  the C++ columnar path and `export.columns_from_events`.
  A new family must also **consume no randomness and move no four-vector**,
  so that switching it on is bit-for-bit invisible everywhere else — the
  structural form of that promise is a `fill(Event&) const` that takes no
  `Rng&` (`rc.hpp`), and the runtime one is a test that compares two
  pipelines' hadronizer-hook draw logs with `==`
  (`tests/test_rc_pipeline.cpp`, T6).
- **Tensor-sector RC band anchors.** `PipelineConfig::rc = PipelineRc::Off` is
  the DEFAULT and is today bit for bit, down to the bytes of an `--rc off`
  npz and HepMC3 file. When it is on, `rc_delta` interpolates log-linearly in
  x between two anchors that are **not the same kind of object**, and the
  distinction is the choice:
  (a) `RC_DELTA_HIGH_X = 0.015` at `RC_X_HIGH = 0.16` is a quoted
  **uncertainty** (JLab E12-13-011 / PR12-13-011; x = 0.16 is that
  experiment's own lower kinematic edge, arXiv:2506.04506 p. 8 — the proposal
  itself is unpublished, so cite it by page or drop the anchor). **There is
  no alternative value to band it against and the record names none**
  (checked 2026-09-06, `docs/open_items/run_2026-09-06/phase_A_numbers.md`
  §A2): "cite it by page or drop it" is an *action*, not a second reading.
  It is unbanded because nothing exists to band it against, not because it is
  better known than (b).
  (b) `RC_DELTA_LOW_X = 0.30` at `RC_X_LOW = 0.01` is the **size of a
  correction this generator does not apply**, taken as a 1σ band — the
  conservative end of Gakh–Shekhovtsova's (hep-ph/0403262, **zero INSPIRE
  citations**) 10–30 %. It is NOT a measured residual. And it is that
  paper's panel value **carried upward in x**, not read at 0.01: the one
  panel its sentence covers spans x = 0.00226–0.00966 at Q² = 0.1, and
  **0.113 is the value it actually reads at x = 0.00966, the x nearest this
  anchor** (0.266 at its own bottom), so the shipped 0.30 errs **wide by
  ×2.65** — the safe direction for a half-width, the wrong one for a quoted
  precision. Both alternatives are **priced, not adopted** (registry row 17,
  `docs/open_items/run_2026-09-06/phase_A_numbers.md` §A2): half-widths on
  A_zz at x = 0.01, Q² = 5 of **1.358472e−04 / 1.204512e−04 / 5.116913e−05**,
  unmoved at x ≥ 0.16, and no clipped-event fraction moves at all.
  (c) `RC_DELTA_LOW_X_OPTIMISTIC = 0.19` IS a measured residual: HERMES's own
  fractional RC systematic at its lowest-x bin (2×10⁻³ on A_zz = −1.06×10⁻²,
  hep-ex/0506018 Table II). **Run 0.19 and 0.30 both; never quote one row
  alone**, exactly as with the FSI σ_XN band.
  (d) The ⁶Li elastic form-factor **normalisations are the MEASURED moments**
  — `LI6_MU_N = +0.8220473 μ_N`, `LI6_QUADRUPOLE_FM2 = −0.0818 fm²` (TUNL
  A = 6) — and never VMC: Wiringa–Schiavilla's Q(⁶Li) = −0.23(9) fm² is 3×
  the measured one, and `OPEN_ITEMS_SOLUTIONS.md` §3–4's (restated in §11.3)
  standing instruction
  forbids deriving a ⁶Li tensor input from those wave functions. Their
  *shapes* are unfitted starting values and are banded by `--rc-fq-scale`
  (0/1/2 — σ^el_T is QUADRATIC in it, so RUN the band, never rescale one row)
  and `--rc-tail-tensor-scale` (0.5/1/2, the η F_m² sector `fq_scale` does
  not span).
- **Where data files live at run time.** `data_dir()` (`cluster.hpp`) is
  `$LIPOLGEN_DATA_DIR` when set and non-empty, else the compiled-in
  `LIPOLGEN_DATA_DIR_DEFAULT`, which CMake sets to `${CMAKE_SOURCE_DIR}/data`.
  Nothing resolves a data path against the working directory.
- **Luminosity shares.** `Optics::lumi_fraction` multiplies COUNTS and never
  cross sections (`PipelineConfig::apply_optics_lumi_fraction`, default true) —
  the same share rule that keeps `Scenario::run_share` out of
  `sigma_per_category_pb()`.
