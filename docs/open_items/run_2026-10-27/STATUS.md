<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-27 — status board

The author's decisions on the priced options of run 2026-10-26
(`DECISIONS.md`: D4 the Jaffe–Manohar Δ convention adopted in the kernel, D5
`PolradFull` default, D6 the VMC cluster waves default, D7 the grid
structure-function backends default), researched, designed, implemented and verified by Opus workers
under the supervising session. Branch `claude/run-2026-10-26`, on top of
`ca9b23a`.

| phase | state | commit | what | record |
|---|---|---|---|---|
| R1 — which Δ normalisation the literature uses (decision 1's question) | DONE | `fbc24f2` | two readers (1989–2016 + the LOI; 2019–2026 + the tree and the letter) and a judge | `phase_R1_delta_convention.md` |
| D5 design — `PolradFull` default | DONE, checked | `fbc24f2` | inventory of every pin, label and document; adversarial check | `plan_D5_polradfull_default.md` |
| D6 design — VMC cluster-wave default | DONE, checked | `fbc24f2` | same; the `LI6_CLUSTER_POLARIZATION` question surfaced for the author | `plan_D6_cluster_wave_default.md` |
| D7 design — grid backend defaults | DONE, checked | `fbc24f2` | same; CI / missing-set failure mode designed | `plan_D7_backend_defaults.md` |
| D4 design — the Jaffe–Manohar convention in the kernel | DONE, checked | `fbc24f2` | inventory of every site the ½ touches: the a₂ term, the finite-γ path, the bindings, the toy's amplitude, the bound row's S2/A7 logic, the port gate, every cos 2φ pin, the documents and their attribution; adversarial check | `plan_D4_jm_convention.md` |
| I4 — the Jaffe–Manohar convention implemented | DONE, verified and fixed (F4) | `9779edd` | per `plan_D4` amended by its adversarial check: `DELTA_COS2PHI_NORM` = 0.5 at the two evaluation sites (every `a_cos2phi` and every Δ line of a₂ exactly half, the whole a₂ on the massless path; everything else bit-identical); T1–T5 and Y1–Y9 re-taken with dated comments, the missed `ds16_over_tree` pin too, N1 and N2 added; the port gate, C++ and Python, reads polligen's coefficient off `xsec.json` against a literal ½ and takes a Jaffe–Manohar re-dump unchanged, at the pre-D4 tolerances; `meta["delta_convention"]` on `Pipeline.generate()` and so every CLI .npz (HepMC3 carries none: an open point) and the run banner, the label tied to the constant by a `static_assert`; the convention paragraph with the one-line rule and the corrected attribution (HJM has no Δ term); row 13 PASS, relabelled, V6 kept in TIGHT; the report regenerated (row 13 only). Verification: 2 medium and 6 low findings, none refuted, all fixed in F4 (the record's §10); six drifted citations in the plan records named for the supervisor; every gate green after F4 | `phase_I4_jm_convention.md` |
| I5 — `PolradFull` default implemented | DONE, verified and fixed (F5) | `cb75ec7` | per `plan_D5` amended by its adversarial check: `RcOptions::tail_model` = `PolradFull` (one line, `rc.hpp`), the CLI default READ from it; the `PolradFull` preconditions bind only where a tail is built and `m_lepton` joins the tail sub-knobs refused off the tail (the two code consequences); `--rc off` and every explicit-model run bit for bit (A/B column hashes); the t-peak pins named (7 C++ cases / 30 assertions failed on the flip alone, each a t-peak pin; B6 found by the step-1 gate), runtime pins per spec — the knob matrix's pin is a sentinel resolved on the inclusive specs only, so no cell that ran stops running (511 / 151 → 517 / 145: six skipped cells now run); every t-peak number in the code and documents labelled and the default's re-measured beside it (×768 → ×1152, the QE shares 22 / 73 / 99.9 → 39 / 99.8 / 100 %, …), two readings that change recorded (the `qe_tensor_scale` stand-in is 14.5× the band at x = 0.30 on the default; `qe_kf_gev` is the smaller QE knob there); row 7 PASS, relabelled "polrad-full (shipped), t-peak (option)", a dated addendum after the hashed block; T6 16.3 s (kept on the default); gates green: C++ 419 / 419, pytest 1469 / 149 / 0, bench 350 / 1 / 1 → the report test 26 / 26 after regeneration, strict 1241 / 95 / 7 / 0 / 24, R5 0 refused, R6 0 broken, `--compare` 0 changes, `--check` OK; no tail-table cache and the HepMC3 provenance gap filed. Verification: 1 medium and 7 low findings, none refuted, all fixed in F5 (the record's §12): `OPEN_ITEMS_SOLUTIONS.md` §9's `w_tail − 1` table and its "five to eight orders" headline labelled `t-peak`, the default's values beside them (×1.02–×1.28 at x = 0.01, Q² ≥ 5; ×4.9 to ×1.2e4 at x ≥ 0.10 and at Q² = 2); `validate()`'s `sp_tensor_scale` clause binds only where a tail is built and `RcModel`'s refusal says "no tail is BUILT" there (no accept/refuse outcome moves; C++ checks appended to the two D5 cases); the remaining unlabelled t-peak numbers (σ^el_T/σ^el_U on both C0 edges, the C0 ratios, cost and spreads, the §9 ladder, whose `qe_kf_gev` and `qe_suppression` swap on the default) labelled with the default's beside them; the `m_lepton` comment in `constants.hpp`, one knob-cell label, the record's reproduction claim, the "+1.9 … +2.5 %" phrasing and the 0.106 s / 0.16 s pair corrected; one stale D5-case comment found by F5. Gates after F5: C++ 419 / 419 (17 244 240 assertions), pytest 1469 / 149 / 0 (per-file counts I5's; the 7 relabelled ids map one to one), bench 351 / 1 / 0, strict 1241 / 95 / 7 / 0 / 24, R5 0 refused, R6 0 broken (1432 citations), `--check` OK, SPDX 133 | `phase_I5_polradfull.md` |
| I6 — VMC cluster waves default implemented | DONE, verified and fixed (F6) | `c52c164` | per `plan_D6` amended by its adversarial check: `CLUSTER_WAVE_DEFAULT` = `VmcAV18` in ONE place (`cluster.hpp`); `PipelineConfig::cluster_wave` reads it, the bindings export it with `cluster_wave_name`, and the CLI default (`None`) and help, the example and `price_defaults.py` READ it; the factories and `BreakupOptions` keep the Hulthén (DD-1: `tagged.json` and every factory pin untouched, the dump refuses a sibling-side flip); the three code consequences: `validate()` refuses a value MOVED from the default off the tags, `make_config` / the CLI an explicit `--cluster-wave` there, and the off-channel knob row reads "not read on channel …", never a bare `vmc`; `meta["cluster_wave"]` on every npz and a banner line on every tagged run (HepMC3 carries none: filed); a missing table names its escape; the example's `--beta` / `--pd` refusals read the knob table (channel × tier × wave) and its banner the true tier; A/B over 28 file pairs: every explicit-wave, inclusive and coherent run column-identical, and the pre-I6 `vmc` run = the post-I6 default in every column (in `meta` only the new key and the knob row's `at_default` / reason); the physics lens reproduces the design (⁶Li A_zz^tag −1.2069 / −0.9533 → −0.5191 / −0.1993, RP 0.0247 → 0.0342, ⁷Li ⟨k⟩ 0.2885 → 0.1861 GeV); the flip alone failed 1 C++ case (4 assertions) and 14 pytest tests — six Hulthén-default pins re-taken, eight citation shifts repaired by `--fix` and two re-read, re-recorded blocks; the +18 assertions attributed to three per-event loops; the knob matrix re-taken with the wave in base AND var (517 / 145 before and after, every id mapped); every default-run number the design, the check and this phase found re-measured on both waves and labelled (D3 FSI, the RC clamp on its published recipe 245 / 520 → 16 / 507, the tag fractions, the ⁷Li FSI ceiling, the T1 virtuality and redraws; one site the design listed, `08_channel_priority.md` ingredient I's 99.5 %, only in F6); the α–d D/S answered on the α–d channel (Hulthén η = +0.070, the wrong sign; VMC −0.048; GK −0.025 ± 0.012); `LI6_CLUSTER_POLARIZATION` / `LI6_B1_RANK2_TRANSFER` KEPT, the inconsistency put to the author first; no verdict moved (Deeps keeps its declared Hulthén verdict, AV18 recorded; labels on rows 10, 11, 16 and the static row; `--compare` 2 MOVED generator strings, 13 pass / 10 fail as before, T6 unchanged, `--check` OK); `price_defaults.png` regenerated; registry rows 9, 10, 11 (and 23) amended, the D2 pointer added, the D6 entry's corrections proposed. Gates: C++ 422 / 422 (17 244 309 assertions), pytest 1472 / 149 / 0, bench 350 / 1 / 1 → the report test 26 / 26 after regeneration, strict 1241 / 95 / 7 / 0 / 24, R5 0 refused, R6 0 broken (1432 citations, 77 records), SPDX 133. Verification: 1 medium and 8 low findings, none refuted, all fixed in F6 (the record's §12): `price_defaults.png`'s fixed 0.1 … 2.6 axis had dropped the ×4.78 `hulthen` A_zz^tag(0.30) price — the axis now follows the drawn points and the figure was regenerated from a table identical to I6's; in `08_channel_priority.md` ingredient I's 99.5 % labelled per wave (64.5 % on the default), ingredient G's lever and item turned to `--cluster-wave hulthen` against the default and "taken by D6", and the D6 log row's inverted ratio corrected (the default ×0.43 / ×0.21 of the opt-in); the missing-table message offers `--cluster-wave hulthen` on the tagged channels only (an inclusive `li6-convolution` run was told to use it and then refused); the T1 redraw fractions corrected to 0.038 / 0.022 / 0.002 % on Hulthén; a `static_assert` beside `CLUSTER_WAVE_DEFAULT` ties the texts that call vmc the default; `validate()` refuses a non-finite `cluster_vmc_mc_sigma` (NaN had passed and failed later as "the plan has no rate"; a finite bound is left to the author); the example's row in `PHYSICS_CHANNELS.md`, the leftover Deeps labels, registry row 11's wording, and this record's moved-citation list (a fourth, plan_D6's test_pipeline.py 77-82) and knob-matrix count (42 renamed + 4 same-id swaps). Gates after F6: C++ 422 / 422 (17 244 315 assertions), pytest 1472 / 149 / 0 (per test identical to I6's), bench 351 / 1 / 0, strict 1241 / 95 / 7 / 0 / 24, R5 0 refused, R6 0 broken (1432 citations, 77 records), `--check` OK, SPDX 133 | `phase_I6_cluster_wave.md` |
| I7 — grid backends default implemented | DONE, verified and fixed (F7) | `7fd441c` | per `plan_D7` amended by its adversarial check: `UNPOL_SF_DEFAULT` = `Ct18Nlo` and `POL_SF_DEFAULT` = `NnpdfPol` in ONE place (`pipeline.hpp`), read by the struct, the bindings (exported), `make_config`, `cli.DEFAULTS` and the examples, two `static_assert`s tying the texts; the grids ATTACHED LAZILY by the LHAPDF tier's filler when a `Pipeline` is built (nothing loaded by `PipelineConfig()` / `make_config()` / the CLI import; `install_lhapdf_sf_filler()`, `attach_shipped_sf_backends`, the core-only path `toy_config()`; `b1_unpol` stays eager); every grid FROZEN at its knot edges (LHAPDF `nearest`, `logcubic`, ForcePositive 0 set on the object, plus a clamp) — an `lhapdf.conf` differing from this store's in the interpolator and ForcePositive (not the extrapolator; the fix stage added an `Extrapolator: error` variant) leaves every default number bit-identical (the same pair moved the pre-D7 ⁷Li rate 4.5e−4; `DECISIONS.md`'s +5 % is a different comparison, continuation against `nearest`) and `meta` records the extrapolation and both below-grid fractions; a pre-flight refuses a missing set by name before any LHAPDF `Config` (exit 1, never 134; LHAPDF's text kept); one kernel-beside-selectors rule, struck slots `Custom`-only, the T2 species rule in C++, the positivity refusal computing its cure; the class-level kernels and every reference JSON untouched, the toy NAMED the pre-D7 run bit for bit (24 runs + the kernel); `--b1-unpol` toy, x_max 1.0, EPPS21 and R1998 not moved. Measured: ⁶Li σ ×0.83862 (591846.17 → 496334.16 pb), coherent ×0.75541, ⁷Li ×0.83604, A_zz ×1.19243, 39.23 % of the rate on the frozen edge (48.37 % coherent); six `helicity-flip` runs, the b₁ band's ×2 row (needs 0.85), `li6-convolution` without `--b1-unpol` and `--inclusive-b1` on the ⁶Li α tag refused at their old settings, the inclusive-deuteron split measured; `make_config` 0.14 ms (37 ms when naming the grid pre-D7), first `Pipeline` 0.25 s, a 5000-event run 0.38 s (0.17 s on the toy); ⟨x⟩ +23 %. The flip alone failed 51 C++ cases and 100 pytest tests; every one re-taken with a dated D7 comment or answered by naming the toy; the knob matrix 627 / 627 cells with their c52c164 outcomes (52 ids relabelled); the I4 / I5 / I6 toy numbers re-measured on D7 and labelled; verdicts moved only where the declared rule names the shipped default (`t1_nmc_f2d_over_f2p` FAIL → PASS; the rule's split recorded for the author); row 7 on the toy Born; nine record citations annotated `(as of c52c164)`, one re-pointed, the census 110 → 119; `build_deps.sh` pins the three sets, CI gains a sets-absent step and a reworked wheel smoke (both re-enacted); USAGE §10 and PACKAGING §2.5 state the C++ API contract and the HepMC3 gap. Gates: C++ 427 / 427 (17 244 377), pytest 1482 / 149 / 0 (737 s), bench 351 / 1 / 0 after the report, strict 1241 / 95 / 7 / 0 / 24, R5 0 refused, R6 0 broken (1432 citations, 78 records), SPDX 136, report 14 pass / 9 fail, `--compare` 34 changes, `--check` OK  Verification: 1 medium and 9 low/info findings, none refuted; F7 (interrupted 2026-10-30, resumed 2026-11-03; the record's §16): a missing set blocks the four T1 rows `blocked (environment)` again (each harness raises LHAPDF's words; make_report end-to-end tested with the sets hidden) and `price_defaults` marks the probes of an environment default; the HEADLINE literals checked against the shipped defaults; the polarised meta/banner keys credited only where g₁ ran; LO and free-nucleon caveats wherever "shipped" is said; python/README's T2 example works on a default config; the t5 shipped-row error RECORDED without a pass bound (1e−15 kept on the toy); USAGE §8 labelled and re-measured, stale wording corrected; build_deps re-verifies a present set; the lhapdf.conf-only advice fixed. Gates after F7: C++ 427 / 427 (17 244 377), pytest 1497 / 149 / 0 (854 s), bench 364 / 1 / 0, strict 1241 / 95 / 7 / 0 / 24, R5 0 refused, R6 0 broken (1432 citations, 78 records), SPDX 136, report 14 pass / 9 fail, `--check` OK | `phase_I7_backends.md` |
| V — adversarial verification of the implemented tree | DONE | (this commit) | four lenses over `7fd441c` — physics, honesty, code, regression: 0 high, 5 medium (honesty: unlabelled pre-D7 continued numbers in OPEN_ITEMS §14, the I7 record's quote and freeze attribution, no consolidated author list; physics: the LO pair's |A₁| > 1 described too narrowly; code: a foreign object under a grid selector ran under the grid's label), 2 low-medium, the rest low / info; none refuted | `phase_V_verification.md` §1–§4 |
| F — fixes | DONE | (this commit) | every finding applied (none rejected): the grid-slot refusal with C++ and Python tests; |A₁| ≤ 1 broken from x ≳ 0.83 (n) / 0.87 (d, p) stated and filed under the LO development item; every pre-D7 continued and toy number found labelled with the frozen / shipped value beside; DECISIONS D4 (I4 §8 item 8's corrections), D6 (VMC's own 1.9 σ, the unsurveyed premise) and D7 ("grid sets … free nucleon", the Q²-only Q₀ freeze, the two lhapdf.conf comparisons) corrected; the T1 cuts printed from `q2_floor` (report regenerated: 3 MOVED strings, 14 / 9 unchanged); a missing `lhapdf.conf` blocked (environment); the knob asymmetry stated in CONVENTIONS and pinned; pythia-flagged and shipped-tail knob cells (641 cells); the D5 texts tied by a test; CI key -v2; the author list below | `phase_V_verification.md` |

Ground rules of the run: the author's four decisions are the only
defaults that move; every number that depends on them is re-taken with a
dated comment naming the decision, or the test passes the old option
explicitly when its purpose is to pin the old physics; tolerances declared
before verdicts; a disagreement recorded, never tuned; the hashed design
block of `t4_djangoh_rad_noRad.py` is never edited; nothing published or
sent; the sibling repository read-only; fetched papers stay outside the
repository.

## For the author

Every open item of the run, each with the record section that carries it
(collected by the fix stage, 2026-11-03).

1. **The ⁶Li cluster constants (D6).** `LI6_CLUSTER_POLARIZATION` 0.811228
   and `LI6_B1_RANK2_TRANSFER` 0.921947, read by the default inclusive run,
   against the VMC wave's 0.887076 / 0.982576 — keep, re-derive, or 0.848;
   and the deuteron species' 0.9325. `phase_I6_cluster_wave.md` §9 item 1;
   `DECISIONS.md` D6.
2. **EPPS21 and R1998 (D7)**: the EPPS21 ⁶Li ratio by default for lithium
   beams (no pipeline path applies it) and R1998 in place of `r_sigma_lt`.
   `phase_I7_backends.md` §13 item 7; `DECISIONS.md` D7.
3. **The LO development item, with the |A₁| violation**: both shipped grids
   are NLO fits evaluated at LO; at LO the pair breaks |A₁| ≤ 1 from
   x ≳ 0.83 (free neutron) and x ≳ 0.87 (deuteron, proton), and `--x-max
   0.95` prevents the refusal, not the violation. `DECISIONS.md` D7;
   `phase_I7_backends.md` §1, §3.5.
4. **Two CLI plans refused at the default**: `--isotope d --plan
   helicity-flip` and `--channel tagged-d-p --plan helicity-flip` exit 1 at
   x_max 1.0 (cure `--x-max 0.95`). `phase_I7_backends.md` §3.5.
5. **The T1 "shipped default" rule, two readings**: I7 re-pointed the T1
   headlines to the shipped default (F₂ᵈ/F₂ᵖ FAIL → PASS); I6 kept Deeps on
   its declared Hulthén. The other reading makes the report 13 / 10.
   `phase_I7_backends.md` §6, §13 item 1; `phase_I6_cluster_wave.md` §9
   item 12.
6. **The t5 TIGHT edge (I4)**: whether it keeps V6, now the pre-D4 reading.
   `phase_I4_jm_convention.md` §8 item 5.
7. **The t5 grid row has no numerical-error bound (I7)**: recorded, not
   bounded; 1e−15 stays on the toy. `phase_I7_backends.md` §13 item 2, §16.
8. **A harness exits 2 when a set is missing (I7)** (README convention 6),
   where the brief said 0; `make_report` exits 0 with `blocked
   (environment)`. `phase_I7_backends.md` §16, left item 2.
9. **CMake `FATAL_ERROR` for a missing `lhapdf-config`** unless
   `-DLIPOLGEN_WITH_LHAPDF=OFF` is explicit (I7). `phase_I7_backends.md` §13
   item 3.
10. **HepMC3 provenance (all four decisions)**: a HepMC3 file records neither
    the Δ convention, the tail model, the cluster wave nor the SF backends;
    adding it changes the sha256-pinned `docs/HEPMC3_CONVENTION.md`.
    `phase_I4_jm_convention.md` §8 item 10; `phase_I5_polradfull.md` §9
    item 2; `phase_I6_cluster_wave.md` §9 item 4; `phase_I7_backends.md` §13
    item 6.
11. **A finite bound on `cluster_vmc_mc_sigma`** (I6). `phase_I6_cluster_wave.md`
    §9 item 11.
12. **The ⁷Li FSI ceiling on the default wave** (I6): 4.41 % / 9.52 % of the
    weight table clipped at 40 mb; band it, widen the grid, or print it.
    `phase_I6_cluster_wave.md` §9 item 3.
13. **A polligen VMC port-gate block** (I6, a candidate).
    `phase_I6_cluster_wave.md` §9 item 2.
14. **I5's two changed readings**: the `qe_tensor_scale` stand-in is 14.5×
    the band at x = 0.30 on the default tail; `qe_kf_gev` is the smaller QE
    knob there. `phase_I5_polradfull.md` §9 item 5.
15. **The sibling sites to change (I4)**: `asymmetries.py`,
    `evgen/polligen/xsec.py`, the letter's Eqs. (1) and (3) and their
    consumers carry the old cos 2φ coefficient; whether the letter means the
    amplitude or Δ to halve is yours. `phase_I4_jm_convention.md` §7.
16. **USAGE §8's D7 throughput** (0.65 / 0.36 / 1.08 M ev/s) has no surviving
    log and was not taken on an idle machine. `phase_I7_backends.md` §16,
    left item 3; `phase_V_verification.md` H6.
17. **The knob asymmetry across D5 / D6 / D7**: an explicit default on a run
    that never reads the knob is refused (D6 `--cluster-wave vmc`), accepted
    and recorded `refused` (D5 `--rc-tail-model polrad-full`), or accepted
    and labelled (D7 `--pol-sf nnpdfpol`) — stated, not unified.
    `docs/CONVENTIONS.md` "Refuse or label"; `phase_V_verification.md` C4.

Also open, smaller: `phase_I4_jm_convention.md` §8 items 1–4, 6, 9 (the sign
of Δ, the exact JM coefficient, the spin-3/2 slot, α_s/8π vs α_s/16π, the
coherent `amp`, Δ on the finite-γ path); `phase_I5_polradfull.md` §9 items 1,
3, 4, 7 (no tail-table cache, platforms, T6's CI margin, the unreproduced
5.18398e−04); `phase_I6_cluster_wave.md` §9 items 6, 7, 9, 10;
`phase_I7_backends.md` §13 items 4, 5, 8, §14 (the registry texts, proposed,
not applied).
