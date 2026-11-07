<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-03 — status board

The 17 open items of run 2026-10-27, decided under the author's delegation
("Go with your recommendations for the 17 items, just record the choice and
why"), each choice and its reason in `DECISIONS.md` (A1–A17, the
supervisor's choices under that delegation), implemented by Opus workers
under the supervising session. Branch `claude/run-2026-10-26`, on top of
`70982c9`.

| phase | items | state | commit | what | record |
|---|---|---|---|---|---|
| A — the choices recorded | A1–A17 | DONE | `8eea2c7` | the author's words quoted; 17 entries, each with its reason, what moves and what does not; records-only entries A2, A8, A14, A15 are closed by `DECISIONS.md` itself (their registry rows in P4) | `DECISIONS.md` |
| P1 — the inclusive polarization factors | A1 | DONE | `3d059d3` | `LI6_CLUSTER_POLARIZATION` 0.811228 → 0.848 (= `LI6_POLARIZATION_VMC_SIX_BODY`), `LI6_B1_RANK2_TRANSFER` 0.921947 → 0.982576 (VMC α–d, pinned at 5e-7), `P_D_DEUTERON` 0.045 → 0.0576 (AV18, pinned at 5e-7); the pre-A1 values named (`_HULTHEN`, `_SCENARIO`) and reached through `li6_ion` / `deuteron_ion` / `default_inclusive_kernel(…, b1_transfer)`, which the port gates use (no reference JSON touched); the Hulthén tag bit for bit; meta keys and a banner line; A/B: tagged / coherent / ⁷Li and the explicit-old-value runs column-identical, inclusive vector ×1.04533, rank-2 ×1.06576, deuteron vector ×0.97973; report 14 / 9 unmoved (one sub-row FAIL → PASS); finding: the default inclusive ⁶Li helicity-flip is now refused at the top cell | `phase_P1_polarization_factors.md` |
| P2 — the positivity bound | A3, A4 | DONE | `df56035` | one `PositivityClippedG1` in every `InclusiveKernel`: g₁ᵖ, g₁ⁿ clipped to ±F₁ on the kernel's own F₁, so `g1_nucleus` and g₂^WW see it; `meta["g1_clipped_rate_fraction"]` (2–5e−7 of the rate on the default helicity runs, NaN where g₁ is not read) and a banner line; toys bit-identical (13 runs, A₁ dump); of the seven top-cell refusals five build, the inclusive d configs 0 / 1 stay (rank-2 on top of the clipped vector); A4: `tagged-d-p` helicity-flip builds, `--isotope d --plan helicity-flip` stays refused; finding: the default ⁶Li helicity-flip / apar builds again, so P1's x_max 0.95 test helpers are removed (the supervisor's correction of 2026-11-06 under A1; no pin moved); report 14 / 9, 0 changes | `phase_P2_positivity.md` |
| P3 — the report's rules | A5, A6, A7 | DONE | `5cd19d2` | convention 8 of `validation/benchmarks/README.md` (the headline on the shipped default; an opt-in object ON in its shipped form; the one exception row 16), line-neutral; headlines re-pointed under A5 and the supervisor's note of 2026-11-06: Deeps AV18 579.12/187 (Hulthén 627.85 recorded) FAIL → FAIL; BONuS AV18 + Glauber 40 mb 3/12, 0/16 (Hulthén 4/12, 0/16 recorded) FAIL → FAIL; deuteron statics AV18 η −1.39σ, **FAIL → PASS by the rule** (the all-waves rule FAIL recorded); HERMES b₁ `miller_shipped` 5.26/6, A_zz 7.31/6 on the shipped F₁ᵈ (cdks_conv_mstw 22.26 recorded) FAIL → FAIL; Boeglin unchanged (AV18 `nk[av18]` 11.54/21 recorded); A6: the edge carried V6 as `PRE_D4_OVER_JM = 2`, every TIGHT saturation halves (toy 0.321/0.381/0.428 → 0.160/0.190/0.214), the ×2 edge recorded, no status moves; A7: grid-row error ≤ `QUAD_REL_TOL` × \|M₂\| (271–448× inside), chosen after the values; report 14 / 9 → 15 / 8, 6 compare changes; findings: the deuteron inclusive kernel's b₁ = 0.3275 × MillerB1, BONuS's CHECKS miss on the AV18 weight (recorded) | `phase_P3_report_rules.md` |
| P4 — infrastructure and provenance | A9, A10, A11, A12, A17 (+ registry rows for A2, A8, A14) | DONE | (next commit) | A17: one rule — an explicit default on an unread knob accepted, row `not-read`, label in `meta` (`cluster_wave vmc` off the tags, `rc_tail_model` and every band/tail knob at its default where band/tail does not apply); label-only, the accepted defaults bit-identical to the unnamed runs; the rows it does not yet reach named in CONVENTIONS; A11: \|σ\| ≤ 3 refused beyond; A12: the clipped cell fraction on every FSI banner (⁷Li 4.41/9.52 %, ⁶Li 0.70/5.25 %, d 1.35 %), registry row 28; A10: `run_provenance` the one list, `lipolgen.<key>` GenRunInfo attributes = the npz meta, HepMC3 files differ from pre only by those lines, T6 C24 + sha256 re-pinned; A9: configure-time FATAL_ERROR, both out-of-tree builds measured; A2/A8/A14 dated notes | `phase_P4_infrastructure.md` |
| P5 — the VMC port-gate block | A13 | TODO | | `li6_alpha_vmc` dumped additively, gates at 1e−8 | `phase_P5_port_gate.md` |
| P6 — the throughput | A16 | TODO | | idle-machine re-measurement with its log | `phase_P6_throughput.md` |
| V — verification | all | TODO | | two adversarial lenses over the implemented tree (physics and honesty; code and regression) | `phase_V_verification.md` |
| F — fixes | all | TODO | | every finding applied or refuted with the reason | `phase_V_verification.md` |
| S — summary | | TODO | | the one page | `SUMMARY.md` |

Ground rules of the run: a default moves only under an entry of
`DECISIONS.md`; every number that depends on a moved default is re-taken
with a dated comment naming the entry, or the test passes the old value
explicitly where its purpose is to pin the old physics; tolerances declared
before verdicts; a disagreement recorded, never tuned; `validation/reference/*.json`
untouched except A13's additive block; `docs/HEPMC3_CONVENTION.md` edited by
A10 only, T6 re-pinned; the hashed design blocks never edited; the sibling
repository read-only; nothing published or sent; every phase committed and
pushed on its own.
