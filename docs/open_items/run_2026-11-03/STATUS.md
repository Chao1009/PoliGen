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
| P1 — the inclusive polarization factors | A1 | DONE | (next commit) | `LI6_CLUSTER_POLARIZATION` 0.811228 → 0.848 (= `LI6_POLARIZATION_VMC_SIX_BODY`), `LI6_B1_RANK2_TRANSFER` 0.921947 → 0.982576 (VMC α–d, pinned at 5e-7), `P_D_DEUTERON` 0.045 → 0.0576 (AV18, pinned at 5e-7); the pre-A1 values named (`_HULTHEN`, `_SCENARIO`) and reached through `li6_ion` / `deuteron_ion` / `default_inclusive_kernel(…, b1_transfer)`, which the port gates use (no reference JSON touched); the Hulthén tag bit for bit; meta keys and a banner line; A/B: tagged / coherent / ⁷Li and the explicit-old-value runs column-identical, inclusive vector ×1.04533, rank-2 ×1.06576, deuteron vector ×0.97973; report 14 / 9 unmoved (one sub-row FAIL → PASS); finding: the default inclusive ⁶Li helicity-flip is now refused at the top cell | `phase_P1_polarization_factors.md` |
| P2 — the positivity bound | A3, A4 | TODO | | g₁ clipped to ±F₁ per nucleon, the clipped share in `meta` and on the banner; the two refused plans re-measured at the default | `phase_P2_positivity.md` |
| P3 — the report's rules | A5, A6, A7 | TODO | | the headline convention written once, Deeps re-pointed; V6 out of TIGHT; the grid row's bound | `phase_P3_report_rules.md` |
| P4 — infrastructure and provenance | A9, A10, A11, A12, A17 (+ registry rows for A2, A8, A14) | TODO | | CMake `FATAL_ERROR`; HepMC3 `GenRunInfo` attributes with the T6 re-pin; the σ bound; the FSI banner line; the one not-read rule | `phase_P4_infrastructure.md` |
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
