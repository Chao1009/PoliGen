<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-17 — status board

The three items run 2026-11-10 left for the author, answered on 2026-11-17
("1. go with your recommendation, is it physics-based? 2. OK. 3. sounds
good. Any gaps need to be closed? If so, continue working on them"): K1 the
author's decision; its details and K2–K4 the supervisor's readings and
choices under those words, each with its reason in `DECISIONS.md`;
implemented by Opus workers under the supervising session. Branch
`claude/run-2026-10-26`, on top of `f104682`.

| phase | items | state | commit | what | record |
|---|---|---|---|---|---|
| A — the decisions recorded | K1–K4 | DONE | `a1f2b4d` | the author's words quoted; the answer to "is it physics-based?"; four entries, each with its reason, what moves and what does not | `DECISIONS.md` |
| P1 — the b₁ ratio reading | K1 | DONE | (this commit) | `Options::b1_reference_q2` (class level 0 = frozen), run default 2.5 and `--b1-ref-q2`; b₁/F₁ Q²-independent to 2.2e−16; `--b1-ref-q2 0` byte-identical to `a1f2b4d` on 48 runs; d confinement 3.87e−7 → 1.66e−9 (cfg 1), ⁶Li 3 → 0 cells outside the domain; HERMES χ²(b₁) 5.26 → 6.24, FAIL → FAIL; report 3 changes (HERMES only) | `phase_P1_b1_ratio.md` |
| P2 — what K1 leaves to re-take | K1 notes, K4 | TODO | | the two rulings of 2026-11-19 (K1 inside `li6-convolution`; the shipped-kernel helper follows the run surface); the documents' rank-2 numbers; the throughput; K2 is closed by P1's measurement | `phase_P2_after_k1.md` |
| V — verification | all | TODO | | two adversarial lenses over the implemented tree | `phase_V_verification.md` |
| F — fixes | all | TODO | | every finding applied or refuted with the reason | `phase_V_verification.md` |
| S — summary | | TODO | | the one page | `SUMMARY.md` |
