<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-10 — status board

The nine items run 2026-11-03 left for the author, decided on 2026-11-10
("1. go with your recommendation. 2. re-dump, for the others go with your
recommendation"): J1 and J2 the author's decisions, J3–J9 the supervisor's
choices under the author's delegation, each with its reason in
`DECISIONS.md`; implemented by Opus workers under the supervising session.
Branch `claude/run-2026-10-26`, on top of `b3f8858`.

| phase | items | state | commit | what | record |
|---|---|---|---|---|---|
| A — the decisions recorded | J1–J9 | DONE | `dd7f4d4` | the author's words quoted; nine entries, each with its reason, what moves and what does not | `DECISIONS.md` |
| P1 — the re-dump | J2 | DONE | `1c98f62` | one dump run in place: 3170 numbers in six JSONs + the stale manifest string, a second run byte-identical; the 9.4e−13 traced to a one-ULP arccos change of the magic-angle input; every gate green at its present tolerance (435/435, 1523/155, bench 368/1), xsec headroom 0.056 → 0.0022 of the 1e−12 | `phase_P1_redump.md` |
| P2 — the deuteron's b₁ slot | J1 | DONE | (next commit) | `MillerB1` unmapped on the inclusive deuteron (b₁, b₂ ×3.0531989383; ⁶Li, ⁷Li, tagged, coherent byte-identical; HERMES: the kernel's own b₁ = the headline, χ²(b₁)/6 14.56 → 5.26); exposure measured: all 24 default deuteron plan × config runs refused at x_max 1.0 by the rank-2 term on the x = 0.955 row (1.6e−9 of the rate; b₁/F₁ up to 6.99 against the pure-state bound 3), 0.95 builds all but helicity-flip c0 (0.87); treatment NOT implemented — recommendation (b), a visible bound −3/2 ≤ b₁/F₁ ≤ 3, for the supervisor's note | `phase_P2_deuteron_b1.md` |
| P2b — the deuteron slot confined to the spin-1 positivity domain | J1, the supervisor's note of 2026-11-12 | TODO | | b₁ ← min(max(b₁, −(3/2)F₁), 3(F₁ − \|g₁\|)) on the inclusive deuteron, visible in `meta` and on the banner; the refusal stays as backstop; ⁶Li untouched | `phase_P2b_spin1_domain.md` |
| P3 — harness rows | J3, J4, J7 | TODO | | the BONuS robustness pin; the HERMES in-grid sub-row; one exit code for an unusable verdict set | `phase_P3_harness_rows.md` |
| P4 — provenance | J5, J6 | TODO | | the clip share's bracket; a HepMC3 file states its provenance or its absence; the extrapolation keys in the one list | `phase_P4_provenance.md` |
| P5 — records and the throughput | J8, J9 | TODO | | the development item brought up to date; the 2 M-event recipe re-measured on the idle machine | `phase_P5_throughput.md` |
| V — verification | all | TODO | | two adversarial lenses over the implemented tree | `phase_V_verification.md` |
| F — fixes | all | TODO | | every finding applied or refuted with the reason | `phase_V_verification.md` |
| S — summary | | TODO | | the one page | `SUMMARY.md` |
