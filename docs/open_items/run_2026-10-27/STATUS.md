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
| R1 — which Δ normalisation the literature uses (decision 1's question) | RUNNING | — | two readers (1989–2016 + the LOI; 2019–2026 + the tree and the letter) and a judge | `phase_R1_delta_convention.md` |
| D5 design — `PolradFull` default | RUNNING | — | inventory of every pin, label and document; adversarial check | `plan_D5_polradfull_default.md` |
| D6 design — VMC cluster-wave default | RUNNING | — | same; the `LI6_CLUSTER_POLARIZATION` question surfaced for the author | `plan_D6_cluster_wave_default.md` |
| D7 design — grid backend defaults | RUNNING | — | same; CI / missing-set failure mode designed | `plan_D7_backend_defaults.md` |
| D4 design — the Jaffe–Manohar convention in the kernel | TODO (after R1) | — | inventory of every site the ½ touches: the a₂ term, the finite-γ path, the bindings, the toy's amplitude, the bound row's S2/A7 logic, the port gate, every cos 2φ pin, the documents and their attribution; adversarial check | `plan_D4_jm_convention.md` |
| I4 — the Jaffe–Manohar convention implemented | DONE, verified and fixed (F4) | (this commit) | per `plan_D4` amended by its adversarial check: `DELTA_COS2PHI_NORM` = 0.5 at the two evaluation sites (every `a_cos2phi` and every Δ line of a₂ exactly half, the whole a₂ on the massless path; everything else bit-identical); T1–T5 and Y1–Y9 re-taken with dated comments, the missed `ds16_over_tree` pin too, N1 and N2 added; the port gate, C++ and Python, reads polligen's coefficient off `xsec.json` against a literal ½ and takes a Jaffe–Manohar re-dump unchanged, at the pre-D4 tolerances; `meta["delta_convention"]` on `Pipeline.generate()` and so every CLI .npz (HepMC3 carries none: an open point) and the run banner, the label tied to the constant by a `static_assert`; the convention paragraph with the one-line rule and the corrected attribution (HJM has no Δ term); row 13 PASS, relabelled, V6 kept in TIGHT; the report regenerated (row 13 only). Verification: 2 medium and 6 low findings, none refuted, all fixed in F4 (the record's §10); six drifted citations in the plan records named for the supervisor; every gate green after F4 | `phase_I4_jm_convention.md` |
| I5 — `PolradFull` default implemented | TODO | — | per `plan_D5`; gates; report regenerated | `phase_I5_polradfull.md` |
| I6 — VMC cluster waves default implemented | TODO | — | per `plan_D6`; gates; report regenerated | `phase_I6_cluster_wave.md` |
| I7 — grid backends default implemented | TODO | — | per `plan_D7`; gates; report regenerated | `phase_I7_backends.md` |
| V — adversarial verification of the implemented tree | TODO | — | four lenses: physics, honesty, code, regression | this file, "Verification" |
| F — fixes | TODO | — | | |

Ground rules of the run: the author's four decisions are the only
defaults that move; every number that depends on them is re-taken with a
dated comment naming the decision, or the test passes the old option
explicitly when its purpose is to pin the old physics; tolerances declared
before verdicts; a disagreement recorded, never tuned; the hashed design
block of `t4_djangoh_rad_noRad.py` is never edited; nothing published or
sent; the sibling repository read-only; fetched papers stay outside the
repository.
