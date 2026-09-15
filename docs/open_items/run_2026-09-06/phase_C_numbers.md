# Phase C numbers — run 2026-09-06, task C1

**The decay-lepton reconstruction efficiency the O5 chain is missing.**
Performed 2026-09-15. The survey this rests on is
`phase_C_survey.md` (read-only survey of `../PolarizedLithiumSim`); this file
is the record of what was *done* in LiPolGen as a result.

`PLAN.md` C1 set the fork explicitly: *"Derive a BOUND from them — stated as a
bound, with its source — and re-run `validation/o5_a2_reach.py` with it as a
fourth ladder row. **If the sibling has nothing usable, say so and stop.**"*

**The sibling has nothing usable. This phase took the second branch.**

---

## C1. The one sentence, with what was looked at

> **Every module of `../PolarizedLithiumSim/fastsim/polli_fastsim/` (11 files,
> 3725 lines), all of `tools/fullsim/` and `tools/analysis/`, `plans/03` and
> `plans/09`, and all 54 entries of `refs/refs_dict.json` plus two reference
> PDFs read directly with `pdftotext`, were surveyed on 2026-09-15, and the
> whole sibling tree was grepped for the five shapes the brief named (a
> per-track efficiency, a J/ψ → ℓℓ reconstruction efficiency, an `eic-smear`
> parameterisation, a fullsim output with a reconstructed J/ψ, an ePIC
> Yellow-Report number the project already assumes).
> `PolarizedLithiumSim` contains **no measurement, no parameterisation and no
> assumed value of a central-detector lepton reconstruction efficiency at
> ePIC.** It contains one per-track efficiency,
> `HadronResponse.eff_track` = 0.95 at `evgen/polligen/hfs.py:251`, which that
> file labels **"stand-in" on the same line that states it**; one *constructed*
> scattered-electron ID profile, `eps_eid` at `evgen/polligen/reco.py:482`,
> whose own docstring says *"No ePIC electron-ID efficiency curve exists in any
> ePIC document"*; and no J/ψ reconstruction, no `eic-smear` parameterisation
> and no EICrecon output of any kind — `eicrecon` is invoked in zero scripts.
> **So the per-lepton reconstruction efficiency in the O5 chain remains
> UNBOUNDED rather than merely unmeasured, no constant was added, no fourth
> ladder row exists, and no number in this tree moved.**

### C1.1 Re-verified here, not taken on the survey's word

Every load-bearing fact above was re-checked in this phase against the sibling
working tree (read-only; nothing in `PolarizedLithiumSim` was modified):

| claim | check | result |
|---|---|---|
| the 0.95 and its label are on one line | `sed -n '251p' evgen/polligen/hfs.py` | `      eff_track      plateau tracking efficiency (0.95, stand-in)` |
| it is the shipped default | `sed -n '293p' evgen/polligen/hfs.py` | `                 eff_track=0.95, eta_cal=3.7, e_min_photon=0.1,` |
| the only lepton-side efficiency is constructed | `sed -n '482,487p' evgen/polligen/reco.py` | `def eps_eid(eta)` … *"a CONSTRUCTED eta profile, not a published curve"*, *"No ePIC electron-ID efficiency curve exists in any ePIC document."* |
| the reference corpus size | `len(refs_dict.json["entries"])` | **54** |
| no reconstruction was ever run | `grep -rl eicrecon --include=*.sh --include=*.py` | **0 files**; `find -name '*.root' -o -name '*.hepmc*'` → **0** |
| the acceptance figures quoted below | `o5_a2_reach.lepton_pair_acceptance` at ⟨W⟩ = 30.2 | 0.4477 (\|η\| < 1), 0.8888 (< 2), **0.9940** (< 3.5) |

## C2. Why each candidate is not a bound

Ten candidates were classified in `phase_C_survey.md` §C-S3. Eight are
placeholders or measurements of a different quantity. The two that are neither
still do not bound ε_pair, and the reasons are the deliverable:

**(a) `HadronResponse.eff_track` = 0.95 — a stand-in with no direction.**
It is the *plateau* of a logistic p_T turn-on (centre 0.2 GeV, width 0.05 GeV)
gating any charged particle inside \|η\| ≤ 3.5 under a **pion mass hypothesis**
(`hfs.py:360-362`, `:372`), applied to **hadronic-final-state particles**, not
to identified leptons. It carries no PID, no material budget, no magnetic
field and no purity term, and it is uniform in η across the whole coverage.
Its source is *stated as none*. **A bound needs a direction and this has
neither an upper nor a lower one**: nothing says a real ePIC per-lepton
efficiency lies above it or below it. It stays what `o5_a2_reach.py` already
called it — a stand-in used to *show what a plausible one costs* — and this
phase did not promote it.

**(b) ATHENA Table 5, `refs/2210.09048.pdf` — a real bound on the wrong
quantity.** *"> 99.8% pion rejection with 95% electron efficiency at
p ≥ 0.1 GeV/c"*, footnote *b* *"Based on simulation for a standalone bECal"*,
Fig. 9 caption *"no other materials are placed in front of the calorimeter and
no magnetic field is involved"*. This **is** an upper bound — on the
calorimetric electron-**ID** efficiency of a barrel EMCal of the ATHENA bECal
design, under the assumption that material in front and a 1.7 T field can only
lower a calorimetric ID efficiency. Five things stand between it and ε_pair,
and the third is fatal for this sample:

1. it is **ATHENA, not ePIC** — a proposal for IP6 that was not built;
2. it is a **standalone calorimeter**, by the paper's own caption;
3. it is the **barrel only** — Table 5 prints no η range at all, and at the O5
   sample's own ⟨W⟩ = 30.2 GeV the J/ψ sits at y = −0.39, where only
   **0.4477** of pairs have *both* leptons inside \|η\| < 1 (0.8888 inside 2,
   0.9940 inside 3.5). A barrel number covers **45 % of the pairs** and is
   silent on the other 55 %;
4. it is an **ID efficiency given a cluster** — no tracking, no vertexing, no
   pair-finding;
5. it is an **electron** number. The paper offers one qualitative sentence on
   muons and no number, and nothing anywhere in the sibling quantifies μ⁺μ⁻.

**Extending it to the endcaps would mean inventing the endcap values**, which
is exactly what `eps_eid` already does and labels as invention. It does not
bound ε_pair and it is not quoted as one anywhere in this tree.

### C2.1 The direction below the band's low rung is now better known

Not a number, and recorded as a direction only. `evgen/polligen/reco.py:496-498`
attributes to ECCE (NIM A 1055 (2023) 168464 §3.5.2) EEMC ~95 %, **BEMC ~70 %**,
FEMC ~90–95 %. `refs/2207.09437` is **not** in `refs/` and there is no ECCE
entry among the 54, so it **cannot be verified here** and is not adopted. Its
direction is nonetheless worth having: *if* any proposal-era electron-ID number
were ever adopted as the per-lepton stand-in, the one that applies where 45 % of
these pairs land is **0.70, not 0.95** — so the 0.95/track currently used to
draw the band's low end is **not even the pessimistic choice among the numbers
the sibling itself names**. The openness below the band's low rung is therefore
*larger* than the current text implies, not smaller. That is a direction, not a
number, and it changes nothing in code.

## C3. What changed, and where

Nothing numeric. Each site gained **one sentence** recording the survey and
naming what was looked at; the existing claim ("no decay-lepton reconstruction
efficiency exists in this tree at all") is **confirmed**, not replaced, and the
band is untouched at both ends.

| site (the ground rules' enumeration) | where the sentence landed |
|---|---|
| `validation/o5_a2_reach.py` | the module docstring's assumption list (a new bullet after the decay-lepton one); `lepton_pair_acceptance`'s docstring; the `EPS_TRACK_SIBLING_STANDIN` comment block; the printed §3b(c); the printed "WHAT IS STILL NOT ESTABLISHED" list |
| `docs/open_items/run_2026-09-03/phase_C_numbers.md` §C2 | §C2.3c(c) (the derivation), §C2.7 (the band bullet), §C2.8 item 3b (the list that claims completeness) |
| `docs/OPEN_ITEMS_SOLUTIONS.md` §11.3 | §11.3b where the claim lives, and §11.3b's "say what is missing" list; plus the §11 registry row at the top of the file |
| `include/lipolgen/cluster_config.hpp` | the `a2_from_geometry` docstring, inside the four-things list |
| `python/lipolgen/configs.py` | the module docstring's four-things list |
| `docs/open_items/run_2026-09-03/mantysaari_collaboration_draft.md` | §4 item 4 (the draft's own live statement) **and** the quoted letter body, so the letter does not claim less than the tree knows |
| `docs/open_items/run_2026-09-03/SUMMARY.md` | "What the run did NOT establish" item 1 |
| `docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md` §B25 | (c), the O5 verdict paragraph. **(d) "No number in the tree moves under any option" still holds.** |

**The count stays FOUR.** The survey confirms the second of the four
unestablished things; it does not add a fifth and does not remove one.

### C3.1 The one consequential edit, and it is bookkeeping

Adding 11 lines to `cluster_config.hpp` and 7 to `configs.py` shifted 12 point
citations in `docs/PHYSICS_CHANNELS.md` (which is under
`check_physics_channels_links.py`). The gate caught all 12 — it was **1236 /
96 / 7 / 12 broken, rc 1** — and `--fix` relocated each onto the declaration it
already named. The resulting diff of `PHYSICS_CHANNELS.md` is **two table rows
and nothing but line numbers**: `cluster_config.hpp` 679→690, 715→726;
`configs.py` 124→131, 133→140, 152→159, 156→163, 163→170, 174→181, 178→185,
179→186, 332→339 (twice). No recorded range's *content* changed, so
`--record-ranges` was **not** run and nothing was re-blessed.

## C4. What this phase deliberately did NOT do

* **No constant was added to `include/lipolgen/coherent.hpp`.** There is
  nothing to name. A constant beside `COHERENT_JPSI_EFF_IR8_LI7` whose
  provenance read "the sibling's own stand-in" would be a placeholder promoted
  to an interface, which is the defect this run's earlier passes kept removing.
* **No fourth ladder row.** The corrected-ε_det ladder in the printed report is
  unchanged: `as shipped` → `+ beam-energy leg` → `+ decay-lepton GEOMETRY` →
  `BAND LOW` → `BAND TOP` ×2 → the three de-squeezed-optics rungs.
* **The self-check was not re-keyed.** There is no row to pin the absence of.
  `_self_check`'s existing lepton assertions are unchanged and still pass —
  `assert 0.0 < ch["a_geom"] <= 1.0`, `assert 0.98 < ch["a_geom"] < 1.0`
  ("geometry is not where the factor is"), and
  `assert 0.0 < ch["down_standin"] < ch["a_geom"]` (the pair stand-in can only
  cost) — and its docstring's "the lepton factor must NOT exceed 1" still
  describes them.
* **The headline band is untouched**, and was re-measured to confirm it:

| rung | before | after |
|---|---|---|
| headline (as shipped) | S = 2.618, 3 σ at 13.1 fb⁻¹/u | **identical** |
| BAND LOW | S = 2.627, 13.0 | **identical** |
| BAND TOP, species lo / hi | S = 2.842 / 3.292, 11.1 / 8.3 | **identical** |
| de-squeezed optics, low / top-lo / top-hi | 0.734 / 0.794 / 0.920 at 167.1 / 142.7 / 106.4 | **identical** |
| `_self_check` | OK | **OK** |

## C5. Two findings the survey surfaced that are NOT this phase's to take

Both are recorded here and **neither was applied**, because `PLAN.md` C1's
second branch is "say so and stop" and neither is the quantity C1 asked for.
They are put in front of the author as decisions.

### C5.1 `HfsModel` names no class in the sibling (survey §C-S5)

`validation/o5_a2_reach.py:437` and `:466` cite
`HfsModel(pt_min_track=0.2)` / `HfsModel(eff_track=0.95)`;
`docs/OPEN_ITEMS_SOLUTIONS.md` §11.3b and
`run_2026-09-03/phase_C_numbers.md` §C2.3c(c) / §C2.8 repeat it.

**Re-verified in this phase:** `grep -rn HfsModel --include=*.py` over the whole
sibling returns **0 hits**; `HadronResponse` returns **26**. The correct
citation is `polligen.hfs.HadronResponse`, `evgen/polligen/hfs.py:251` for the
label and `:293` for the default.

**The substance of the claim is exactly right** — that file does label the 0.95
a "stand-in" in as many words, on the line that states it, and this phase read
that line. The defect is that the citation is **unresolvable as written**: a
reader who greps the sibling for `HfsModel` finds nothing and cannot confirm the
label. It is a naming defect in a citation, not a factual one, and fixing it
moves no number. **Left untouched**, deliberately: C1's branch says change
nothing else, and the sentences this phase added cite `evgen/polligen/hfs.py:251`
by file and line, which does resolve.

### C5.2 The 17.75 % carries no detector efficiency either (survey §C-S6)

This is a **bound**, it is real, and it is on the **recoil** leg, not the
leptons — so it is not what C1 asked for and it was not applied. Re-verified in
this phase by `pdftotext` on `refs/2511.05638.pdf`, §IV last sentence
(p. 3→4), verbatim:

> *"The current simulation only accounts for the acceptance effect and does not
> incorporate the efficiencies of the detector. Additionally, we did not account
> for the efficiency and acceptance of the reconstructed distribution, which is
> to showcase our detection performance at its detector level form."*

> **BOUND.** `COHERENT_JPSI_EFF_IR8_LI7` = 0.1775 is an **UPPER bound** on the
> far-forward intact-⁷Li tagging efficiency × acceptance at the IR-8 secondary
> focus.
> **Source:** W. Chang, E.-C. Aschenauer, A. Jentsch, A. Kumar, Z. Tu, Z. Yin,
> PRD 113 (2026) 032018 (`refs/2511.05638.pdf`), §IV last sentence for the
> exclusion, §V for the number.
> **The assumption that makes it a bound:** that a Roman-Pot detector
> efficiency, a reconstruction efficiency and a reconstruction acceptance are
> each ≤ 1, so restoring any of them can only **lower** the product. Nothing
> quantifies how far — the paper does not, and this tree does not.

**Why this matters and why it is an author decision.** 0.1775 is the single
largest multiplier in the O5 chain. The tree currently describes it as an
"efficiency"; by the paper's own §IV it is a **pure geometric acceptance**
carrying no detector efficiency on the recoil leg either. **The direction is
DOWN, and it does not get absorbed the way the lepton factor did.** The chain's
whole two-omission balance — ×1.1224 UP (beam energy) against ×0.897 DOWN
(decay leptons) = 1.007, which is the stated reason neither may be quoted alone
— has room for exactly those two. This bound is a *third* factor, on the same
side as the leptons: it **compounds with the decay-lepton omission and breaks
that cancellation** rather than being cancelled by the beam-energy leg.
Adopting it would touch every site the ground rules enumerate and would change
how the band's assumption list reads, which is exactly why this phase did not
take it unasked. `phase_C_survey.md` §C-S6 carries the full statement.

*(This sharpens rather than contradicts the existing assumption list, which
already says the chain has no decay-lepton acceptance or efficiency. True — and
the recoil leg has no detector efficiency either. Both omissions point down.)*

## C6. Suites, gates, tallies

Baseline measured on the clean tree before any edit; "after" measured on the
rebuilt tree.

| | baseline | after |
|---|---|---|
| `build/lipolgen_tests` | 408 cases / 408 passed / 0 failed / 0 skipped, **17 241 975** assertions, rc 0 | 408 / 408 / 0 / 0, **17 241 975**, rc 0 |
| `python -m pytest python/tests -q` | **1023 passed, 151 skipped**, 296.6 s, rc 0 | **1023 passed, 151 skipped**, rc 0 |
| `python3 validation/check_spdx_headers.py` | 102 files, all carry the header, rc 0 | 102 files, rc 0 |
| `python3 validation/check_physics_channels_links.py` | 1236 / 96 ranges / 7 external / **0 broken** strict; 19 / 115 / 0 / 0 `[SPIN32]`; 0 / 0 / 8 / 0 `[PYTHIA_BRIDGE]`; rc 0 | identical, rc 0 **after `--fix` relocated the 12 drifted anchors** (§C3.1) |
| `python3 validation/o5_a2_reach.py` | — | rc 0, `self-check: OK`, every rung bit-identical (§C4) |

**Nothing regressed. No number moved. The tree gained one sentence at each of
the eight sites the ground rules name, and twelve line numbers in
`PHYSICS_CHANNELS.md` followed the lines they already pointed at.**
