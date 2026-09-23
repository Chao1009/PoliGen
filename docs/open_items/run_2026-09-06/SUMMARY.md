# Run 2026-09-06 — one page

What the follow-on run **priced**, **closed** and **fixed**, phase by phase, in
the numbers a reader would cite — **each with the window it was measured in** —
and, at the end, what it did **not** establish. Every number here is quoted from
the phase record named beside it, or re-measured at the close-out and said to
be. `STATUS.md` is the board, `PLAN.md` the brief. **The decision registry is
still the 2026-09-03 run's**: `../run_2026-09-03/STATUS.md`'s table, stated in
full in `../run_2026-09-03/AUTHOR_DECISIONS.md`. Cite a decision by its ROW
NUMBER.

**The one sentence.** This run **priced decisions and took none**: the three
costs the 2026-09-03 run left in the word "unpriced" were each *implemented
opt-in at a no-op default and run*, an exact radiative tail and a sourced
Q(⁷Li) were added behind selectors, a bound that does not exist was reported as
the negative it is, the documentation gate was taught three rules it had been
missing — and **one shipped output moved, because it was wrong**: the tagged
sector's S–D interference phase.

**Start → end.** Baseline `5af0427`, clean: 401 doctest cases / 17 240 286
assertions / **1 skipped**, 926 pytest + 112 skipped, gate 1206 / 96 / 7 / 0.
End, **re-measured 2026-09-16 at `86cd1e6` plus the applied findings and this
close-out** (the brief forbids `git commit`, so phases A–D are the hashes on
the board and the rest is the working tree):

* `build/lipolgen_tests` **408 cases / 408 passed / 0 failed / 0 skipped**,
  **17 241 975 assertions**. The suite's one **unconditional** skip — RC test
  T9, the `PolradFull` Rosenbluth limit — is **gone** since `af3f415`, which
  is why the tally reads 0 rather than 1. **The remaining skips are still
  conditional and still real**: a case that is nothing but an optional tier or
  an optional data tree is skipped at REGISTRATION time, so this 0 also says
  that **both optional data trees were present on this machine** — with the
  MSTW grid or `data/vmc` absent the same tree would report a nonzero count,
  and registry row 25 is open on the three in-case guards that would not.
* `python -m pytest python/tests -q` **1085 passed / 151 skipped / 0 failed**
  (264 s).
* `validation/check_physics_channels_links.py` strict: **1240 refs / 95 ranges
  / 7 external / 0 broken / 24 allow-listed** on `PHYSICS_CHANNELS.md`,
  **19 / 115 / 0 / 0 broken / 10 allow-listed** on
  `theory/SPIN32_FINITE_GAMMA.md`, **0 / 0 / 8 external / 0 broken** on
  `PYTHIA_BRIDGE.md`. `--audit-ranges` (R5): **209 range citations on 179
  blocks, 0 REFUSED**. `--records-only` (R6): **747 citations in 43 dated run
  records, 0 broken** — 43 because this page is one of them. `check_spdx_headers.py` **103/103**.
* `--fix` and `--record-ranges` were **run, not assumed**, and both are
  no-ops: **0 fixed** in all three strict documents, and the range sidecar is
  **byte-identical** across `--record-ranges`
  (md5 `132bbaeb1e042c06ee1c22dd619e9b8b` before and after).
* **`validation/reference/`, by sha256 against the project baseline
  `a94fd6e`:** `b1_default_li6`, `beams`, `bookkeeping`, `coherent`,
  `spectator`, `spin`, `xsec` **identical**; `tagged.json` and
  `_manifest.json` **differ, by registry row 26 and nothing else**. The
  brief's "byte-identical to `a94fd6e`" must be read that way for this run,
  and §B27(d) says so in the registry. (2026-09-23: `_manifest.json` is byte-identical to `a94fd6e` again; `../run_2026-09-23/phase_A_port_gate.md` §4. The same day `b1_default_li6.json` gained a `provenance` label — its whole-file sha256 moved `d7bd8ce6…` → `8373db8e…`, every numeric block byte-identical; `../run_2026-09-23/phase_B3_chain_rc.md` §4.)

---

## A — the three unpriced costs, priced; the ⁶Li tables regenerated as a band (`cdd8591`, `phase_A_numbers.md`, `phase_A_li6_tables.md`)

**Row 3 — the third R option, `--r-source`** (§A1). The option the registry
called *"the better third option, not taken here"* is now a real hook:
`PipelineConfig::r_source` threads ONE R into `Li6ConvolutionOptions::r_func`
**and** `InclusiveKernel::Options::r_func`. `unset` is the default and does not
enter the branch — bit for bit **by construction**; `sigma-lt` is measured
bit-identical (162 values over x = 0.05/0.10/0.30 × Q² = 2.5/5 ×
y = 0.1/0.5/0.9, plus σ_pb, the three per-category cross sections and every
column of a 2000-event run). Under `r1998`, at **y = 0.5, Q² = 2.5**, the
tensor weight K/D_φ — and A_zz, which is exactly −(2/3) of it at θ_m = 0 —
moves **−3.8357 % / +26.3988 % / +3.3518 %** at x = 0.05 / 0.10 / 0.30,
against option (ii)'s **−5.5148 % / +24.6099 % / +2.6629 %**. Four things the
numbers say that an estimate would have missed: sharing the hook does **not**
shrink the shift (it is *larger* at x = 0.10); the **cos 2φ amplitude separates
the two options exactly** (Δ carries no R, so (ii) leaves it bit-identical and
(iii) moves it purely through D_φ, −8.33 % / −6.73 % / −3.14 %); **(iii)'s
shift is y-dependent and (ii)'s is not** (at x = 0.05, Q² = 2.5 it runs
−5.4705 % at y = 0.1 → +2.3678 % at y = 0.9, **straddling zero** — so a (iii)
number quoted without its y is not a statement about a run); and **(iii) moves
the unpolarised rate**, σ_pb **−0.684681 %** (591846.161405217 →
587793.902124821 pb, ⁶Li inclusive, config 1, x_max 0.95, `tensor-thirds`
(0.7, 0.6), 2000 events, seed 7), which (ii) cannot. **The row does not close
by choosing (iii)**: it becomes *(iii) with which R?*

**Row 17 — the RC band's low-x anchor** (§A2). `RC_DELTA_LOW_X` = 0.30 is the
source panel's value **carried upward in x**; **0.113 is what that panel reads
at x = 0.00966**, the x nearest `RC_X_LOW` = 0.01, so the shipped anchor errs
**wide by ×2.65** — the safe direction for a half-width and the wrong one for a
quoted precision. Band half-width on `A_zz` at δ_low = **0.30 / 0.266 /
0.113** (⁶Li, config 1, inclusive, `--plan tensor-thirds --pzz 0.6`, θ_S = 0,
Q² = 5 GeV², P_zz = +1, **this generator's own ⁶Li b₁**):
**1.358472e−04 / 1.204512e−04 / 5.116913e−05** at x = 0.010,
**1.459067e−04 / 1.308566e−04 / 6.313127e−05** at 0.063,
**9.680016e−05 / 8.798804e−05 / 4.833349e−05** at 0.100, and
**1.920691e−05 on all three** at 0.160. The anchor is carried in full only at
x ≤ 0.01 (×0.8867 / ×0.3767), **compressed** between the anchors and
**identically zero** at x ≥ 0.16 — so the registry's own "≈ ×0.89 / ≈ ×0.38"
is exact at the anchor and **overstates the reduction everywhere else**. The
band still leads the whole radiative tail at x = 0.01 by **×768 / ×681 /
×289**, so the RC budget's ordering survives all three. **The HIGH anchor has
no alternative to run**: the record names none.

**Row 20's D13 — `--pzz-mode`** (§A3). Measured at the standard configuration
(inclusive, beam config 1, 100 000 events, seed 20260713,
`--plan helicity-flip --pz 0.7 --pe 0.7`, `--b1-model miller`) with
`--pzz 0.5`, a value inside **both** domains. **On ⁷Li, the isotope the row is
about, the typed fill moves NO observable** — σ agrees to **1 ulp**
(1.970e−16 relative), A_∥ to **6.2e−14 σ_stat** — while still moving
**0.1060 %** of the sample (106 events) into different (x, Q²) cells, so the
file is not bit-identical. On ⁶Li σ moves **−0.0023527 %** and A_∥
**−4.27e−6 σ_stat**. **What moves on both is the recorded alignment**, the
divisor of every tensor estimator: T 0.4 → 0.5 and P_zz 0.409403 → 0.5, i.e.
δ(A_zz) and δ(cos 2φ) at fixed N move **−20 %** and **−18.1195 %**.
`--pzz 0.6` at `--pz 0.7` is **refused** at J = 3/2 with its edge named, never
clamped.

**A4 — the ⁶Li tables are regenerated, as a BAND** (`phase_A_li6_tables.md`).
Same observables, same (x, Q²) points, same object, over
`--b1-unpol {toy, ct18nlo, mstw}` × `--unpol-sf {toy, ct18nlo, mstw}`; the
`toy` column reproduces the 2026-09-02 tables to every printed digit, which is
what makes the other two a band rather than a different measurement. **The
3 × 3 matrix is 7 legal cells, not 9.** At x = 0.10, Q² = 2.5, y = 0.5, D_φ
spans a factor **1.12** (8.422582 … 9.450533) while **A_zz spans 2.23**
(+1.100752e−05 … +2.458518e−05), and that span is driven **×2.22 by the
numerator**, along `--b1-unpol`.

**Nothing in phase A moved a default**, and `--r-source`, `--pzz-mode` and
`--rc-sp-tensor-scale` took the knob-provenance table from 66 rows to **69 on
a default run** (printed by the CLI banner itself) and the matrix test to
**12 specs × 81 knob variants = 627 cells**, both re-counted at the close-out.

## CW — the tagged sector's S–D interference sign, unplanned (`a7b3d18`, `53f6951`, `phase_CW_numbers.md`)

**The one default that moved this run, and it moved because it was wrong.**
`TaggedModel::build_amp2` summed the partial waves without the i^L phase; a
partial-wave amplitude needs φ_L = i^L ψ_L, whose observable relative part is
`(-1)^floor(L/2)`. Three independent derivations
(`../../benchmarking/07_cw_sign_investigation.md`) agree and the run calls it
**certain**: the Cosyn–Weiss Eq. (6.12) identity now holds to
**8.881784e−16** over every cell of the 280 × 96 grid on both wave functions,
against **2.740499e+00** before, and CW TABLE II's three cell-centre rows
reproduce to **4.36e−06 / 2.77e−06 / 3.64e−07**.

**What moved:** A_zz^tag(k = 0.20 GeV), acceptance-weighted at the YR
high-acceptance optics, **+0.8450 → −1.2069** (Hulthén β = 0.30) and
**+0.4518 → −0.5191** (VMC AV18) — and the corrected VMC column is, digit for
digit, the `vmc-flipD` *control* column the 2026-09-03 run published. At the
CLI's own default fill the ⁶Li Hulthén tag-fraction prediction goes
**0.027030 → 0.020396** (−24.5 %), and the `tensor-thirds` categories
0.028127 → 0.018403 (azz±) and 0.017775 → 0.037222 (azz0).

**What did NOT move, measured:** σ_pb to ≤ 2e−16 relative; every inclusive,
coherent and ⁷Li tagged cell byte-identical on two seeds; the spin-blind ⁶Li
α-tag rate to **≤ 4.0e−16** relative (1–2 ulp, and the VMC channel
bit-identical); P_D, a norm; the four dilutions to **≤ 1.4e−6 absolute**,
because they are angle-integrated and ∫Θ₀Θ₂ dc = 0 kills the cross term;
`LI6_CLUSTER_POLARIZATION`, a closed form in P_D; and the **+2.069 %** embedded-
deuteron defect of registry row 11, which re-measures as 1.020687624664 against
1.020687500612.

**A tag fraction at a tensor-polarised fill is NOT invariant**, and three
documents had said it was: it is the same acceptance-weighted integral as
A_zz^tag. The uniform-M mix is unmoved to 15 digits; the equal-thirds average
moves **+3.1e−6**, a per-M norm-residual effect (Σ n_M k² differs by 5e−5
between M = 0 and ±1) — **not** a summation-order effect, which is what the
figure first published here said.

**The reference re-pin was deliberate and is registry row 26.**
`polligen`'s own `tagged._amp2_table` carries the same missing phase, so it
cannot be the source for those two blocks; `validation/reference/tagged.json`'s
`channels/deuteron/model` and `channels/li6_alpha/model` are dumped from **this**
library instead (8 fields each; the ⁷Li block came through **byte-identical**),
and `tagged.json` is no longer an external port gate for them. (Until 2026-09-23; see below.) The priced
alternative — keep `polligen`'s blocks under an `xfail` naming the bug — was priced before the sibling fixed its own phase (2026-09-15 19:50, commit `1066555`; the fixed polligen agrees with the re-pin at 2e-13, re-measured 2026-09-16), and the polligen port gate WAS restored on 2026-09-23 — every block of `tagged.json` re-dumped from `polligen` at the sibling's `c0f86a8`, the two re-pinned blocks moving by at most 1.33e−13 (li6_alpha) / 2.28e−13 (deuteron) relative, `../run_2026-09-23/phase_A_port_gate.md`; the `xfail` alternative was a
gate that is green only because it is told to expect failure.

## B — the s/p tensor bound, the exact tail, and Q(⁷Li) (`af3f415`, `phase_B_numbers.md`)

**B1 — the tensor fraction of the s-/p-peaks is bounded, not computed, and at
the standard points the bound is EMPTY.** `--rc-sp-tensor-scale` (default 0)
prices the **elastic** s/p column, and on ⁶Li `--config 1`,
s = 3980.0000000000005 GeV², the production grid (101 × 77 nodes, 3051 accepted
cells), `ho` C0 edge, it recovers **0.0 %** of the tensor-fraction collapse
×0.66139 / ×0.0031829 / ×6.6076e−05 at x = 0.01 / 0.10 / 0.30, Q² = 5 —
bit-identical there, because ⁶Li's coherent form factor is 45+ decades down at
the s/p vertex — and **at most 21.72 %** of it anywhere on the grid. The piece
that carries the collapse at x ≥ 0.10 is the **quasi-elastic** s/p column, and
it is bounded by **nothing**. `sp_tensor_scale` defaults to 0, every default is
bit-identical, and **no registry row was added**: its price is recorded into
the existing §B15.

**B2 — `RcTailModel::PolradFull`, POLRAD Eq. (18) + Appendix B + Eq. (A.4).**
The exact tail exists, is opt-in, and **is not inside the band the two t-peak
models were being run as**: over the sampler's own 3051 accepted cells it lies
between them on **1725 (56.5 %)** and **outside on 1326 (43.5 %)**, and above
both in the Q² ≥ 20 GeV², y ≤ 0.9 window. The σ-weighted mean-shift ratio is
**0.428** — a ratio of means over cells, **not** an event-weighted number and
**not** a mean of per-cell fractions, which is 6.483, and three sites had
called it event-weighted. So the two t-peak models are a **price range, not a
confidence interval**. It also corrects a shipped sentence: *"POLRAD supplies
no tensor s/p peak"* is true of **Eq. (38)** and **false of the paper**. And it
leaves the quasi-elastic half of B1's gap exactly where B1 left it — a nucleon
has no tensor structure function, so the quasi-elastic tail is tensor-blind at
all three peaks on all three models. **T9's unconditional skip is gone**: the
Rosenbluth limit it guards is now a real gate.

**B3 — Q(⁷Li) is sourced, committed and gated, and it PAYS registry row 20's
D11.** `LI7_QUADRUPOLE_FM2` = **−4.06 fm²** at `include/lipolgen/rc.hpp:726`,
beside `LI6_QUADRUPOLE_FM2`, from TUNL's A = 5, 6, 7 evaluation (Tilley *et
al.*, NPA 708 (2002) 3) — whose A = 7 half prints `Q = −40.6 ± 0.8 mb` and
whose A = 6 half prints the `Q = −0.818(17) mb` that **is**
`LI6_QUADRUPOLE_FM2`: one document, one sign convention, one unit rule, which
is why −4.06 was taken over the −4.00(3) the note had quoted from memory. The
row's own predicted price is **confirmed to six digits**: the A = 7 α–t gate
ratio is **0.871265** against −4.00 and **0.858389** against −4.06 (12.87 % vs
14.16 % low), which does not change the verdict "13–14 % low". §4.2's gate is
real code now (`li7_alpha_t_quadrupole`, reproducing the note's offline recipe
**bit for bit**, pinned by doctest T13 and pytest G8) — and it is a **reported
ratio, not a b₁ verdict**: b₁(⁷Li) is still unimplemented and still blocked on
D2.

## C — the decay-lepton reconstruction efficiency: a negative, stated as one (`aefd535`, `phase_C_numbers.md`, `phase_C_survey.md`)

**No defensible bound exists**, anywhere in the sibling `PolarizedLithiumSim`.
Read end to end for one: every module of `fastsim/polli_fastsim/`, all of
`tools/fullsim/` and `tools/analysis/`, `plans/03`, `plans/09` and all 54
`refs/` entries. What is there is one per-track stand-in that the file itself
labels a stand-in (`polligen.hfs.HadronResponse`, `evgen/polligen/hfs.py:251`),
one **constructed** scattered-electron ID profile
(`evgen/polligen/reco.py:482`), no J/ψ reconstruction and no EICrecon output at
all. The one published lepton-identifying number in the corpus is **ATHENA's**
standalone-barrel **electron** ID efficiency — not ePIC, not reconstruction, no
muon counterpart, and quoted for |η| < 1, where only **0.4477** of pairs have
**both** leptons inside: it covers 45 % of the pairs and is silent on the other
55 %. **It does not bound ε_pair.** So the factor is **unbounded rather than
merely unmeasured**, and that sentence is now at all eight sites the O5 record
enumerates.

**A bound WAS found — on a different leg, and it points down.** By Chang *et
al.*'s own §IV, `COHERENT_JPSI_EFF_IR8_LI7` = **0.1775** — the single largest
multiplier in the O5 chain — is a **pure geometric acceptance** and therefore
an **UPPER bound** on the far-forward intact-recoil tagging efficiency ×
acceptance, with no Roman-Pot detector efficiency on the recoil leg either.
Direction **DOWN**, unquantified. It is a **third** factor on the leptons'
side, so it **compounds with the decay-lepton omission and breaks** the
×1.1224-UP / ×0.897-DOWN = 1.007 balance that is the stated reason neither of
those two may be quoted alone. **Nothing was folded in**: the tree states the
bound at the constant's home, in `o5_a2_reach.py` and in
`OPEN_ITEMS_SOLUTIONS.md` §11.3, and the price is recorded into registry
**row 21 / §B25**. The O5 band is untouched at both ends, re-verified by
running the script at the close-out: headline **S = 2.618 σ, 3 σ at
13.1 fb⁻¹/u**; LOW **2.627** at 13.0; TOP **2.842 … 3.292** at 11.1 … 8.3;
the de-squeezed optics rung **0.734 σ at 167.1 fb⁻¹/u**; `self-check: OK`.

## D — the documentation gate's residue class (`ba773c1`, `phase_D_numbers.md`)

Three rules, **each chosen by a measurement whose recipe is a committed script**
(`validation/record_rule_sweep.py`). **R5**: `--record-ranges` refuses a block
the citing sentence does not evidence, judged **per citation** — 89 of 180
recorded blocks were refusable, 22 wrong blocks were re-pointed by content, and
103 per-citation exemptions carry pins. **B2**: an unnamed point citation in a
strict document must quote something on the cited line — 30 of 43 flagged, 7
stale ones re-pointed. **R6**: a relaxed gate over `docs/open_items/`, which
repaired **266 broken citations**; its point rule (a ≥ 6-character token, or
≥ 5 digits, absent from both adjacent lines) was chosen because it takes the
±1-drift acceptance from **29 % to 12.7 % / 8.2 % with zero honest breaks**,
while no range rule reaches an acceptable price (+1 still passes 89.0 %).
Every attack from three review rounds is now caught. The residues are
**recorded, not hidden**: anchor-read point pins are checked for presence and
not uniqueness, and a neighbour's pin text counts as evidence.

## E — close-out (working tree, 2026-09-15/16)

Every tally in *Start → end* above was re-measured here, not carried forward;
`--fix` and `--record-ranges` were run and are no-ops. **Eight further
instances of the run's recurring defect — a claim stated more strongly than
measured, or corrected at one site and not another — were found across
nineteen sites**, each fixed where it stood:

1. **"All 25 decisions."** The registry has been **27** since 2026-09-15.
   Corrected at `README.md`, at `../run_2026-09-03/SUMMARY.md`'s header and at
   its item 9, which now counts the shapes at HEAD: **eight applied** (five
   *confirm or revert* — 1, 11, 18, 19, 26 — two *confirm* — 4, 12 — and one
   *applied in part* — 6), one decided on the record (7), three carried
   forward (22, 23, 24) and **fifteen open**.
2. **"Row 17 still carries the word 'unpriced'."** It has not since
   2026-09-06: §B13(d), the index's shape column and `STATUS.md` row 17 all
   said PRICED while `AUTHOR_DECISIONS.md`'s Part B item 8 said the opposite,
   and `../run_2026-09-03/SUMMARY.md` item 8 still listed all three costs as
   unpriced. Both corrected: **all three are priced, and none is decided.**
3. **"Nothing in that repository was run, and no number from it is quoted"**
   (registry row 27 and §B28(d)). The first half is true; the second is not —
   the sibling's `evgen/README.md:113-115` numbers and its published figure
   are quoted, by reading, and §11 of `phase_CW_numbers.md` itemises eight
   such sites. The one **measured** number is this tree's own reproduction of
   that curve value, **−0.4786 pre-fix → +0.9279 after**. Corrected at both
   registry sites and in this board's CANDIDATE section, whose "was
   read-modified" also read as though the files had not been read.
4. **"Neither was applied" / "Left untouched, deliberately"** (phase C §C5 and
   §C5.1) — while **the same commit `aefd535` applied C5.1's citation fix at
   three sites**, and this board's own C row said so two sentences earlier.
   `HfsModel` names no class in the sibling (0 grep hits against 26 for
   `HadronResponse`); the three live sites were repaired in that commit and
   the two remaining **dated-record** sites at the close-out. C5.1 is a
   citation repair and **not an author decision at all**; C5.2's bound is the
   only decision that section carries.
5. **"All nine `validation/reference/*.json` are byte-identical to
   `a94fd6e`"** (`../run_2026-09-03/SUMMARY.md`) — true of that run, and not
   of the tree since 2026-09-06, which the same page's item 6b describes
   without correcting the sentence. Re-measured by sha256 and boxed there.
6. **Two stale suite tallies presented in the present tense**: `README.md`
   still named the doctest skip "RC test T9, a code path v0 does not
   implement" — gone since `af3f415`, and the tree now measures **0 skipped** —
   and `docs/PACKAGING.md` still said "both trees now read 401 / 17 240 286".
   Both dated and corrected.
7. **"101 headers" as the present cost of changing the licence** (§B21(d) and
   the index table's shape for row 24), where the gate measures **103** and
   the state board and `README.md` already said 103. The *historical*
   sentences — "E2 stamped the identifier onto 101 files" — are right and were
   left alone; only the two that price the decision **today** were corrected.
8. **A decision assigned to a registry row the row did not carry.** This
   board's C row says *"whether to fold a lower value in is registry §B25's"*,
   and §B25 said nothing about the recoil-leg bound — so phase C's one real
   decision had a pointer and no home. §B25(c) and (d) now carry it, with what
   is measured (the band is untouched at both ends) and what is not (there is
   **no** measured price for folding a lower ε in, because nobody has a number
   to fold — which is exactly what an upper bound with no lower one means).
   The board's "New registry rows" section, which had been silent about phase
   C altogether, now says why C added no row.

**And the close-out caused, then repaired, an instance of the very residue
class phase D built R6 for.** Correcting `README.md` and
`OPEN_ITEMS_SOLUTIONS.md` moved lines below them; the strict gate caught **one
broken citation** and R6 **eighteen**, and every one was re-pointed **by
content** — the old target line was reconstructed and located in the new file —
not by adding the same offset to a number. Both gates read 0 broken afterwards.
Had those rules not been built in phase D, nineteen citations would have gone
stale silently in the same pass that was tidying the run.

---

## What this run did NOT establish

Carried forward from `../run_2026-09-03/SUMMARY.md`, each with what this run
changed — and the items this run opened itself.

1. **No decay-lepton reconstruction efficiency exists, and it is now
   UNBOUNDED rather than unmeasured.** The sibling was surveyed end to end and
   supplies none (phase C). Only the geometry is bounded, A_geom = 0.994 at
   10 × 99.5. **And the recoil leg is now bounded only from ABOVE**
   (`COHERENT_JPSI_EFF_IR8_LI7` = 0.1775, a pure geometric acceptance), which
   points the same way as the lepton omission and breaks the chain's stated
   ×1.1224 / ×0.897 cancellation. `Optics::lumi_fraction` is still the single
   correction that on its own restores a NO (0.734 … 0.920 σ,
   106.4 … 167.1 fb⁻¹/u), and the far-forward working point is still unchosen.
2. **Whether O5's band TOP crosses 3 σ.** Unchanged by this run and
   re-verified bit-identical after phase C: 2.842 … 3.292 σ, so it straddles.
   The low edge, 2.627 σ at 13.0 fb⁻¹/u, holds.
3. **b₁ beyond A = 2 is not validated.** The gate is A = 2; the α–d step is
   untested, the ±100 % band stays, and the *shipped* unpolarised default is
   outside the window the lift was granted on. **What changed:** the ⁶Li
   tables are now a **band** over the two axes (A4) rather than a single
   `ToyF2` column, and which axis moves what is measured rather than assumed.
4. **b₁(⁷Li) is not implemented** and is still blocked on D2, deliberately:
   its sign flips with the unpolarised backend. **What changed:** the second
   half of this item is CLOSED — Q(⁷Li) is sourced, in the tree and gated
   (B3), and the A = 7 gate is a real pinned test reporting **0.858389**,
   14 % low. That gate validates the α–t **wave function**, not b₁.
5. **The exact Mo–Tsai elastic tail was never obtained** — `[MT69]` is not in
   this tree and no number from it is quoted. **What changed:** POLRAD's own
   exact tail is implemented (`--rc-tail-model polrad-full`) and checked
   POLRAD-internally and against the leading-log fallback — **not** against an
   external exact tail — and it lies **outside** the two t-peak models on
   **43.5 %** of the sampler's cells, which makes those two a price range and
   not a bracket. The **elastic** s/p tensor fraction is now bounded and the
   bound is **empty at the standard points** and ≤ 21.72 % anywhere; the
   **quasi-elastic** s/p column is bounded by nothing and **is not computed by
   the exact model either** — a nucleon has no tensor structure function. The
   polarised quasi-elastic tail is still *priced, not computed*; no RC
   calculation exists for a tagged tensor asymmetry; and everything the band is
   calibrated on is still **deuteron**.
6. **O3 in the form asked** — the α–d split from GFMC configurations — is
   still unanswerable here: this tree has no GFMC configurations. Unchanged.
7. **Nothing was pushed, published or sent.** `ci.yml` has still never been
   parsed by GitHub, no wheel was uploaded, and the Mäntysaari-group letter is
   still drafted and unsent — with its factor list now carrying one entry
   bounded above and three not bounded at all (row 21 / §B25).
8. **The sibling repository was not touched.** `PolarizedLithiumSim` carries
   the same inverted S–D phase at `tagged.py:247` (until its own commit `1066555` of 2026-09-15 fixed it and regenerated the figure — folded A_zz +0.491 → −0.843), in a **published figure**
   and in five of its own test pins. It was **read and never modified**; the
   one execution is recorded in `AUTHOR_DECISIONS.md` §B28(d) — on 2026-09-16,
   after the sibling's own fix `1066555`, its `polligen.tagged.TaggedModel`
   was run read-only to produce that section's 1.33e-13 / 2.28e-13
   comparison, which cannot be obtained by reading. This item said "never
   executed" until 2026-09-16. Whether its numbers are regenerated is registry
   **row 27**.
9. **Twenty-seven decisions sit in the registry and this run took none of
   them.** It priced three that had been left unpriced (rows 3, 17, 20/D13),
   paid one (row 20's D11), moved the evidence under three more (rows 9, 10,
   21), added two (rows 26, 27) — and **applied** exactly one thing that
   changes a shipped output, the S–D phase, which is row 26's *applied —
   confirm or revert*. **Pricing a cost is not taking it.**
10. **Three rules were added to the documentation gate and three residues were
    left, named** (phase D §D4): anchor-read point pins are checked for
    presence and not uniqueness; a neighbour's pin text counts as evidence;
    and no range rule at an acceptable price catches a +1 drift (89.0 % still
    pass). R5 and R6 are gates on **citations**, not on physics.
