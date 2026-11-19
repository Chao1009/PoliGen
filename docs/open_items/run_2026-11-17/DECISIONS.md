<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Decisions of 2026-11-17

Taken on the three items run 2026-11-10 left for the author
(`../run_2026-11-10/SUMMARY.md` "For the author"; the same list on
`../run_2026-11-10/STATUS.md`). The author's words, in full:

> 1. go with your recommendation, is it physics-based? 2. OK. 3. sounds
> good. Any gaps need to be closed? If so, continue working on them

So: **K1 is the author's decision** — on item 1 the author took the
recommendation as it stood in the summary, and asked a question that is
answered under K1. **Every detail under K1, and K2–K4, are the SUPERVISOR's
readings and choices** under those words — to be labelled so wherever they
are quoted, never as the author's words. The numbering follows the summary's
list (K1 is its item 1, …); K4 is the author's closing question. Each entry
gives the choice (bold), why, what changes in the tree and what does not.
The phases that implement them, their records and gate tallies are on the
board, `STATUS.md`.

Ground rules carried over from the previous runs: a default moves only under
an entry below; every number that depends on a moved default is re-taken with
a dated comment naming the entry (`K<n> (2026-11-17)`), or the test passes
the old value explicitly where its purpose is to pin the old physics;
tolerances and model choices are declared before verdicts (where one is
chosen after the values were seen, the record says so); a benchmark that
disagrees is recorded, never tuned; `validation/reference/*.json` are not
touched; `docs/HEPMC3_CONVENTION.md` is not edited; the hashed design blocks
of `t4_djangoh_rad_noRad.py` and `t2_boeglin_nk.py` and the first 160 lines
of `t2_hermes_b1_table2.py` are never edited; the four author decisions of
2026-10-27 (D4–D7), the seventeen entries of 2026-11-03 (A1–A17) and the
nine of 2026-11-10 (J1–J9, with their dated notes) are not re-opened; the
sibling repository `PolarizedLithiumSim` is read-only; nothing is published
or sent.

**K1 — the digitized deuteron b₁ tables are read as the ratio b₁/F₁ at a
reference scale: b₁(x, Q²) = b₁ᵗᵃᵇˡᵉ(x) × F₁(x, Q²) / F₁(x, Q₀²) with
Q₀² = 2.5 GeV², so that A_zz ∝ b₁/F₁ does not depend on Q² at fixed x.**
*The author's decision, on the recommendation of
`../run_2026-11-10/SUMMARY.md` item 1, which read: "read the toy as the
ratio b₁/F₁ at a reference scale (a scale-independent A_zz) — as its own
decision: it moves every rank-2 number of ⁶Li at every Q²." Why: `MillerB1`
is a Q²-independent table used against an evolving F₁, so b₁/F₁ grows with
Q² at fixed x (0.54 → 6.99 at x = 0.955 from Q² = 2.5 to 1897 GeV² on the
deuteron, against the pure-state bound 3) — the cause of every high-x rank-2
violation run 2026-11-10 met.*

*The author's question, "is it physics-based?" — the supervisor's answer.
Yes in what it removes; an approximation, the literature's own, in what it
puts in its place.*

1. *What it removes is unphysical.* At leading order b₁ is a charge-weighted
   combination of the same quark densities F₁ averages — the densities of
   the target in its helicity states, the ones the J1 note's positivity
   domain is written in. Their Q² evolution does not depend on the state of
   the target, and it keeps each of them non-negative. A b₁ that stands
   still while F₁ falls by a factor 13 therefore leaves the domain for no
   reason but the table's missing scale (the supervisor's statement of the
   textbook argument, not re-read in a source for this entry).
2. *What it puts in its place is what both published deuteron descriptions
   do.* In the standard convolution (Cosyn, Dong, Kumano, Sargsian, PRD 95,
   074036, Eq. (16), as read in
   `../run_2026-09-03/phase_A_miller_normalisation.md`) b₁(x, Q²) =
   ∫ dy/y δ_T f(y) F₁ᴺ(x/y, Q²): all of the Q² dependence enters through F₁.
   Kumano's parametrisation of the HERMES data (PRD 82, 017501; the text of
   arXiv:1005.4524 read for this entry) writes the tensor-polarised
   distributions as one weight function times the unpolarised ones,
   δ_T q(x) = δ_T w(x) q(x), "ignoring the scaling violation", with "Q² =
   2.5 GeV², which is about the average Q² of the HERMES measurements".
   CDKS show their curve "at Q² = 2.5 GeV² for comparison with the HERMES
   data". K1 is that weight-function form applied to the structure
   functions themselves.
3. *It is the assumption polarised DIS makes for A₁* (g₁ = A₁ F₁ with A₁
   taken independent of Q²), here for A_zz.

*What it is not: a derivation.* b₁/F₁ is exactly scale-independent only if
b₁ and F₁ evolve alike in shape, which QCD does not promise (b₁ changes
sign, and a tensor-polarised gluon mixes in). The exact treatment is the
evolution of tensor-polarised parton densities, which needs a parton
decomposition neither digitized curve provides: Miller's (PRC 89, 045203;
the text of arXiv:1311.4561 read for this entry) is a pion-exchange term
computed with pion distributions "at Q² = 1.17 GeV²" plus a hidden-colour
term whose strength is chosen to reproduce the HERMES point at x = 0.452
(⟨Q²⟩ = 4.69 GeV²) — the curve has no single scale of its own. The
reference scale is therefore a choice, and its price is measured (below).

*The details the recommendation left open — the SUPERVISOR's choices, all
declared here before any number was taken on the new reading:*

- **(a) Q₀² = 2.5 GeV².** The scale at which both descriptions above meet
  the HERMES data; inside the shipped CT18NLO grid (floor 1.677 GeV²);
  between Miller's 1.17 and 4.69. One named constant.
- **(b) The F₁ is the kernel's own, per nucleon, at both scales** — the
  number the kernel divides by in A_zz is the number it multiplied by. At
  Q² = Q₀² every rank-2 number is exactly what it is today.
- **(c) Scope: every place a digitized deuteron table is the rank-2 input**
  — `MillerB1` (the default) and `CdksB1` (`--b1-model cdks`; its table is
  the Q² = 2.5 GeV² curve by name) on the inclusive ⁶Li and deuteron
  kernels the run surface builds, and the opt-in `--inclusive-b1` slot of
  tagged ⁶Li α. Not `li6-convolution`, which convolves F₁ᴺ at the running
  scale already (to be confirmed in the code; if it is not so, recorded, not
  changed here).
- **(d) On the run surface only.** A kernel built at class level and the
  factories keep the frozen table unless asked, as D6, D7 and the J1 note
  did — so the polligen port gates and `b1_default_li6.json` are untouched
  bit for bit.
- **(e) One knob, `--b1-ref-q2`** (default 2.5; `0` selects the frozen
  table, the pre-K1 reading, which stays constructible and selectable),
  labelled in `meta`, the banner and the knob table under A17's rules; no
  accepted value is silently read as another.
- **(f) b₂ follows b₁** by the same factor. The J1 note's confinement of
  the deuteron stays and acts on the final b₁.
- **(g) Where F₁(x, Q₀²) ≤ 0 the scaled b₁ is zero** (the only point of the
  domain there); whether a run ever reaches it is measured.

*What changes:* every rank-2 number of the default inclusive ⁶Li and
deuteron runs away from Q² = 2.5 GeV² (b₁, b₂, A_zz, the tensor weights) —
measured, with the price of the reference scale (Q₀² at the grid floor
and at 4.69 GeV²) beside it; the share of the deuteron's confinement; the
HERMES row, whose headline is by convention 8 the b₁ a default deuteron run
carries at each bin's ⟨Q²⟩, with the frozen table recorded as a sub-row and
the status whatever its declared rule gives; the pins and documents that
state the old reading, each re-taken with `K1 (2026-11-17)`. *What does
not:* ⁷Li, the coherent channel, the tagged channels at their defaults, the
vector sector, every unpolarised number; the class-level kernel, the
reference tables and their gates; `li6-convolution`.

**K2 — ⁶Li is left as K1 leaves it; its three joint-domain cells are
re-measured after K1.** *The supervisor's reading of the author's "OK" on
item 2, which read: "Item 1 would settle both; otherwise extend the
confinement to ⁶Li." If the measurement after K1 finds no ⁶Li cell of a
default window outside the spin-1 domain, the item is closed by K1. If it
finds one, extending the J1 note's confinement to ⁶Li is the fallback this
"OK" covers; it would be its own phase with its own measurement, labelled
the supervisor's choice.*

**K3 — the smaller items of the summary's item 3 stay as they are
recorded.** *The supervisor's reading of the author's "sounds good":
`--b1-model cdks` on the deuteron stays refused; `pol_sf` keeps its
configuration-level rule; the empty-map corner of `set_run_attributes`
stands; the T6 row's 30 s ceiling is not moved (a re-run settles a loaded
machine, and a looser ceiling would hide a real slow-down). None is
scheduled.*

**K4 — the gaps: what K1 itself opens is closed in this run; the rest is
recorded with the reason it is not.** *The supervisor's assessment, under
the author's "Any gaps need to be closed? If so, continue working on
them". Closed here, because K1 opens them or makes them stale:* the same
reading on the second digitized table and on the opt-in tagged slot
(K1 (c)); the HERMES row re-pointed by its own rule; the rank-2 numbers
the documents quote (the A_zz price of the channel ranking, USAGE,
PHYSICS_CHANNELS), re-taken with a date; the throughput, if K1 costs an F₁
evaluation per call (J9's recipe). *Not closed — none blocks a default run,
and each needs the author or an input that does not exist:*

1. *The leading-order evaluation of the NLO sets (J8)* — the largest
   physics gap the report still shows: the shipped CT18NLO, evaluated with
   leading-order formulas, sits 8.6 % above the NMC F₂ᵈ (the two absolute
   F₂ᵈ rows FAIL), and the Bjorken integral of the shipped pair is 12 %
   high. Closing it replaces D7's shipped pair of evaluations, which is the
   author's decision; recommended as the first item of the next run.
2. *Miller's normalisation* (per deuteron or per nucleon, a factor 2 on the
   default A_zz; `docs/OPEN_ITEMS_SOLUTIONS.md` §10) — the author's
   recorded decision; only a numerical reproduction of his Eq. (20) closes
   it.
3. *b₁(⁶Li) beyond the α–d picture* — no calculation exists for any A > 2.
4. *The recorded disagreements with data on shipped defaults* (BONuS,
   Deeps, the ⁶Li charge form factor's first zero) — recorded, not tuned.

## Dated notes

**K1 (c) note (2026-11-19) — `li6-convolution` is only partly at the running
scale; K1 is extended to the digitized column inside it.** *The supervisor's
ruling after phase P1 (`phase_P1_b1_ratio.md` §5 item 1), declared before
any number is taken on it. (c) left the convolution out on a premise it
marked "to be confirmed in the code", and the code does not confirm it: of
the four terms only the two α–d D-wave terms fold F₁ at the running Q²; the
embedded-deuteron S term and the CG D-wave term fold the digitized CDKS
Q² = 2.5 GeV² column, which has no Q² dependence — and that part is the
larger one ((s + cg)/total = 1.68 / −0.83 / 1.15 / 0.72 / 0.91 at (x, Q²) =
(0.1, 2.5) / (0.1, 100) / (0.3, 10) / (0.5, 100) / (0.7, 1000) on the
default object, P1's measurement). It is a digitized deuteron table used as
a rank-2 input, which is what the author's K1 is about, so it gets the same
reading rather than a second convention inside one backend. P1 changed
nothing there, as (c) said; phase P2 does, under these details:*

- *Under the integral the embedded deuteron's column is read as
  b₁ᵈ(z, Q²) = b₁ᵈ,ᵗᵃᵇˡᵉ(z) × F₁ᵈ(z, Q²)/F₁ᵈ(z, Q₀²) at the folded argument
  z, with the convolution's OWN per-nucleon deuteron F₁ — the one its α–d
  D-wave term folds (the `--b1-unpol` backend) — at both scales. Not the
  kernel's F₁ at x as in (b): the factor belongs under the integral, where
  the only deuteron F₁ is the convolution's, and all four terms then take
  their Q² dependence from one F₁. Zero where F₁ᵈ(z, Q₀²) ≤ 0, as (g).*
- *An option of the convolution object, default 0 = the frozen column: the
  class level, the factories at their defaults, the reference tables and
  their gates stay bit for bit, as (d). At Q² = Q₀² every number is the
  frozen one — the published band tables are at Q² = 2.5 GeV² and do not
  move.*
- *The run surface hands `--b1-ref-q2` to it, so the knob is READ on
  `--b1-model li6-convolution` (0 = the frozen column), with (e)'s rule
  against the Q² range of the `--b1-unpol` grid. The kernel-level factor of
  (b) is not applied on this model as well: one reading, once.*
- *This does not make b₁/F₁ of this backend independent of Q² (the folded F₁
  and the kernel's F₁ at x differ), and is not meant to.*
- *Measured in P2: the four terms and the total before → after; the band
  rows on runs; which `--x-max` builds; `--b1-ref-q2 0` against P1's tree
  byte for byte.*

*The backend stays opt-in under its mandatory 100 % band; the ruling is
listed for the author in the summary, and `--b1-ref-q2 0` undoes it.*

**K1 (d) note (2026-11-19) — `lipolgen.shipped_inclusive_kernel` is the run
surface, so its default follows it.** *The supervisor's ruling. P1 read
(d)'s "the factories" as covering this Python helper and left its new
keyword at 0. But the helper is by its own definition (D7) "the inclusive
kernel a default run builds", and after P1 it is that only when asked.
Phase P2 makes its `b1_reference_q2` default `B1_REFERENCE_Q2` (0 still
selects the frozen table); `default_inclusive_kernel`, the class-level
factory the port gates and `b1_default_li6.json` stand on, keeps 0.
Benchmark rows built on the helper then carry the run's b₁, as convention 8
asks; a row that moves is re-taken by its own rule and listed, and a test
that pins the frozen reading passes 0 explicitly.*

**K2 note (2026-11-19) — closed by K1.** *Phase P1's measurement M6: over
the three default windows ⁶Li has no cell outside the spin-1 domain (3 → 0;
max b₁/F₁ 2.290 → 0.178, the joint maximum 1.046 → 0.342). By K2's own
wording the item is closed; the confinement is not extended to ⁶Li.*

**K4 note (2026-11-19) — one more gap, recorded and not closed.** *5. The
HepMC3 provenance of K1: the new `b1_reference_q2` key is in the npz `meta`
only, because `docs/HEPMC3_CONVENTION.md` lists the `lipolgen.*` run
attributes and is not edited; a HepMC3 file therefore does not say which
reading made it (`lipolgen.b1_model` is `miller` under both). Amending the
convention is the author's.*
