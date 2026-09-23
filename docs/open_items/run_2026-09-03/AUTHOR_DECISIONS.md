# Run 2026-09-03 — the author-decision batch

Drafted in Phase F by the decisions reviewer, from the working tree at
`dffe94e` (10 commits ahead of `origin/master` = `a94fd6e`, tree clean).
**Nothing in this file is decided.** Every section states one question, its
options, what the run's evidence says, the measured consequence of each option
on shipped or published numbers, and what changes in the tree once the author
answers. Numbers are quoted from the run's own records with the file that
measured them; where a cost was *not* measured the section says so instead of
estimating one.

State of the tree this file was drafted against, re-measured here (not copied):

| what | measured 2026-09-05 at `dffe94e` | matches `STATUS.md` phase E row |
|---|---|---|
| `build/lipolgen_tests` | **400 cases / 17 240 262 assertions / 1 skipped / 0 failed** | yes |
| `python -m pytest python/tests -q` | **925 passed, 112 skipped** (144.64 s) | yes |
| `validation/check_physics_channels_links.py` (strict, `env.sh` sourced) | **1145 references / 97 ranges / 7 external / 0 broken / 6 allow-listed**; `SPIN32_FINITE_GAMMA.md` 19 / 115 / 0 broken; `PYTHIA_BRIDGE.md` 8 external | yes |
| `validation/check_spdx_headers.py` | **101 files**, all carry the identifier | yes |

*(Phase F then applied its 48 review findings and wrote the close-out
documents on the same day, and three of those four numbers moved — by construction, not by regression: **401 doctest
cases / 17 240 286 assertions / 1 skipped / 0 failed** (one case added, `item
5v`, splitting gate condition 3 out of checklist item 5, plus the `REQUIRE`s of
the cases whose data guard became a `doctest::skip` decorator), **926 pytest
passed / 112 skipped** (one test added, for `--lumi` alone), and
**1206 references checked (strict) / 96 ranges / 7 external / 0 broken / 6
allow-listed** (re-measured 2026-09-16 from `git archive bd775bc`; the
**1205 / 97** published on the day was the pre-commit working tree's) — 27
citations added when the findings were applied, 33 more
when the close-out closed the consistency lens's remaining gaps, all
re-measured here at the end of the day. SPDX is unchanged at 101/101.
Nothing in this file's evidence depends on those four numbers.)*

**Ordering.** Sections are ordered by how much a published number depends on the
answer: first the decisions that scale or sign every default tensor number,
then those that move an opt-in path or a pinned document, then those that move
nothing numerical, then the outward-facing ones (a name, a push, a wheel, a
letter). Each heading carries the decision's registry id, because the tree
kept **four** registries when this file was drafted — as of 2026-09-05
`STATUS.md`'s decision table is the ONE registry, its rows 16–25 carry what it
was missing, and §10/§15.5/the top table point at its row numbers; the four
were:
`STATUS.md` decision-table rows 1–15; `OPEN_ITEMS_SOLUTIONS.md` §10's close
conditions 4 and 6; `OPEN_ITEMS_SOLUTIONS.md` §15.5's D1–D13 (the ⁷Li rank-2
list, mirrored as D1–D12 in `phase_D_li7_rank2.md` §7); and the three
"author to confirm" rows (3, 4, 13) of `OPEN_ITEMS_SOLUTIONS.md`'s top table,
which predate this run but were built on by it. Decisions the run **applied**
(a default or a behaviour changed, evidence called certain) are marked
*applied — confirm or revert*; nothing else in the tree moved a default.

Part B at the end lists what this review found wrong with the bookkeeping
itself — decisions missing from the table, a dependency corrected at one site
and not another, and a decision whose framing does not match the tree.

## The index — every registry row, and where it is stated in full

Two files, one job each. **`STATUS.md`'s decision table is the REGISTRY**: it
holds the row numbers, and a decision is cited by its row number and never by
its position. **This file is where each of those rows is STATED IN FULL** —
question, options, what the run's evidence says, the measured consequence of
each option, and what changes in the tree once it is answered. Nothing is
decided in either.

The mapping is total in both directions, and was checked row by row on
2026-09-05:

| `STATUS.md` row | stated in full at | shape |
|---|---|---|
| 1 `CdksB1` stops halving the CDKS column | §B12 | applied — confirm or revert |
| 2 Miller's ×0.5 stays | §B2 | open, *likely not certain* |
| 3 `Li6ConvolutionOptions` keeps `r_sigma_lt` | §B7 | open (third option implemented opt-in and **priced 2026-09-06**; the row becomes *(iii) with which R?*) |
| 4 the ⁶Li publication ban is lifted, the band is not | §B4 | applied — confirm |
| 5 the A = 2 gate's wave default stays `kFdeutFile` | §B17 | open |
| 6 the published `A_zz(Born)` column is wrong ×3.253983 | §B11 | **applied in part — confirm** (the corrected column is published beside the original, 2026-09-06, which executes option (ii) in all but the deletion); **the deletion of the original is still open** |
| 7 the Gakh–Shekhovtsova δ(x) shape is rejected | §B14 | decided on the record |
| 8 `eps_b0` stays −0.08 | §B5 | open (cost recorded) |
| 9 the tagged deuteron control keeps Hulthén | §B9 | open — **evidence moved 2026-09-06**, see the box in §B9(c): the Cosyn–Weiss gate now runs on AV18, and `tagged.json`'s spin-1 model blocks are re-pinned from this library. No decision taken **on THIS row** — the phase fix and the re-pin themselves are **row 26 / §B27**, added 2026-09-15 |
| 10 `LI6_CLUSTER_POLARIZATION` stays 0.811228 | §B6 | open (band mandatory) — **checked 2026-09-06 against the S–D interference fix: does not inherit it.** 0.811228 = (1 − 1.5 P_D^{αd})(1 − 1.5 P_D^{d}) is a closed form in two D-state probabilities, and P_D is a norm, phase-blind |
| 11 the embedded deuteron follows `--cluster-wave` | §B10 | applied — confirm or revert. **Re-measured 2026-09-06**: §B10(d)'s +2.069 % is unmoved by the S–D interference fix (1.020687624664 vs 1.020687500612) |
| 12 the coherent \|t\| ceiling stays 0.2, on the anchor range | §B16 | applied — confirm |
| 13 the author's name | §B22 | open — a placeholder ships |
| 14 push `ci.yml`, or not (coupled to pushing the run) | §B23 | open |
| 15 publish an `auditwheel` wheel, or not | §B24 | open |
| 16 the α–d source default stays `FitRescaled` | §B8 | open |
| 17 the RC band's low-x anchor | §B13 | open (alternatives **priced 2026-09-06**) |
| 18 G3a's counting window, and the floor this run changed | §B4 | applied — confirm or revert |
| 19 `PDF:PomSet` = 11 is refused | §B19 | applied — confirm or revert |
| 20 the thirteen ⁷Li rank-2 items of §15.5 | §B3 (D2), §B18 (D13), §B20 (the other eleven) | open (**D13 priced and D11 PAID, both 2026-09-06**; D1, D3–D10 and D12 untouched, and the row stays open on them and on D2) |
| 21 send the Mäntysaari-group letter | §B25 | open — **evidence moved 2026-09-15** (§B25(c)): the decay-lepton caveat surveyed and still UNBOUNDED, and the recoil leg's `COHERENT_JPSI_EFF_IR8_LI7` = 0.1775 now an UPPER bound, direction-only. No option's number moves; only (i)'s text does |
| 22 `TENSOR_LL_SIGN` = −1 | §B1 | carried forward — confirm only |
| 23 the ⁶Li effective polarization vs the naive 1/3 | §B6 | carried forward — confirm only |
| 24 the licence, GPL-3.0-or-later | §B21 | carried forward — confirm only; **103 files stamped** (101 when E2 stamped them, re-measured 2026-09-16) |
| 25 the MSTW rows still tallied as PASSED when the grid is absent | §B26 | open — phase F's own |
| 26 the tagged S–D interference phase, and the `tagged.json` re-pin | §B27 | **applied — confirm or revert** (added 2026-09-15; the one shipped output that moved in `91e48b9..HEAD`) |
| 27 `PolarizedLithiumSim`'s published A_zz^tag numbers carry the same inverted sign | §B28 | open — about another repository; nothing in this tree changes either way (added 2026-09-15) |

**The one section with no registry row is §B15**, and deliberately: the run
does not call the promotion of its no-op opt-ins (`c0_shape`, `tail_model`,
`qe_tensor_scale`, `a_transfer_frac`, `cluster_vmc_mc_sigma`, `pol_sf`) an
author decision — it is listed so the batch is complete, and every one of those
knobs moves nothing at the CLI default whichever way it goes. If the author
wants it decided, it becomes a new row; nobody here made that call.
*(This sentence read "it becomes row 26" until 2026-09-15, when rows 26 and
27 were added for the tagged S–D phase and the sibling's published numbers.)*

**Three entries are not one-to-one, and they are rows 20, 23 and 18.** Row 20
is three sections, because §15.5's thirteen items split into one that gates
the b₁ ban, one that would move a shipped fill, and eleven of design. §B6
states rows 10 and 23 together, because they are the same constant seen from
two sides. §B4 states rows 4 and 18 together, because the counting window and
its floor are the conditions the ban lift was read through. Every other row is
exactly one section, and no other section is shared.

---

## Part A — the batch

### B1. The tensor-sector sign, `TENSOR_LL_SIGN` = −1 — *carried forward, confirm only* (`STATUS.md` decision **row 22**; `OPEN_ITEMS_SOLUTIONS.md` top table row 3)

**(a)** Is A_zz = −(2/3) b₁/F₁ (Cosyn et al. Eq. 27, HERMES Eq. 6, POLRAD
Eqs. 9/10) the convention of this program, or the repository's earlier
+(2/3) b₁/F₁ transcription?

**(b)** −1 (the literature; in the tree since 2026-08-29) or +1 (the earlier
private convention, reachable by setting the one constant back).

**(c)** `OPEN_ITEMS_SOLUTIONS.md` row 3: "decided by literature: −1 … author
to confirm". `include/lipolgen/constants.hpp:17` calls the same choice "author
decision, plans/08 D1" as if already taken. The two sites disagree about
whether this is open; the evidence (four independent sources give the same
sign) favours −1.

**(d)** Not a run measurement — it is exact: +1 flips the sign of every tensor
asymmetry and every tensor weight in every channel, and every rtol-1e−12
reference carrying a tensor entry (`validation/reference/{xsec,b1_default_li6,
tagged,spin,bookkeeping}.json` all mention b₁/A_zz/tensor; the `b1_default_li6` one is a self-pin of LiPolGen's own C++, not a polligen reference). Nothing this run
published depends on the sign being *re*-decided; everything depends on it
staying.

**(e)** Confirm: `OPEN_ITEMS_SOLUTIONS.md` row 3's "author to confirm" is
closed and the row says so. Revert: `constants.hpp:40` one constant, every
tensor reference JSON re-dumped, every tensor-sector document.

### B2. Miller's b₁ table stays halved, `B1_MILLER_TABLE_TO_PER_NUCLEON` = 0.5 (`STATUS.md` row 2; `OPEN_ITEMS_SOLUTIONS.md` §10 condition 6)

**(a)** Is Miller's digitized b₁ᵈ (PRC 89:045203) per deuteron, so that the ×0.5
to per nucleon is right, or per nucleon, so that it should not be halved?

**(b)** Keep 0.5 (per deuteron — the status quo) or set the constant to 1
(per nucleon).

**(c)** `phase_A_miller_normalisation.md` §§4, 7, 8. For per deuteron: his
Eq. (1) densities are "in a target hadron", Eq. (5) is a light-cone correlator
in the normalised deuteron state, Eq. (6)'s ½ is consumed by the spinless
pion's quark-spin average — no 1/A in Eqs. (1)/(5)/(6)/(20). For per nucleon:
his Table I transcribes HERMES's per-nucleon values unrescaled, Fig. 5
overlays them, and P₆q is tuned to one of them. The paper is self-inconsistent
by exactly this factor. The run's verdict: **likely per deuteron, not
certain**. What would settle it: Miller, or his Eq. (20) at x = 0.012
(10.5 × 10⁻² ⇒ per nucleon; 5.25 × 10⁻² ⇒ per deuteron), which needs a pion
PDF set the tree does not have.

> **Priced against data 2026-09-23** (`../run_2026-09-23/phase_B1_spin_deuteron.md` §3, `validation/benchmarks/t2_hermes_b1_table2.py`, no tuning; re-run 2026-09-23 by the Documents stage with the same output): against HERMES Table II (six bins, stat ⊕ syst), keep 0.5 gives χ²(b₁ᵈ) = 5.262/6 (4 bins within 1σ), drop it 5.297/6 (5 within); χ²(A_zzᵈ, F₁ rebuilt with MSTW + R1990) = 5.550 vs 4.514 (R1998: 5.410 vs 4.587); b₁ = 0 gives 21.80. At x = 0.012 the two options read 5.71 and 11.43 (×10⁻²) against HERMES 11.20 ± 6.17. **HERMES does not decide this row** (Δχ² = +0.035 on b₁), and Miller's P₆q was tuned to a HERMES point, so neither χ² is an independent test.
>
> | configuration | χ²(b₁ᵈ)/6 | bins within 1σ | χ²(A_zzᵈ)/6 (R1990; R1998) | bins within |
> |---|---|---|---|---|
> | `MillerB1` shipped (×0.5 — keep) | 5.262 | 4/6 | 5.550; 5.410 | 4/6 |
> | Miller ×1 (drop the 0.5) | 5.297 | 5/6 | 4.514; 4.587 | 5/6 |
> | `DeuteronConvolutionB1` + `MstwSF` (the A = 2 gate's configuration) | 22.264 | 2/6 | 22.009; 22.008 | 2/6 |
> | b₁ = 0 (baseline, not a model) | 21.799 | 2/6 | 21.586; 21.586 | 2/6 |
>
> The status label does not move: *likely, not certain*, keep 0.5.

**(d)** Keep: nothing moves. Drop the 0.5: **every `MillerB1` number doubles,
exactly** — the DEFAULT inclusive tensor rate, A_zz on the default `--b1-model
miller` path, `validation/reference/b1_default_li6.json` (rtol 1e−12, T9),
the `toy_b1` pins in `tests/test_sf.cpp`, and B11's column (which would then
be wrong by a further ×2). No run measurement exists for the doubled state
because it is a linear factor.

**(e)** Keep: close the row. Drop: `constants.hpp:108`, re-dump
`b1_default_li6.json` (`validation/dump_b1_default_li6.py`), re-pin T9 and
`tests/test_sf.cpp`, and every document table quoting a default ⁶Li A_zz.
Sites that will need the pointer changed: `constants.hpp:107`,
`sf.hpp:398`, `CONVENTIONS.md:222`, `USAGE.md:199`, `README.md:188-189` (re-pointed 2026-09-16; +10 from README's own close-out edits).

### B3. The unpolarised nucleon input — `ToyF2` or MSTW2008 LO — settled once, for both isotopes (`STATUS.md` decision **row 20**; `OPEN_ITEMS_SOLUTIONS.md` §15.5 **D2**)

**(a)** Should the shipped default nucleon input for the b₁ convolution
(`PipelineConfig::b1_unpol`, Phase A) and for every kernel
(`PipelineConfig::unpol_sf`, Phase D) stay `toy`, or become the input the
A = 2 gate passes with?

**(b)** (i) `toy` everywhere (status quo); (ii) `mstw` for `b1_unpol` only
(the convolution's input, which is what the gate tests); (iii) `mstw` or
`ct18nlo` for `unpol_sf` too (the kernel's F₂/F₁ everywhere).

**(c)** The A = 2 gate passes at 0.843243 **only** with MSTW2008 LO; the
shipped `ToyF2` gives 0.440, outside [0.5, 2] (`phase_A_numbers.md` §0, B4
below). `OPEN_ITEMS_SOLUTIONS.md` §15.2: "the unblocking decision is D2: one
unpolarised backend, decided once for both isotopes" — the ⁷Li α+t
convolution's **sign** flips with the backend at four of six x points
(`phase_D_li7_rank2.md` §5.4). The run kept `toy` under its ground rule
(default bit-for-bit), not on evidence; the evidence favours the gate's
input. A non-toy default has a build cost the run states: `--b1-unpol mstw`
"needs the optional PYTHIA tier and is refused, never silently downgraded,
without it" (`README.md:188-189`, re-pointed 2026-09-16), so a core-only build could not run the
default.

**(d)** Measured by Phase D on the ⁶Li inclusive channel
(`OPEN_ITEMS_SOLUTIONS.md:2303-2325`, `USAGE.md:697-707`): under `ct18nlo` σ
×0.79847, A_zz ×1.25239, A_∥ ×0.82519; `mstw` A_zz ×1.26044. ⁷Li inclusive σ
×0.795936 (`ct18nlo`) / ×0.791318 (`mstw`); tagged-⁷Li-α ×0.792543 /
×0.788565 (`OPEN_ITEMS_SOLUTIONS.md` §15.3). On b₁ itself mstw/toy = 1.848 /
1.276 / 0.817 at x = 0.10 / 0.30 / 0.50 (`README.md:194`, re-pointed 2026-09-16). Option (i)
leaves the default's ⁶Li convolution numbers **outside the lift** (B4).

**(e)** (i): nothing. (ii): `pipeline.hpp:857` default, `cli.py`'s
`--b1-unpol` default, the `validate()` collision rule at `pipeline.cpp:1033`,
`meta["b1_unpol"]` expectations in `python/tests`; no reference JSON (the
Miller default does not read it). (iii): `pipeline.hpp:878`, every kernel's
F₂ — `b1_default_li6.json`, `xsec.json` and the other rtol gates move, the
PYTHIA/LHAPDF tier becomes a default requirement, and the `TensorSF`
provenance tables in `USAGE.md` §2a/§14 are republished.

### B4. The ⁶Li publication ban: what the lift covers, and G3a's counting window (`STATUS.md` decision rows **4** and **18**; `OPEN_ITEMS_SOLUTIONS.md` §10 "G3a's stated limitation"; the floor change — *applied*)

**(a)** Does the author accept the lift **as scoped** — a pass for the
configuration MSTW2008 LO + CDKS Eq. (21), not for the shipped default — and
G3a as a *qualified* pass inside a counting window (0, 1.0] that excludes the
third sign change the computed curve has over the reference's own domain?

**(b)** Ceiling: (i) keep (0, 1.0] as a stated scope choice and quote
"exactly two sign changes in (0, 1.0], and a third at x ≈ 1.22 that the window
excludes"; (ii) widen to the reference's domain (0, 1.59] — then G3a fails on
every configuration and the criterion's "three parts, all must hold" does not
hold; (iii) restate the ban's trigger as G3b alone (design §5.4's Escalation
clause reads only the peak ratio). Floor: confirm or revert the change
[0.02, 1.0] → (0, 1.0] the run made in Phase A (`tests/test_b1_nuclear.cpp`,
design §5.4 amended).

**(c)** `phase_A_numbers.md` §0: G3b = 0.843243 (AV18 + MSTW + Eq. (21)),
1.000338 with CD-Bonn as well; 0.440 on `ToyF2`; 0.520 at Eq. (17).
`OPEN_ITEMS_SOLUTIONS.md:925-981`: at the AV18 third crossing the digitized
x·b₁ is 7.0–8.5 % of its own peak (13.5–14.5 % at κ = 1), "a real sign
disagreement in the several-per-cent-of-peak region, not noise"; only the
CD-Bonn rows (0.8–1.0 % of peak) are in the digitisation floor. The run's own
words on the floor: "the phase that re-opened the window argument moved the
edge that cost it nothing and stopped at the edge that would have cost it the
clause." The evidence favours (i)+(iii) for the *lift* (Escalation reads G3b)
and leaves the *criterion* honestly failed under (ii).

**(d)** No shipped number depends on the window: the gate is
`DeuteronConvolutionB1`, refused as a beam species (`pipeline.cpp:778-786`).
What depends on it is every sentence "the A = 2 gate passes" and the lift
itself, in the **17 files** that carry it — the measured count
(`git grep -il 'publication ban'` returns the same 17 at `ac22331` and at
`HEAD`, re-measured 2026-09-05; `phase_A_gate_mechanics.md` §4 enumerates
them), now stated with that basis in `STATUS.md`'s phase A row and
`phase_A_numbers.md`, with `PLAN.md` A7's "fourteen" labelled the pre-survey
estimate and `phase_A_numbers.md`'s "sixteen" the mid-task one (Part B
item 9). Under (ii) item 10's close and the lift are re-argued from G3b alone.

**(e)** (i): close the row. (ii): `tests/test_b1_nuclear.cpp` G3a window,
design §5.4, `OPEN_ITEMS_SOLUTIONS.md` §10, and the lift restated in all
17 files.
(iii): design §5.4's criterion text only.

### B5. `CoherentScenario::eps_b0` stays −0.08, a deuteron-sized number (`STATUS.md` row 8; amended by row 12)

**(a)** Should the coherent tensor slope modulation keep the default that
implies Q_charge(⁶Li) = −0.9345 fm² — 11.42× the measured −0.0818 fm² — or
move into the ⁶Li band?

**(b)** (i) keep −0.08 (status quo, pinned); (ii) the α+d model −0.0527;
(iii) GFMC −0.0171; (iv) the measured −0.0070; (v) keep the constant but
publish nothing from a single row (already the documented rule).

**(c)** `phase_C_numbers.md` §C4, `coherent.hpp:455-481`,
`design_G_cluster_config.md:1500-1510`: the default is a deuteron number
(set from the deuteron's own +0.105 at B_d = 33.1, opposite sign, never
rescaled); the honest ⁶Li band at B = 50 is −(0.0070 … 0.0527), a factor-7.5
band whose 1σ on η spans a sign change. Kept only because
`validation/reference/coherent.json` pins it. `eps_b0` and `slope_b` are not
independent (band the product δ = eps_b0·B).

**(d)** (i): every *generated* coherent tensor number — `a2_deformation`,
`cos2phi_coefficient`, the sampled azimuth — is 11.4× the measured-quadrupole
expectation; a₂(±1, 0.3) = +0.300, larger than the deuteron's own −0.28.
(ii)–(iv): the positivity edge moves 0.245 → 0.372 / 1.146 / 2.80 GeV²
(`t_positivity_edge`, `phase_D_numbers.md` §D5; 2.7990 on the derived
−0.0070024), `coherent.json` moves at rtol 1e−12, every pinned c₂ and azimuth
moves. The O5 J/ψ verdict does **not** move — it was computed on
`a2_from_quadrupole`, not on `eps_b0` (`phase_C_numbers.md:1298-1304`).
**Follow-on that waits on this row:** un-truncating the φ sampler (Σ_m p_m
e^{−ΔB_m|t| cos 2φ} instead of 1 + c₂ cos 2φ) — physically right, moves every
pinned coherent azimuth, filed by row 12 as the successor to this decision.

**(e)** (i)/(v): close the row. (ii)–(iv): `coherent.hpp:481` one line,
re-dump `coherent.json`, T10b/T23, the `USAGE.md` §"ΔB, eps_b0" tables,
`PHYSICS_CHANNELS.md:806`, `OPEN_ITEMS_SOLUTIONS.md` §11 (which still carries
the pre-row-12 dependency sentence — Part B, finding 2).

### B6. `LI6_CLUSTER_POLARIZATION` stays 0.811228; the `--cluster-wave vmc` drift is documented, not closed (`STATUS.md` decision rows **10** and **23**; carries forward `OPEN_ITEMS_SOLUTIONS.md` top table row 4)

**(a)** Which effective ⁶Li polarization does the inclusive channel ship, given
that the cluster product 1 − 1.5 P_D(α–d) × 1 − 1.5 P_D(d) spans 0.811 … 0.905
on the wave functions in the tree and the ab-initio six-body number is 0.848?

**(b)** (i) 0.811228 (Hulthén product; status quo, −4.3 % below 0.848);
(ii) 0.905427 (VMC α–d P_D; +6.8 % above); (iii) 0.887076 (VMC α–d + AV18
deuteron; +4.6 % above); (iv) 0.848 itself (`LI6_POLARIZATION_VMC_SIX_BODY`,
bypassing the product); (v) an opt-in `PipelineConfig` override (priced in
`phase_C_numbers.md` §C5.5: moves every inclusive g₁ by +11.6 % away from
0.848).

**(c)** `phase_C_numbers.md` §C5.5, `beams.hpp:50-76`, `CONVENTIONS.md:80-107`:
the drift between a tagged `vmc` row and the inclusive row of one programme is
**+11.61 %** vector / **+6.58 %** rank-2; on the Hulthén default the
"same wave function in two experiments" claim holds to 1.22e−5. The formula
carries the error, not the choice of P_D (`P_D_LI6` = 0.0867 was a fitted
effective depolarisation, not an α–d D-state probability). The run's reading:
substitution is the wrong fix; the evidence favours keeping (i) with the band
0.81 … 0.91 mandatory on every quoted inclusive ⁶Li polarization.
Carried forward: `OPEN_ITEMS_SOLUTIONS.md` row 4 still reads "author to
confirm" on the 2026-08-29 closure (whole-nucleus P_p = P_n = 0.85 ± 0.03
for ⁶Li, 0.87 / −0.03 for ⁷Li) and does not point at this row.

> **(c)'s LAST NUMBER HAS MOVED, 2026-09-06 — the leg is unchanged in
> direction and TIGHTER.** The "same wave function in two experiments" residual
> is now **7.84e−6**, not 1.22e−5: the tagged `vector_dilution()` re-measures
> **0.8699431789** against the closed form 0.869950 (6.82e−6 absolute) after
> the tagged S–D interference sign fix
> (`../../benchmarking/07_cw_sign_investigation.md`;
> `../run_2026-09-06/phase_CW_numbers.md`). The dilution is angle-integrated,
> so the fix moves only the 96-cell midpoint-quadrature residual of
> `∫Θ₀Θ₂ dc = 0` — **+11.61 % / +6.58 %**, the two numbers the decision turns
> on, are unchanged: re-measured **0.116126** vector (0.116131 pre-fix) and
> **0.065762** rank-2 (unmoved at six decimals), both inside `T26`'s 1e−4 pins
> on 0.116131 / 0.065762.
> **No decision leg changed and the row stays open**; what changed is one
> measured residual. Live sites carrying the new pair: `CONVENTIONS.md:93`,
> `PHYSICS_CHANNELS.md`, `beams.hpp`, `tagged.hpp`.

**(d)** (i): nothing moves; the 11.61 % split stands between runs. (ii)–(iv):
every inclusive ⁶Li g₁ / A_∥ scales by new/0.811228 (×1.1161, ×1.0935,
×1.0453); `validation/reference/beams.json` and `xsec.json` (both carry
0.811228) move at rtol 1e−12; `tests/test_beams.cpp:106` (anchor read 2026-09-15 "0.81123, 5e-6").

**(e)** (i): close the row and link `OPEN_ITEMS_SOLUTIONS.md` row 4 to it.
Otherwise: `beams.hpp:134`, re-dump `beams.json`/`xsec.json`, T-pins in
`tests/test_beams.cpp`, `test_tagged.cpp`, `test_pipeline.cpp`,
`python/tests/test_module.py`, the `CONVENTIONS.md`/`PHYSICS_CHANNELS.md:83,248`
sentences.

### B7. `Li6ConvolutionOptions` keeps R = `r_sigma_lt`, the gate keeps `r1998` (`STATUS.md` row 3; `OPEN_ITEMS_SOLUTIONS.md` §10 condition 4)

**(a)** Should the opt-in ⁶Li convolution backend use the R the A = 2 gate was
validated with (`r1998`), the kernel's own R (`r_sigma_lt`, status quo), or
one R hook threaded through `default_inclusive_kernel` into both?

**(b)** (i) `r_sigma_lt` (status quo); (ii) `r1998`; (iii) one shared hook.

**(c)** `USAGE.md:282-310`, `OPEN_ITEMS_SOLUTIONS.md:1091-1118`,
`b1_nuclear.hpp:487-507`: the shipped observable is a ratio K/D_φ whose
denominator carries `InclusiveKernel`'s R; a different numerator R does not
cancel — (1 + r1998)/(1 + r_sigma_lt) = 1.088 at x = 0.1, 1.039 at x = 0.3
(Q² = 2.5). The evidence favours (i) inside the pipeline and (ii) on the gate,
and the run itself names (iii) "the better third option, not taken here".

**(d)** (ii) measured on the opt-in `--b1-model li6-convolution` x·b₁:
−5.5 % / **+24.6 %** / +5.4 % / +2.7 % / −1.9 % at x = 0.05 / 0.10 / 0.20 /
0.30 / 0.50, all in the two orbital terms (term (1) bit-identical); on the
gate the reverse swap is +2.7 % on G3b (0.8432 → 0.8662). The default
`miller` path and every reference JSON are untouched by any option.

**(iii) was priced on 2026-09-06** (`../run_2026-09-06/phase_A_numbers.md`
§A1; it was unpriced when this section was drafted). It is implemented
**opt-in** as `PipelineConfig::r_source` / `--r-source {unset, sigma-lt,
r1998}` (`RSource`), which threads ONE R object into
`Li6ConvolutionOptions::r_func` *and* `InclusiveKernel::Options::r_func`.
`unset` is the default and does not enter the wiring branch, so today is bit
for bit **by construction**; `sigma-lt` installs the shared object and is
measured bit-identical to it (162 values over x = 0.05/0.10/0.30 ×
Q² = 2.5/5 × y = 0.1/0.5/0.9, plus σ_pb, all three per-category cross
sections and every column of a 2000-event run), which is why its
`knob_provenance` row reads `not-read`; `r1998` moves both halves. It is
refused off `--b1-model li6-convolution`, which is already inclusive-only and
⁶Li-only. **Six things the numbers say, and they do not all point the same
way:**

1. **The tensor weight, at y = 0.5, Q² = 2.5.** K/D_φ — and A_zz, which is
   exactly −(2/3) of it at θ_m = 0, so this is one number and not two —
   moves **−3.8357 % / +26.3988 % / +3.3518 %** at x = 0.05 / 0.10 / 0.30
   under (iii), against (ii)'s **−5.5148 % / +24.6099 % / +2.6629 %**. The
   gap is **1.7 / 1.8 / 0.7 points**.
2. **Sharing the hook does NOT shrink the shift.** At x = 0.10 (iii) is
   *larger* than (ii). This section's own argument — that a numerator R "does
   not cancel" — is right about the mechanism and would be wrong if read as
   "one hook restores the cancellation": it restores part of it, and the ⁶Li
   orbital terms are not proportional to F₁ anyway (they convolve F₁ᵈ against
   a density integrating to zero and see the *slope* of R).
3. **The cos 2φ amplitude separates the two options exactly.** Its numerator
   is Δ = `toy_delta_gluon`, which carries no R at all, so **(ii) leaves it
   bit-identical** and (iii) moves it **purely through D_φ**:
   −8.33 % / −6.73 % / −3.14 % at the same three x, Q² = 2.5. An extraction
   that reads Δ off the cos 2φ amplitude is untouched by (ii) and is not by
   (iii).
4. **(iii)'s shift is y-DEPENDENT and (ii)'s is not** (agreeing to 2 ulp,
   because D_φ does not move and K ∝ b₁ at the default b₂ = 2x·b₁). At
   x = 0.05, Q² = 2.5 the (iii) shift runs **−5.4705 % at y = 0.1 →
   −3.8357 % at y = 0.5 → +2.3678 % at y = 0.9**, straddling zero; at
   x = 0.30, Q² = 5 it runs +2.2010 % → +1.6954 % → −0.0289 %. **A (iii)
   number quoted without its y is not a statement about a run** — a cost (ii)
   does not carry.
5. **(iii) moves the unpolarised RATE and (ii) cannot.** σ_pb
   **−0.684681 %** (591846.161405217 → 587793.902124821 pb) on ⁶Li inclusive,
   config 1, x_max 0.95, `tensor-thirds` (0.7, 0.6), 2000 events, seed 7 —
   because the kernel's `r_func` reaches F₁, F_L, D(y) and `ToyG1`, which is
   the whole point of sharing it. `--b1-unpol`, by contrast, leaves the
   spin-blind cell cross section bit for bit.
6. **The mismatch factor is a Q² = 2.5 statement.** At x = 0.30, Q² = 5,
   `r1998` = 0.130957 is *below* `r_sigma_lt` = 0.163636 and
   (1 + r1998)/(1 + r_sigma_lt) = **0.971916**, so the sign of the whole
   effect flips against this section's 1.088 / 1.039.

So the row does not close by choosing (iii): **(iii) still requires naming an
R**, and the question becomes *(iii) with `sigma-lt`, which reproduces today,
or (iii) with `r1998`, which costs the six lines above?*

**(e)** (i): close the row — and note that the tree now carries the *option*
(iii) whether or not the row takes it, at its no-op default. (ii):
`Li6ConvolutionOptions::r_func` (null ⇒ `r_sigma_lt`; the gate's null ⇒
`r1998`), the T-pins on the convolution, the `USAGE.md` §2a cost table.
**(iii): the wiring is DONE** (2026-09-06) — `RSource`,
`PipelineConfig::r_source`, the `default_inclusive_kernel` parameter, the
`validate()` refusal, the `knob_provenance` row, `meta["r_source"]`, the CLI
flag and banner line, `lipolgen.R_SOURCE` / `make_config(r_source=...)`, and
the gates (`tests/test_b1_nuclear.cpp` T18; `test_b1_model.py` P10;
three `test_knob_provenance.py` matrix cells). Taking (iii) as the DEFAULT is
what remains, and it is a one-line change of `PipelineConfig::r_source`'s
initialiser plus the four pointer sites, which have been rewritten to say
"opt-in" rather than "hard-locked" (`CONVENTIONS.md`, `USAGE.md` §2a,
`README.md`, `PHYSICS_CHANNELS.md` §3) — `phase_D_li7_rank2.md:430`, which
notes the choice "matters more" for ⁷Li, is a dated record and is left as
written.

### B8. The α–d source default stays `FitRescaled` (`STATUS.md` decision **row 16**, added 2026-09-05; `phase_C_numbers.md` §C3.5; `USAGE.md` §9; `OPEN_ITEMS_SOLUTIONS.md` §11 O2)

**(a)** Should `ClusterConfigSampler` draw the α–d radial functions from the
smoothed fit (`li6.adr.fit`, rescaled to the raw normalisation) or from the
raw VMC overlap (`li6.ad`, `OverlapRaw`)?

**(b)** (i) `FitRescaled` (status quo); (ii) `OverlapRaw`.

**(c)** `phase_C_numbers.md` §C3: the raw R₂ node near 1.07–1.1 fm is real
(3.3 σ) and irrelevant (r < 1.5 fm is −0.09 % of q_int); against GK's
η = −0.025(12) the two are −2.0 σ and −2.4 σ — no measurement prefers one;
the raw block is rough at ≈1 σ per point everywhere; the fit is biased −2.7 σ
(on MC error) in the 6–8 fm window in the direction that makes the D-wave
excess look smaller. The evidence favours (i) as the numerical object and
records its bias.

**(d)** (ii) moves Q_charge by +8.7 % and η by +11.5 %, inside a band already
a factor 7.5 wide; **no `validation/reference/*.json` depends on
`cluster_config`** (checked by the run), so the rtol gates do not move. What
moves is pinned documentation: design G §8/§2.7, T17, T22, T22b, T22c, the
`USAGE.md` table, the `PHYSICS_CHANNELS.md` row, the pytest anchors.

**(e)** (i): add the row to the table and close it. (ii):
`cluster_config.hpp:214` one line plus the pins above.

### B9. The tagged deuteron control keeps Hulthén at P_D = 0.045; AV18 is opt-in (`STATUS.md` row 9)

**(a)** Should the control channel's deuteron follow the same
`--cluster-wave vmc` switch that gives AV18 P_D = 0.057600, by default?

**(b)** (i) Hulthén default, AV18 opt-in (status quo); (ii) AV18 default.

**(c)** `phase_C_numbers.md` §C5.4: `fdeut.av18`'s u(k), w(k) are the p–n
waves (the earlier "no VMC d → p+n table" was false); the control's job is to
be the Cosyn–Weiss tagged limit of the analytic family the ⁶Li channel is
built from; `fdeut.av18` prints no MC errors, so the AV18 row carries no band;
`validation/reference/tagged.json` pins Hulthén. The evidence favours (i).

> **(c) HAS MOVED, 2026-09-06 — two of its four legs changed, and they pull
> opposite ways.** (1) *"The control's job is to be the Cosyn–Weiss tagged
> limit"*: the Cosyn–Weiss gate **now runs on AV18**, not on Hulthén
> (`../../benchmarking/07_cw_sign_investigation.md` §8). CW's TABLE II is
> quoted for AV18, and the Hulthén pair's f₂/f₀ **never reaches √2** anywhere
> on the grid, so it cannot carry CW's landmarks — the old gate's agreement
> with "k = 0.30 GeV" was a coincidence of β = 0.30. That leg now argues
> **for** (ii), or at least says the control is already AV18 where it matters.
> (2) *"`tagged.json` pins Hulthén"*: still true, but the two spin-1 `model`
> blocks were **deliberately re-pinned on 2026-09-06** from this library
> rather than from `polligen` (whose `_amp2_table` carries the S–D sign bug),
> so "the reference pins it" is now a statement about a file this project
> generates, not an external port gate. That leg is **weaker** than it was. (2026-09-23: back to an external port gate — `tagged.json` is again `polligen`'s in every block. The row stays open.)
> Legs (3) *no MC errors* and (4) *the analytic family* are untouched. **The
> row stays open**; what changed is the evidence, not a decision.

**(d)** (ii): P_D +28 %, **−2.03 %** / **−1.18 %** on the channel's vector /
tensor dilutions; the relative S–D sign does not flip; `tagged.json` moves at
rtol 1e−12. *(**Re-measured 2026-09-06** after the S–D interference fix: the
two percentages are **unchanged at the three digits row 9 prints** —
−2.03 % vector and −1.18 % tensor. At six digits the **vector** leg moves in
the last place, −2.026820 % → −2.026832 %, and the tensor leg not at all,
−1.182136 % → −1.182144 % (−1.18214 % at five decimals either way) — because both
dilutions are angle-integrated and the S–D cross term drops out of them by
L-orthogonality. And *"the relative S–D sign does not flip"* is still true of
the **stored** ψ₂ = +W, which is what it was about; what changed is that the
**amplitude** now applies φ_L = i^L ψ_L, which it did not before.)*

**(e)** (i): close the row. (ii): the deuteron-channel default in
`tagged.hpp` (`DEUTERON_AV18()`) / `breakup.hpp` (`BreakupOptions::source`),
re-dump `tagged.json`, T25.

### B10. The embedded deuteron follows `--cluster-wave`, and the flag is refused where it is never read — *applied, confirm or revert* (`STATUS.md` row 11)

**(a)** Confirm the run's choice among follow / refuse / document for a
`--cluster-wave vmc` α-tag run whose embedded deuteron used to stay Hulthén.

**(b)** (i) follow (applied); (ii) refuse the combination; (iii) document only.

**(c)** `phase_C_numbers.md` §C5.5b: §C5.4's own rationale forbids two
wave-function families in one channel; `DEUTERON_AV18()` and
`BreakupOptions::source` both take the flag, so it means one deuteron
everywhere it is read.

**(d)** The defect was **+2.069 %** on every polarized tagged-α observable
under the flag (exact, 1.020687209533 at every (x, Q²) — *re-measured
2026-09-06 after the S–D interference fix as the ratio of the two channels'
`vector_dilution`: **1.020687624664** against **1.020687500612** before, i.e.
+2.069 % on both, unmoved in the seventh figure, because the dilutions are
angle-integrated*); the repair moved the
opt-in path by **−2.027 %**; the Hulthén default is bit for bit (0 reference
JSONs, 0 pre-existing assertions moved; +45 assertions added). A `vmc` run's
whole-nucleus reading is 0.887076 (0.886169 ± 0.0011 end to end on 400 k
events). `validate()` now refuses `cluster_wave` on `inclusive`/`coherent`.

**(e)** Confirm: close the row. Revert: `tagged.cpp`/`breakup.cpp` back to the
unconditional `DEUTERON()`, T27 and the two pytests removed, the five
`CONVENTIONS.md`/`USAGE.md`/`PHYSICS_CHANNELS.md` sites restated.

### B11. The published `A_zz(Born)` column is wrong by ×3.253983 — **republished BESIDE the original on 2026-09-06** (`STATUS.md` row 6; `../run_2026-09-02/phase_C_numbers.md` §8.1c-corr / §8.2-corr / §8.3-corr)

**(a)** Republish `run_2026-09-02/phase_C_numbers.md` §8.1c/§8.2's A_zz(Born)
column and everything that inherits it, or leave the note?

**(b)** (i) leave, with row 6's note; (ii) republish.

**(c)** `phase_B_numbers.md` §B3.2: the column reproduces
`azz(toy_b1(x, q2, f1, B1Mode::Digitized), …)` — the deuteron Miller table —
whereas the shipped kernel has used `Li6B1(MillerB1)` since before the RC
commit. Row 6 gives the file as `phase_C_numbers.md` without its directory;
it is the **2026-09-02** run's file (this run's `phase_C_numbers.md` has no
§8.1c).

**(d)** Measured at Q² = 5: published −1.473482e−03 / −4.975091e−03 /
−1.585230e−04 against this generator's −4.528242e−04 / −1.528923e−03 /
−4.871661e−05 at x = 0.01 / 0.10 / 0.30. The τ, w_hi − 1, ΔA_zz and "band
half-width 4.4204e−04" numbers inherit it; σ^el_T, σ^q_U, r_U do not. No code
or reference depends on the column. **Coupled to B2**: if the Miller 0.5 is
dropped the corrected column doubles again.

**REPUBLISHED BESIDE THE ORIGINAL ON 2026-09-06.** The thing that takes no
decision away was done: the original tables **stay**, labelled with what they
were computed with, and the recomputation is a dated table next to each —
`../run_2026-09-02/phase_C_numbers.md` **§8.1c-corr** (all eight (x, Q²)
points × both C0 edges), **§8.2-corr** (the five-row band, with τ, w_hi − 1
and the half-width) and **§8.3-corr** (the six run-level rows that inherit the
column). `OPEN_ITEMS_SOLUTIONS.md` §9 carries the same pairs at the same date.
So option (i) "leave with the note" and option (ii) "republish" are **both
still open**, and (ii) now costs nothing but a deletion.

**The correction is one constant on `A_zz` and is not one constant on anything
derived from it.** `published ÷ shipped = 3.253983147` at every point, exactly
`1 / (LI6_B1_RANK2_TRANSFER × LI6_B1_PER_NUCLEON)` = `1/(0.921947 × 2/6)` —
the two factors `Li6B1` applies and `toy_b1` does not — and it is exact
because the shipped `b2` is the Callan–Gross `2·x·b₁` (measured
`tables().b2 == 2*x*tables().b1` as doubles), which makes `azz` linear in b₁.
`ΔA_zz = [2 r_T − A_zz r_U]/(1 + r_U)` is **not** proportional to it.

**Four things the recomputation established that the note did not.**
1. **The sign flip at x = 0.01 is six entries, not one.** `ΔA_zz` changes sign
   on **both** C0 edges at **all three** Q² of §8.1c; the 2026-09-04 note
   recorded `(0.01, 5)` `ho` alone. At x = 0.10 and 0.30 the sign survives and
   the magnitude falls by ×0.32 and ×0.41 (`ho`) — not by 1/3.253983, because
   only part of ΔA_zz carries `A_zz`.
2. **Which ΔA_zz spreads move is decided by one rule.** A knob that leaves
   `r_U` alone has an `A_zz`-independent spread: `fq_scale` (8.040830e−07) and
   `tail_tensor_scale` (3.615e−08 against the published 3.608e−08, +0.21 %, i.e. the third significant figure) keeps its spread, while
   `qe_suppression` 1.680e−07 → **5.167e−08** and `qe_kf_gev` 1.878e−07 →
   **5.776e−08** shrink by ×0.31 and `c0_shape` 2.877e−08 → **3.312e−08**
   *widens* by 15 %.
3. **The RC budget's ORDER changes at second place** — and for a bookkeeping
   reason: §9's ladder compared `qe_tensor_scale`, which was corrected when it
   was added, against four rows that were not. All six on one footing:
   `fq_scale` 8.041e−07 ≫ **`qe_tensor_scale` 1.162e−07** > `qe_kf_gev`
   5.776e−08 > `qe_suppression` 5.167e−08 > `tail_tensor_scale` 3.615e−08 >
   `c0_shape` 3.312e−08. `qe_tensor_scale` rises from fourth to second because
   the three below it shrank.
4. **Every conclusion that does not contain `A_zz` is unmoved, and was
   re-measured rather than assumed.** §9's conclusions 1 (the σ^el_T sign
   change in x) and 2 (the σ^q_U/σ^el_U ladder) re-measure identical to every
   printed digit. Conclusion 3 — "the tail is not the leading RC systematic" —
   depends on the column in both of its numbers and **stands**: ×767.6
   (1.358472e−04 ÷ 1.769797e−07), the 768 already on record. §8.3's
   *"dominant … by a factor 300"* does not reproduce as any well-defined ratio
   even on its own column (549.7 or 201.6, depending which two numbers are
   meant) and is restated there as **168.9×** the largest knob spread on the
   half-width and **61.9×** on the 0.19 ↔ 0.30 spread — the band still
   dominates on both readings.

**One reproduction exception, recorded rather than smoothed.** Every `ho` row
of §8.1c and every row of §8.2 reproduces to every printed digit. **Three of
the eight `vmc-ft` ΔA_zz entries do not**: +7.9026e−08 / +3.2291e−07 /
+9.7195e−07 / +2.5575e−08 published against +7.902804e−08 / +3.229021e−07 /
+9.720112e−07 / +2.557553e−08 re-measured — 2.4 × 10⁻⁵ to 6.3 × 10⁻⁵ relative,
in the fifth significant figure. Deterministic in this tree (two processes,
bit-identical; explicit `--rc-c0-shape ho` bit-identical to the default), so it
is a difference against the 2026-09-04 build and **its cause is not
established**. It is four to five orders of magnitude below every conclusion
drawn from the column, and the `ho` edge — which carries all of them — is
exact.

**Cost of each option, now that both exist.** (ii) costs **zero shipped
numbers**: no source file changed, no `validation/reference/*.json` carries an
`A_zz` or an `rc_*` field, no assertion moves, and the recomputation is
documentation only. **B2 still gates the VALUE, not the publication**: if the
Miller 0.5 is dropped the corrected column doubles again, so a reader must
take the corrected column with B2's answer — which is exactly why both columns
are printed side by side rather than one replacing the other.

**(e)** (i): delete nothing — the correction tables already stand beside the
originals and are labelled as corrections. (ii): delete the superseded columns
from `phase_C_numbers.md` §8.1c/§8.2/§8.3 and `OPEN_ITEMS_SOLUTIONS.md` §9,
after B2 is answered.

### B12. `CdksB1` stops halving the CDKS column, `B1_CDKS_TABLE_TO_PER_NUCLEON` = 1 — *applied, confirm or revert* (`STATUS.md` row 1; §10 condition 5)

**(a)** Confirm that the digitized CDKS b₁ column is per nucleon and is no
longer converted.

**(b)** (i) 1.0 (applied); (ii) revert to 0.5.

**(c)** CDKS Eq. (10) carries an explicit 1/A; the text under their Eq. (16)
says "b₁ is defined by the one per nucleon"; f(y) is normalised to one
nucleon; F₁ᴺ = (F₁ᵖ + F₁ⁿ)/2; HERMES's published b₁ᵈ is per nucleon by their
Eq. (5), confirmed by inverting Table II in all six bins (mean ratio 0.946
against 0.473 for per deuteron). The run calls this **certain**.

**(d)** `--b1-model cdks` and `b1_convolution` doubled: `tests/test_sf.cpp`
pins −2.816198975086724e−04 → −5.632397950173448e−04 (x = 0.3) and
5.8517162680014353e−05 → 1.1703432536002871e−04 (x = 0.05). No reference JSON
moved (the only b₁ reference is `"b1_model": "miller"`); the default is
untouched.

**(e)** Confirm: close the row. Revert: `constants.hpp:120` and the two pins.

### B13. The RC band's low-x anchor: `RC_DELTA_LOW_X` = 0.30 at `RC_X_LOW` = 0.01, and what it is anchored to — **priced 2026-09-06** (`STATUS.md` decision **row 17**, added 2026-09-05; `rc.hpp` `RC_DELTA_LOW_X`; `phase_B_numbers.md` §B6; `../run_2026-09-06/phase_A_numbers.md` §A2)

**(a)** Keep the tensor-RC band's low-x anchor at δ = 0.30, x = 0.01, now that
the run has shown its basis is "the Q² = 0.1 panel's lowest-x value, carried
upward in x to 0.01" and not "the conservative end of a Q² spread"?

**(b)** (i) keep 0.30 at 0.01 (status quo, basis restated); (ii) re-site to
the panel's own value at its top, |δ| = 0.113 at x = 0.00966; (iii) 0.266 at
x = 0.00226 (the value that rounds to the paper's "30 %").

**(c)** `rc.hpp:277-288`, `phase_B_numbers.md` §B6: 0.01 lies above the only
panel the 10 %/30 % sentence covers (x = 0.00226–0.00966 at Q² = 0.1); all
three readings are conservative in magnitude (0.30 > 0.266 > 0.113) and
magnitude is the only property w = 1 ∓ δτ uses. "Re-siting an anchor after
seeing why it sat where it did is the author's call."

**(d)** **PRICED 2026-09-06** — the alternatives were run, not estimated
(`../run_2026-09-06/phase_A_numbers.md` §A2;
`python/tests/test_rc_low_x_anchor.py` pins the table). The shipped output is
still unaffected under every option, because the band is opt-in
(`--rc tensor-band`), `--rc off` is byte-identical, and **no
`validation/reference/*.json` carries an `rc_*` field at all** — so all three
options are free of the rtol-1e−12 gate.

**The three-row table.** Band half-width on `A_zz` = δ(x)·|A_zz|, at §9's own
configuration (⁶Li, `--config 1`, `--channel inclusive`,
`--plan tensor-thirds --pzz 0.6`, θ_S = 0, Q² = 5 GeV², P_zz = +1) and with
**this generator's own ⁶Li b₁** — `A_zz` = −4.528242e−04 / −1.316774e−03 /
−1.528923e−03 / −1.280460e−03 at x = 0.010 / 0.063 / 0.100 / 0.160, i.e. the
value corrected in `phase_B_numbers.md` §B3.2 and **not** §9's own
×3.253983 deuteron-b₁ column:

| δ_low | x = 0.010 | x = 0.063 | x = 0.100 | x = 0.160 |
|---|---|---|---|---|
| **(i) 0.30 — shipped** | **1.358472e−04** | **1.459067e−04** | **9.680016e−05** | 1.920691e−05 |
| (iii) 0.266 — the panel at x = 0.00226 | 1.204512e−04 | 1.308566e−04 | 8.798804e−05 | 1.920691e−05 |
| **(ii) 0.113 — the panel AT x = 0.00966** | **5.116913e−05** | **6.313127e−05** | **4.833349e−05** | 1.920691e−05 |
| ratio, (iii) ÷ (i) | ×0.886667 | ×0.896851 | ×0.908966 | ×1 |
| ratio, (ii) ÷ (i) | ×0.376667 | ×0.432682 | ×0.499312 | ×1 |

**0.113 is the value the panel actually reads at the x nearest the anchor**
(x = 0.00966; 0.01 lies above the panel's top), so the shipped 0.30 errs
**wide by ×2.65** — the safe direction for a half-width and the wrong one for
a quoted precision.

**(c)'s own arithmetic was right at one x and wrong as a rule.** The
"≈ ×0.89 / ≈ ×0.38" above reproduces **exactly** at x = 0.01 (0.886667,
0.376667 — the ratio of the anchors), because δ = δ_low is returned exactly at
x ≤ `RC_X_LOW`. It **overstates the reduction everywhere else**: the
log-linear interpolation carries only part of the move (×0.4327 and ×0.4993 at
x = 0.063 and 0.100 under (ii), against the ×0.3767 at the anchor), and at
x ≥ `RC_X_HIGH` = 0.16 the E12-13-011 anchor pins δ = 0.015 and the half-width
**does not move at all**. The whole decision therefore lives in x < 0.16.

**Four further facts the decision needs, each measured.**
1. **The ordering of the RC budget survives every option.** At x = 0.01,
   Q² = 5 the band leads the **whole** radiative tail (ΔA_zz = −1.769797e−07)
   by **×768 / ×681 / ×289** at δ_low = 0.30 / 0.266 / 0.113, and at the other
   three x by ×17 000–×45 500. Choosing (ii) does not make the tail the
   leading RC systematic.
2. **The band's peak does not move**, it sharpens. It is at x = 0.063 on all
   three anchors — 1.459067e−04 / 1.308566e−04 / 6.313127e−05 — and the
   0.063 : 0.010 ratio rises 1.074 → 1.086 → 1.234 as the anchor falls.
3. **What moves in a file**: exactly **two of the 51 columns**
   (`rc_tensor_lo`, `rc_tensor_hi`) and **two of the 60 `meta` keys**
   (59 at `cdd8591`, before B1 added `rc_sp_tensor_scale`)
   (`rc_delta_low_x` and `knob_provenance`, which records it) on a 2000-event
   ⁶Li inclusive run at seed 99 — no kinematics, no `weight`, no `rc_tail`.
   Mean |w_hi − 1| **6.153874e−05 → 5.484602e−05 → 2.472882e−05**.
4. **On the tagged band, where the clamp bites, the anchor changes the clipped
   fraction by NOTHING** — not "little": the **same 245** of 20 000 ⁶Li
   tagged-α events (1.2250 %, seed 1) and the **same 520** of 20 000 ⁷Li
   (2.6 %, seed 11) clip at all three, and every non-band column is
   byte-identical. `RcModel::clamp_tau` clips |τ| against `band_tau_max` and
   δ never enters it. What the anchor sets **at** the clip is the width there:
   |τ| = 1 exactly, so those events' edges are exactly 1 ∓ δ_low —
   **[0.700, 1.300] / [0.734, 1.266] / [0.887, 1.113]** — and the mean
   |w_hi − 1| over the ⁶Li 20 k falls 3.384715e−02 → 3.006184e−02 →
   1.302792e−02. (Two bookkeeping notes from the same runs: the ⁶Li
   **1.2250 %** is **seed 1** and sits mid-scatter in a ten-seed scan spanning
   1.0950–1.3650 %, mean 1.2485 %, sd 0.0923 — the ⁷Li 2.6 % is seed 11, and
   seed 1 there gives 2.51 %.)

   > **The ⁶Li figures in this item were RE-MEASURED 2026-09-15.** As first
   > written they were **124 = 0.6200 %**, a ten-seed span 0.4850–0.6200 %
   > (mean 0.5280 %, sd 0.0411) and mean |w_hi − 1| 2.587195e−02 →
   > 2.297940e−02 → 9.962921e−03 — all measured in Phase A at `cdd8591`,
   > **before** the tagged S–D interference fix `a7b3d18`, which moves the
   > spin-1 (M, cos θ_k) draw `RcModel::tagged_tau` reads. The ⁷Li figures are
   > J = 3/2 and reproduce unchanged. **The structural claim this item makes —
   > that the anchor moves the clipped fraction by nothing — is unaffected and
   > was re-verified on the fixed build at all three anchors.** The price this
   > item carries (the ×2.65 width at the clip) is likewise unchanged.

**The HIGH anchor: the record names no alternative, and this run does not
invent one.** `RC_DELTA_HIGH_X` = 0.015 at `RC_X_HIGH` = 0.16 is
E12-13-011's 1.5 %, which lives in the **unpublished proposal only** — the
published companion arXiv:2506.04506 contains no radiative-correction
discussion at all. Grepped across the tree: there is **no second reading, no
`_OPTIMISTIC` partner and no band edge** for it anywhere; the only alternative
any document states is an *action* — "cite it by page for the 1.5 %, or drop
the anchor" (`design_C_tensor_rc.md`, `CONVENTIONS.md` (a)) — which is not a
value, so there is nothing to run. That is stated rather than estimated, and
it is why the x ≥ 0.16 column above is a single number.

**(e)** (i): add the row to the table and close it. (ii)/(iii): `rc.hpp` at
`RC_DELTA_LOW_X` and `RC_X_LOW`, T3's `==` anchors (which are **symbolic** —
`rc_delta(RC_X_LOW) == RC_DELTA_LOW_X` — so they follow the constant and none
of them is a hard-coded 0.30), the pinned table and its four assertions in
`python/tests/test_rc_low_x_anchor.py`, the band tables in `USAGE.md`,
`OPEN_ITEMS_SOLUTIONS.md` §9, `PHYSICS_CHANNELS.md`, `CONVENTIONS.md` (b),
`design_C_tensor_rc.md` Q7, and the `--rc-delta-low-x` help + banner reason.
**No reference JSON moves and no existing assertion breaks**: every gate on
the anchor is symbolic (`rc_delta(RC_X_LOW) == RC_DELTA_LOW_X`,
`|w − 1| ≤ RC_DELTA_LOW_X`, `worst_lo ≈ 1 − RC_DELTA_LOW_X`), and the only
literal `0.30`s in the suites are a bad-anchor `CHECK_THROWS` argument
(`tests/test_rc.cpp`) and an explicitly-passed `rc_delta_low_x=0.30`
construction (`python/tests/test_rc.py`), neither of which pins the default.
So the blast radius is documentation plus the one pinned table above.

### B14. The Gakh–Shekhovtsova δ(x) shape is rejected on the record (`STATUS.md` row 7)

**(a)** Confirm that the log-linear interpolation between two anchors stays,
and the paper's own δ(x) shape is not digitised and adopted.

**(b)** (i) keep log-linear (applied); (ii) adopt the shape.

**(c)** `phase_B_numbers.md` §B6.3, from the paper's arXiv EPS:
δ ≡ (Δσ_RC − Δσ_Born)/Δσ_Born is singular across the b₁ zero-crossing the RC
itself moves (x₀ 0.20118 → 0.18467 at Q² = 4), spanning −1.691 … +7.810 there
and −0.887 … +4.091 at Q² = 10; not monotone (−0.400 at x = 0.85 against
δ_high = 0.015); no panel covers x = 0.01, Q² = 5; it is a ratio against
*their* Born.

**(d)** (ii) would emit negative weights from w = 1 ± δτ and falsify T3 by
construction; no shipped number depends on (i) (the band is opt-in). By-product
already flagged: the "10–30 %" magnitude is the Q² = 0.1 panel's alone,
quoted at Q² = 5.

**(e)** (i): close the row. (ii): a digitisation task plus a new `RcModel`
shape, T3 rewritten.

### B15. Opt-in physics left at its no-op default — promotion (phase-B items B1/B2/B3/B6 of `phase_B_numbers.md`; `phase_C_numbers.md` C5.2; `OPEN_ITEMS_SOLUTIONS.md` §14; **not claimed as author decisions by the run — listed so the batch is complete**)

**(a)** Should any of the terms the run added behind a knob, each argued to be
physically present, become the default?

**(b)** Per knob: `RcOptions::c0_shape` `Ho` → `VmcFt`; `tail_model` `TPeak`
→ `TPeakPlusLL` **or, since 2026-09-06, `TPeak` → `PolradFull`** (a third
value of the SAME knob, not a new one: POLRAD Eq. (18) + Appendix B +
Eq. (A.4), the exact τ_A tail — see (d)); `qe_tensor_scale` 0 → 1;
`a_transfer_frac` 0 → 0.5/1;
`cluster_vmc_mc_sigma` 0 → ±1; `pol_sf` `toy` → `nnpdfpol`. **Added
2026-09-06: `RcOptions::sp_tensor_scale` 0 → 1**, the direct sibling of
`qe_tensor_scale` one level down (the tensor fraction of the leading-log
s-/p-peaks), which is refused unless `tail_model = TPeakPlusLL` and therefore
cannot become a default without (b)'s second entry going first.

**(c)** The run's ground rule (new physics opt-in, default bit for bit) is the
only reason each sits at zero; each header argues the term is real.

**(d)** Measured by the run: C0 shape — the sign of (1/6)σ^el_T/σ^el_U at
x = 0.10 is **indeterminate** across the band (+1.5595e−04 is one edge, not a
result). s-/p-peaks — the "< 1 % at Q² ≥ 20" acceptance is met
event-weighted (+0.61 %) and fails per cell (331 of 1356 cells, worst ×6444
at x = 0.79, y = 0.0088). `qe_tensor_scale = 1` — Δ(ΔA_zz) 0.086 % / 0.0003 %
/ 0.25 % of the band at x = 0.01 / 0.10 / 0.30. `a_transfer_frac` — A_zz band
half-widths 1.358472e−04 / 1.518818e−04 / 1.921170e−04 at f = 0 / 0.5 / 1.
`cluster_vmc_mc_sigma` — 0.02 % on the tagged tensor observable. `pol_sf` —
read on no channel at the CLI default plan. The first four are inside `--rc`,
itself opt-in, so they move nothing at the CLI default whichever way they go.
**`sp_tensor_scale = 1`, priced 2026-09-06** (`../run_2026-09-06/phase_B_numbers.md`
§B1) — Δ(ΔA_zz) is **exactly 0** at x = 0.01 / 0.10 / 0.30, Q² = 5 for scale
0.5 **and** 1, because the coherent s/p vertex sits at Q′² = 4.37 / 4.94 /
4.98 GeV² where ⁶Li's form factor is 45–51 decades down, so scale 1 is
bit-identical to 0 there; it is non-zero on **77 of 3051** accepted cells
(x ≤ 7.94e−03, y ≥ 0.366), reaching **582.9 %** of the band half-width at
x = 4.169e−04, y = 0.9692, and it recovers at most **21.72 %** of the
tensor-fraction collapse ×0.66139 / ×0.0031829 / ×6.6076e−05 that
`tail_model = TPeakPlusLL` causes at those same three x. So this entry is a
price tag whose price is **zero where the run is quoted** — and the piece that
carries that collapse, the **quasi-elastic s/p column**, has no knob at all.

**`tail_model = PolradFull`, priced 2026-09-06** (`../run_2026-09-06/phase_B_numbers.md`
§B2). This is the entry that changed most, because there is now a **computed**
third answer where (b) previously offered only a lower bound and a stated
model. Measured: **it is not inside the pair** — over the sampler's 3051
accepted cells it lies between `TPeak` and `TPeakPlusLL` on **1725 (56.5 %)**
and outside on **1326 (43.5 %)**, covers a **median 0.6555** of the gap **over
the 3027 cells whose gap is nonzero** (0.6620 over all 3051; on the other 24
`TPeakPlusLL` equals `TPeak` exactly) and the **ratio of the σ-weighted mean
shifts** is **0.428** — a ratio of means, *not* a σ-weighted mean of the
per-cell fractions, which is 6.483 — and in the Q² ≥ 20, y ≤ 0.9 window it is
**above both**. Whole-run mean `rc_tail` (200 k, seed 1234) **at the CLI's
default fill P_z = 0.7**: **1.021778527** (`TPeak`,
the default) → **1.040274188** (`TPeakPlusLL`) → **1.029702912**
(`PolradFull`); at **P_z = 0** (the plan both test suites use, and the plan
this run's own §B2 published) **1.021836305 → 1.040368933 → 1.029775347** on
the same build. Δ(ΔA_zz) against `TPeak` is **0.0435 % / 0.3073 / 3.254 %**
of the band half-width at x = 0.01 / 0.10 / 0.30 — about **half** the
`TPeakPlusLL` edge's. It also **retires the reason (b)'s `sp_tensor_scale`
entry exists**: it computes the s-/p-peaks' Eq. (A.4) tensor content instead of
bounding it, and scored against it the bound was not even one-sided
(`r_T/r_U` × the t-peak: ×0.79395 / ×0.0020032 / ×0.00022591 computed against
×0.66139 / ×0.0031829 / ×6.6076e−05 bounded). **What it does NOT change**: the
quasi-elastic tensor tail is still exactly zero on all three models, because
Eq. (A.5) has `Im₅…₈ ≡ 0` for a spin-½ target — so this entry does not close
the gap the paragraph above names, and `qe_tensor_scale` is still the only
stand-in for it. **Not validated against Mo–Tsai**, which is not in this tree;
checked POLRAD-internally and against the leading-log fallback only. It costs
~10 s of tail-table build and needs `n_eta ≥ 64`.

> **CORRECTION, 2026-09-15 (residue pass) — this paragraph published an
> event-weighted triple with no P_z, and a mislabelled 0.428.** Both fixed
> above, re-measured on the working-tree build: the triple is
> **1.021778527 / 1.040274188 / 1.029702912** at `tensor_thirds_plan(0.7, 0.6)`
> (the CLI default) and **1.021836305 / 1.040368933 / 1.029775347** at
> `tensor_thirds_plan(0.0, 0.6)`; **0.428** is the RATIO of the σ-weighted mean
> shifts (0.427958 unclipped, 0.427368 at the shipped `tail_max = 10`), not a
> σ-weighted mean of the per-cell coverage fractions, which is **6.483**; and
> the **median 0.6555** is over the **3027** cells with a nonzero
> `TPeak → TPeakPlusLL` gap (**0.6620** over all 3051). The cell-σ-weighted
> census is P_z-free — measured identical on both plans. The decision this
> paragraph asks for is unchanged: no default moved.

> **2026-09-23 — the ⁶Li C0 zero priced against data (`OPEN_ITEMS_SOLUTIONS.md` open item Q1; filed here because §B15's `c0_shape` is its nearest relative — NO registry row names `LI6_FF_HO_*` or Q1, measured by `grep`, and none is created here).** `validation/benchmarks/t3_li6_charge_ff_fb.py` (`../run_2026-09-23/phase_B2_nuclear.md` §3): the shipped `HoSpin1FF` C0 zero q₀ = **3.0998 fm⁻¹** lies **+0.2713 fm⁻¹ (+9.6 %)** above the band [2.6944 (UVa Fourier–Bessel zero), 2.8284 (Li *et al.* 1971 \|F_L\|² minimum)] fm⁻¹ — a recorded FAIL; T11's [2.9, 3.3] window does not reach the band (gap 0.072 fm⁻¹) and was not moved. Inside the model the C2 fill-in cannot close it (C2 shares the C0 monopole, so the model's \|F_L\|² minimum is q₀ to 1.7e−11); only Coulomb distortion, not computed, remains. **What adopting the FB shape would move** (a scratch `Spin1ElasticFF` subclass, never in the tree; F_q scaled with F_c, F_m shipped; every `ho` entry reproduces §B1.5 to every printed digit), Q² = 5: (1/6)σ^el_T/σ^el_U at x = 0.01 / 0.03 / 0.10 / 0.30 goes **−5.094e−04 / −7.358e−04 / +1.560e−04 / +1.329e−02** (`ho`) → **−5.613e−04 / −8.263e−04 / −8.167e−05 / +1.574e−02** (FB), against `vmc-ft`'s −5.470e−04 / −8.043e−04 / −4.070e−05 / +1.556e−02 — **the data-derived shape lies OUTSIDE the [ho, vmc-ft] band at every one of these points**, and it lands on the NEGATIVE side of the x = 0.10 sign question (it does not settle it: the FB row's provenance is UNVERIFIED). σ^el_U moves ×1.019 / ×1.038 / ×1.152 / ×3.532. ΔA_zz from the whole tail moves −1.770e−07 → −2.241e−07 (x = 0.01), +2.127e−09 → +2.033e−09 (0.10), +5.130e−12 → +1.093e−11 (0.30) — at most 2.2e−07 absolute, three orders below the 4.4e−04 band half-width, so the tail stays *not the leading RC systematic*. **As an option for this section's `c0_shape` (b):** a third C0 edge `fb` (the UVa row, R = 6.0 fm) — costs a new `Spin1ElasticFF` class, a vendored-row provenance still unverified, and every `--rc` tensor-fraction number at x ≥ 0.1 moving outside today's band; nothing at the CLI default moves (`--rc` is opt-in). **Not applied; no label moved.**

**(e)** Per knob: its default line in `rc.hpp:431/1070/1079/1096/1189`
(those five line numbers are as of the 2026-09-03 tree and have since shifted
— `rc.hpp` grew from **1581 lines at `cdd8591`** to 1788 in the working tree of that day over the 2026-09-06 B3 and B1 edits, and stands at **2154 as committed at `af3f415`** (`wc -l`, re-measured 2026-09-15; "HEAD to working tree" named no commit and went false the same week) —
and they are left as recorded rather than silently renumbered),
`pipeline.hpp`, the `USAGE.md` §"knobs" table, and the shipped npz for `--rc`
runs. `sp_tensor_scale`'s own default line is `rc.hpp:1452`, its numerator
term `rc.cpp:1524` and its refusal `rc.cpp:1010`, all as of 2026-09-06
**before** the B2 edits of the same day, which grew `rc.cpp` by ~470 lines and
moved every one of them — left as recorded rather than silently renumbered,
the same convention this paragraph already uses. `PolradFull`'s own entry
points are `polrad_full_sigma_el` / `polrad_full_sigma_qe_u` /
`polrad_im_el_spin1` in `rc.hpp` and `rc.cpp`, and its `meta` value is
`rc_tail_model = "polrad-full"` on the row §B15 already owns.

### B16. The coherent |t| ceiling stays 0.2 GeV², justified by the anchor range, not by positivity — *applied, confirm* (`STATUS.md` row 12)

**(a)** Confirm the reversed dependency: the ceiling's primary reason is the
digitised anchor range (|t| ≤ 0.30, knob-independent); the positivity edge is
the contingent second reason and is derived, not the derivation.

**(b)** (i) as applied; (ii) restore "positivity-derived", which would license
extrapolating a linear-in-|t| fit to 2.80 GeV² at the measured quadrupole.

**(c)** `phase_D_numbers.md` §D5, `phase_D_small_items.md` §D5.5: the two
reasons agree only at the shipped `eps_b0`; edge = 2(1/|P_zz| − amp)/(|eps_b0|·B)
= 0.2450 at −0.08, 0.3719 / 1.1462 / 2.80 at the band (2.7990 on the derived
−0.0070024; "quote 2.80, never 2.8000").

**(d)** Nothing numerical: on 200 000 default coherent events ⟨|t|⟩ =
0.019963, max |t| = 0.197575, and `sample_t` renormalises on [0, t_max], so
the truncation redistributes exp(−B t_max) = 4.5e−5 of the rate.

**(e)** (i): close the row and fix the un-amended site in
`OPEN_ITEMS_SOLUTIONS.md` §11 (Part B, finding 2). (ii): `coherent.hpp:60-90`
prose only.

### B17. The A = 2 gate's wave-function default stays `kFdeutFile` (AV18) although CD-Bonn + MSTW reproduces CDKS on every landmark (`STATUS.md` row 5)

**(a)** Should `DeuteronConvolutionB1::Options::wave` default to the tabulated
AV18 file (status quo) or to the analytic CD-Bonn, CDKS's own wave function?

**(b)** (i) `kFdeutFile`; (ii) `kCdBonn`.

**(c)** `phase_A_numbers.md` §8, `phase_A_cdbonn.md`: CD-Bonn + MSTW + Eq. (21)
gives 1.000338 with zeros 0.0641 / 0.4570 — a residual below the figure's
digitisation error, specific to CD-Bonn *and* MSTW together, and "not
three-digit agreement". Flipping is a one-line change plus a re-pin "and
belongs to a task that says so".

**(d)** No shipped number: the object is the gate, refused as a beam species
(`pipeline.cpp:778-786`). (ii) moves every pinned gate number in layer 3 and
checklist item 4 (the verdict row 0.843243 → 1.000338) and every published row
of `phase_A_numbers.md` §§0–7, `USAGE.md:365-370`, `README.md:152-164`.

**(e)** (i): close the row. (ii): `b1_nuclear.hpp:659`, `tests/test_b1_nuclear.cpp`
pins, the three documents.

### B18. `--pzz` on `helicity-flip` — not read at the default fill, honoured at the opt-in `--pzz-mode typed`; **priced 2026-09-06** (`STATUS.md` decision **row 20**, added 2026-09-05; `OPEN_ITEMS_SOLUTIONS.md` §15.5 **D13**; `phase_D_numbers.md` §D1.5 F5)

**(a)** Honour `--pzz` on the `helicity-flip` plan, refuse it there, or leave
it documented?

**(b)** (i) leave documented (applied: the banner prints the fill's own
moments, "--pzz is not read by this plan"); (ii) refuse `--pzz` with that
plan; (iii) honour it — **built opt-in and priced 2026-09-06** as
`--pzz-mode {ladder,typed}`, default `ladder`, which is (i) bit for bit.

**(c)** `phase_D_numbers.md:238-256`: `--plan helicity-flip --pzz 0.6` at
P_z = 0.7 fills T = 0.4 from the max-entropy ladder; 0.6 is outside the
plan's domain by 0.02 (edge T = 0.58 at P_z = 0.7, where p(−½) = 0). The run
moved no fill "because that would change shipped numbers".

**(d)** (iii) moves a shipped fill (T 0.4 → the typed value, or a refusal at
0.6). (i)/(ii) move nothing.

**PRICED 2026-09-06** (`../run_2026-09-06/phase_A_numbers.md` §A3), at the
standard configuration — inclusive, beam config 1, 100 000 events, seed
20260713, `--plan helicity-flip --pz 0.7 --pe 0.7`, the shipped
`--b1-model miller` — with `--pzz 0.5`, a value inside **both** domains, so
the ladder and the typed fill can be compared at all:

* **On ⁷Li, the isotope this row is about, (iii) moves NO observable.** The
  rank-2 sector of ⁷Li is identically zero and both fills honour `--pz`, so
  the only fill moments the kernel sums — Σ p_m = 1 and ⟨J_z⟩/J = P_z — are
  the same numbers: σ agrees to **1 ulp (1.970 × 10⁻¹⁶ relative)**,
  σ[apar+] to 3.939 × 10⁻¹⁶ and σ[apar−] **exactly**, and A_∥ from the rate
  moves +4.019 × 10⁻¹⁶ absolute = **6.2 × 10⁻¹⁴ of its own statistical
  error**. The whole difference is arithmetic. It is *not* nothing, though:
  **0.1060 % of the 100 000 events (106) land in a different (x, Q²) cell**,
  because cell selection is a discrete function of weights that differ in the
  last bit — so the file is not bit-identical and a reference JSON built
  under one mode would not reproduce under the other.
* **On ⁶Li (J = 1) the rank-2 sector is live and (iii) does move the rate.**
  σ 591783.2520093301 → 591769.3290332475 pb = **−0.0023527 %**; per category
  −0.0023541 % (apar+) and −0.0023514 % (apar−); A_∥ from the rate
  −0.00117171560617921 → −0.00117174317396206 = **−4.27 × 10⁻⁶ σ_stat**, so a
  1σ A_∥ shift would need N ≈ 5.5 × 10¹⁵ events. 1.2010 % of the events change
  cell. **The size is the b₁ model's, not the window's**: at `--x-max 0.95`
  the same run gives −0.0023527 %, while `--b1-model li6-convolution` there
  gives **−0.00010719 %** and `cdks` **−0.00006617 %** — a factor 22 to 36
  below `miller`.
* **What moves on BOTH isotopes is the recorded alignment**, and it is the
  divisor of every tensor estimator: T 0.4 → 0.5 (+25 %) at J = 3/2 and
  P_zz 0.409403 → 0.5 (+22.1292 %) at J = 1, i.e. δ(A_zz) and δ(cos 2φ) at
  fixed N move **−20 %** and **−18.1195 %**. The vector moment does not move
  at all — both fills honour `--pz` — so A_∥'s own divisor P_e P_z is
  untouched, which is why the A_∥ effect above is a rate effect and not a
  normalisation one.
* **The refusal, shown**: `--plan helicity-flip --pzz 0.6` — the registry's
  own line — at `--pz 0.7` is REFUSED at J = 3/2 with the edge named
  ("0.26 <= t <= 0.58 at this pz", p(−1/2) = −0.005), and ACCEPTED at J = 1,
  where 0.6 is inside the wider spin-1 domain 0.1 ≤ P_zz ≤ 1 and costs
  σ −0.0049496 %, A_∥ −8.99 × 10⁻⁶ σ_stat, 2.541 % of events re-celled.
  Nothing is clamped.

**So the row's own sentence was right and is now quantified**: on ⁷Li (iii)
"moves a shipped fill" and no observable — 1 ulp — while it moves 106 events
of 100 000; on ⁶Li it moves the rate by 2.4 × 10⁻⁵ of itself. What (iii)
really buys is the ALIGNMENT ITSELF, which is the thing a tensor programme
divides by, and what it really costs is that `--pzz 0.6` at `--pz 0.7` stops
running on ⁷Li instead of quietly filling 0.4.

**(e)** (i): add the row to the table and close it. (ii): `make_plan` /
`cli.py` validation. (iii): DONE opt-in 2026-09-06 — `--pzz-mode`
(`cli.py:99`), `PZZ_MODES` and `make_plan(pzz_mode=...)`
(`python/lipolgen/__init__.py:401`, `:769`, the ONE place the string becomes
`HelicityFlipOptions::use_explicit_pzz`), `KnobRunContext::pzz_mode`
(`pipeline.hpp:1283`) with its `knob_provenance` row (`pipeline.cpp:2407`),
the spin-1 refusal taught to name its edges (`spin.cpp` `spin1_populations`),
`tests/test_bookkeeping.cpp` and `tests/test_spin.cpp` (2 cases, 40
assertions), `python/tests/test_pzz_mode.py` (25 tests) and two
`test_knob_provenance.py` matrix variants (24 cells). **Which references
carry a `helicity-flip` fill, measured 2026-09-06**: exactly one of the nine,
`validation/reference/bookkeeping.json`, and it pins **five** of them at
rtol 1e-12 — `helicity_flip_j12`, `helicity_flip_j1_maxent_anchor`,
`helicity_flip_j32_maxent_anchor` and `helicity_flip_j1_tilted_offset` from
the LADDER branch and `helicity_flip_j1_explicit_pzz` already from
`use_explicit_pzz = true`; the other eight reference files carry no fill at
all (`grep -c "apar\|helicity"` returns 0 on each, 26 on that one and 2 on
`_manifest.json`, which only names two `SpinCategory` rows). The opt-in
default leaves all nine byte-identical (`git diff --stat
validation/reference/` empty) and every gate case passing. **Making `typed`
the default would move those four ladder rows**, and with them every
published `helicity-flip` number; this run did not make it the default.

### B19. `PDF:PomSet` = 11 is refused — *applied, confirm* (`STATUS.md` decision **row 19**, added 2026-09-05; `src/pythia/pythia_bridge.cpp`; `STATUS.md` phase D row)

**(a)** Confirm that a previously accepted input is now an error.

**(b)** (i) refuse (applied); (ii) accept with a warning; (iii) accept
silently (the pre-run state).

**(c)** Measured, 20 000 coherent events at ⁶Li config 1: `PomHISASD` returns
densities only after `setXPom`, which the bridge never calls —
`n_pom_flavour_fallback` = 20000/20000, so `meta["pom_set"] = 11` would record
a knob that did not run. Same rule as `validate()`'s other refusals.

**(d)** No default moves (`pom_set` = 6). The D4 band is stated over the
twelve genuine DPDF fits about set 6: ⟨n_ch⟩ −2.6 % / +9.9 %, kaons ×2.8; the
T0 columns are bit-identical across the fourteen sets that run.

**(e)** Confirm: close. (ii)/(iii): `pythia_bridge.cpp:853-861`, the
`pythia_bridge.hpp:252` and `cli.py:485/1126` sentences, one pytest.

### B20. The deferred ⁷Li b₁ design — D1, D3–D12 (`STATUS.md` decision **row 20**, added 2026-09-05; `OPEN_ITEMS_SOLUTIONS.md` §15.5; `phase_D_li7_rank2.md` §7)

**(a)** None of these moves a shipped number today — ⁷Li's rank-2 sector is
exactly zero and now loud — but each must be answered before any ⁷Li tensor
number exists. D2 is B3 above; D13 is B18.

| id | question | options | what the run measured | what changes when decided |
|---|---|---|---|---|
| D1 | which decomposition | α + t (S = 1.008, purely orbital) / ⁶Li + n (S = 0.682) / refuse | they disagree in sign and by up to 8×; non-orthogonal, so summing double counts (§5.5) | which `B1Model` gets written |
| D3 | normalisation target | S_αt = 1.0084 / 1 | 0.83 % — a provenance decision; S > 1 means the ⁶Li "non-α–d 18 % gets b₁ = 0" rule has no ⁷Li form | one normalisation line, stated |
| D4 | F₁ of the triton | isoscalar shortcut / true Z = 1, N = 2 (`triton_sf.hpp` exists) | 0.6 … 8.6 %; copying ⁶Li's `f1_alpha_ = f1_d_` would be a silent error | the convolution's F₁ input |
| D5 | κ: CDKS Eq. (17) or (21) | 17 / 21 | +0.35 … +69 % on ⁷Li against −1 … +7 % on ⁶Li; the ⁶Li default's justification inverts | `finite_q_delta` on the ⁷Li object — and the same question the ⁶Li backend answers with Eq. (17) |
| D6 | which A_zz | A_T = −b₁/F₁ / A_zz^{(3/2)} = −(2/3) b₁/F₁ | a factor 3/2; the spin-1 `azz()`'s 2/3 has no J = 3/2 counterpart | `asymmetries.cpp`, `SPIN32_FINITE_GAMMA.md` §6.4 test 10 |
| D7 | `ClusterPartialWave` refuses odd L | teach the type / keep working around | the i^L phase is unobservable for a single wave | `b1_nuclear.cpp:244,267` |
| D8 | `LightConeDensities` interface | L-generic alignment slot / explicit coefficient / keep the A = 2 shape | reusing φ₂ for L = 1 is "numerically exact, semantically a lie" | an interface change that would also derive ⁶Li's 1.5 and 6/√2 |
| D9 | b₂_32 | inherit 2x·b₁ / state it | none | one stated line |
| D10 | Δ_32 (cos 2φ) | none / reuse ⁶Li's toy 1e−2 | 3·Q_NN = ±3 at J = 3/2 gives a larger amplitude for the same Δ | `delta_32_func` |
| D11 | Q(⁷Li) constant — **PAID 2026-09-06** | −4.00(3) / −4.06 fm², with a source → **−4.06 fm² adopted** | not in the tree; ratio 0.871 → 0.858 — **confirmed: 0.871265 → 0.858389** | done: `LI7_QUADRUPOLE_FM2` = −4.06 beside `LI6_QUADRUPOLE_FM2` in `include/lipolgen/rc.hpp:726`, and §4.2's gate is `li7_alpha_t_quadrupole` + doctest T13 + pytest G8 |
| D12 | the J = 3/2 tensor run plan | two-state T = ±1 contrast / a thirds pattern | the pure-alignment fill is rank-3 clean by construction | a new plan in `spin.hpp`, `cli.py` |

**(e)** Deciding any of them changes nothing shipped until step 3 of
`phase_D_li7_rank2.md` §8 is taken; D11 is the one the run recommends
committing first (a gate before a b₁).

**CLOSED FOR ROW 20's D11, 2026-09-06** (task B3;
`../run_2026-09-06/phase_B_numbers.md` §B3). The constant is sourced,
committed and bound: **`LI7_QUADRUPOLE_FM2` = −4.06 fm²** at
`include/lipolgen/rc.hpp:726`, beside `LI6_QUADRUPOLE_FM2`, from TUNL's
A = 5, 6, 7 evaluation (Tilley *et al.*, NPA 708 (2002) 3), whose A = 7 half
prints `Q = −40.6 ± 0.8 mb (1988DI1B)` and whose A = 6 half prints the
`Q = −0.818(17) mb (1998CE04)` that IS `LI6_QUADRUPOLE_FM2` — **one document,
one sign convention, one unit rule**, which is what the brief's "match
conventions" required and is why the −4.06 option was taken rather than the
−4.00(3) the note quoted from memory. This row's own predicted price is
**confirmed to six digits**: the A = 7 α–t gate ratio is **0.871265** against
−4.00 and **0.858389** against −4.06 (12.87 % vs 14.16 % low), a 1.5 % swap
that does not change the verdict "13–14 % low". §4.2's gate is now real code —
`li7_alpha_t_quadrupole` in `b1_nuclear.hpp`, reproducing the note's offline
numpy recipe **bit for bit** — pinned by `tests/test_b1_nuclear.cpp` T13
(26 assertions) and `python/tests/test_li7_rank2.py` G8 (4 tests). It is a
**reported ratio and not a b₁ verdict**: b₁(⁷Li) is still unimplemented and
still blocked on D2, and G8 re-asserts ⁷Li's exactly-zero rank-2 sector in the
very test that runs the gate. **No new registry row was added.** D1, D3–D10
and D12 are untouched.

### B21. The licence, GPL-3.0-or-later — *carried forward, confirm only; E2 built on it* (`STATUS.md` decision **row 24**; `OPEN_ITEMS_SOLUTIONS.md` top table row 13)

**(a)** Confirm GPL-3.0-or-later as the project licence — which this run
stamped on 101 source files (`SPDX-License-Identifier`, measured by
`validation/check_spdx_headers.py`) and wrote into `CITATION.cff:13`, with
`LICENSE` already at `a94fd6e`.

**(b)** (i) GPL-3.0-or-later (in the tree); (ii) GPL-3.0-only. Anything more
permissive is excluded by the row's own note (forced by HepMC3/LHAPDF).

**(c)** `OPEN_ITEMS_SOLUTIONS.md` row 13 still reads "author to confirm"; no
site says it was confirmed. B24's licence statement (the wheel is a combined
work under GPL-3.0-or-later) rests on it.

**(d)** No number. Changing the identifier is **103** headers (101 when E2
stamped them; re-measured 2026-09-16 by the gate itself) + `CITATION.cff` +
`README.md`'s licence section + `ci.yml:33`/`PACKAGING.md`'s combined-work
paragraphs, all gated by `test_spdx_headers.py` / `test_release_metadata.py`.

**(e)** (i): close row 13 with a date. (ii): the files above and the two gates'
expected string.

### B22. The author's name (`STATUS.md` row 13)

**(a)** Whose name goes in `CITATION.cff:20` and `AUTHORS:4`, where the run
left `<AUTHOR NAME — to be filled by the author>`?

**(b)** A name; or leave the placeholder (the gate accepts either, but not a
half-filled pair).

**(c)** `phase_E_numbers.md` §E2.4: the git identity here is a project name
and the configured email is not a name, so nothing could stand in; no
copyright line was added to any source file.

**(d)** No number.

**(e)** Two files, kept in sync by `python/tests/test_release_metadata.py`,
which checks "either still the placeholder everywhere or replaced everywhere".

### B23. Push the CI workflow — and the tree makes this inseparable from pushing the run (`STATUS.md` row 14)

**(a)** Should `.github/workflows/ci.yml` and `.github/scripts/build_deps.sh`
go live on GitHub?

**(b)** The row frames this as "push them, or not". The tree does not offer
that choice as framed: both files are **inside commit `7f68339`**, the
workflow's trigger is `on: push: branches: [master]` plus `pull_request` and
`workflow_dispatch` (`ci.yml:37-41`), and `master` is 10 commits ahead of
`origin/master` (= `a94fd6e`) with `PLAN.md`'s rule "the user pushes". So the
real options are: (i) push the run as it stands — GitHub runs the four-job
workflow **on that push**, on its runners, fetching from the six upstream
URLs, and that first run is the one nothing has verified; (ii) add one commit
before pushing that deletes `.github/` or narrows `on:` to
`workflow_dispatch`, then push — the run goes up with the workflow inert until
the author triggers it; (iii) do not push the run.

**(c)** `phase_E_numbers.md` §E4.2–E4.5, `ci.yml:1-34`, `PACKAGING.md` §1:
exercised locally end to end 2026-09-05 (Ubuntu 22.04.5 / glibc 2.35 / gcc
11.4.0 / 8 cores) — `build_deps.sh` 337.28 s (more than half of it Python
bindings nothing here imports), job 2 configure 3.08 s, build 95.93 s / 0
warnings, C++ 400 / 17 240 262 / 1 skipped bit-for-bit, pytest 923 / 114
(no `.git`, no sibling `PolarizedLithiumSim`), cache 431 MiB / 113.5 MiB
compressed — all of it measured at `7f68339`, before phase F added one doctest
case and one pytest. On GitHub, where `actions/checkout` supplies the `.git`
and the sibling is still absent, the same job on the current tree would be
expected to read **925 / 113**. Not established: `act`/docker are not
installed, so nothing parsed the file as GitHub parses it — cache round-trip,
the image's apt set, 4-core timings (extrapolated 8–12 min cold),
expression validity (`${{ env.* }}` in four places). The workflow has no
upload step. The row's "no Actions run exists for this repository" was not
verified by this review (no network).

**(d)** No number in the tree moves under any option. Cost of (i): a first
unverified run on GitHub's runners at push time. Cost of (ii): one commit.

**(e)** (i): nothing. (ii): `git rm -r .github` or the `on:` block, one
commit, and `PACKAGING.md` §1 / `phase_E_numbers.md` §E4.10 restated.

### B24. Publish an `auditwheel`-repaired wheel — and what to do about the data trees and the baked-in paths (`STATUS.md` row 15)

**(a)** Publish the `manylinux_2_35_x86_64` wheel the run built, repaired,
verified and discarded?

**(b)** (i) do not publish; (ii) publish as built; (iii) publish with the
PYTHIA `xmldoc`/`pdfdata` and LHAPDF set trees inside (~300 MiB, and a
separate redistribution question for the LHAPDF grids); (iv) rebuild inside an
older manylinux image for a wider tag. (ii)–(iv) each also need an answer on
the builder's home directory surviving as a literal string in the artefact.

**(c)** `phase_E_numbers.md` §E4.6–E4.8, `PACKAGING.md` §2–3: `pip wheel .`
87.77 s → 1 779 767 B; `auditwheel repair` 5.35 s → 7 232 971 B (×4.06),
vendoring `libpythia8` 14.0 MB (GPL-2.0-or-later), `libLHAPDF` 1.24 MB
(GPL-3.0), `libHepMC3.so.4` 1.18 MB (GPL-3.0). In a `bwrap` namespace with the
deps prefix replaced by an empty tmpfs: `import lipolgen` with all tiers,
the pure-C++ generator and the vendored VMC tables work (⟨x⟩ =
0.027023337580082865 bit-identical); `LhapdfSF` and `--hadronize` **fail**
(their compiled-in defaults are this machine's absolute paths); with the two
trees at a different path and the three env vars exported the T2 HepMC3
output is byte-identical (sha256 `68ed6f06…e405ac`). `libLiPolGenCore.so`
keeps `RUNPATH $ORIGIN:/home/cpeng/Projects/polli/deps/install/lib`, and four
more `…/deps/install/…` strings sit in the vendored `.so`s. Licence: the
combination ships under GPL-3.0-or-later (B21).

**(d)** No number in the tree. Under (ii) a user without the two data trees
cannot run `LhapdfSF` or `--hadronize`.

**(e)** (i): nothing. (ii)–(iv): a release process outside the tree; if the
paths are to be scrubbed, `CMakeLists.txt`'s RPATH for the core library and a
rebuild.

### B25. Send the Mäntysaari-group letter — all of it, (a)/(b) only, or wait (`STATUS.md` decision **row 21**, added 2026-09-05; `mantysaari_collaboration_draft.md`; `OPEN_ITEMS_SOLUTIONS.md` §11.2/§13; `PLAN.md` C6)

**(a)** Send the reconciled draft, and in what form?

**(b)** (i) send (a), (b) and (c) with (c) on photoproduction and its four
unestablished factors inside the ask; (ii) send (a) and (b) only, which are
channel-independent; (iii) wait for further work.

**(c)** `phase_C_numbers.md` §C2 (four passes), §C6: the O5 verdict is
**MARGINAL, as a band** — S = 2.63 σ at the band's low edge and 2.84 … 3.29 σ
at its top, 3 σ at 8.3 … 13.0 fb⁻¹/u, inside {1, 10, 100} at both ends, with
ε_det the dominant unquantified factor and `Optics::lumi_fraction` the single
correction that on its own restores a NO (0.73–0.92 σ, 106–167 fb⁻¹/u).
`a2_from_quadrupole` is a closed form, not a Good–Walker amplitude. Nothing
has been sent (the draft says so at its top and bottom). **Surveyed 2026-09-15
(run 2026-09-06, task C1):** the sibling `../PolarizedLithiumSim` was read end
to end for the decay-lepton reconstruction efficiency the ask's fourth caveat
names — every module of `fastsim/polli_fastsim/`, all of `tools/fullsim/` and
`tools/analysis/`, `plans/03`, `plans/09` and all 54 `refs/` entries, plus
whole-tree greps — and supplies none, so that caveat stays **unbounded rather
than bounded** and the draft's wording of it is unchanged
(`run_2026-09-06/phase_C_numbers.md` §C1, `phase_C_survey.md`).

**And the ask's factor list gained a BOUND on a different leg, 2026-09-15**
(`run_2026-09-06/phase_C_numbers.md` §C5.2; `phase_C_survey.md` §C-S6). By
Chang *et al.*'s own §IV — *"only accounts for the acceptance effect and does
not incorporate the efficiencies of the detector … nor the efficiency and
acceptance of the reconstructed distribution"* — `COHERENT_JPSI_EFF_IR8_LI7`
= **0.1775**, the single largest multiplier in the O5 chain, is a **pure
geometric acceptance** and therefore an **UPPER bound** on the far-forward
intact-recoil tagging efficiency × acceptance. Direction **DOWN**,
**unquantified**: every omitted factor is ≤ 1 and nothing says how far. It is
a **third** factor on the leptons' side, so it **compounds with the
decay-lepton omission and breaks** the ×1.1224-UP / ×0.897-DOWN = 1.007
balance that is the stated reason neither of those two may be quoted alone.
**The tree STATES the bound and folds NO value in** — applied 2026-09-15 at
the constant's home (`coherent.hpp:250-260`), in `o5_a2_reach.py` and in
`OPEN_ITEMS_SOLUTIONS.md` §11.3 — so *"send the ask with its four
unestablished factors inside it"* now means sending a list in which one entry
is bounded ABOVE and three are not bounded at all.

**(d)** No number in the tree moves under any of (i)–(iii) — the letter is
the artefact and the band was measured before it. **Including the recoil-leg
bound: it is a direction, not a value**, so the O5 band is untouched at both
ends and was re-verified bit-identical after it was recorded (2.618 σ at
13.1 fb⁻¹/u headline; LOW 2.627 at 13.0; TOP 2.842 … 3.292 at 11.1 … 8.3;
`o5_a2_reach.py` `self-check: OK`, every rung bit-identical). **What has no
measured price and is stated as such:** folding any lower ε in. Nobody has a
number to fold — that is the point of the bound — so the cost of doing it
cannot be quoted, only its direction (down, on top of the decay leptons).
Option (i) is the only one whose TEXT changes: the ask's fourth caveat reads
differently once one of its factors is bounded above.

**(e)** Nothing in the tree; the draft is the artefact. The bound's own
sentences are already in the tree at the three sites above and are not
contingent on this row.

> **2026-09-23 — two wording corrections to the draft, requested by the benchmark run and NOT applied** (the draft is outward-facing; it is the author's). (1) **G-17** (`../run_2026-09-23/phase_B2_nuclear.md` §8): the draft's "anchored on the *measured* asymptotic α–d D/S ratio η = −0.025 ± 0.006 ± 0.010, George & Knutson" → "anchored on the asymptotic α–d D/S ratio η = −0.025 ± 0.006 ± 0.010 from George & Knutson's restricted phase-shift analysis" — the paper's title says it is a restricted phase-shift analysis, not a measurement of d + α tensor analysing powers, and η is exactly linear in the quadrupole dial, so it is a consistency band on that dial. The record names the site at the draft's lines 102-104 (anchor read 2026-09-23 "anchored on the *measured*"); the letter body repeats the claim at its line 531 (anchor read 2026-09-23 "dial anchored on the measured asymptotic D/S ratio"), which the record did not list. (2) `mantysaari_collaboration_draft.md:67` (`../run_2026-09-23/phase_B3_chain_rc.md` §6.2): "An unpolarized coherent-rate baseline that IS citable" → "An unpolarized exclusive-VM coherent-rate scale that IS citable" — 85.2 % (J/ψ), 95.4 % (φ), 97.2 % (ρ) of eSTARlight's rate is at Q² < 0.1 GeV², where `coherent.hpp` generates nothing. Outside the letter, every live site a grep finds carries both corrections as of 2026-09-23 — the last of them only after the run's fix stage (George–Knutson: `tests/test_cluster_config.cpp` comments at T22b/T23, `docs/PHYSICS_CHANNELS.md` §13's Good–Walker item, `OPEN_ITEMS_SOLUTIONS.md` §11.2 and §11.3b, `python/lipolgen/configs.py`; eSTARlight: `docs/PHYSICS_CHANNELS.md` §13's eSTARlight item, `docs/benchmarking/00_in_tree_checks.md` §8, `docs/benchmarking/01_generators.md` §4, `docs/references/00_corpus.md`, `docs/references/01_generators-chain.md`, `docs/open_items/physics_literature.md`, the `estarlight_li6_coherent` docstring in `python/bindings.cpp`). Not changed: the dated records, which carry an adjacent dated note for George–Knutson (`run_2026-09-02/design_G_cluster_config.md`, `run_2026-09-02/phase_G_numbers.md`, `run_2026-09-03/phase_C_numbers.md`, `run_2026-09-03/PLAN.md`) and none for eSTARlight (`run_2026-09-02/estarlight_li6.md`'s title, `run_2026-09-02/STATUS.md`, left as written on their date); the TEST_CASE name "the eSTARlight 6Li baseline" in `tests/test_coherent.cpp` (its comment carries the scope); and `06_critic.md`'s D-4 row label, which names the entry it critiques. Nothing else in the letter moves, and nothing has been sent. The row's options (i)–(iii) are unchanged; a send under (i) or (ii) would carry the two sentences as they stand unless the author applies them.

---

### B26. The MSTW rows still tallied as PASSED when the grid is absent — split T16/T17 out, or leave them (`STATUS.md` decision **row 25**, added 2026-09-05)

**(a)** `tests/test_b1_nuclear.cpp` T16 (the `--b1-unpol mstw` reach row) and
T17 (the third-crossing row of the verdict configuration) still guard their
MSTW block with an in-case `if (!mstw_grid_present()) { MESSAGE(...) }`, so
with the grid absent doctest tallies each case as **passed** while that row was
never evaluated. Split each into its own registration-time-skipped case, as
T1v and (2026-09-05) `item 5v` are — or leave them?

**(b)** Phase F applied the split to the rows the ⁶Li ban lift RESTS on and to
the cases that are nothing but MSTW: `item 5v` (gate condition 3, the 1.000338
CD-Bonn + MSTW peak ratio, decision row 4's second ratio) is now its own
`doctest::skip(!mstw_grid_present())` case, and so are the four
`tests/test_mstw_sf.cpp` cases (which previously reported "3 passed,
0 assertions" with the grid absent). T16 and T17 were left, because each sits
in a case whose ToyF2 rows are shipped-default rows that must still run and
whose split would duplicate a lambda; each prints what it did not measure.

**(c)** Measured with the grid absent
(`LIPOLGEN_PYTHIA8_PDFDATA=/nonexistent/pdfdata`): before the split,
`-tc="*item 5*,*T17*"` reported `5 passed | 0 failed`; after it,
`-tc="*item 5v*,MstwSF*"` reports `0 passed | 402 skipped`. With the grid
present nothing changes: the cases run, and no pin moved.

**(d)** No shipped number is involved either way. Leaving T16/T17 costs the
tally's honesty on those two rows only; splitting them costs two more test
cases and a duplicated helper.

**(e)** If split: two new `TEST_CASE`s in `tests/test_b1_nuclear.cpp`, and
`docs/USAGE.md` §2a's sentence — which now names exactly which rows are
registration-skipped and which are tallied as passed — becomes "all of them".

### B27. The tagged S–D interference phase, and the `tagged.json` re-pin that followed it — *applied, confirm or revert* (`STATUS.md` decision **row 26**, added 2026-09-15; `src/core/tagged.cpp` `build_amp2`; `../run_2026-09-06/STATUS.md` CW row; `../../benchmarking/07_cw_sign_investigation.md`)

**(a)** Two questions that travel together. **(a1)** Confirm that
`TaggedModel::build_amp2` must apply the partial-wave phase φ_L = i^L ψ_L
(written `(-1)^floor(L/2)` on the real radial amplitudes), which it did not
before `a7b3d18` — so **A_zz^tag flips sign on every spin-1 tagged channel at
the shipped defaults**. **(a2)** Confirm that
`validation/reference/tagged.json`'s two spin-1 `model` blocks stay re-pinned
**from this library**, rather than from `polligen`, whose own
`tagged._amp2_table` carries the same missing phase.

**(b)** (a1): (i) keep the phase (applied); (ii) revert to the unphased sum.
(a2): (i) keep the re-pin from this library (applied) — `tagged.json` stops
being an external port gate for those two blocks; (ii) restore `polligen`'s
blocks as a **known-wrong port gate under an `xfail`** that names the bug — an alternative OVERTAKEN on 2026-09-15: the sibling fixed its own phase (commit `1066555`), so the polligen port gate can simply be RESTORED with no `xfail`, the fixed polligen agreeing with the re-pin at 2e-13 (re-measured 2026-09-16) —
keeping the external gate at the cost of a permanently failing pin.

> **2026-09-23 — (a2) overtaken in the tree:** the carry-through was removed and `tagged.json` is dumped from `polligen` in every block again (provenance `"polligen"`, the fix commit `1066555` named); the two blocks moved by ≤ 1.33e−13 / 2.28e−13 relative. (a2)(i)'s "`tagged.json` stops being an external port gate" no longer holds, and (a1)(ii) now also costs a failing `polligen` port gate. The row stays *applied — confirm or revert* for (a1). `../run_2026-09-23/phase_A_port_gate.md`.

**(c)** The run calls (a1) **certain**, on three independent derivations
(`07_cw_sign_investigation.md`). The Cosyn–Weiss Eq. (6.12) identity
A_zz^wf = [(2 f₀ + f₂/√2)(f₂/√2)/(f₀² + f₂²)](1 − 3 cos²θ_k) holds to
**8.882e−16** on every cell of the 280 × 96 grid on both wave functions, and
CW TABLE II's three cell-centre rows reproduce to residuals **4.36e−06 /
2.77e−06 / 3.64e−07** (`tests/test_tagged.cpp:723` (anchor read 2026-09-15 "-1.9371243623"),
`tests/test_tagged.cpp:724` (anchor read 2026-09-15 "+0.9993127695"),
`tests/test_tagged.cpp:726` (anchor read 2026-09-15 "+0.9673403636"), all
pinned at atol 1e−4). A regression guard at
`tests/test_tagged.cpp:741` (anchor read 2026-09-15 "ik30] < -1.5") pins the sign
that separates the two amplitudes, so the pre-fix behaviour cannot return
silently. This is a **bug fix against a published formula**, not a modelling
choice — which is the argument for treating it differently from the other
applied rows, and the reason a row exists at all is that it is nonetheless the
one shipped output that moved in `91e48b9..HEAD`.

**(d)** **Measured cost of the applied option** (re-measured 2026-09-15 on
`86cd1e6`). A_zz^tag(k = 0.20 GeV): **+0.8450 → −1.2069** (Hulthén),
**+0.4518 → −0.5191** (VMC). At the CLI's own default fill (`--pz 0.7`
max-entropy ladder) the ⁶Li Hulthén tag-fraction prediction goes
**0.02703 → 0.02040**, −24.5 %, and the tensor-thirds categories go
0.02813 → 0.01840 (azz±) and 0.01778 → 0.03722 (azz0). **What does NOT move:**
σ_pb to ≤ 2e−16 relative; every inclusive, coherent and ⁷Li tagged cell
byte-identical on two seeds; P_D (a norm, phase-blind); the four dilutions to
≤ 1.4e−6 absolute (angle-integrated, so ∫Θ₀Θ₂ dc = 0 kills the cross term and
only the 96-cell quadrature residual survives); `LI6_CLUSTER_POLARIZATION`
(§B6, checked). **Reference files:** `b1_default_li6` (a self-pin of LiPolGen's own C++, not a polligen reference), `beams`, `bookkeeping`,
`coherent`, `spectator`, `spin`, `xsec` all sha256-identical to `a94fd6e`;
`tagged.json` differs **only** in `channels/deuteron/model` and
`channels/li6_alpha/model` (8 fields each — e.g. deuteron `vector_dilution`
0.932494769105 → 0.932496109312) plus two added provenance keys, with the
`li7_alpha` block and every non-`model` part identical; `_manifest.json`
differs in the one `tagged.json` provenance line. (Since 2026-09-23 `_manifest.json` is byte-identical to `a94fd6e` again and `tagged.json` differs from it only in the two spin-1 model blocks and the two provenance keys; the same day `b1_default_li6.json` gained a `provenance` label, whole-file sha256 `d7bd8ce6…` → `8373db8e…` with every numeric block byte-identical — `../run_2026-09-23/phase_B3_chain_rc.md` §4.) **The brief's expectation
`validation/reference/` byte-identical to `a94fd6e` must therefore be read as
byte-identical EXCEPT those two `model` blocks and the manifest line, by this
row.** **Cost of (a1)(ii):** restores a sign the run calls refuted at blocker
severity, re-breaks the CW gate, and moves every number above back.
**Cost of (a2)(ii):** an `xfail` on a reference pin, i.e. a gate that is green
only because it is told to expect failure.

**(e)** Confirm both: close the row; the provenance note already in
`validation/README.md` and `_manifest.json` stands. Revert (a1): one line in
`src/core/tagged.cpp` `build_amp2`, the five re-pinned test sites listed in
`../run_2026-09-06/phase_CW_numbers.md` §9, and a re-dump of `tagged.json`.
Revert (a2) alone: re-run `validation/dump_polligen_reference.py` for those two
blocks and mark the resulting pin `xfail` with the bug named. (Moot since 2026-09-23: the restored `polligen` blocks carry the fix, so there is nothing to `xfail`.)

### B28. `PolarizedLithiumSim`'s published A_zz^tag numbers carry the same inverted sign (`STATUS.md` decision **row 27**, added 2026-09-15; `../run_2026-09-06/STATUS.md` CANDIDATE section)

> **ANSWERED IN THE SIBLING (2026-09-15 19:50, its commit `1066555`, author Chao Peng — not this run):** `evgen/polligen/tagged.py` now applies the relative phase (−1)^(L//2), the AV18 deuteron control gate was added, money plot 4 was regenerated (90° curve −0.482 → +0.922; folded +0.491 → −0.843; `evgen/README.md:117` now reads "the folded A_zz reads −0.84"), and the sibling's commit records agreement with LiPolGen's re-pinned reference tables at 2e-13 — re-measured here 2026-09-16 against the sibling at fe1e58e: worst relative difference 1.33e-13 (li6_alpha) / 2.28e-13 (deuteron) between the fixed polligen's model blocks and LiPolGen's re-pinned `tagged.json`, inside the C++ gate's rtol 1e-12. The "+0.49 … −0.48" text quoted below exists only at the pre-fix commit 2a27972. On 2026-09-23 that agreement was used: `tagged.json`'s port gate against `polligen` was restored in every block (`../run_2026-09-23/phase_A_port_gate.md`, re-measured at the sibling's `c0f86a8`: the same 1.33e−13 / 2.28e−13).

**(a)** *(as it stood 2026-09-15 before the sibling's fix)* The sibling repository's `tagged.py:247` built its amplitude the way
this library did before `a7b3d18`, so its published A_zz^tag numbers — and the
figure made from them — carry the inverted S–D interference sign. Should they
be regenerated?

**(b)** (i) fix `tagged.py:247` and regenerate the numbers and the figure;
(ii) fix the code and flag the existing figure as superseded, without
regenerating; (iii) record the finding only — the status quo until 2026-09-15 19:50, when the sibling took option (i) itself. **Row 27 is therefore CONFIRM ONLY: the question is answered in the other repository.**

**(c)** The defect is **inherited, not independent**: it is the same missing
φ_L = i^L ψ_L, and §B27(c)'s three derivations settle it. This is a decision
rather than a cost inside row 26 because it is about **another repository**
with its own published output, and because nothing in this tree changes
whichever way it goes.

**(d)** **Not measured here, with ONE exception, stated as such —
precisely.** On 2026-09-16 — and only then, after the sibling's own fix
`1066555` — `evgen/polligen/tagged.py`'s `TaggedModel` WAS executed,
read-only, to produce the box above: the 1.33e-13 (li6_alpha) / 2.28e-13
(deuteron) comparison against LiPolGen's re-pinned `tagged.json`. That
comparison cannot be made by reading, and the sibling ships no dumped model
JSON, so the blocks had to be computed. **No other script of that tree was
run, and nothing in it was written by this project.** This paragraph said
"nothing ... was EXECUTED ... and no number was produced from it" until
2026-09-16, which the box eleven lines above already contradicted. Other
numbers ARE quoted **from** it by reading: its `evgen/README.md:113-115` says *"the folded A_zz reads +0.49
and −0.07 at k ≈ 0.33 GeV/c where the 90° curve says −0.48"*, and
`evgen/money_tagged_azz_6Li.png` is named as the figure that carries the same
sign. The one number that was **measured** is this tree's own reproduction of
that curve value — **−0.4786 pre-fix, +0.9279 after** the fix
(`../run_2026-09-06/phase_CW_numbers.md` §11.2, row 5) — which is a LiPolGen
measurement standing beside the sibling's published figure, not a re-run of
it. What is otherwise known is the code path and that it matches the pre-fix
one. *(This paragraph read "no number from it is quoted" until 2026-09-15,
which `STATUS.md` row 27 repeated and §11 of the phase record contradicted at
eight itemised sites; both are corrected.)* (i) costs
a regenerated figure and whatever cites it; (ii) costs a flag and leaves the
numbers wrong; (iii) costs a published figure that this project has established
carries a wrong sign.

**(e)** (i) or (ii): a change in a repository this run did not touch, plus
whatever announcement its publication history requires. (iii): close the row
with the finding on the record.

---

## Part B — what this review found wrong with the decision bookkeeping

Measured against the whole diff `a94fd6e..dffe94e`, in the working tree.
Nothing below was fixed **by this review**; phase F then applied the
bookkeeping fixes on 2026-09-05, and each item below says what was done. The
decisions themselves are still open — applying a fix here never means taking
one.

1. **Row 14's framing does not match the tree.** `ci.yml` is committed in
   `7f68339` with `on: push: branches: [master]`; the run's ten commits are
   unpushed and the user pushes. "Push the workflow, or not" is therefore
   "push the run, or first remove/disable the workflow" — B23 restates it.
   The file's own header (`ci.yml:3-5`) knows this ("none of the run's
   commits … have ever been pushed"); the row did not say it.
   **APPLIED 2026-09-05:** `STATUS.md` row 14 now states the coupling and the
   three options; the decision itself is untouched.

2. **The row-12 dependency reversal was corrected at two sites and not at a
   third, in a live document.** Row 8 and `coherent.hpp:468-481` now say the
   ceiling does *not* follow from `eps_b0`; `OPEN_ITEMS_SOLUTIONS.md:1540-1555`
   (§11's O4 entry) still says "`COHERENT_T_MAX_DEFAULT` = 0.2 is a consequence
   of that oversized `eps_b0` (the positivity edge would sit at |t| = 2.8 GeV²
   …)" — the pre-reversal sentence, with the "2.8" the run elsewhere replaced
   by "2.80" — while §11.6 of the same file
   (`OPEN_ITEMS_SOLUTIONS.md:2069-2082`) states the reverse.
   `phase_C_numbers.md` §C4.5 and `phase_D_small_items.md` §D5.1 carried the old
   sentence too; those are dated records, but neither was marked superseded.
   **APPLIED 2026-09-05:** the §11 O4 sentence is rewritten (the edge is the
   consequence, the ceiling is not), both records carry a dated "superseded by
   D5" note, and the last bare "2.8" is now "2.80".

3. **Decisions the tree calls the author's that the table does not carry**
   (each is a section above): the α–d source default (B8: `phase_C_numbers.md`
   §C3.5 "an author decision", `USAGE.md:2576`, `OPEN_ITEMS_SOLUTIONS.md` §11 O2, "`FitRescaled` **stays the
   default by author decision**");
   the letter (B25: "the author's decision alone", four sites); the RC anchor
   (B13: `rc.hpp:288` "the author's call"); G3a's ceiling (B4:
   `OPEN_ITEMS_SOLUTIONS.md:970` "the design's and the author's call") and
   the floor change the run *did* make; the thirteen ⁷Li items (B3, B18, B20:
   `OPEN_ITEMS_SOLUTIONS.md` §15.5 is titled "The author decisions this
   needs"); the PomSet-11 refusal (B19, applied); and the three pre-existing
   "author to confirm" rows (B1, B6's carry-forward, B21), one of which (the
   licence) the run stamped onto 101 files.
   **APPLIED 2026-09-05:** all of them are now `STATUS.md` decision rows
   16–24, each pointing at its section here.

4. **Four registries, and the table is pointed at by only five of its own
   rows.** Sites in code and docs cite `STATUS.md` decision rows **7, 8, 12,
   13, 14/15** only (row 7 `phase_B_numbers.md:1664`; row 8 `coherent.hpp:84`,
   `USAGE.md:1566`, `OPEN_ITEMS_SOLUTIONS.md:2144`; row 12
   `phase_D_numbers.md:1432`; row 13 `test_release_metadata.py:7`,
   `check_spdx_headers.py:14-16`, `README.md:290`; rows 14/15 `README.md:268` (anchor read 2026-09-16 "activates it")
   and `README.md:256` (anchor read 2026-09-16 "PACKAGING.md")). Rows 1–6 and 9–11 were pointed at by none: their
   dependents cite `OPEN_ITEMS_SOLUTIONS.md` §10 (conditions 4/6),
   `phase_A_miller_normalisation.md`, or `phase_C_numbers.md` §C5 instead
   (`constants.hpp:107`, `CONVENTIONS.md:106/233/266`, `USAGE.md:207/292`,
   `beams.hpp:71`). Since 2026-09-05 `README.md:167-177` cites rows 1, 2 and 3
   by number, so those three now have a pointing site and rows 4–6 and 9–11
   still have none. Each decision is stated consistently across its sites
   (numbers checked: 24.6 %, 0.057600, 11.61 %, 2.069 /
   2.027 %, 0.9345 / 11.42, 0.245 / 2.80, 0.440 / 0.843243 / 1.000338), but
   not **once**.
   **APPLIED 2026-09-05:** `STATUS.md`'s table is named the ONE registry in its
   own preamble, and §10 conditions 4/6, §15.5 and the top-table rows 4 and 13
   now point at its row numbers.

5. **"D2" names three different decisions in the same document set**: the
   2026-09-03 sweep's D2 (SF injection, `OPEN_ITEMS_SOLUTIONS.md` row 14 and
   `STATUS.md`'s phase D row), §15.5's D2 (the ⁷Li unpolarised backend,
   `PLAN.md` "decisions D2 and D11"), and plans/08's D2
   (`OPEN_ITEMS_SOLUTIONS.md` row 3 "unblocks D2", a file not in this tree).
   Likewise "D1".
   **APPLIED 2026-09-05:** `PLAN.md` C-phase text, `OPEN_ITEMS_SOLUTIONS.md`
   row 3 and §15.5's preamble now qualify every "D<n>" with its registry.

6. **Row 6's file pointer is ambiguous.** "`phase_C_numbers.md` §8.1c/§8.2"
   is `run_2026-09-02/phase_C_numbers.md`; the file of the same name beside
   the table has no §8.1c. **APPLIED 2026-09-05:** row 6 now spells the run
   directory out.

7. **The tensor sign's status differs by site**: `constants.hpp:17` "author
   decision, plans/08 D1" (taken); `OPEN_ITEMS_SOLUTIONS.md` row 3 and
   `docs/surveys/needs_survey.md:65,402` "author to confirm" / "unresolved
   author decision" (open). `needs_survey.md:65` also cites
   `fastsim/polli_fastsim/asymmetries.py:41` with `TENSOR_LL_SIGN = +1.0` — a
   sibling repository, not this tree, where the constant is −1.0.
   **APPLIED 2026-09-05:** `needs_survey.md` says so at both sites and
   `OPEN_ITEMS_SOLUTIONS.md` row 3 records the split status; whether to CLOSE
   the row is B1 and is untouched.

8. **Costs not measured where the table implies a priced choice**: row 3's
   third option (one R hook) is named "the better" and unpriced (B7); the
   RC-anchor alternatives (B13) and honouring `--pzz` on `helicity-flip`
   (B18) have no measured effect on any observable.
   **APPLIED 2026-09-05:** row 3 says "unpriced" in the row itself, and the
   two new rows 17 and 20 carry the word for the other two.
   **CLOSED FOR ROW 3, 2026-09-06:** the third option is implemented opt-in
   (`--r-source`) and priced in `../run_2026-09-06/phase_A_numbers.md` §A1;
   §B7(d) and `STATUS.md` row 3 carry the numbers, and the word "unpriced" is
   gone from both.
   **CLOSED FOR ROW 20's D13, 2026-09-06:** honouring `--pzz` is implemented
   opt-in (`--pzz-mode typed`) and priced in §A3 of the same file; §B18(d) and
   `STATUS.md` row 20 carry the numbers, and the measured answer is that on
   ⁷Li it moves no observable (σ to 1 ulp) while moving 0.106 % of a
   100 000-event sample into other cells, and on ⁶Li it moves the rate by
   −0.0023527 %.
   **CLOSED FOR ROW 17, 2026-09-06:** the RC-anchor alternatives were RUN, not
   estimated — §B13(d)'s three-row band table
   (`../run_2026-09-06/phase_A_numbers.md` §A2,
   `python/tests/test_rc_low_x_anchor.py`) — and the word is gone from row 17
   as well.  *(This line read "Row 17 still carries the word" until
   2026-09-15; §B13(d), the index table's shape column and `STATUS.md` row 17
   had all said PRICED since 2026-09-06, so the tree disagreed with itself at
   this one site.  **All three of the run's unpriced costs are now priced**,
   and `SUMMARY.md` item 8 carries the same correction.)*

9. **The count of ban-lift sites is stated three ways**: `PLAN.md` A7
   "fourteen files", `phase_A_numbers.md` status note "the sixteen files",
   `STATUS.md` phase A row "across 18 files". They may be three moments in
   time; none said so. **APPLIED 2026-09-05:** the measured count is **17**
   (`git grep -il 'publication ban'` at `ac22331` and at `HEAD`;
   `phase_A_gate_mechanics.md` §4 enumerates them), it is stated with that
   basis at all three sites, and the two older counts are labelled
   pre-survey / mid-survey.

10. **Row 5 sat out of order** (after row 12) in the table — cosmetic, but
    the batch above renumbers, so cite `STATUS.md` rows by number, not
    position. **APPLIED 2026-09-05:** row 5 now sits between rows 4 and 6, and
    the table's preamble says to cite by row number.

11. **Unverifiable here, stated as such**: row 14's "no Actions run exists for
    this repository" (network); every "measured on GitHub" expectation in
    `ci.yml`.

12. **This file's own `path:line` anchors are checked by nothing, and they
    went stale the same day they were written.** It was written against
    `dffe94e`; the close-out then rewrote `README.md`, `docs/USAGE.md`,
    `docs/PHYSICS_CHANNELS.md` and `include/lipolgen/rc.hpp`, and its anchors
    landed on blank lines, on comments, on a table header and in the wrong
    paragraph — the same drift `validation/check_physics_channels_links.py`
    exists to catch. **APPLIED 2026-09-05:** every `path:line` in this file
    was re-read against the tree; **38 anchors moved**, and three were added
    (`README.md:256` (anchor read 2026-09-16 "PACKAGING.md"), `README.md:268` (anchor read 2026-09-16 "activates it"), `README.md:313` (anchor read 2026-09-16 "and CITATION") — the rows 15/14/13 citations item 4's
    list was missing). **NOT applied, and an option rather than an
    oversight:** putting this file into the gate's `EXTRA_DOCS`. Measured with
    the gate's own `REF`/`RANGE` regexes over this file with THIS ITEM
    EXCLUDED (writing the item changes the count), it would see **11 of 81
    `path:line` strings** — 6 points and 5 ranges, 9 of them `README.md` —
    because the gate's `PATH` matches only a repository-relative path (or bare
    `README.md`) and this file writes most of its anchors as bare basenames
    (`rc.hpp:270`, `USAGE.md:207`). Real coverage therefore means first
    rewriting every anchor in repository-relative form, and adding the
    document also moves every sentence in the tree that counts the gate's
    documents: `README.md:112` "over three documents", `SUMMARY.md:106` "the
    strict gate covers **three** documents", and the "both extra documents" of
    `docs/PACKAGING.md:164` and `phase_E_numbers.md:707`, plus §E1's own
    record of why those two and not others. Both are the author's call, not a
    correction pass's.

Everything else the review checked holds: no default other than
`B1_CDKS_TABLE_TO_PER_NUCLEON` (row 1) moved between `a94fd6e` and `dffe94e`
(`constants.hpp`, `PipelineConfig`, the CLI defaults, `finite_q_delta`,
`COHERENT_T_MAX_DEFAULT`, `eps_b0`, `LI6_CLUSTER_POLARIZATION`,
`alpha_d_source` all checked against `a94fd6e`); the two behaviour changes
(`validate()` refusing `cluster_wave` where unread, PomSet 11 refused) are
refusals of inputs, not changes to what a default run produces; the suites and
gates re-measured above match the phase E row exactly.
