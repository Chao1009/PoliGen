<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-03 — summary (one page)

Supervisor, 2026-11-10, on `claude/run-2026-10-26` from `70982c9` (the end of
run 2026-10-27) to `34215c4` plus this file's own commit, the last of the run. The author's instruction, verbatim: "Go with your
recommendations for the 17 items, just record the choice and why." Every
choice below is therefore the SUPERVISOR's, under that delegation, recorded
with its reason in `DECISIONS.md` (A1–A17 and the dated notes under A1, A5,
A6, A13, A17); the author's words are quoted there and nowhere paraphrased
as a decision. Commits: A `8eea2c7` · P1 `3d059d3` · P2 `df56035` · P3
`5cd19d2` · P4 `4b981b4` · P4b `b591880` · P5 `888ca3b` · P6 `5372f2d` ·
V/F `34215c4` · S this file's commit. Phase records `phase_P*.md`, `phase_V_verification.md`;
the board `STATUS.md`.

## What moved, and what it measured

| decision | what moved | measured consequence |
|---|---|---|
| A1 (P1) | `LI6_CLUSTER_POLARIZATION` 0.811228 → 0.848, `LI6_B1_RANK2_TRANSFER` 0.921947 → 0.982576, `P_D_DEUTERON` 0.045 → 0.0576 (the AV18+UX family; the old values keep names, the Hulthén tag is bit for bit the old physics) | ⁶Li vector ×1.04533, rank-2 ×1.06576, deuteron vector ×0.97973; everything else column-identical |
| A3/A4 (P2) | \|g₁\| ≤ F₁ per nucleon on every kernel (`PositivityClippedG1`), one meta key, one banner line; refusals stay for what the clip does not reach | the clip binds on 2–5 × 10⁻⁷ of the rate; five of seven default top-cell refusals build, the two inclusive deuteron ones stay (rank-2 on the toy b₁); toys bit-identical by column diff |
| A5/A6/A7 (P3) | convention 8: a row's headline is the shipped default configuration; Deeps, BONuS, deuteron statics and HERMES b₁ re-pointed, the declared configurations recorded as sub-rows; V6 out of the TIGHT edge; the grid row's error bound `QUAD_REL_TOL × \|M₂\|` | report 14 / 9 / 0 → **15 / 8 / 0** — the one status move (deuteron statics FAIL → PASS, η −1.39σ on AV18) is the rule's and labelled; TIGHT saturations halve (the harness carried V6 as ×2) |
| A17 (P4, P4b) | an explicit default on a knob the run never reads is accepted and labelled `not-read`; a moved value stays refused; no built run has a `refused` row | the five accepted defaults bit-identical to the unnamed runs; only labels moved |
| A9, A10, A11, A12 (P4) | missing `lhapdf-config` a configure-time error unless `-DLIPOLGEN_WITH_LHAPDF=OFF`; `lipolgen.<key>` GenRunInfo attributes = the npz meta (T6 clause C24, sha re-pinned); \|σ_VMC\| ≤ 3; the FSI clipped-cell fraction on every FSI banner, registry row 28 | HepMC3 files differ from before only by the attribute lines; FSI clip 4.41 / 9.52 % (⁷Li), 0.70 / 5.25 % (⁶Li), 1.35 % (d) of cells |
| A13 (P5) | the `li6_alpha_vmc` block in the polligen port gate at a declared 1e−8, landed by splice | residuals 5.389e−9 (vector), 3.196e−9 (tensor); the in-place dump is no longer byte-reproducible on this machine (below) |
| A16 (P6) | throughput re-measured on the idle machine, raw log kept | 0.540 / 0.331 / 1.071 M ev/s by the 2026-10-30 recipe (0.61 / 0.33 / 1.05 at 2 M events); this run's code costs ×0.98 / ×0.96 / ×0.98 |
| A2, A8, A14, A15 | records only | — |

Gates at the end: C++ 435 / 435 cases, pytest 1523 passed / 155 skipped,
bench subset 368 / 1 / 0, citations 1243 / 0 broken (records 1432 / 0),
SPDX 136, doc-link 112; `make_report --check` OK at 15 / 8 / 0.

## Verification (V) and fixes (F)
Two read-only Opus lenses over `70982c9..5372f2d` (physics and honesty; code
and regression) found 11 + 8 findings, none a blocker; 27 adversarial
verdicts confirmed 18 and disputed one. Phase F (`34215c4`) applied 15 (stale sentences after A17, A1 and A3 dated
line-neutrally; the HERMES and EPPS21 headlines named as convention 8's
opt-in clause; the A11 refusal text printing the refused value exactly, with
a test; 71 re-pin comments dated; two board hashes and a registry title;
the dump script's shadowing wrapper and two header leftovers fixed line for
line; two cited range blocks re-recorded after their rows were re-read),
applied one in part (the labelled pre-A1 history stays), recorded one as
adequate (both sides of A11 pinned in C++ and Python) and left two to the
supervisor's notes. No physics default, constant, tolerance, registry label
or reference number moved; no check loosened. Three supervisor's
clarifications of 2026-11-10 (under A1, A5, A13) rule on the findings that
touched the notes.

## For the author — exposed by the run, not decided
1. **The bare deuteron reads the ⁶Li rank-2 transfer.** The shipped inclusive
   deuteron kernel's b₁ is (2/6) × `LI6_B1_RANK2_TRANSFER` × MillerB1 =
   0.32753 × MillerB1 — the ⁶Li α–d mapping applied to a deuteron. HERMES
   χ²(b₁)/6 is 14.56 on it against 5.26 on MillerB1 (the headline, named as
   convention 8's opt-in clause); it drives the two surviving deuteron
   top-cell refusals. *Recommendation:* give the deuteron its own b₁ slot
   (MillerB1 unmapped) and then re-open the toy b₁'s high-x behaviour, which
   a ~3× larger deuteron rank-2 term at the top cell will expose.
2. **The reference dump is no longer byte-reproducible here.** Run on
   2026-11-06, `dump_polligen_reference.py` moves 3170 committed numbers of
   six JSONs at the rounding level (≤ 2 ULP; `xsec.json` up to 9.4e−13
   relative among values above 1e−6, against the port gates' rtol 1e−12) with
   the sibling unchanged; the 2026-10-26 libc/libpython upgrade is the
   likely, unproven cause. The committed tables still pass every gate. One
   stale manifest string (the sibling's tagged.py line range) awaits the
   next real dump.
3. **BONuS on the AV18 weight:** the header CHECKs miss (|dv| 5.5e−4 /
   8.3e−4 vs 5e−4; bin averaging −3.2 % vs 2 %); recorded, no verdict flips.
4. **HERMES low-x bins** sit on the frozen CT18NLO edge (F₁ᵈ 12.40 vs the
   implied 7.04 at x = 0.012).
5. **The meta clip share counts cell centres only** (`g1_clipped_rate_fraction`).
6. **A10:** three C++ examples without a `Pipeline` write no `lipolgen.*`
   attributes; the C++ examples lack the two `*_sf_extrapolation` keys.
7. **A8:** `t5_delta_moment_bound` exits 0 (`EnvironmentBlocked`) for a
   present-but-wrong set version.
8. **The LO development item** (|A₁| > 1 from x ≳ 0.83 at LO) stands from run
   2026-10-27, as does EPPS21 as an option the inclusive rate does not apply.
9. **Throughput recipe:** 200 k-event runs read the inclusive rate 10 % low
   (clock ramp); use 2 M events next time.
