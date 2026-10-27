<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — one page

**What this run is.** The author asked for a review of the benchmark report,
a ranking of the channels by physics weight
(`../../benchmarking/08_channel_priority.md`) and the unfinished benchmark
tasks planned and executed (`PLAN.md`). On a machine that, unlike the
2026-09-26/27 one, has the three LHAPDF sets, the pinned PYTHIA 8.317,
`gfortran`, the network and the institution's APS licence, ten phases ran
(A, B1, B2, B3, C1, C2, C4, C5, D, P; C3 deferred with design) and the
integration stage (E) folded them in. The board, the files and the gate
tallies are `STATUS.md`; the priced options are `DECISIONS_PRICED.md`; the
regenerated report's `--compare` is `phase_E_report_compare.md`. **No physics
default, constant, tolerance or registry label moved, and no file under
`validation/reference/`; one vendored benchmark reference was corrected as
PLAN.md D2 allowed (Q_d 0.2859(3) → 0.2860(15) fm² in
`validation/benchmarks/data/deuteron_static_measured.json`, a dated
transcription correction, the old value kept, recorded only, no verdict
moved). No decision was taken, nothing was sent. Phases 0, 0r and B3 were
committed one by one, A–E together as `6b1b63e`, and phase F — four
adversarial lenses and the fixes their findings needed (`STATUS.md`
"Verification"), plus the §7 map extension — as `1c5cf68`; all pushed.**

## What it measured

- **The report, regenerated with every row run** (`docs/benchmarking/REPORT.md`,
  the M2 baseline): **23 harnesses — 13 pass, 10 fail (recorded), 0 blocked**
  (2026-09-27: 18 — 6 / 8 / 3, 1 blocked (environment)); 41 changes against
  the 2026-09-27 report, every one traced to a phase; regenerated again in
  phase F with **0 changes** against `6b1b63e`'s (the date, revision,
  runtimes, the plan-row legend sentence and the figure footers moved, no
  number); one figure per row that
  draws one (19 of 23, `make_report.py --plots`); `--check` passes.
- **The LHAPDF sub-rows, measured for the first time:** CT18NLO fails absolute
  NMC F₂ᵈ (3401.8/146, 8.6 % high; `LhapdfSF`'s massless leading-order charm
  and bottom are 3.6 % of the model F₂ᵈ on average, 7.0–7.6 % at x < 0.05,
  and without them it is 4.6 % high and x-dependent — phase F, recorded, not
  a verdict; the rest is unexplained, the NLO coefficient functions one candidate) while F₂ᵈ/F₂ᵖ passes (215.5/189); NNPDFpol1.1 fails world g₁ᵈ
  (COMPASS 45.0/15, HERMES 52.6/15, E143 51.0/28) and its Bjorken integral
  sits at 98 % of g_A/6 without C_Bj; NMC F₂(⁶Li)/F₂(D) passes again.
- **The three BLOCKED rows:** row 6 (Rand–Frosch–Yearian) **PASS** — ⁶Li
  F_M² against the tree's unfitted F_mag shape, 13.35/6 with one free scale
  (p = 0.038; s = 1 recorded at p = 0.0145); row 7 (DJANGOH 4.6.21 built,
  the elastic radiative tail isolated) **PASS** — `polrad-full` within 2.6 %
  of DJANGOH's exact proton tail in four Q² bins, the shipped `t-peak` 39 % /
  10 % low below Q² = 100 GeV² (a lower bound, as documented); row 10
  (BONuS) **FAIL, recorded** — the Glauber weight's backward variation
  exceeds the data's flat band in 24 of 28 DIS blocks, while the data fall
  toward 90° as the weight does.
- **Five new rows:** R against the world's R (row 12, **PASS** for the
  shipped `r_sigma_lt` on NMC R, 12.54/10, which however fails the SLAC sets);
  the gluonic Soffer-type bound on Δ's n = 2 moment (row 13, **PASS**: the
  shipped toy's generated amplitude uses 13–19 % of the room — DS16 Eq. (5)
  and the LOI's Eq. (6) make the kernel's cos 2φ coefficient twice theirs per
  unit Δ; 6–9 % under the first-pass identification A7, declared before the
  number and shown inconsistent with DS16's definition by the harness's
  second pass); HERMES F₂ᵈ (row 14, **FAIL** for every
  backend; HERMES's b₁ᵈ is per nucleon against HERMES's own F₂ᵈ);
  Wiringa–Schiavilla's ⁶Li form factors (row 15, **PASS one-sided**: the
  shipped C0 zero 3.0998 ≥ 2.8143 fm⁻¹); Boeglin's d(e,e′p)n (row 16,
  **PASS** for all three waves, only under A = 0.20, the upper end of the
  paper's 10–20 % FSI estimate, standing for FSI and the CC1 off-shell
  prescription together — at A = 0.10 the Hulthén fails, 51.57/21,
  p = 2.2 × 10⁻⁴, phase F, post hoc; a shape test, whose free scale absorbs
  the Hulthén's ×1.34 magnitude over CD-Bonn in the window).
- **The default-input price table of the ranked observables**
  (`validation/price_defaults.py`, N = 200 000): R (`--r-source r1998`,
  priced on the opt-in `li6-convolution` path; refused on the shipped `miller`
  b₁) moves claim 1's cos 2φ
  amplitude −15.2 / −14.7 / −7.8 / +4.2 % over the four sweet-spot super-bins
  x = 0.011 / 0.028 / 0.071 / 0.141 (−15.6 / −15.1 / −7.1 / +4.6 % at the cell
  centres, row 12 — twice the "−3 … −8 %" the ranking carried) at fixed Δ/F₁,
  the shipped toy's convention (its Δ ∝ F₁ moves with R); per unit Δ — the
  letter's extraction, `moment_B`, any calculation — R moves a₂ by
  ≤ 4.4 × 10⁻⁵ and F₂ moves it 1:1 (toy → CT18NLO ×1.36 / 1.19 / 0.96 / 0.90;
  phase F); A_zz
  ×1.25–1.26 under either grid backend;
  A_∥ ×0.51–0.53 under NNPDFpol1.1; A_zz^tag ×0.43 / ×0.21 under
  `--cluster-wave vmc` (the VMC AV18+UX α–d cluster wave — mostly its D state,
  P_D 0.0867 → 0.0194 — and the embedded AV18 deuteron); c₂ ×0.27 / ×0.16 at
  the measured quadrupole.
- **The least-verified inputs, read against their papers:** the g₁ᵈ D-state
  convention holds (all four sets per nucleon, uncorrected — the verdict
  column stands); Q_d's reference corrected to Bishop–Cheung's own
  0.2860 ± 0.0015 fm² (recorded only) — both it and the literature's
  0.2859(3) are superseded by Puchalski, Komasa, Pachucki, PRL 125 (2020)
  253001, Q_d = 0.285699(15)(18) fm², against which the recorded IA
  distances move by about +0.1 point (phase F; re-vendoring is the author's
  call); η_d and α_s(M_Z) = 0.1180(9) confirmed.
- **Gates:** C++ 416/416 cases, 17 242 712 assertions (phase E; no C++ file
  changed since). On the final tree (phase F): pytest **1456 passed / 155
  skipped / 0 failed**; with `LIPOLGEN_BENCH=1 -k bench` **351 passed / 1
  skipped**; the three citation-gate modes exit 0 (the records: 838
  citations, 0 broken); SPDX 133 files; `make_report.py --check` OK — the
  details in `STATUS.md` (phases E and F).
- **The verification (phase F):** four adversarial lenses, 38 findings
  (physics 8, one major; honesty 14; code 9; regression 7); none refuted; 33
  applied, 4 in part, 1 (the §7 map, E-2 / E-4) applied afterwards by
  the supervising session in the phase-F commit. The major one: "F₂ cancels in the cos 2φ amplitude and R is the
  lever" holds only for the shipped toy's Δ ∝ F₁ — the letter's extraction
  rides on F₂ (rank 1's C4 1 → 2, score 28 → 29, no rank moved).

## What it did NOT establish

- **Nothing validates the tensor observable itself** (`../../benchmarking/BENCHMARK_PLAN.md`
  §6, U-1 … U-12, unchanged): no polarized e + Li datum exists.
- Row 12's PASS is not evidence that `r_sigma_lt` is right (it fails
  Whitlow and E140, and NMC's correlated reading); row 13's bound says the
  shipped Δ is *allowed*, never that it is right — and whether the letter's
  scenario is allowed depends on a factor-2 normalisation question (below);
  row 15 is one-sided (the figure does not draw the C0 zero); row 16 passes
  only at the upper end of the paper's FSI estimate and as a shape (its free
  scale absorbs the magnitude); row 6 needs a free scale and does not reach F_m's
  second lobe; row 7 is one beam, one shared dipole, the proton only, with an
  unexplained −1.34 % in its top Q² bin; row 10's FAIL refutes the rule's
  premise, not the weight's direction.
- No ePIC chain leg ran; Abbott's deuteron A(Q²) / t₂₀ is deferred with design
  (no generation path reads a deuteron elastic form factor).
- No default was switched: the prices are inputs to the registry at M2.

## For the author

1. **The decisions priced** (`DECISIONS_PRICED.md`, none taken): the default R
   (claim 1: ×0.844–1.046 on the amplitude at fixed Δ/F₁; registry row 3,
   `Li6ConvolutionOptions`, is the separate convolution-path question); the
   default unpolarized and polarized backends (both fail their absolute
   data); D2 / row 9 (`--cluster-wave vmc` — the VMC AV18+UX α–d cluster wave
   and the embedded AV18 deuteron — moves A_zz^tag ×0.43 / ×0.21, which no
   p–n anchor tests; Hulthén's η is +17.2σ off);
   `eps_b0` (row 8, c₂ ×0.27); row 2's premise confirmed on HERMES's own F₂ᵈ;
   and the candidates — the M2 clause, a deuteron `TabulatedSpin1FF`, shipping
   the Wiringa–Schiavilla table, `TPeak` vs `PolradFull`, the RC tail table's
   +2 % interpolation bias.
2. **The factor 2 in Δ's normalisation** (`phase_C2_delta_bound.md` §7, §9):
   the cos 2φ kernel term carries Δ with twice the coefficient that DS16's
   definition (Eqs. (3)–(5)) and the LOI's own Eq. (6) give. The letter's bag
   scenario is therefore either 2× its stated amplitude, or 2× the bag value
   at 1.35–1.48 of the gluonic Soffer bound. Decide which. A copy of
   Jaffe–Manohar PLB 223 settles which printed formula is theirs. Beside it:
   the shipped Δ has the opposite sign of the letter's `moment_A` and is 4–7×
   smaller (`phase_A_baseline.md` §4.3).
3. **The documents to file in `PolarizedLithiumSim/refs/`** — fetched
   2026-10-26 into the session scratch directory, which is not permanent:
   the BONuS Supplemental Material (`BONuSprcSupMat2.txt`, sha256
   `18c7eb384c59c854e9e7926a1f96ec6840a9cc05674142bdfefee837b3f65f0e`),
   Rand–Frosch–Yearian (`7cbb2d21…31d5bbe2`), Suelzle–Yearian–Crannell,
   Rodning–Knutson, Bishop–Cheung, Close–Kumano, Sather–Schmidt, Moniz
   *et al.*, SMC PRD 58 112001, COMPASS PLB 769 (SCOAP3), PDG 2026 Table 1.1
   and the QCD review, and the DJANGOH 4.6.21 source and manual (each sha256
   in `DECISIONS_PRICED.md` §3 item 7 and the phase records).
4. **The licence question:** the BONuS supplemental tables and RFY's Tables I,
   II, VII, VIII, XIV are vendored as numbers (facts quoted with attribution)
   from APS-copyrighted documents read through the institution's licence, and
   so are two single values in `deuteron_static_measured.json`, η_d =
   0.0256(4) (PRC 41 898) and Q_d = 0.2860(15) fm² (PRA 20 381); the Licence
   lines of the WS98 and Boeglin CSVs say nothing of the papers is reproduced
   while their headers quote a caption and short sentences with attribution
   (to reword at their next revision); the
   DJANGOH cards vendored verbatim come from repositories that state no
   licence. Confirm, or replace them by a pointer and a hash.
5. **Report upstream, if wanted (nothing was sent):** PYTHIA 8.317's
   `Pythia8::loadExtraPython` execs every file of a `python/` directory in the
   current working directory at import (`phase_0_environment.md`; worked
   around once, in `_report.import_pythia8()`); HEPData ins894309's Table 2
   cross-section column repeats the proton's; COMPASS PLB 769 Table A.1 prints
   four high-x g₁ᵈ values that differ from HEPData's, the latter being the
   self-consistent ones (`phase_D_inputs.md` §1.3).
6. **Questions only the author can answer:** SMC in the g₁ᵈ verdict now that
   its paper is identified; NMC R's correlated reading as row 12's verdict (it
   would flip it to FAIL); HERMES F₂ᵈ's normalisation treated as one
   correlated component (row 14); WS98's C2 plotting convention (a 3/√2, ask
   R. B. Wiringa) and its 0.7116 normalisation of RFY's points; whether the
   tagged ⁷Li rate should read the struck polarization the Ciofi–Simula
   breakup assigns; re-declaring row 10's tolerance in a future run.
7. **Left for the author in code comments** (no C++ file changed this run):
   `include/lipolgen/rc.hpp`'s `RcTailModel::PolradFull` point (3) and the first
   bullet of "FIVE THINGS THIS IS NOT" (partly stale since row 7);
   `include/lipolgen/b1_nuclear.hpp`'s Q_d attribution; `tests/test_cluster.cpp`'s
   CD-Bonn η comment; `tests/test_rc.cpp`'s "NO 6Li table is shipped" comment.
   The exact texts are in `STATUS.md` ("What the run did not apply").
8. **Housekeeping:** phase D's `pip install --user pdg` also put SQLAlchemy
   2.0.55 and greenlet 3.5.6 into `~/.local`; the pytest gate's `polligen`
   fixture leaves gitignored `__pycache__/*.pyc` in the sibling repository.
9. **Phase F's items in files this run did not touch** (`STATUS.md`, "What
   the run did not apply"): E-3 (T5 for the Δ sector's inputs, a proposal;
   E-2 / E-4 were applied in the phase-F commit: `fsi.*` now re-runs the
   deuteron-control T2 rows and `rc.*` the T3 form-factor rows); the census in `python/tests/test_doc_link_gate.py` before
   `phase_C1_r_world.md`'s remaining citations can carry `(anchor read)` /
   `(as of)` annotations; matplotlib in CI (every `plot()` test skips
   there); re-vendoring Q_d at the 2020 value, and the two CSVs' Licence
   lines, at those files' next revision.
