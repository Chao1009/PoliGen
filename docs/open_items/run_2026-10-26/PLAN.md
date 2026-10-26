<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — plan: the channel ranking, and the unfinished benchmark work it orders

The author asked for three things: a review of the benchmark report
(`../../benchmarking/REPORT.md`, 2026-09-27), a ranking of the channels by
physics importance so that benchmarking and development follow it
(`../../benchmarking/08_channel_priority.md`, written first and reviewed by
two adversarial readers before this plan was frozen), and the unfinished
tasks planned and executed — the execution by Opus 5.5 workers, the planning,
integration and verification by the supervising session. Progress is on
`STATUS.md` beside this file; the one page is `SUMMARY.md`.

## Baseline (measured at the start, `0cea6c0` = `master`, tree clean)

A fresh checkout on a machine that had **no dependency stack and no build**
(`../deps/install` absent; Python 3.10.12; numpy 2.2.6; scipy 1.15.3; gcc
11.4.0; 4 cores). The 2026-09-26 report was made on a machine with **no
LHAPDF set and PYTHIA 8.312**, and says it is to be regenerated where the
sets are installed. So phase A builds the stack the project pins
(`.github/scripts/build_deps.sh`: HepMC3 3.3.0, LHAPDF 6.5.5 with CT18NLO,
NNPDFpol11_100 and EPPS21nlo_CT18Anlo_Li6, PYTHIA 8.317) and measures the
baseline on it. Three things this machine has that the 2026-09-26 one had
not: the network (HEPData's JSON API answers the standard Python client;
arXiv, INSPIRE and GitHub are reachable), **institutional access to APS**
(the BONuS Supplemental Material and Phys. Rev. 144 (1966) 859 download with
curl's own identity — the two documents rows 10 and 6 were BLOCKED on), and
`gfortran` (a DJANGOH build is possible — row 7). Elsevier (ADNDT) does not
serve a PDF to a non-browser client and is not pursued.

## Ground rules

Those of `../run_2026-09-03/PLAN.md` and `../run_2026-09-06/PLAN.md`,
unchanged: **no default moves; every decision stays in the registry
(`../run_2026-09-03/STATUS.md`, cite by row number); a benchmark that
disagrees with the tree is recorded, never tuned; every tolerance is declared
in the harness header before the verdict and the header says when the number
was previewed first; a knob that did not run is labelled, never recorded as
run; a count that does not reproduce from committed state is not published;
nothing is pushed, published or sent.** Three rules of this run's own:

1. **No client-identity spoofing, ever.** A document that a plain `curl` or
   the standard Python `requests` client cannot fetch is recorded as not
   reachable and named for the author, as the 2026-09-23 run did; nothing is
   fetched by pretending to be a browser. The APS PDFs above came back to
   curl's own identity through the institution's licence and are kept
   **outside the repository** (`PolarizedLithiumSim/refs/` is the author's
   corpus; the session scratch directory holds them until then). Only
   **numbers** are vendored, with provenance, under
   `validation/benchmarks/data/` — the publishing convention of every
   harness already in the tree.
2. **A worker owns the files it is given and no other.** Shared files
   (`BENCHMARK_PLAN.md`, `validation/benchmarks/README.md`, `README.md`,
   `REFERENCES.md`, the report) are edited by the integration stage alone,
   from the workers' recorded results.
3. **Every new harness is a complete row**: the file, its vendored data with
   a sha256 re-checked on every read, a `python/tests/test_bench_*.py` with
   the file guards (point count, first and last value, one flipped byte
   refused) and the exit-code tests, a `REPORT |` line, `tiers.py` coverage,
   and a header that states the three normalisation rules for that row.

## Figures — one per row, so the report grows a plot at each step

The author asked (2026-10-26) for a benchmarking plot at the end of each
step where one applies, concise, so that the report is built up
incrementally. The rule: **every harness whose comparison can be drawn
implements `plot(report, out_dir)` through `validation/benchmarks/_plots.py`
(one style: the reference as ink points with errors, each generator
configuration in its fixed palette slot, a pull or ratio panel where the
tolerance is a χ², the status in the title) and writes
`docs/benchmarking/figures/<harness>.png`; `make_report.py --plots` calls it
after `run()` and links the figure in the row's Detail.** Figures are a view
of REPORT.json, never a record (PNG bytes are not in `--check`). Applies to:
every data comparison (T1–T4), the price table (A3: σ ratios per channel and
backend), the positivity saturation (C2), the Deeps/BONuS shapes, the η/Q_d
and r_ch comparisons (a dot-and-band plot), the DJANGOH per-bin ratios, the
PYTHIA closure ratios. Does not apply to: EPIOS/EST (identities at 1e−14),
weighted/unweighted (one pull), T6 (a clause table). A new worker **P**
retrofits the figures of the 18 existing harnesses and wires `--plots`.

## Phase A — the baseline the plan has been asking for since 2026-09-26

- **A0** (done by the supervising session before the workers started;
  `phase_0_environment.md`): the stack built, LiPolGen built, the C++ suite
  clean (416/416), the first report with the LHAPDF sets, and one
  environment defect found and worked around — PYTHIA 8.317's Python
  module executes every file of a `python/` directory in the current
  working directory at import, so `import pythia8` fails from the
  repository root; `_report.py` `import_pythia8()` is the one workaround,
  used by T4, T6 and their test. The T4 full-row test still pins numbers
  taken on 2026-09-26 with the unpinned PYTHIA 8.312; on the pinned 8.317
  the row passes (worst 0.076 against 0.20) and the pin is re-taken by the
  integration stage with a dated note.
- **A1** Build the stack into `../deps/install` with the project's own
  script; configure and build LiPolGen; run `build/lipolgen_tests` and
  `python -m pytest python/tests -q`; record the tallies against the
  2026-09-27 ones (416 cases / 403 passed / 13 failed — all 13 "Info file not
  found for PDF set"; 39 failed / 1209 passed / 163 skipped — all 39 the same
  cause). **Expected here: 0 of either kind of failure.** Any other failure
  is a finding.
- **A2** `LIPOLGEN_BENCH=1 make_report.py --all --compare docs/benchmarking/REPORT.json`:
  the first numbers for row 4 on this tree's own build, the CT18NLO and
  NNPDFpol11 sub-rows of the three T1 harnesses, and
  `t5_sum_rules[bjorken[nnpdfpol11]]`; PYTHIA 8.317 for T4 and T6. Every
  status change and headline move is a line of `phase_A_baseline.md`.
- **A3** The per-channel **rate and default-input price table** that
  `08_channel_priority.md` §5.1 asks for: σ per channel at ⁶Li / ⁷Li config 1
  under the shipped defaults and under each opt-in backend
  (`--unpol-sf {toy,mstw,ct18nlo}`, `--pol-sf {toy,nnpdfpol}`,
  `--r-source`, `--cluster-wave`), from `meta` of short runs — numbers, not
  decisions.

## Phase B — the three rows that were BLOCKED on a document

- **B1 — row 10, BONuS** (`t2_bonus_spectator_shape.py`): the Supplemental
  Material of Tkachenko *et al.*, PRC 89 (2014) 045206, is 100 blocks of
  R_D/S (data over BONuS's own spectator simulation) against cos θ_pq at
  3 Q² × 5 W* × 4 p_s (0.078–0.135 GeV/c) for the 4 and 5 GeV beams, plus
  the F₂ⁿ/F₂ᵈ and F₂ⁿ/F₂ᵖ tables. Vendor the R_D/S blocks with provenance and
  the sha256; wire the harness as the comparison its own header anticipates:
  in the backward PWIA window the ratio's flatness bounds what any FSI or
  off-shell effect may do, and LiPolGen's `GlauberFsiWeight` variation across
  that window on the deuteron control at the BONuS p_s is held to it (a
  must-not-contradict rule, declared before the number); the forward rise is
  recorded against the weight's direction. State what R_D/S can and cannot
  say (it is a ratio to a simulation, light-cone-corrected and
  acceptance-folded). Tests. The row's status moves from BLOCKED to a
  measurement.
- **B2 — row 6, Rand–Frosch–Yearian** (`t3_li_magnetization_rfy.py`): from
  the fetched Phys. Rev. 144 (1966) 859 — Table I (⁶Li) and Table II (⁷Li)
  are the measured magnetic form factors against q, Tables VII/VIII the
  fitted analyses — transcribe the lithium rows into
  `rfy1966_li_magnetization.csv` in the header's own format (both ⁶Li fits
  and every ⁷Li row), reading the PDF pages as images where the OCR text is
  ambiguous and recording page and table for every number. Run the harness
  unchanged; add the measured F_M(q) points as a second, recorded sub-row
  against the tree's `(1 − q²/q_z²) exp(−q²b²/4)` shape. Tests.
- **B3 — row 7, DJANGOH with the elastic tail on**
  (`t4_djangoh_rad_noRad.py`): build DJANGOH 4.6.x from the public source
  (`github.com/spiesber/DJANGOH`, self-contained Fortran plus LHAPDF, which
  the deps prefix provides), run e⁻ 18 × p 275 NC Rad = 1 and Rad = 0 with
  IEL2, IEL31–IEL33 ≠ 0 **and** WMIN ≤ M_p, at a sample size whose
  statistical error on each four-bin ratio is below the ratio's distance
  from 1, and tabulate σ per Q² bin. Vendor the table with the full input
  cards, the DJANGOH version and commit, the LHAPDF set, the seeds and the
  logs' sha256. Then re-state the comparison: with the elastic tail on, the
  two sides share the O(α) elastic term and the harness can define a
  tolerance for the first time — declare it in the header before the
  number, from the term-by-term accounting the header already carries
  (DJANGOH adds the inelastic lepton-vertex pieces LiPolGen does not), or
  record why no tolerance is definable even now. Keep the ePIC-sample rows as
  the recorded history.

## Phase C — new anchors for the top-ranked sectors (ranks 1–2, ingredients A, C, F)

- **C1 — R data** (`t1_r_world.py`, new): H1 F_L (`10.17182/hepdata.45340`,
  the only R inside the EIC x window), Whitlow R1990 (`…29555`) and SLAC
  E140 (`…22468`) from HEPData's JSON API; `r_sigma_lt` and `r1998` (clipped
  and unclipped) against them, p(χ², ndf) ≥ 0.01 per data set on the DIS
  cut, stat ⊕ sys diagonal, full covariance recorded where the record gives
  it. The one-page bias table ingredient A asks for: (1 + r1998)/(1 +
  r_sigma_lt) across the window. In the same worker, because the vendoring
  pattern is identical: HERMES F₂ᵈ (`…66147`) as a sub-row of `t1_nmc_f2d`'s
  family and the final COMPASS g₁ᵈ (`…78374`) as a sub-row of
  `t1_g1d_world` (both `02` §10 ranks 2 and 4), each with its publishing
  convention stated.
- **C2 — positivity and the Soffer bound** (`t5_positivity_bounds.py`,
  new): the spin-1 positivity bounds on the b₁, b₂ (and the b₃/b₄ scenario
  shapes) against F₁, F₂ and the gluonic Soffer bound on `toy_delta_gluon`
  over the generator's (x, Q²) grid — closed forms, zero free parameters,
  the first gate of any kind on the Δ sector (`04_theory.md` T-13, T-18;
  the bound's source stated and its transcription checked against the
  paper). A violation is a recorded FAIL at the cell, never a clip.
- **C3 — the deuteron form factors in the RC control**
  (`t2_deuteron_ff_abbott.py`, new): Abbott *et al.* t₂₀ / A(Q²)
  (`10.17182/hepdata.40433`, 5 tables; `…31372`) against the shape-blind HO
  stand-in `HoSpin1FF` for the deuteron (`02` §5.3): a recorded shape
  comparison at the measured Q², and the elastic-tail fraction's sensitivity
  to it on the deuteron control. The Abbott EPJ A 7 (2000) 421
  parametrizations are in the sibling corpus (`refs/nucl-ex_0002003.pdf`)
  if a closed-form world fit is wanted beside the points.
- **C4 — the ⁶Li elastic form factors at the theory level**
  (`t3_li6_ff_wiringa_schiavilla.py`, new, lower priority): digitize
  Wiringa–Schiavilla's F_C0 / F_C2 / F_M1 (`refs/nucl-th_9807037.pdf`, a
  vector figure; the sibling's `tools/digitize_figure.py` with its
  calibration check) into the empty `data/ff/li6_elastic.csv` slot's format
  and compare `HoSpin1FF`'s C0 zero and slope and the F_M shape of B2 to it —
  closes `04` T-20 at the theory level; the data level stays rows 5 and 6.
- **C5 — d(e,e′p)n momentum distribution** (`t2_boeglin_nk.py`, new, lower
  priority): Boeglin *et al.* PRL 107 (2011) 262501
  (`refs/1106.0275.pdf`) reduced cross sections in the 35°–45° θ_nq window,
  digitized with a stated error, against the three deuteron waves' n(k) out
  to 0.55 GeV/c — T25's model-vs-model made model-vs-data.

## Phase D — the "least verified inputs" read against their papers

- **D1** The T1 references: whether the published g₁ᵈ of COMPASS
  (`refs/1612.00620.pdf`), HERMES (`refs/hep-ex_0609039.pdf`) and E143
  (`refs/hep-ph_9802357.pdf`) is D-state corrected, which R NMC used, whether
  NMC's normalisation is inside "sys", E143's beam energies and the SMC
  paper — each resolved by reading the paper and written into the harness
  header in place of UNVERIFIED, with the page. Where a convention changes
  the comparison (the D-state factor is 7 % on the isoscalar combination),
  the harness grows the corrected sub-row beside the existing one; the
  verdict column is not re-chosen after the fact.
- **D2** η_d and Q_d: Rodning–Knutson PRC 41 (1990) 898 and Bishop–Cheung
  PRA 20 (1979) 381 are fetched. Read the values the papers actually
  conclude — the first grep of Bishop–Cheung shows **0.2860 ± 0.0015 fm²**
  in its text where the tree transcribes 0.2859(3) — and record what each
  number is (a measurement, a theoretical value, a re-evaluation), with the
  page; correct the vendored `deuteron_static_measured.json` only if the
  transcription is wrong, and then as a dated correction with the old value
  kept in the provenance. Same for α_s(M_Z) in `t5_sum_rules` (PDG 2026).
  Also Suelzle 1967 (`PhysRev.162.992`, fetched) for the ⁷Li HO parameters
  `03` §7 rank 6 names.

## Phase E — documents, report, decisions

- **E1** `make_report.py --all` then `--check`; `REPORT.md` / `REPORT.json`
  regenerated on this machine (the M2 baseline moves — that is the point;
  `--compare` against the committed file lists every move).
- **E2** `BENCHMARK_PLAN.md` §3 (M1 status note dated here) and §4 (rows 6,
  7, 10; the new rows numbered 12–); `validation/benchmarks/README.md`;
  `README.md`'s benchmark line; `REFERENCES.md` items 8, 9, 59 and the
  fetched 1a items; `PHYSICS_CHANNELS.md` only where a sentence became stale.
- **E3** `DECISIONS_PRICED.md` here: the registry rows this run makes
  decidable from numbers (the default unpolarized backend; `r_sigma_lt`
  vs `r1998`; D2's Hulthén at M2; row 8's `eps_b0`; row 2 on HERMES) — each
  with its options and measured cost, **none taken**.
- **E4** `STATUS.md`, `SUMMARY.md`, this file's status notes; the
  `08_channel_priority.md` §7 log entry if any ranking changed.

## Phase F — verification

Four adversarial lenses over the uncommitted tree (physics; honesty — every
number traced to a window, every verdict to a declared tolerance; code —
tests, exit codes, sha256 guards, the gates; regression — the 18 existing
harnesses and every gate bit-identical where no input changed), then a fix
stage, then the report regenerated once more. The gates: `./build/lipolgen_tests`,
`python -m pytest python/tests -q`, `LIPOLGEN_BENCH=1 python -m pytest
python/tests -q -k bench`, `python3 validation/check_physics_channels_links.py`
(strict), `--records-only`, `--audit-ranges`, `validation/check_spdx_headers.py`,
`make_report.py --check`.

## Not in this run, and why

- **The 27 registry decisions.** Priced where this run produces a price
  (E3); taken by nobody here.
- **The ePIC chain legs** (abconv, npsim, EICrecon, Rivet): a container stack
  this machine does not have mounted; they validate the record, not the
  physics (`08_channel_priority.md` §5.4); M3.
- **The Mäntysaari-group letter and the dipole-model run** (row 21): outward.
- **The sibling repository** (`PolarizedLithiumSim`): read, never edited.
- **Default switches** suggested by the ranking (the unpolarized backend,
  `eps_b0`, the deuteron wave): the price is this run's deliverable; the
  switch is the author's at M2.
