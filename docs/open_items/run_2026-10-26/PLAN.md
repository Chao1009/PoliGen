<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — plan: the channel ranking, and the unfinished benchmark work it orders

The author asked for three things: a review of the benchmark report
(`../../benchmarking/REPORT.md`, 2026-09-27), a ranking of the channels by
physics importance so that benchmarking and development follow it
(`../../benchmarking/08_channel_priority.md`), and the unfinished tasks
planned and executed — the execution by Opus 5.5 workers, the planning,
integration and verification by the supervising session; every substantial
step committed and pushed; a concise benchmarking figure at the end of each
step where one applies. The ranking and this plan were reviewed by two
adversarial readers before the workers started (49 findings; `STATUS.md`
phase 0r) and revised. Progress is on `STATUS.md`; the one page is
`SUMMARY.md`.

## Baseline (measured at the start, `0cea6c0` = `master`, tree clean)

A fresh checkout on a machine that had **no dependency stack and no build**
(`../deps/install` absent; Python 3.10.12; numpy 2.2.6; scipy 1.15.3; gcc
11.4.0; 4 cores). The 2026-09-26 report was made on a machine with **no
LHAPDF set and PYTHIA 8.312**, and says it is to be regenerated where the
sets are installed. Phase 0 built the stack the project pins
(`.github/scripts/build_deps.sh`: HepMC3 3.3.0, LHAPDF 6.5.5 with CT18NLO,
NNPDFpol11_100 and EPPS21nlo_CT18Anlo_Li6, PYTHIA 8.317) and measured the
baseline on it (`phase_0_environment.md`). Three things this machine has
that the 2026-09-26 one had not: the network (HEPData's JSON API answers
the standard Python client; arXiv, INSPIRE and GitHub are reachable),
**institutional access to APS** (the BONuS Supplemental Material and Phys.
Rev. 144 (1966) 859 download with curl's own identity — the two documents
rows 10 and 6 were BLOCKED on), and `gfortran` (a DJANGOH build is possible
— row 7, which was blocked on a run, not a document). Elsevier (ADNDT) does
not serve a PDF to a non-browser client and is not pursued.

## Ground rules

Those of `../run_2026-09-03/PLAN.md` and `../run_2026-09-06/PLAN.md`,
unchanged: **no default moves; every decision stays in the registry
(`../run_2026-09-03/STATUS.md`, cite by row number); a benchmark that
disagrees with the tree is recorded, never tuned; every tolerance is declared
in the harness header before the verdict and the header says when the number
was previewed first; a knob that did not run is labelled, never recorded as
run; a count that does not reproduce from committed state is not published;
nothing is published or sent.** Rules of this run's own:

1. **No client-identity spoofing, ever.** A document that a plain `curl` or
   the standard Python `requests` client cannot fetch is recorded as not
   reachable and named for the author, as the 2026-09-23 run did. The APS
   documents came back to curl's own identity through the institution's
   licence. **The author's request of 2026-10-26 ("plan for the unfinished
   tasks, then execute") is read as the release of what the 2026-09-27
   decisions file left "to the author" — the blocked rows 6, 7, 10 and the
   documents they wait for.** The documents are kept **outside the
   repository** (the session scratch directory; filing them into the
   author's corpus `PolarizedLithiumSim/refs/`, and whether the licence
   allows the supplemental tables to be vendored, are listed for the author
   in `SUMMARY.md`); each one's DOI and sha256 is recorded in the vendored
   file's provenance. Only **numbers** are vendored, with provenance, under
   `validation/benchmarks/data/` — the publishing convention of every
   harness already in the tree.
2. **A worker owns the files it is given and no other.** Shared files
   (`BENCHMARK_PLAN.md`, `validation/benchmarks/README.md`, `README.md`,
   `REFERENCES.md`, the report, `make_report.py`) are edited by the
   integration stage (or the figures worker, for its own code paths) alone,
   from the workers' recorded results.
3. **Every new or re-wired harness is a complete row, and every newly
   vendored file is guarded**: the file, its sha256 re-checked on every read
   (exit 2 on a mismatch), a `python/tests/test_bench_*.py` with the file
   guards (point count, first and last value, one flipped byte refused) and
   the exit-code tests, a `REPORT |` line, `tiers.py` coverage, a header that
   states the three normalisation rules for that row, and the README's
   expensive-row rule honoured (> 30 s, or any generation /
   generator-vs-generator sample, behind `LIPOLGEN_BENCH=1`).
4. **A tolerance declared after its number was seen is not a tolerance.**
   Where this run has already previewed a number (the BONuS R_D/S tables were
   read by the supervising session; r_mag(⁶Li) = 2.9469 fm is printed in the
   committed report; the T4 worst distance 0.076; the CT18NLO and NNPDFpol11
   sub-rows), the harness header says so, and the rule is written in terms
   of what the reference states, before the comparison is computed.
5. **Cut order.** Workers run in rank order of `08_channel_priority.md`
   (A, C1, C2, B1, B2, P, D, C4/C5), with B3 in its own stream; if the run
   is cut short, what is finished is integrated and the rest is listed.

## Figures — one per row, so the report grows a plot at each step

**Every harness whose comparison can be drawn implements `plot(report,
out_dir)` through `validation/benchmarks/_plots.py` (one style: the
reference as ink points with errors, each generator configuration in its
fixed palette slot, a pull or ratio panel where the tolerance is a χ², the
status in the title) and writes `docs/benchmarking/figures/<harness>.png`;
`make_report.py --plots` calls it after `run()` and links the figure in the
row's Detail.** Figures are a view of REPORT.json, never a record (PNG bytes
are not in `--check`). Applies to: every data comparison (T1–T4), the price
table (A3), the Δ bound's saturation (C2), the Deeps/BONuS shapes, the η/Q_d
and r_ch comparisons (dot-and-band), the DJANGOH per-bin ratios, the PYTHIA
closure ratios. Does not apply to: EPIOS/EST (identities at 1e−14),
weighted/unweighted (one pull), T6 (a clause table). Worker **P** retrofits
the figures of the existing rows and wires `--plots`.

## Phase A — the baseline the plan has been asking for since 2026-09-26

- **A0** Done by the supervising session before the workers started
  (`phase_0_environment.md`): the stack built, LiPolGen built, the C++
  suite clean (416/416; 17 242 712 assertions), pytest 1269 passed / 158
  skipped / 0 failed, the first report with the LHAPDF sets, and one
  environment defect found and worked around — PYTHIA 8.317's Python module
  executes every file of a `python/` directory in the current working
  directory at import, so `import pythia8` fails from the repository root;
  `_report.py` `import_pythia8()` is the one workaround, used by T4, T6 and
  their test. The T4 full-row pin was re-taken on the pinned 8.317 (worst
  distance 0.076, W2 n_ch_eta[−2,0), a low-count bin; 0.023 on the unpinned
  8.312 of 2026-09-26; the row passes at 0.20) — recorded in
  `phase_0_environment.md` with the T4 sub-rows copied from the harness.
- **A1** The tallies and the `--compare` lines, written up by worker A as
  `phase_A_baseline.md`, with what the CT18NLO / NNPDFpol11 verdicts mean for
  ingredients A and H.
- **A2** (folded into A1) every status change and headline move, one line
  each.
- **A3** The **default-input price table of the ranked observables**, not of
  σ alone: per channel, under each knob **the knob-provenance table reports
  as read on that channel and plan** (refusals and `not read` labels
  recorded as such, never as numbers): inclusive ⁶Li — the cos 2φ amplitude
  at the four sweet-spot cells and the run-level A_zz under `--unpol-sf
  {toy,mstw,ct18nlo}` (tensor-thirds / transverse-tensor plans; `--pol-sf`
  is not read there — say so; `--r-source` is accepted only with `--b1-model
  li6-convolution`, so price it on that path and label it); inclusive ⁷Li —
  A_∥ under `--pol-sf {toy,nnpdfpol}` and `--unpol-sf` on `helicity-flip
  --pe 0.7`; tagged ⁷Li α — A_∥ and the spectator spectrum's mean k under
  `--cluster-wave`, `--triton-sf`, `--pol-sf`; tagged ⁶Li α — A_zz^tag(k)
  at k = 0.20 / 0.30 GeV under `--cluster-wave`; coherent ⁶Li — the
  generated c₂ under `eps_b0` at the shipped −0.08 and the measured-quadrupole
  −0.0070 and under `--coherent-f0` 0.02 / 0.08 (`eps_b0` through the
  config/API, it has no CLI flag); plus σ per channel as the rate line.
  Numbers from `lg.run` meta and the returned columns, N = 200 000, seed
  20261026; one figure (`price_defaults.png`). No default moves.

## Phase B — the three BLOCKED rows (6, 10 on a document; 7 on a run)

- **B1 — row 10, BONuS** (`t2_bonus_spectator_shape.py`): the Supplemental
  Material of Tkachenko *et al.*, PRC 89 (2014) 045206, is 100 blocks of
  R_D/S (data over BONuS's own spectator simulation) against cos θ_pq: at
  the 4 GeV beam 3 Q² (0.93, 1.66, 3.38 GeV²) × 5 W* (1.17 … 2.44 GeV) × 4 p_s
  (0.078–0.135 GeV/c) = 60 blocks, at 5 GeV 2 Q² × 5 W × 4 p_s = 40; plus the
  F₂ⁿ/F₂ᵈ and F₂ⁿ/F₂ᵖ tables. The existing header's "6 W* × 5 Q² × 4 p_s" is
  wrong and is corrected. Vendor the blocks with provenance and the sha256.
  **The comparison is a new design, argued as such**: read the paper for
  what the denominator simulation contains; **declare a DIS window before
  reading any block** (W* ≥ 2 GeV and Q² ≥ 1 GeV², above the generator's
  own scenario floor; the resonance-region blocks recorded as out of model);
  then, in the backward PWIA window (cos θ_pq ≤ −0.3), the ratio's
  consistency with flat bounds any FSI/off-shell effect, and LiPolGen's
  `GlauberFsiWeight` (a library evaluation on the deuteron control at the
  BONuS p_s, not a generation) is held to that flatness band — the rule and
  its tolerance written in the header before the number, with the note that
  the supervising session previewed the tables; the forward rise recorded
  against the weight's direction. State what R_D/S can and cannot say.
  Tests (a new test file; the old `test_bonus_*` skips re-pointed by the
  integrator). The row's status moves from BLOCKED to a measurement.
- **B2 — row 6, Rand–Frosch–Yearian** (`t3_li_magnetization_rfy.py`): the
  paper (fetched) does not print magnetization radii — Table VII (⁶Li) and
  VIII (⁷Li) give harmonic-well parameters (a₀ = 2.16 ± 0.04 F free, 1.78 F
  from charge scattering, "any value in 1.6 < a₀ < 2.2 is possible") and
  Tables I and II the measured F_M points against q. **The primary
  comparison is therefore the measured F_M(q) points of Tables I and II
  against the tree's F_mag shape** — transcribed from the page images,
  vendored as `rfy1966_li_magnetic_ff.csv` with a sha256 check added to the
  harness — with a tolerance declared first (a shape χ² with one free scale,
  the paper's own errors; the ⁷Li half "not implemented in the tree",
  recorded). The radius comparison the harness was written for is kept as a
  **recorded** sub-row only if the a₀ → ⟨r²⟩_mag formula, the centre-of-mass
  and nucleon-size choices and which RFY error enters are all written in the
  header before the derivation, with the note that the tree's 2.9469 fm was
  previewed in the committed report; otherwise a₀ is vendored verbatim and
  the radius rule retired with a dated note.
- **B3 — row 7, DJANGOH, the elastic radiative tail isolated**
  (`t4_djangoh_rad_noRad.py`): build DJANGOH 4.6.x from the public source
  (`github.com/spiesber/DJANGOH` at `e356b84`, self-contained Fortran plus
  LHAPDF from the deps prefix). **Configuration and tolerance in the header
  before any run**: IEL2 stays 0 (the manual defines it as the
  non-radiative elastic peak, not a tail); the elastic radiative tail is
  isolated either by IEL31–IEL33 on with the inelastic channels INC31–INC34
  off, or by an IEL-on / IEL-off pair at identical cuts whose difference over
  the Born is the tail; WMIN ≤ M_p so the tail's vertex is inside the
  acceptance, and the generator's W_l² ≥ 9 GeV² applied on the DJANGOH side
  by cut or by selection so that the two windows match; the PDF set and the
  integration points fixed in advance (control first: reproduce the ePIC
  samples' eight numbers with IEL = 0, WMIN = 3). The compared quantity is
  DJANGOH's elastic-tail-only ratio per Q² bin against LiPolGen's t-peak and
  polrad-full tails on the proton, with the tolerance taken from the
  header's own POLRAD Eq. (44) vs HERACLES accounting (the leading-log vs
  exact difference, ~5–10 %, or tighter if the manual states DJANGOH's
  elastic tail is exact Mo–Tsai). The full Rad/noRad ratio with the tail on
  is recorded history, never the comparison (it is dominated by inelastic
  lepton radiation LiPolGen does not have). If no isolating configuration
  exists, the row stays BLOCKED with the new reason — a documented negative
  is a result.

## Phase C — new anchors for the top-ranked sectors (ranks 1 and 4; ingredients A, C)

- **C1 — R data** (`t1_r_world.py`, new): **NMC's R** (inside the NMC F₂
  record, `10.17182/hepdata.32752`, 0.002 < x < 0.12 — the one data set that
  covers the letter's sweet spots) as the primary sub-row; H1 F_L
  (`…45340`, 0.00024 < x < 0.0036) as the low-x extrapolation check; Whitlow
  R1990 (`…29555`) and SLAC E140 (`…22468`) labelled as `r1998`'s own fit
  input (a closure, not a test). `r_sigma_lt` and `r1998` (clipped and
  unclipped) at every point; p(χ², ndf) ≥ 0.01 per data set on the DIS cut,
  stat ⊕ sys diagonal, full covariance recorded where given; the bias table
  (1 + r1998)/(1 + r_sigma_lt) at the four sweet-spot cells and across the
  window. In the same worker: HERMES F₂ᵈ (`…66147`) as `t1_hermes_f2d.py`,
  the datum that makes "HERMES's b₁ is per nucleon" checkable (`02` §10
  rank 2). The final COMPASS g₁ᵈ is **already** the `compass` sub-row of
  `t1_g1d_world` (ins1501480 table 1 via nnpdf-data) — not duplicated; a
  cross-check that the HEPData record equals the vendored table is optional.
- **C2 — the one bound on Δ** (`t5_delta_moment_bound.py`, new): the
  gluonic Soffer-type bound on the n = 2 moment (Detmold–Shanahan PRD 94
  (2016) 014507 Eq. (36), `refs/1606.04505.pdf`, with the Jaffe–Manohar
  n = 2 relation between Δ's moment and the gluon transversity moment),
  applied to the shipped `toy_delta_gluon` at scale 10⁻² (the pipeline's
  default, not the function's 10⁻³) and to the letter's `moment_A` /
  `moment_B` models (sibling `polli_fastsim/delta_models`, read only), with
  **every assumption declared in the header before the number**: the gluon
  PDF set used for f₁^g (CT18NLO at the deps prefix), the ⁶Li gluon-helicity
  model (NNPDFpol11 Δg per nucleon times the effective polarization, or
  Δg = 0 as a band edge), f₁LL^g = 0 (or a band). The result is **a bound
  under assumptions**: PASS/FAIL per Q², the saturation fraction recorded.
  Not in this row: T-13 (spin-3/2, trivially satisfied on a zero sector),
  b₃/b₄ (higher-twist slots), and any "positivity" without a cited spin-1
  source — if a spin-1 helicity-positivity source is found on disk (Cosyn
  *et al.* 2410.12764 or HJM), the b₁ bounds go in as a second, sourced
  sub-row; otherwise "bound not sourced" is written.
- **C3 — Abbott t₂₀ / A(Q²)** — **deferred with design.** No generation
  path reads a deuteron elastic form factor (`rc_tail` ≡ 1 on tagged
  channels; the inclusive deuteron `HoSpin1FF::for_ion(DEUTERON())` throws),
  so a comparison would move nothing and `HoSpin1FF`'s t₂₀ fails by
  construction (one monopole for F_c and F_q). The item is a development
  one: a deuteron `TabulatedSpin1FF` from the Abbott EPJ A 7 (2000) 421
  parametrization, which would also unblock the inclusive-deuteron RC run —
  filed in `DECISIONS_PRICED.md` as a registry candidate.
- **C4 — the ⁶Li elastic form factors at the theory level**
  (`t3_li6_ff_wiringa_schiavilla.py`, new): digitize Wiringa–Schiavilla's
  F_C0 / F_C2 / F_M1 (`refs/nucl-th_9807037.pdf`, a vector figure; the
  sibling's `tools/digitize_figure.py`) and vendor the curves **only under
  `validation/benchmarks/data/`** in the `data/ff/li6_elastic.csv` slot's
  columns — `data/ff/` stays empty (`tests/test_rc.cpp` asserts nothing is
  shipped there; shipping it is a registry item). Compare `HoSpin1FF`'s C0
  zero and slope and the F_M shape to the curves; must-not-contradict
  tolerance declared first; closes `04` T-20 at the theory level only.
- **C5 — d(e,e′p)n momentum distribution** (`t2_boeglin_nk.py`, new, last):
  Boeglin *et al.* PRL 107 (2011) 262501 (`refs/1106.0275.pdf`) reduced
  cross sections in the 35°–45° θ_nq window, if the figure is a vector
  drawing; **the paper's own 10–20 % FSI and the de Forest CC1 off-shell
  prescription declared as an allowance band before the verdict**, or the
  comparison made as ratios to the paper's own PWIA curve — against the
  three deuteron waves' n(k) out to 0.55 GeV/c.

## Phase D — the "least verified inputs" read against their papers

- **D1** The T1 references: **the harness header already makes the
  verdict conditional on the D-state convention** ("IF that holds, g1d_tree
  (with the factor) is the like-for-like comparison, which is why it carries
  the verdict") and already records the other reading as `dis_no_dstate`.
  So: declare, before reading, that the verdict follows the convention the
  papers state; read COMPASS (`refs/1612.00620.pdf`), HERMES
  (`refs/hep-ex_0609039.pdf`), E143 (`refs/hep-ph_9802357.pdf`, and whether
  hep-ex/9705012 is the right identifier), NMC's R and normalisation, SMC;
  write each resolution into the header with the page in place of
  UNVERIFIED; if the papers put the verdict on the `dis_no_dstate` column,
  swap the verdict column with a dated note (no new sub-row); record that
  phase 0 previewed NNPDFpol11 failing on the current column.
- **D2** η_d and Q_d: Rodning–Knutson PRC 41 (1990) 898 and Bishop–Cheung
  PRA 20 (1979) 381 are fetched. Read what each paper CONCLUDES (not merely
  quotes) and what it is — the first grep of Bishop–Cheung shows
  0.2860 ± 0.0015 fm² in its text where the tree transcribes 0.2859(3);
  correct `deuteron_static_measured.json` only if the transcription is
  wrong, as a dated correction with the old value kept; state whether any
  verdict moves. Also α_s(M_Z) in `t5_sum_rules` (PDG 2026), the Suelzle
  1967 ⁷Li HO parameters and the Moniz 1971 ⁶Li k_F row (recorded only).

## Phase P — figures for the existing rows, and `--plots`

Worker P: `plot()` for every existing harness where a figure applies, the
`make_report.py --plots` wiring and its tests, the figures written once;
`phase_P_figures.md`.

## Phase E — documents, report, decisions

- **E1** `make_report.py --all --plots` then `--check`; `REPORT.md` /
  `REPORT.json` regenerated (the M2 baseline moves — that is the point;
  `--compare` lists every move).
- **E2** `BENCHMARK_PLAN.md` §3: a **dated status note only** (no change to
  the adopted milestone rule — the proposed clause goes to
  `DECISIONS_PRICED.md`); §4 rows 6, 7, 10 updated and the new rows
  numbered 12–; `validation/benchmarks/README.md`; `README.md`'s benchmark
  line; `REFERENCES.md` items 8, 9, 59 and the fetched 1a items (obtained,
  kept outside the repository); `PHYSICS_CHANNELS.md` only where a sentence
  became stale.
- **E3** `DECISIONS_PRICED.md`: the registry rows this run prices (the
  default unpolarized backend; `r_sigma_lt` vs `r1998`; D2 at M2; row 8's
  `eps_b0`; row 2 on HERMES) and the candidates it files (the M2 clause; the
  deuteron `TabulatedSpin1FF`; shipping the WS table; the refs/ filing and
  the licence question) — options and measured cost, **none taken**.
- **E4** `STATUS.md`, `SUMMARY.md` (with the items for the author), the
  `08_channel_priority.md` §7 log line.

## Phase F — verification

Four adversarial lenses over the uncommitted tree (physics; honesty — every
number traced to a window, every verdict to a declared tolerance; code —
tests, exit codes, sha256 guards, the gates, the figures; regression — the
18 existing harnesses and every gate bit-identical where no input changed),
then a fix stage, then the report regenerated once more and the run
committed and pushed.

## Not in this run, and why

- **The 27 registry decisions.** Priced where this run produces a price
  (E3); taken by nobody here.
- **The ePIC chain legs**: no container stack mounted; they validate the
  record, not the physics; M3.
- **The Mäntysaari-group letter and the dipole-model run** (row 21): outward.
- **The sibling repository**: read, never edited.
- **Default switches** suggested by the ranking: the price is this run's
  deliverable; the switch is the author's at M2.
- **Abbott's deuteron form factors** (C3): deferred with design, above.
