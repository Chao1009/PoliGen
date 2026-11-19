<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — phase P: a figure for each existing row that can draw one, and `make_report.py --plots`

Worker P, 2026-10-26 (measured 23:20–23:47 UTC on the run's 4-core machine,
matplotlib 3.10.9, HEAD `2d89a64` plus the uncommitted work of this run).
Plan: `PLAN.md`, the section "Figures — one per row" and phase P. **Nothing
in the tree moved**: no physics default, constant, registry label,
tolerance, reference-JSON number, `run()`, printed line or verdict changed
(measured, §5). A figure is a view of a row of `REPORT.json`, never the
record: every number a figure draws is in the row it illustrates, and
`make_report.py --check` does not look at figures.

## 1. What was done

1. **`make_report.py --plots`** (§4): after `run()` each child calls the
   harness's own `plot(report, out_dir=<report dir>/figures)` when the
   module defines one, puts the written path, relative to the report
   directory, into the row as `figure`, and REPORT.md links it as an image
   under the row's heading in *Detail*; the header says how many rows carry
   one. `--check` ignores every figure field the way it ignores the date and
   the runtimes. Without `--plots` the report is byte for byte what it was.
2. **`plot(report, out_dir=None)` and `--plot`** on 11 of the 18 harnesses
   that existed before this run, through `validation/benchmarks/_plots.py`
   (the one style: the reference as ink points with errors, each
   configuration in its fixed palette slot, a pull panel under a χ²
   comparison, the status in the title). `run()` untouched; `main()` (or the
   `__main__` block) runs `run()` exactly as before and, only when `--plot`
   is on the command line, then `plot()` and one extra line
   `figure: <path>`; the exit rule (0 ran, 2 broken) unchanged.
3. **One test per figure** in the harness's existing test file,
   `test_<harness>_plot_writes_png` (a non-empty PNG with the PNG signature
   in a tmp dir), and two `--plots` tests in `test_bench_report.py`.
4. **The 11 figures written once** into `docs/benchmarking/figures/` by
   running each harness with `--plot` (§7), each looked at and its label
   collisions fixed before that run.

## 2. The rows that got a figure (11), what each draws, and the channel it sits under

The status of every row is the harness's, unchanged; the numbers are the
ones the row printed when the figure was drawn (23:43 UTC), in the windows
the harness declares. Channel ranks and ingredients are
`../../benchmarking/08_channel_priority.md` §2 and §4.

| harness | the figure | numbers drawn, with their window | status | sits under (08) |
|---|---|---|---|---|
| `t1_nmc_f2d` | F₂ᵈ per nucleon vs x, NMC on the DIS cut (ink; the 3 points outside it hollow), ToyF2 / MSTW2008 LO / CT18NLO (LO formula) at each point's own (x, Q²) as ticks; pull panel (th − data)/σ | χ²/n on Q² ≥ 1, W² ≥ 4 GeV² (Q² 1.25–65 GeV² on the cut), stat ⊕ sys diagonal: ToyF2 8307.9/155 (p = 0), MSTW 416.5/155, CT18NLO 3401.8/146 (its grid floor 1.677 GeV²); tolerance p ≥ 0.01 | FAIL | ingredient A, F₂ (ranks 1 rate, 4) |
| `t1_nmc_f2d_over_f2p` | F₂ᵈ/F₂ᵖ vs x, the same layout | ToyF2 2161.4/211 (p = 0), MSTW 247.2/211 (p = 0.044), CT18NLO 215.5/189 (p = 0.090); Q² 1.11–99.03 GeV² on the cut | FAIL | ingredient A |
| `t1_g1d_world` | g₁ᵈ per nucleon vs x, COMPASS / HERMES / E143 (ink, one marker each) and SMC (recorded only, hollow), the four backends as ticks; pulls on the three verdict sets | ToyG1: COMPASS 107.8/15, HERMES 69.4/15, E143 54.9/28 (Q² 1.1–60.8 GeV² on the cut); the other backends' χ² are in the row (NNPDFpol1.1 45.0/15, 52.6/15, 51.0/28) | FAIL | ingredient H, g₁ (ranks 3, 5) |
| `t2_hermes_b1_table2` | x·b₁ᵈ per nucleon in HERMES's six bins (ink, stat ⊕ syst), the five model configurations as markers side by side, b₁ = 0 as the zero line; pull panel with the per-bin tolerance \|pull\| ≤ 1 shaded | χ²(b₁)/6 and bins within: CDKS conv. + MSTW 22.26 (2/6, the row's headline), CDKS conv. ToyF2 22.67 (2/6), Miller × 0.5 shipped 5.26 (4/6), Miller × 1 5.30 (5/6), CDKS digitized 21.78 (2/6), b₁ = 0 21.80 (2/6) | FAIL | ingredient D, registry row 2 (rank 4) |
| `t2_deuteron_static` | dot-and-band: η of the three shipped waves against η_d = 0.0256(4) with the 2σ tolerance shaded; Q_d (impulse approximation) beside it against 0.2860 ± 0.0015 fm² (phase D's correction, 2026-10-26; the figure regenerated) | η pulls CD-Bonn −0.07σ, AV18 −1.39σ, Hulthén (tagged default) +17.20σ; Q_d −5.42 / −5.71 / +1.52 % (recorded, not asserted) | FAIL | ingredient G (ranks 6, 7) |
| `t2_deeps_spectator_tail` | backward-window data / fitted model per p_s bin (0.30–0.53 GeV/c), colour = wave, line style = prescription (LC, instant form, LC × Glauber) | Hulthén LC χ²/ndf 627.85/187 (p = 4.9 × 10⁻⁴⁹), AV18 579.12/187, CD-Bonn 278.81/187, cos θ_pq ≤ −0.3, one free scale per (Q², W*, cos) group; Hulthén LC runs 0.79 / 0.95 / 1.27 / 1.71 / 1.41 over the five p_s bins | FAIL | ingredients G, I (tagged, ranks 6, 7) |
| `t3_li6_charge_ff_fb` | \|F_C0(q)\|/Z of the shipped `HoSpin1FF` and of the UVa Fourier–Bessel density on a log axis, the band [q_FB, q_Li71] shaded, the tree's zero marked | q₀ = 3.0998 fm⁻¹, +0.2713 fm⁻¹ above the band [2.6944, 2.8284] fm⁻¹ (T11's [2.9, 3.3] not moved) | FAIL | ingredient F, the C0 zero (rank 4; ≤ 2.2 × 10⁻⁷ on ΔA_zz per 08 §5) |
| `t3_li6_rms_devries` | dot-and-band: the three de Vries Table I radii and Angeli–Marinova (closure, hollow) against the tree's r_ch line, the PASS window shaded | r_ch = 2.5710 fm; de Vries pulls +0.62 / +0.22 / +0.01σ, window [2.46, 2.64] fm; Angeli −0.46σ (the tree's own input) | PASS | ingredient F (an input's consistency, not a prediction) |
| `t3_nmc_li6_over_d` | NMC F₂(⁶Li)/F₂(D) vs x (ink; the 9 points below the LHAPDF grids hollow), EPPS21 / isoscalar CT18NLO at each on-grid point's own (x, Q²), the valence window shaded; pull panel (data − tree)/σ | x ≥ 0.30: 4.627/4 (p = 0.328); on grid (Q² ≥ 1.69 GeV²) 26.94/15 (p = 0.029); tolerance p ≥ 0.01 on both | PASS | the EPPS21 baseline whose valence depletion (`EMC_VALENCE_DEPLETION_EPPS21`, `include/lipolgen/sf.hpp`) scales the polarized-EMC curves (rank 5) |
| `t4_pythia_ep_closure` | bridge / stock per observable in W1 and W2 with statistical errors, the stock reference the line at 1, D-1's ±20 % band | worst \|ratio − 1\| = 0.076 (W2 n_ch_eta[−2,0)) over 15 compared sub-rows; W2 n_ch_eta[−inf,−2) not compared (5 stock particles < 100), drawn hollow | PASS (opt-in) | ingredient J, T2 hadronization |
| `t5_sum_rules` | the Bjorken running integral ∫_{x_min}^1 (g₁ᵖ − g₁ⁿ) dx divided by g_A/6·C_Bj at each Q² (the reference is the line at 1, ±10 % shaded), ToyG1 and NNPDFpol1.1 at Q² = 2, 5, 10 GeV²; the BC / g₂^WW closure in the title line | BC closure PASS 45/45; ToyG1 diverges (Q² = 5: 0.6430 at x_min = 10⁻⁶ against 0.1854 ± 0.0023); NNPDFpol1.1 +19.3 / +12.4 / +9.6 % at Q² = 2 / 5 / 10, its last decade still +2.2 % of the reference (not converged by the declared 1 % rule — what `phase_A_baseline.md` §3 already corrected in `phase_0_environment.md`) | FAIL | ingredient H |

What the figures show at a glance that the one-line rows do not (read off
the PNGs, i.e. off each row's own theory column and points; no verdict and
no number of a row changes by it): ToyF2 lies above every NMC F₂ᵈ point
below x = 0.03 and misses in both directions above it, while CT18NLO's
excess is not confined to small x but runs across the whole range; in g₁ᵈ
the toy's large negative values (down to −0.96) at COMPASS's lowest x drive
COMPASS's χ², while HERMES and E143, whose points start near x = 0.03, fail
above it; the CDKS configurations sit below HERMES's positive b₁ᵈ in the
three lowest-x bins, which only the two Miller curves follow; the Hulthén η
is the only η outside the 2σ band; the default wave's Deeps shape miss is
largest at p_s = 0.46 GeV/c; the tree's C0 zero lies where its own |F_C0|/Z
has already fallen to ~10⁻⁵, more than two orders below the FB curve there.

## 3. The rows that did not get a figure (7), and why

| harness | why no `plot()` |
|---|---|
| `t5_epios_source_modes` | an identity: max moment round-trip residual 0 over 8 modes against 1e-12 — a number at zero has no comparison to draw (plan, "Figures": "Does not apply to: EPIOS/EST") |
| `t5_est_identity` | an identity at the 1e-14 level (1.138e-14 against 2.31e-14) — same reason |
| `t5_weighted_unweighted` | one pull (+0.359 against \|pull\| ≤ 3) — the plan's "weighted/unweighted (one pull)" |
| `t6_hepmc3_convention` | a clause table (23/23) — the plan's "T6 (a clause table)" |
| `t2_bonus_spectator_shape` | owned by worker B1 this run (being re-wired): **not touched**; it already defines `plot()` in the working tree (its figure is B1's) |
| `t3_li_magnetization_rfy` | owned by worker B2 this run: **not touched**; it already defines `plot()` in the working tree (B2's) |
| `t4_djangoh_rad_noRad` | owned by worker B3: **not touched**; `plot()` committed in `372b915` |

So, counted over the 21 harnesses on disk at 23:46 UTC (the 18 plus this run's
`t1_hermes_f2d`, `t1_r_world`, `t5_delta_moment_bound`): **17 define
`plot()`, 4 do not by design**; a full `make_report.py --all --plots` should
therefore say "17 of 21 rows carry one" if every row runs. Measured through
the generator (§5): 11 of 15 on the rows of this phase, 6 of 6 on the six
other-owned rows, into scratch directories.

The five files the task listed as owned by other workers this run but open
to a `plot()` once their working-tree copy was unmodified —
`t1_g1d_world.py`, `t1_nmc_f2d.py`, `t1_nmc_f2d_over_f2p.py`,
`t2_deuteron_static.py`, `t5_sum_rules.py` — were all unmodified (`git status
--porcelain` empty) at the moment of insertion, checked twice by the
insertion script (before reading and again before writing); each got its
`plot()` directly. Nothing is left for the integrator to add for them.

## 4. The `make_report.py` change (`--plots` code paths only)

- New option `--plots`; new hidden child option `--child-report-dir`; new
  constants `FIGURES_DIR = "figures"`, `FIGURE_KEYS = ("figure",
  "figure_error")`, `PLOTS_FLAG = " --plots"`.
- `draw_figure(mod, raw, row, report_dir, buf)`: calls `plot(raw,
  out_dir=<report_dir>/figures)` with stdout captured; a module without
  `plot()` gets nothing and no complaint; `plot()` returning None means
  nothing to draw (a BLOCKED row); a `plot()` that raises, or returns a path
  that is not a non-empty file, costs the figure only — `figure_error` with
  the reason (`_report.classify_exception`'s text), the row's status and
  numbers untouched, and `generate()` echoes `FIGURE   <name>: not drawn --
  <reason>` to stderr. Called only for a row whose status is one a harness
  returns (pass / fail / blocked), never for an error, environment or
  opt-in row.
- `run_harness`, `run_harness_isolated`, `child_main`, `generate` carry the
  report directory through (default None = no figures). `generate` records
  `generated.figures_dir = "figures"` and the ` --plots` token in
  `generated.command` when it draws.
- `render_md`: when `figures_dir` is in the header, a line
  "**figures:** N of M rows carry one, linked under the row's heading in
  *Detail* (… a view of REPORT.json, never compared by `--check`)"; under
  each `### <name>` with a figure, `![<name>](figures/<name>.png)` — a
  relative link; a `figure_error` as a "**figure:** not drawn — …" bullet.
- `--check`: `overlay_volatile` takes every figure field from the committed
  report exactly as it takes the date and the runtimes — `figures_dir`, each
  row's `figure` / `figure_error` (present or absent as committed), and the
  ` --plots` token of the command line — so a report made with `--plots`
  checks clean against a regeneration made without; a changed or missing PNG
  is no difference; `--check` never draws and refuses `--plots` (exit 2);
  `--compare --plots` needs `--out-dir` (exit 2 otherwise: the figures go
  beside a written report). Said in the module docstring (FIGURES) and in
  the `--check` option text.
- Without `--plots`: no new key anywhere and no new line; measured
  byte-identical to HEAD's generator (§5).

## 5. Verification (commands in §7)

- **No printed line and no row moved**: each of the 11 harnesses run twice
  in one process, HEAD's source and the working tree's (a scratch script,
  not committed: both sources loaded as modules with the harness's own
  `__file__`, `run(verbose=True)` each, stdout compared with every line
  matching `\d+\.\d+ s\b` masked, and `_report.jsonable(_report.
  normalise_row(run(), NAME))` compared): printed lines identical (102 / 168 / 270 / 38 / 12 / 28 /
  15 / 13 / 29 / 80 / 24 lines for `t1_g1d_world`, `t1_nmc_f2d`,
  `t1_nmc_f2d_over_f2p`, `t2_deeps_spectator_tail`, `t2_deuteron_static`,
  `t2_hermes_b1_table2`, `t3_li6_charge_ff_fb`, `t3_li6_rms_devries`,
  `t3_nmc_li6_over_d`, `t5_sum_rules`, `t4_pythia_ep_closure`) and the
  normalised REPORT.json row
  (`_report.normalise_row` → `jsonable`) identical, for all 11.
- **The generator without `--plots` is byte-identical to HEAD's**: HEAD's
  `make_report.py` and the new one, `generate(opt_in=True)` on the test
  suite's 12 fake harnesses (the 9 of `FAKES` and the 3 `PLOT_FAKES`),
  runtimes zeroed: REPORT.json and REPORT.md identical.
- **pytest** (only the files touched): `test_bench_t1_nucleon.py`,
  `test_bench_static.py`, `test_bench_t5_closures.py`,
  `test_bench_deeps.py`, `test_bench_report.py`, `test_bench_t4_t6.py`:
  **138 passed, 3 skipped, 0 failed** (29.3 s);
  `test_bench_spin_deuteron.py` + `test_bench_nuclear.py`: **36 passed, 1
  failed** — the failure is `test_bench_spin_deuteron.py::
  test_bonus_row_is_blocked_by_name`, which asserts row 10 BLOCKED; worker
  B1's re-wiring makes it a measurement (FAIL), and the plan leaves the old
  `test_bonus_*` skips to the integrator — not caused by this phase. The 11
  new plot tests alone: **11 passed** (15.2 s).
- **The report rows render**: `make_report.py --all --plots --only <the 15
  rows of this phase> --out-dir <scratch>` → 7 pass, 8 fail, 0 error; 11
  figures, header "11 of 15 rows carry one"; the four rows without `plot()`
  have no image line. The six rows owned by other workers through the same
  generator (scratch): 6 of 6 drawn, no `figure_error`.
- `validation/check_spdx_headers.py`: 131 files, all carry the line.

## 6. What this does not do

It validates nothing: a figure redraws a row. Three figures depart from a
plain redraw of the row's numbers, each said in its `plot()` docstring:
`t3_li6_charge_ff_fb` evaluates the two form-factor curves on a q grid for
drawing (the FB curve from the report's own vendored coefficients, the tree
curve from the shipped `HoSpin1FF`; the zeros, the band and the distance
drawn are the report's); `t2_hermes_b1_table2` draws x·b₁ rather than b₁ so
that all six bins are readable (x is the bin's, so each pull is b₁'s);
`t2_deeps_spectator_tail` draws the weighted means `data_over_model_by_ps`
WITHOUT error bars, because the row does not carry their errors and `run()`
was not to change. The T1 figures plot every point at its own Q², so a
column of points at one x is several Q² (colour is for configurations, as
the task asked). PNG bytes are not reproducible across matplotlib builds
and are not in `--check`; the footer of each PNG carries the date and the
revision it was drawn at (`2d89a64+`, dirty).

## 7. Files touched, and the exact commands

Code: `validation/benchmarks/make_report.py` (the `--plots` paths and the
docstring); `plot()` + `--plot` in `validation/benchmarks/t1_g1d_world.py`,
`t1_nmc_f2d.py`, `t1_nmc_f2d_over_f2p.py`, `t2_deeps_spectator_tail.py`,
`t2_deuteron_static.py`, `t2_hermes_b1_table2.py`, `t3_li6_charge_ff_fb.py`,
`t3_li6_rms_devries.py`, `t3_nmc_li6_over_d.py`, `t4_pythia_ep_closure.py`,
`t5_sum_rules.py`. Tests (appended, one per harness):
`python/tests/test_bench_t1_nucleon.py` (3), `test_bench_static.py` (2),
`test_bench_t5_closures.py` (1), `test_bench_spin_deuteron.py` (1),
`test_bench_deeps.py` (1), `test_bench_nuclear.py` (2),
`test_bench_t4_t6.py` (1), and `test_bench_report.py` (2 `--plots` tests
with their fakes). Figures: `docs/benchmarking/figures/t1_g1d_world.png`,
`t1_nmc_f2d.png`, `t1_nmc_f2d_over_f2p.png`, `t2_deeps_spectator_tail.png`,
`t2_deuteron_static.png`, `t2_hermes_b1_table2.png`,
`t3_li6_charge_ff_fb.png`, `t3_li6_rms_devries.png`,
`t3_nmc_li6_over_d.png`, `t4_pythia_ep_closure.png`, `t5_sum_rules.png`.
This record.

    cd /home/cpeng/Projects/polli/PoliGen && source env.sh
    for h in t1_g1d_world t1_nmc_f2d t1_nmc_f2d_over_f2p t2_deeps_spectator_tail \
             t2_deuteron_static t2_hermes_b1_table2 t3_li6_charge_ff_fb \
             t3_li6_rms_devries t3_nmc_li6_over_d t4_pythia_ep_closure t5_sum_rules; do
        python3 validation/benchmarks/$h.py --plot      # exit 0 each; prints "figure: <path>"
    done
    python3 -m pytest -p no:cacheprovider -q python/tests/test_bench_t1_nucleon.py \
        python/tests/test_bench_static.py python/tests/test_bench_t5_closures.py \
        python/tests/test_bench_deeps.py python/tests/test_bench_report.py \
        python/tests/test_bench_t4_t6.py                 # 138 passed, 3 skipped
    python3 -m pytest -p no:cacheprovider -q python/tests/test_bench_spin_deuteron.py \
        python/tests/test_bench_nuclear.py               # 36 passed, 1 failed (B1's row 10, above)
    python3 validation/benchmarks/make_report.py --all --plots --only \
        t1_g1d_world,t1_nmc_f2d,t1_nmc_f2d_over_f2p,t2_deeps_spectator_tail,t2_deuteron_static,t2_hermes_b1_table2,t3_li6_charge_ff_fb,t3_li6_rms_devries,t3_nmc_li6_over_d,t4_pythia_ep_closure,t5_sum_rules,t5_epios_source_modes,t5_est_identity,t5_weighted_unweighted,t6_hepmc3_convention \
        --out-dir <scratch>/report_P                     # 11 of 15 rows carry one
    python3 validation/check_spdx_headers.py

For the integration stage (E1): `make_report.py --all --plots`, then
`make_report.py --check` — the check passes on the report the first command
writes (§4), and the regenerated figures replace the ones written here.
