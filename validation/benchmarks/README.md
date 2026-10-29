<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# validation/benchmarks — the automated benchmark harnesses

`docs/benchmarking/BENCHMARK_PLAN.md` §3: *a benchmark that is not automated is
not a benchmark.* Every wired row of the plan is one file here. This page states
the conventions once and lists every harness; the evidence for each number is in
the run record named in the last column, and the numbers themselves are
regenerated into `docs/benchmarking/REPORT.md`. (Merged 2026-09-23 from this
file and the nuclear group's `README_nuclear.md`, which it replaces; the nine
harnesses of 2026-09-26 and the report generator added on that run; the
five of 2026-10-26 (`BENCHMARK_PLAN.md` §4 rows 12–16) and one figure per row
(`make_report.py --plots`, `_plots.py`) added on that run — **23 harnesses** in
all, 19 of them drawing a figure.)

## Conventions

A harness is `<tier>_<name>.py` (`t1_` nucleon, `t2_` deuteron, `t3_` nuclear
inputs, `t4_` generator-vs-generator, `t5_` MC closures, `t6_` chain), and each:

1. carries the three normalisation rules of BENCHMARK_PLAN.md §8 in its header,
   stated **for that row** on both sides (per nucleon / per nucleus / per
   deuteron; the A-scaling assumption; the isoscalar denominator), or says
   explicitly why a rule does not apply to it;
2. reads its reference from a **vendored** table under `data/`, whose file
   carries its provenance (URL, date fetched, licence, the paper, the
   publishing convention). **Since 2026-09-26 every harness that reads a
   vendored file records a sha256 of it and re-checks it on every read**; an
   altered file raises `RuntimeError` and the harness exits 2. The digest
   covers every byte of the file in every harness **except row 4**
   (`t3_nmc_li6_over_d`), which hashes only the body after its provenance
   marker: its header (source, licence, the SPDX line) is pinned only by
   `REPORT.json`'s `data_sha256` through `make_report.py --check`. Row 5
   (`t3_li6_charge_ff_fb`) hashed only its ⁶Li row until 2026-09-27 and now
   hashes every byte too (the row digest kept). The nnpdf-data and
   static-observable files also carry the sha256 of their canonical body
   (`_body_sha256`), checked on its own. Three references are
   not files under `data/`: `t5_sum_rules` embeds its PDG constants in the
   harness under a sha256 it re-checks, `t6_hepmc3_convention`'s reference is
   `docs/HEPMC3_CONVENTION.md` (its test pins the document's sha256 and each
   clause's line), and `t4_pythia_ep_closure`'s is a stock PYTHIA run made in
   the same process. No harness and no test needs the network;
3. runs the generator or the library function at the reference's kinematics;
4. prints one `REPORT | name | generator value | reference | tolerance | STATUS`
   line (`PASS` / `FAIL` / `BLOCKED`) and returns the row from `run()` — a
   dict with at least `name`, `generator`, `reference`, `tolerance`, `status`.
   Sub-rows print either as `REPORT | name[config] | …` (`t5_sum_rules`,
   which then prints the row's own line; row 3 prints only its six
   configuration lines) or as `SUBROW | …` (T1, `t2_deuteron_static`); the
   DJANGOH harness (row 7) prints `SUBROW | …` lines and its `REPORT |` line since 2026-10-26 (a `ROW |` line before). **`docs/benchmarking/REPORT.md` and
   `REPORT.json` exist since 2026-09-26** and are generated from the dicts
   `run()` returns, never from the printed lines and never by hand:

       cd LiPolGen && source env.sh
       python3 validation/benchmarks/make_report.py --all      # regenerate both
       python3 validation/benchmarks/make_report.py --all --plots   # ... and one figure per row (docs/benchmarking/figures/)
       python3 validation/benchmarks/make_report.py --check    # still true?
       python3 validation/benchmarks/make_report.py --all --compare docs/benchmarking/REPORT.json

   `make_report.py` runs every `t[1-6]_*.py` in its own child process (a crash
   costs that row only), records a missing LHAPDF set as `blocked
   (environment)` and any other exception as `error`, and writes the
   environment (PYTHIA, LHAPDF and the sets present, numpy, HepMC3) into the
   header. `--check` exits 1 unless the committed files are what a
   regeneration gives (date, revision, runtimes and environment aside);
   `--compare OLD.json` prints every status change or headline-number move and
   exits 1 on any (the M2 rule; a move is relative, with no absolute floor
   since 2026-09-27, so a 1e-14 headline that doubles counts). `REPORT.json` is the M2 baseline;
5. is exercised by a pytest under `python/tests/test_bench_*.py` (rows 1 and 9
   also by `tests/test_bench_spin.cpp`), with a guard on each vendored file
   (point count, first and last value, one flipped byte refused);
6. **exits 0** whenever it ran and printed its row — `PASS`, recorded `FAIL` and
   `BLOCKED` alike, since the row carries the verdict — and **exits 2** only when
   the harness itself is broken (a missing or altered vendored file, a missing
   dependency, any exception). A harness is a measurement, not a CI gate: the
   pytests are the gate. (Uniform since 2026-09-23's fix stage.) A
   **configuration** that needs a missing LHAPDF set is a `BLOCKED` sub-row
   with reason `environment: LHAPDF set <name> not installed` (T1,
   `t5_sum_rules`), never dropped and never faked; a harness whose **only**
   configuration needs one (row 4; row 13 since phase F of run 2026-10-26) raises, exits 2, and is `blocked
   (environment)` in the report. T4 and T6 report a build tier that is off
   (`HAVE_PYTHIA8`, `HAVE_HEPMC3`) or a missing `pythia8` / `pyhepmc` module as
   a `BLOCKED` row with an `environment:` reason; a missing or broken
   `lipolgen` itself is the harness broken, exit 2, as everywhere else (until
   2026-09-27 T4 and T6 turned it into a `BLOCKED` row with exit 0);
7. **draws its comparison, where one can be drawn** (since 2026-10-26):
   `plot(report, out_dir=None)` through `_plots.py` — the one style: the
   reference as ink points with errors, each configuration in its fixed
   palette slot, a pull or ratio panel under a χ² comparison, the status in
   the title — and `--plot` on its own command line (exit rule unchanged);
   `make_report.py --plots` calls it after `run()` and links
   `docs/benchmarking/figures/<harness>.png` under the row's heading in
   `REPORT.md`'s Detail. A figure is a view of `REPORT.json`, never the
   record: `--check` takes every figure field from the committed report as it
   takes the date and the runtimes (PNG bytes depend on the matplotlib
   build). No figure, by design: `t5_epios_source_modes` and
   `t5_est_identity` (identities), `t5_weighted_unweighted` (one pull),
   `t6_hepmc3_convention` (a clause table). CI (`.github/workflows/ci.yml`)
   installs no matplotlib (nor scipy), so every `plot()` test skips there:
   the figures are exercised only where matplotlib is installed (the
   integration runs); `test_bench_report.py` checks, without matplotlib,
   that every figure the committed report links exists (since phase F).

**A benchmark that disagrees with the tree is RECORDED, never tuned.** No
harness moves a default, a registry row or a reference JSON; a `FAIL` is a
measurement with its distance, priced for the registry row it bears on. Every
tolerance is stated in the harness header independently of the verdict, and
where the number was seen before the tolerance was chosen (row 4, the three
T1 rows, CD-Bonn's η in `t2_deuteron_static`, `t3_li6_rms_devries`, the BC
preview in `t5_sum_rules`; since 2026-10-26 also row 10's BONuS tables, read
by the supervising session before the design), the header says so; where a
number was measured first under a tolerance already in place (since
2026-10-26: the T1 CT18NLO / NNPDFpol11 sub-rows, before phase D re-read the
D-state convention, and T4's 0.076 on the pinned PYTHIA 8.317 under the
2026-09-26 tolerance — both in phase 0), the header records it as measured
(T4's since phase F);
rows 6, 7, 15 and 16 (2026-10-26) also hash their tolerance in a design block
written before the first run, and refuse (exit 2) an edited block.

**Blocked rows.** A row whose input is neither on disk nor open access is
`BLOCKED`: its harness exists, and its test skips with a message naming the
exact document to download (citation and link as in
`docs/references/REFERENCES.md` §1; rows 6 and 10 until 2026-10-26, when the
run of that date obtained both documents). A row can also be blocked because the
comparison it was planned as does not exist (row 7, until 2026-10-26), or because the tree as
built cannot define it (`t5_sum_rules[fsi_unitarity]`).

**The opt-in variable.** An expensive row — one that costs more than 30 s, or
one the plan counts as expensive by kind (a full generation, a chain leg, a
generator-vs-generator sample) whatever it costs — runs only with

    LIPOLGEN_BENCH=1 python -m pytest python/tests -q

and is otherwise reported as skipped with that reason. Since 2026-09-26 two
things read it: `python/tests/conftest.py` registers the pytest marker
`bench_expensive` and skips every test carrying it, with a reason naming the
variable, unless `LIPOLGEN_BENCH=1`; and `make_report.py` runs an opt-in row
(a module with `EXPENSIVE = True`, or one in its own `EXPENSIVE` set) only
with the variable or `--all`, else writes it as `skipped (opt-in)`. Both skip
reasons say "> 30 s, or a generator-vs-generator sample" (until 2026-09-27
they said "costs more than 30 s", untrue of T4). Its users
today: `t4_pythia_ep_closure` — its full sample (`test_bench_t4_t6.py::test_t4_full_row_passes_and_is_pinned`)
takes 9.2–10.2 s, under the 30 s line, and is opt-in because the plan counts a
generator-vs-generator sample as expensive, not by its cost — and
`test_bench_report.py`'s `--check` of the committed report (≈ 45 s). Row 7
(22–26 s, measured 2026-10-26 on a shared 4-core machine; its `--diagnostics` table-vs-quadrature check, ~110 s, is opt-in as `test_bench_djangoh_eltail.py::test_diagnostics_table_vs_quadrature`) is not opt-in, so every other row runs in the ordinary
`python -m pytest python/tests -q`.

**Environment.**

    cd LiPolGen && source env.sh     # build/python on PYTHONPATH, LHAPDF sets
    python3 validation/benchmarks/t5_epios_source_modes.py
    python -m pytest python/tests -q -k bench -rs

Row 4 needs LHAPDF with `EPPS21nlo_CT18Anlo_Li6` and `CT18NLO` installed
(`env.sh` points `LHAPDF_DATA_PATH` at them); its tests skip, naming the sets,
if they are not. The T1 CT18NLO / NNPDFpol11 sub-rows and
`t5_sum_rules[bjorken[nnpdfpol11]]` need `CT18NLO` and `NNPDFpol11_100`. T4
and T6 need the PYTHIA 8 and HepMC3 tiers and the `pythia8` and `pyhepmc`
Python modules; T4 records the PYTHIA version it ran (8.312 on 2026-09-26; the
project pins 8.317, run since 2026-10-26 — its Python module imports only from a
working directory with no `python/` entry, so T4, T6 and their tests import it
through `_report.import_pythia8()`,
`docs/open_items/run_2026-10-26/phase_0_environment.md`).

**Tools here that are not harnesses** (their names do not match `t[1-6]_*.py`):
`make_report.py` and `_report.py` (the report, above); `_plots.py` (the one
figure style every `plot()` draws with, convention 7); `tiers.py`
(BENCHMARK_PLAN.md §7 as data: `--since REV` lists the harnesses a diff must
re-run, `make_report.py --since REV` runs them); `vendor_nnpdf_commondata.py`
(dev-time: rebuilds the six `nnpdf_*.json` files from the nnpdf-data 4.1.5
wheel, sha256 `6ed209ad…3490269a`, checking every file it reads against the
wheel's RECORD; `--check` regenerates in memory and exits 1 on any byte
difference). Beside them, one directory up: `../price_defaults.py` (run
2026-10-26, phase A3) — the default-input price table of the ranked
observables, a recipe and not a harness (no reference, no tolerance, no
verdict): each channel's own observable under every knob the knob-provenance
table reports READ, as the ratio to the shipped default at config 1,
N = 200 000, seed 20261026 (~190 s); `--plot` draws
`docs/benchmarking/figures/price_defaults.png` through `_plots.py`; tested by
`python/tests/test_price_defaults.py`.

## The harnesses (status and runtime measured by run 2026-09-26; rows run again, re-wired or added by run 2026-10-26 say so)

Runtime is the wall time of `python3 validation/benchmarks/<file>`,
interpreter start-up included, as the range of three runs measured on
2026-09-27 (UTC) by the integration stage of run 2026-09-26 with
`time.perf_counter` around a subprocess (`/usr/bin/time` is not installed on
that machine: PYTHIA 8.312, no LHAPDF set). Row 4 cannot run there; its
status and runtime are 2026-09-23's, from the development machine. Rows run
again, re-wired or added by run 2026-10-26 (a shared 4-core machine with the
three LHAPDF sets and PYTHIA 8.317; `docs/open_items/run_2026-10-26/`) carry
that run's numbers and runtimes and say so. The status column is the row's own
verdict; `REPORT.md` has every sub-row.

### The plan's first ten (BENCHMARK_PLAN.md §4)

| row | harness | reference, vendored under `data/` | status | the number, with its tolerance | runtime | tests | record |
|---|---|---|---|---|---|---|---|
| 1 | `t5_epios_source_modes.py` | `epios2026_table2_deuteron_source_modes.json` — EPIOS, PRC 113 (2026) 060501, Table II (COSY/ANKE deuteron source; the P_zz column is ideal, only P_z is measured) | **PASS** | max moment round-trip residual **0** over 8 modes, tolerance 1e-12. Recorded, not asserted: the tree's EST (`ladder`) fill matches the ideal P_zz at mode 0 only, **1 of the 6 modes it computes**; at modes 5, 6 (\|P_z\| = 1) it throws and only the closed-form limit matches ("3 of 8" until 2026-09-26 counted those two) | 0.11–0.12 s | `test_bench_spin_deuteron.py::test_epios_*` (5); doctest *"bench row 1: …"* | `phase_B1_spin_deuteron.md` §1 |
| 3 | `t2_hermes_b1_table2.py` | `hermes2005_b1d_table2.json` — HERMES, PRL 95 (2005) 242001, Table II (per nucleon); `whitlow1990_r1990.json` — R1990 (the scan's b₁ = 0.635 read as 0.0635) | **FAIL, recorded** | every bin inside stat ⊕ syst: met by no configuration. χ²(b₁ᵈ)/6 = 22.26 (convolution + MSTW, the A = 2 gate's), 5.26 (shipped Miller ×0.5), 5.30 (Miller ×1), 21.80 (b₁ = 0) — HERMES does not decide registry row 2. χ²(b₁) does not depend on the F₂ backend (pinned without MSTW too) | 4.5–4.9 s | `test_bench_spin_deuteron.py::test_hermes_*` (7), `::test_r1990_*` (2) | `phase_B1_spin_deuteron.md` §3 |
| 4 | `t3_nmc_li6_over_d.py` | `hepdata_ins394050_table1.csv` — NMC F₂(⁶Li)/F₂(D), HEPData ins394050 Table 1, CC0 | **PASS** (2026-09-23; again 2026-10-26 with the sets installed; blocked (environment) where they are missing) | p(χ², ndf) ≥ 0.01 on both sets: **4.6269/4** on x ≥ 0.30 (p = 0.328), **26.941/15** on the 15 on-grid points (p = 0.029); 9 points below the EPPS21 grid reported, not compared | 1.3–1.4 s (2026-09-23) | `test_bench_nuclear.py::test_nmc_*` (4; the table guard runs without the sets) | `phase_B2_nuclear.md` §2 |
| 5 | `t3_li6_charge_ff_fb.py` | `uva_ncd_fb_data_li6.dat` — the ⁶Li row of the UVa NCD archive's `FB_data.dat` (no licence stated; provenance of the row UNVERIFIED) | **FAIL, recorded** | the C0 zero must lie in [2.6944 (FB zero), 2.8284 (Li *et al.* 1971 minimum)] fm⁻¹: `HoSpin1FF`'s is **3.0998**, **+0.2713 fm⁻¹** above; T11's [2.9, 3.3] not moved (open item Q1). A C0 with no zero below 5 fm⁻¹ is a recorded FAIL, not a crash | 0.11–0.14 s | `test_bench_nuclear.py::test_fb_*` (4) | `phase_B2_nuclear.md` §3 |
| 6 | `t3_li_magnetization_rfy.py` | `rfy1966_li_magnetic_ff.csv` — Rand–Frosch–Yearian, Phys. Rev. 144 (1966) 859, Tables I (⁶Li, 7 points) and II (⁷Li, 24 points), F_M² at 180°, read from the page images; `rfy1966_li_magnetization.csv` — Tables VII, VIII, XIV, a₀ and the coupling parameters verbatim (the PDF, sha256 `7cbb2d21…d5bbe2`, obtained 2026-10-26 through the institution's APS licence and kept outside the repository) | **PASS** (BLOCKED until 2026-10-26) | ⁶Li F_M² against `HoSpin1FF`'s unfitted F_mag² (q_z = 1.30 fm⁻¹, b = 1.85 fm), one free scale, p(χ², 6) ≥ 0.01 declared in a hashed design block: **13.35/6, s = 1.419 ± 0.206, p = 0.0378**; s = 1 recorded 17.49/7 (p = 0.0145); ⁷Li not compared (`HoSpin1FF.for_ion(li7())` throws); the paper's own fits re-evaluated as a transcription closure (⁶Li 8/8; ⁷Li 4/5 — the fifth, Table VIII (C) a₀ = 1.66, misses by 0.31, and a post-hoc refit holds 5/5); the 2026-09-23 radius rule retired (RFY print no radius; ADNDT 14's lithium 'magnetization radii' are √(5/2)·a₀). Cross-check (row 15): Wiringa–Schiavilla draw these seven points at 0.7116 × RFY's normalisation, the free scale's factor (1/0.7116 = 1.405) | 0.4–0.5 s | `test_bench_rfy.py` (16) | `run_2026-10-26/phase_B2_rfy.md` |
| 7 | `t4_djangoh_rad_noRad.py` | `djangoh_eltail_ep18x275.json` — private DJANGOH 4.6.21 runs (`spiesber/DJANGOH` at `e356b84`; 40 runs, every card verbatim, each log's sha256) of the isolated elastic tail and its Born, e 18 × p 275, four Q² bins; `djangoh_rad_norad_ep18x275.json` — the ePIC table, recorded history | **PASS** | R = σ_tail(LiPolGen)/σ_tail(DJANGOH): `polrad-full` 1.0199 / 1.0261 / 1.0219 / 1.0098 within \|R − 1\| ≤ 0.03 + 2δ_b (`polrad-full` is the shipped default since 2026-10-27, D5); `t-peak` (shipped until then, an option since) 0.6118 / 0.8968 / 1.0112 / 1.0067 within its one-sided ceiling R ≤ 1.03 + 2δ_b only (39 % / 10 % low at Q² < 100 GeV²); δ_b 3.7–4.8 × 10⁻⁴; control 8/8 ePIC numbers reproduced; the tail table's +1.9 … +2.5 % recorded | 22–26 s | `test_bench_djangoh_eltail.py` (19 since the B3 verification, one opt-in), `test_bench_chain_rc.py` (4) | `run_2026-10-26/phase_B3_djangoh.md` |
| 9 | `t5_est_identity.py` | `est_identity_sources.json` — the EST closed form; COMPASS ⁶LiD (Koivuniemi *et al.*, SPIN 2004) | **PASS** | max \|`spin_temperature_pzz(1, P_z)` − (2 − √(4 − 3P_z²))\| = **1.138e-14** over 19 999 interior points, tolerance **2.31e-14** (bisection bound 2.13e-14; the planned 1e-14 does not hold); COMPASS numbers at the printed digit | 0.14–0.16 s | `test_bench_spin_deuteron.py::test_est_*` (6); doctest *"bench row 9: …"* | `phase_B1_spin_deuteron.md` §2 |
| 10 | `t2_bonus_spectator_shape.py` | `bonus2014_tkachenko_supplemental.json` — CLAS BONuS R_D/S (data over the collaboration's own PWIA spectator simulation, no FSI) vs cos θ_pq, Tkachenko *et al.*, PRC 89 (2014) 045206 Supplemental Material (APS; fetched 2026-10-26 through the institution's licence; numbers only, the original's sha256 recorded inside): 100 blocks (4 GeV 3 Q² × 5 W* × 4 p_s, 5 GeV 2 × 5 × 4; 22 missing rows) and the DIS-region F₂ⁿ/F₂ᵈ, F₂ⁿ/F₂ᵖ table; `clas_db_deuteron_query_2026-09-23.json`, `clas_deeps_klimenko2006_f2P.json` (the 2026-09-23 history, still read for its census) | **FAIL, recorded** (BLOCKED until 2026-10-26) | DIS window W* ≥ 2, Q² ≥ 1 GeV² declared first (28 blocks with backward bins: 12 at 4 GeV, 16 at 5 GeV; 68 out of model). Tolerance: `GlauberFsiWeight` (Hulthén, 40 mb) over its backward mean inside the data's 1σ flat band (stat, PDG scale factor) at every cos θ_pq ≤ −0.3 bin of every block: **4/12, 0/16**; worst v/δ 2.62 / 3.44. The data fall toward 90° as the weight does — χ²_flat 739.4/84 against 226.3/84 for the weight's shape (20 mb 403.3). Recorded: AV18 3/12, 0/16; 20 mb 8/12, 13/16; stat ⊕ syst 8/12, 8/16; the forward rise; the instant-form/light-cone factor (same direction, 40–66 %); ToyF2 F₂ⁿ/F₂ᵖ 40.5/10 (8.30/9 with one scale) | 2.3–2.4 s | `test_bench_bonus.py` (26); `test_bench_spin_deuteron.py::test_bonus_*` (3 history guards; `test_bonus_row_is_blocked_by_name` retired 2026-10-26) | `run_2026-10-26/phase_B1_bonus.md` |

Rows 2 and 8 of the plan have no harness: row 2 was resolved in code on
2026-09-06 (`tests/test_tagged.cpp`'s Cosyn–Weiss gate), and row 8 is a
relabelling (George–Knutson's η is a consistency band on the quadrupole dial,
not a benchmark). Records are under `docs/open_items/run_2026-09-23/`.

### By tier, added 2026-09-26 and 2026-10-26 (BENCHMARK_PLAN.md §4, "Then, by tier", and rows 11–16)

The "row" is the survey row the harness wires (`02_data_nucleon_deuteron.md`,
`03_data_nuclear.md`, `05_chain.md`); a plan-row number (11–16) where the
harness is a §4 row. None was tuned; the record is
`docs/open_items/run_2026-09-26/SUMMARY.md` (the 2026-10-26 rows:
`docs/open_items/run_2026-10-26/`, one phase record each) and the harness
header.

| row | harness | reference, vendored under `data/` | status | the number, with its tolerance | runtime | tests |
|---|---|---|---|---|---|---|
| T1 F-1 | `t1_nmc_f2d.py` | `nnpdf_NMC_NC_NOTFIXED_D.json` — NMC F₂ᵈ per nucleon, NPB 483 (1997) 3, HEPData ins424154 tables 17–32 (CC0), via the PyPI nnpdf-data 4.1.5 wheel (GPL-3.0-or-later); 158 points; read 2026-10-26: 158/158 = the paper's Table 6; NMC's 2.5 % normalisation is outside 'sys' and not in the file; F₂ᵈ extracted with NMC's own R below x = 0.12 and R1990 above (`phase_D_inputs.md` §1.5) | **FAIL, recorded** | p(χ², ndf) ≥ 0.01 on the DIS cut Q² ≥ 1, W² ≥ 4 GeV² (155 points), stat ⊕ sys diagonal: shipped ToyF2 **8307.9/155** (p = 0). Sub-rows: MSTW 416.5/155 (p = 8.2e-26); CT18NLO at LO 3401.8/146 on Q² ≥ 1.677 (FAIL, mean model/data 1.086; measured 2026-10-26, the set installed — blocked (environment) before). Full covariance recorded beside it (21100.7; MSTW 858.5) | 0.78–0.85 s | `test_bench_t1_nucleon.py` (42 for the three T1 rows until 2026-10-26, 50 since — phase D's pins of the LHAPDF sub-rows and the papers read, phase P's figures: `test_f2d_*` 7, `test_dp_*` 6, `test_g1d_*` 11, the rest shared — file guards, blocked LHAPDF configs by name, exit codes, the no-scipy χ² fallback; one needs `$LIPOLGEN_NNPDF_WHEEL`) |
| T1 F-2 | `t1_nmc_f2d_over_f2p.py` | `nnpdf_NMC_NC_NOTFIXED.json` — NMC F₂ᵈ/F₂ᵖ, NPB 487 (1997) 3, HEPData ins426595 tables 2–21, same wheel; 260 points, both the `hepdata` and the five-component `legacy` uncertainties; read 2026-10-26: 260/260 = the paper's Table 3; the five legacy components are its VX, SM, RC, E, E′ columns; no R enters the ratio (R_d = R_p); its normalisation < 0.1 % | **FAIL, recorded** | as F-1, 211 points in the cut: ToyF2 **2161.4/211** (p = 0). Sub-rows: MSTW 247.2/211 (p = 0.044, within its own tolerance; 260.7 through the Q²-frozen `f2n_over_f2p` hook); CT18NLO 215.5/189 on Q² ≥ 1.677 (p = 0.090, within its own tolerance; 239.5 through the hook; measured 2026-10-26, blocked (environment) before) | 0.78–0.88 s | (above) |
| T1 D-2, D-1, D-6 | `t1_g1d_world.py` | `nnpdf_COMPASS15_NC_NOTFIXED_MUD.json` (15), `nnpdf_HERMES_NC_7GEV_ED.json` (15), `nnpdf_E143_NC_NOTFIXED_ED.json` (28), `nnpdf_SMC_NC_NOTFIXED_MUD.json` (13, recorded only) — g₁ᵈ per nucleon, same wheel; the papers, read 2026-10-26, publish g₁ᵈ per nucleon and NOT D-state corrected (COMPASS PLB 769 p. 39, HERMES Eq. (23), E143 Eq. (56), SMC Eq. (13)), so the verdict column stands; E143 is hep-ph/9802357 (the metadata's hep-ex/9705012 is E154), SMC is PRD 58 (1998) 112001 Table IX; four high-x COMPASS points differ from the printed Table A.1, HEPData's being the self-consistent ones (`phase_D_inputs.md` §1) | **FAIL, recorded** | p ≥ 0.01 for each of COMPASS, HERMES, E143 on the cut: ToyG1 as g₁_nucleus(d)/2 = 0.9325·(g₁ᵖ + g₁ⁿ)/2 gives **107.8/15, 69.4/15, 54.9/28**. Sub-rows: MSTW-F₁ 53.2 / 47.1 / 50.2 (FAIL); measured 2026-10-26, the sets installed (blocked (environment) before): ToyG1 on CT18NLO's F₁ 20.2/11, 64.3/11, 59.0/24 on Q² ≥ 1.677 and NNPDFpol11_100's g₁ 45.0/15, 52.6/15, 51.0/28 (both FAIL; COMPASS alone passes on CT18NLO, p = 0.043); the no-D-state variant recorded (worse for NNPDFpol11 on all three) | 0.70–0.81 s | (above) |
| 12 (T1 R-1, R-2, R-3, R-5) | `t1_r_world.py` | `hepdata_ins424154_nmc_r.json` (NMC R, NPB 483 (1997) 3 Table 7, 10 points), `hepdata_ins786161_h1_fl.json` (H1 F_L, PLB 665 (2008) 139, 8 + 25), `hepdata_ins296980_whitlow_r.json` (Whitlow *et al.*, PLB 250 (1990) 193, 89), `hepdata_ins360765_e140_r.json` (E140, PRD 49 (1994) 5641, 11) — HEPData CC0, each checked against the paper's table where it prints one | **PASS** (added 2026-10-26) | p(χ², ndf) ≥ 0.01 per set on Q² ≥ 1, W² ≥ 4, stat ⊕ sys diagonal; status = `r_sigma_lt` on NMC: **12.54/10** (p = 0.25). Sub-rows: `r1998` 6.20/10; H1 8.82/8 and 8.18/8; Whitlow 209.76/68 (FAIL) and 63.05/68; E140 29.66/9 (FAIL) and 12.71/9 (Whitlow and E140 are R1998 fit inputs: a closure for `r1998`); NMC sys correlated (recorded) 60.22/10 and 21.28/10. (1 + r1998)/(1 + r_sigma_lt) at the four sweet-spot cells x = 0.011 / 0.028 / 0.071 / 0.141: 1.1848 / 1.1780 / 1.0762 / 0.9563 | 1.6–1.8 s | `test_bench_t1_r_world.py` (51, shared with row 14) |
| 14 (T1 F-7) | `t1_hermes_f2d.py` | `hepdata_ins894309_hermes_f2.json` — HERMES F₂ᵖ, F₂ᵈ per nucleon, JHEP 05 (2011) 126 Tables 5, 6 (81 + 81; the record's Table 2 σ column is the proton's, not used) | **FAIL, recorded** (added 2026-10-26) | p ≥ 0.01 on Q² ≥ 1, W² ≥ 4, per-bin errors diagonal + the 7.5 % normalisation fully correlated: ToyF2 **5303.5/68**. Sub-rows: MSTW 318.9/68, CT18NLO at LO 1788.2/45; the b₁ inversion reads per nucleon (0.929–1.039) | 0.5–0.6 s | (above) |
| T2 §5.2 | `t2_deuteron_static.py` | `deuteron_static_measured.json` — η_d = 0.0256 ± 0.0004 (Rodning–Knutson, PRC 41 (1990) 898, Eq. (10): a measurement), Q_d = +0.2860 ± 0.0015 fm² (Bishop–Cheung, PRA 20 (1979) 381: a calculation confirming Reid–Vaida's semi-empirical value); read from the papers 2026-10-26 (`phase_D_inputs.md` §2), Q_d corrected from the literature's 0.2859(3), which the paper does not print | **FAIL, recorded** | \|η − η_d\| ≤ 2σ for every shipped deuteron: CD-Bonn 0.025571 (−0.07σ), AV18 0.025045 (−1.39σ), **Hulthén** (the tagged channel's default, β = 0.30 GeV, P_D = 0.045) **0.032478, +17.2σ**. Q_d in impulse approximation recorded, not asserted: 0.2705 / 0.2697 / 0.2903 fm² (−5.42 / −5.71 / +1.52 % of 0.2860; −5.4 / −5.7 / +1.6 % of 0.2859 until 2026-10-26) | 1.3–1.5 s | `test_bench_static.py::test_deuteron_static_*` (9) |
| 11 | `t2_deeps_spectator_tail.py` | `clas_deeps_klimenko2006_f2P.json` (the row-10 file's Deeps table; the 60 blocks at Q² = 1.8 / 2.8, p_s the paper's bin averages 0.30–0.53 GeV/c; 604 of 643 rows, stat = 0 rows dropped) | **FAIL, recorded** | backward-window (cos θ_pq ≤ −0.3) shape, one free scale per (Q², W*, cos) group: χ²/ndf = **627.85/187** for the shipped Hulthén deuteron (light-cone map coded in the harness); AV18 579.12, CD-Bonn 278.81 recorded; tolerance p ≥ 0.01, declared before any χ²; transverse-window Glauber check recorded, not asserted | ≈ 2 s | `test_bench_deeps.py` (16 since 2026-10-26, one of them the figure) | `../../docs/open_items/run_2026-09-27/DECISIONS.md` D3 |
| 16 (T2 §5.4) | `t2_boeglin_nk.py` | `boeglin2011_deuteron_pmiss.csv` — Boeglin *et al.* PRL 107 (2011) 262501 Fig. 2, read back from the vector PDF (35 points with bars at θ_nq = 35 / 45 / 75° and the MS CD-Bonn/Paris PWIA/FSI curves at each) | **PASS** (Hulthén default, under A = 0.20, the upper end of the paper's 10–20 % FSI estimate — the Hulthén fails at 0.10, phase F; added 2026-10-26) | n(k = p_m), one free scale, 22 points at 35 + 45°, bars ⊕ A = 20 %, p ≥ 0.01: Hulthén **21.98/21** (p = 0.40), AV18 11.54, CD-Bonn 8.89; without A 185.9 / 66.5 / 60.3 (all fail); with the MS FSI factor 95.6 / 76.7 / 24.4 (only CD-Bonn) | 0.42–0.54 s | `test_bench_boeglin.py` (13) |
| T3 | `t3_li6_rms_devries.py` | `devries1987_li6_rms.json` — de Vries *et al.*, ADNDT 36 (1987) Table I, 2.54(5) / 2.56(5) / 2.57(10) fm; Angeli–Marinova 2.589(39) fm (closure, recorded); a survey transcription, not re-read | **PASS** (a consistency check of an input, not a prediction: the point radius is a tree input derived from Angeli–Marinova) | r_ch of the shipped `HoSpin1FF` (nucleon-folded, no Darwin–Foldy term) from the Richardson-extrapolated C0 slope: **2.57100 fm** (± 1.4e-9 numerical), within 2σ of each de Vries row ([2.46, 2.64] fm; pulls +0.62 / +0.22 / +0.01σ). The tolerance was chosen after the 2.5710 preview was known; the header says so | 0.10–0.13 s | `test_bench_static.py::test_devries_*` (7) |
| 15 (T3 T-20) | `t3_li6_ff_wiringa_schiavilla.py` | `ws1998_li6_formfactors.csv` — Wiringa–Schiavilla PRL 81 (1998) 4317 Fig. 1, read back from the vector PDF (C0/C2 lines IA and IA+MEC, F_L²/F_T² symbols with Monte Carlo bars, the Stanford/Amsterdam/Saskatoon data; plus 68 rows in `data/ff/li6_elastic.csv`'s columns — `data/ff/` not created); also reads `rfy1966_li_magnetic_ff.csv` for a recorded closure | **PASS** (one-sided; added 2026-10-26) | shipped C0 zero q₀ = **3.0998** fm⁻¹ ≥ q_last − (δ + 0.1) = 2.8143, q_last = 2.9190 the last IA+MEC C0 vertex WS98 draws; the upper side untested (the figure does not draw the zero). Recorded: vmc-ft no zero; F_M1 vs fm χ² 88.9/19; the μ and Q closures fail (q quantization; a 3/√2 in WS98's C2) | 0.40–0.50 s | `test_bench_li6_ff_theory.py` (18) |
| T4 D-1 | `t4_pythia_ep_closure.py` | none vendored: stock PYTHIA `WeakBosonExchange:ff2ff(t:gmZ)` at D-1's settings (`tests/test_pythia.cpp`), same library and default PDF (pSet 13) as the bridge, run in process | **PASS** (opt-in) | \|bridge/stock − 1\| ≤ 0.20 per observable per window (D-1's), at matched (x, Q²), e 10 × p 100 GeV, 20 000 events per side per window: worst **0.076** on the pinned PYTHIA 8.317 (W2 n_ch_eta[−2,0), a low-count bin; ⟨n_ch⟩ 1.007 / 0.997; 2026-10-26), **0.023** on 8.312 (⟨n_ch⟩ 1.003 / 0.998; 2026-09-26), over 15 sub-rows | 9.2–10.2 s (14 s on 2026-10-26's machine) | `test_bench_t4_t6.py::test_t4_*` (10, two of them T4 / T6 exit 2 without `lipolgen`, one the figure; the full row behind `LIPOLGEN_BENCH`, its means re-pinned on 8.317 on 2026-10-26) |
| T5 §1.2 | `t5_sum_rules.py` | embedded, sha256-checked: g_A, M_Z, m_b, m_c from PDG 2026 (`pdg.sqlite` in the PyPI `pdg-2026.0` wheel, CC BY 4.0); α_s(M_Z) = 0.1180(9) not in it — verified 2026-10-26 against PDG 2026 Table 1.1 and the QCD review's Eq. (9.25) (`phase_D_inputs.md`; until then from memory) | **FAIL, recorded** (any sub-row fails) | BC / g₂^WW implementation closure **PASS**: 45/45 cells, worst \|resid\|/bound 0.667. Bjorken on ToyG1 **FAIL**: at Q² = 5 the truncated ∫(g₁ᵖ − g₁ⁿ) = 0.6430 against g_A/6·C_Bj = 0.1854 ± 0.0023, not converged (g₁ⁿ ~ x^−1.22 is not integrable); tolerance converged and within 10 %. Bjorken NNPDFpol11 **FAIL** (measured 2026-10-26, the set installed; blocked (environment) before): 0.2087 / 0.2084 / 0.2083 at Q² = 2 / 5 / 10, +19.3 / +12.4 / +9.6 % from g_A/6·C_Bj and not converged (the last decade adds 2.16 / 2.20 / 2.25 % of the reference against ≤ 1 %) — about g_A/6 without C_Bj, the leading-order g₁ formula's signature; FSI unitarity blocked (not definable in the tree as built) | 1.4–1.6 s | `test_bench_t5_closures.py::test_sr_*` (13; the NNPDFpol11_100 one skips, naming the set, where it is missing) and `::test_t5_sum_rules_plot_writes_png` (the figure, since 2026-10-26) |
| T5 §1.1 | `t5_weighted_unweighted.py` | none: the same physics point by the other route (⁶Li inclusive config 1, x ≥ 0.1, y ≥ 0.5, helicity-flip `apar+`) | **PASS** | \|pull\| ≤ 3: unweighted 89.1934 ± 0.0512 pb vs Mode W 89.1674 ± 0.0516 pb, pull **+0.359** (≈ 3.0e6 events per route, seed 20260926; a no-op Mode W would pull +18.3). The window replaced a first design that could not fail; the header records both | 2.2–2.5 s | `test_bench_t5_closures.py::test_wu_*` (5) |
| 13 (T5 C2) | `t5_delta_moment_bound.py` | `delta_moment_bound_sources.json` — the bound's six sources (DS16 Eqs. (3)–(5), (9)–(14), (32), (36); SS90 (2.5), (2.11), (2.12), (3.5); LOI18 (6)–(9); CvDM17 (34)–(36); CRSZ25 (2); KS20 (22), (67), (68)) transcribed with each PDF's sha256 (no PDF vendored), the CT18NLO and NNPDFpol11_100 set identities (`.info` and member-file sha256, checked against the installed sets: another version is a `BLOCKED` row with an `environment:` reason; a missing set exits 2 and is `blocked (environment)` since phase F), the letter's `moment_A` / `moment_B` configuration, the grids | **PASS** (a bound under declared assumptions; added 2026-10-26) | shipped `toy_delta_gluon` at 10⁻²: \|M₂\|/B = **0.0644 / 0.0803 / 0.0941** at Q² = 2 / 5 / 10 GeV², margin ≥ −1e-12 required — the Jaffe–Manohar normalisation, the kernel's since D4 (2026-10-27), under which A7 is the consistent identification (V6, kept as declared and relabelled the pre-D4 kernel / the letter as printed, ×2: 0.129 / 0.161 / 0.188; tight edge, V6 kept, 0.272 / 0.337 / 0.394); helicity-state rates min W = 0.9494 (0.9493997990 since D4, the Δ term halved), A_T in [−0.101, 0.051] against [−2, 1]. Recorded: the letter's `moment_A` / `moment_B` (dilution 1/3) 0.715 / 0.699 / 0.702 and 0.738 / 0.695 / 0.674 centrally (at face value), 1.43 / 1.40 / 1.40 and 1.48 / 1.39 / 1.35 under V6 (as printed); at dilution 1 2.0–2.2 (4.0–4.4 under V6) | 2.8–3.1 s | `test_bench_t5_delta_bound.py` (21; 23 since phase F — a missing set is now blocked (environment), exit 2, and the LHAPDF search order is LHAPDF's) |
| T6 | `t6_hepmc3_convention.py` | `docs/HEPMC3_CONVENTION.md` (sha256 `56d3a91a…` since D1, 23 clauses, each tied to its line) — a document, not data | **PASS** (FAIL 21/23 until 2026-09-27) | every clause on every event, the doc's 1e-9 — for four-momentum read per component relative to the summed beam energy (200–700× looser on p_x, p_y than `tests/test_hepmc.cpp`'s doctest rule; the verdict is the same under either): **23 of 23** hold on 8 `lipolgen-run` files (800 events: inclusive ⁶Li/⁷Li, tagged ⁶Li-α/⁷Li-α/d-p, coherent at T2; inclusive at T0; one `--rc tensor-band`). C13 (four-momentum) and C14 (charge) failed on the inclusive files, which leave out the (A−1) remnant by design (`docs/USAGE.md` §2); decision D1 (`../../docs/open_items/run_2026-09-27/DECISIONS.md`) scoped the doc: inclusive records balance per nucleon | 3.7–4.2 s until 2026-10-27; 16.2–16.4 s since D5 (2026-10-29, unloaded; 3 runs) (the `--rc tensor-band` file builds the default `polrad-full` table) | `test_bench_t4_t6.py::test_t6_*` (9) |
