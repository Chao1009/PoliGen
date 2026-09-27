<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# validation/benchmarks — the automated benchmark harnesses

`docs/benchmarking/BENCHMARK_PLAN.md` §3: *a benchmark that is not automated is
not a benchmark.* Every wired row of the plan is one file here. This page states
the conventions once and lists every harness; the evidence for each number is in
the run record named in the last column, and the numbers themselves are
regenerated into `docs/benchmarking/REPORT.md`. (Merged 2026-09-23 from this
file and the nuclear group's `README_nuclear.md`, which it replaces; the nine
harnesses of 2026-09-26 and the report generator added on that run.)

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
   DJANGOH harness prints its row on a `ROW |` line. **`docs/benchmarking/REPORT.md` and
   `REPORT.json` exist since 2026-09-26** and are generated from the dicts
   `run()` returns, never from the printed lines and never by hand:

       cd LiPolGen && source env.sh
       python3 validation/benchmarks/make_report.py --all      # regenerate both
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
   configuration needs one (row 4) raises, exits 2, and is `blocked
   (environment)` in the report. T4 and T6 report a build tier that is off
   (`HAVE_PYTHIA8`, `HAVE_HEPMC3`) or a missing `pythia8` / `pyhepmc` module as
   a `BLOCKED` row with an `environment:` reason; a missing or broken
   `lipolgen` itself is the harness broken, exit 2, as everywhere else (until
   2026-09-27 T4 and T6 turned it into a `BLOCKED` row with exit 0).

**A benchmark that disagrees with the tree is RECORDED, never tuned.** No
harness moves a default, a registry row or a reference JSON; a `FAIL` is a
measurement with its distance, priced for the registry row it bears on. Every
tolerance is stated in the harness header independently of the verdict, and
where the number was seen before the tolerance was chosen (row 4, the three
T1 rows, CD-Bonn's η in `t2_deuteron_static`, `t3_li6_rms_devries`, the BC
preview in `t5_sum_rules`), the header says so.

**Blocked rows.** A row whose input is neither on disk nor open access is
`BLOCKED`: its harness exists, and its test skips with a message naming the
exact document to download (citation and link as in
`docs/references/REFERENCES.md` §1). A row can also be blocked because the
comparison it was planned as does not exist (row 7), or because the tree as
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
(16.2–16.5 s) is not opt-in, so every other row runs in the ordinary
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
project pins 8.317).

**Tools here that are not harnesses** (their names do not match `t[1-6]_*.py`):
`make_report.py` and `_report.py` (the report, above); `tiers.py`
(BENCHMARK_PLAN.md §7 as data: `--since REV` lists the harnesses a diff must
re-run, `make_report.py --since REV` runs them); `vendor_nnpdf_commondata.py`
(dev-time: rebuilds the six `nnpdf_*.json` files from the nnpdf-data 4.1.5
wheel, sha256 `6ed209ad…3490269a`, checking every file it reads against the
wheel's RECORD; `--check` regenerates in memory and exits 1 on any byte
difference).

## The harnesses (status and runtime measured by run 2026-09-26)

Runtime is the wall time of `python3 validation/benchmarks/<file>`,
interpreter start-up included, as the range of three runs measured on
2026-09-27 (UTC) by the integration stage of run 2026-09-26 with
`time.perf_counter` around a subprocess (`/usr/bin/time` is not installed on
that machine: PYTHIA 8.312, no LHAPDF set). Row 4 cannot run there; its
status and runtime are 2026-09-23's, from the development machine. The status
column is the row's own verdict; `REPORT.md` has every sub-row.

### The plan's first ten (BENCHMARK_PLAN.md §4)

| row | harness | reference, vendored under `data/` | status | the number, with its tolerance | runtime | tests | record |
|---|---|---|---|---|---|---|---|
| 1 | `t5_epios_source_modes.py` | `epios2026_table2_deuteron_source_modes.json` — EPIOS, PRC 113 (2026) 060501, Table II (COSY/ANKE deuteron source; the P_zz column is ideal, only P_z is measured) | **PASS** | max moment round-trip residual **0** over 8 modes, tolerance 1e-12. Recorded, not asserted: the tree's EST (`ladder`) fill matches the ideal P_zz at mode 0 only, **1 of the 6 modes it computes**; at modes 5, 6 (\|P_z\| = 1) it throws and only the closed-form limit matches ("3 of 8" until 2026-09-26 counted those two) | 0.11–0.12 s | `test_bench_spin_deuteron.py::test_epios_*` (5); doctest *"bench row 1: …"* | `phase_B1_spin_deuteron.md` §1 |
| 3 | `t2_hermes_b1_table2.py` | `hermes2005_b1d_table2.json` — HERMES, PRL 95 (2005) 242001, Table II (per nucleon); `whitlow1990_r1990.json` — R1990 (the scan's b₁ = 0.635 read as 0.0635) | **FAIL, recorded** | every bin inside stat ⊕ syst: met by no configuration. χ²(b₁ᵈ)/6 = 22.26 (convolution + MSTW, the A = 2 gate's), 5.26 (shipped Miller ×0.5), 5.30 (Miller ×1), 21.80 (b₁ = 0) — HERMES does not decide registry row 2. χ²(b₁) does not depend on the F₂ backend (pinned without MSTW too) | 4.5–4.9 s | `test_bench_spin_deuteron.py::test_hermes_*` (7), `::test_r1990_*` (2) | `phase_B1_spin_deuteron.md` §3 |
| 4 | `t3_nmc_li6_over_d.py` | `hepdata_ins394050_table1.csv` — NMC F₂(⁶Li)/F₂(D), HEPData ins394050 Table 1, CC0 | **PASS** (2026-09-23; blocked (environment) where the sets are missing) | p(χ², ndf) ≥ 0.01 on both sets: **4.6269/4** on x ≥ 0.30 (p = 0.328), **26.941/15** on the 15 on-grid points (p = 0.029); 9 points below the EPPS21 grid reported, not compared | 1.3–1.4 s (2026-09-23) | `test_bench_nuclear.py::test_nmc_*` (4; the table guard runs without the sets) | `phase_B2_nuclear.md` §2 |
| 5 | `t3_li6_charge_ff_fb.py` | `uva_ncd_fb_data_li6.dat` — the ⁶Li row of the UVa NCD archive's `FB_data.dat` (no licence stated; provenance of the row UNVERIFIED) | **FAIL, recorded** | the C0 zero must lie in [2.6944 (FB zero), 2.8284 (Li *et al.* 1971 minimum)] fm⁻¹: `HoSpin1FF`'s is **3.0998**, **+0.2713 fm⁻¹** above; T11's [2.9, 3.3] not moved (open item Q1). A C0 with no zero below 5 fm⁻¹ is a recorded FAIL, not a crash | 0.11–0.14 s | `test_bench_nuclear.py::test_fb_*` (4) | `phase_B2_nuclear.md` §3 |
| 6 | `t3_li_magnetization_rfy.py` | *none on disk* — awaits `rfy1966_li_magnetization.csv` (format in the header) | **BLOCKED** | needs Rand–Frosch–Yearian, Phys. Rev. 144 (1966) 859, or de Jager *et al.*, ADNDT 14 (1974) 479 Table V — REFERENCES.md §1 #8, #9. Tree side measured, uncompared: r_mag(⁶Li) = 2.9469 fm | 0.10–0.11 s | `test_bench_nuclear.py::test_rfy_*` (2; one skips naming both) | `phase_B2_nuclear.md` §4 |
| 7 | `t4_djangoh_rad_noRad.py` | `djangoh_rad_norad_ep18x275.json` — `eic/InclusiveDjangohSamples` at `9869d9a`, four Q² bins, e 18 × p 275 GeV | **BLOCKED — not a benchmark as planned** | none definable: all eight DJANGOH logs have IEL2 = IEL31 = IEL32 = IEL33 = 0 (elastic radiative tail OFF), and `rc_tail` is the elastic (+ nuclear QE) tail only, so the two share no O(α) term. Both sides printed: DJANGOH 1.0725 / 1.0894 / 1.1606 / 1.2280, LiPolGen proton t-peak 1.0162 / 1.0132 / 1.0150 / 1.0207. The samples' W_h ≥ 3 GeV cut also excludes the elastic tail (W_h = M_p) whatever IEL is — DJANGOH's `HSPRLG` zeroes IEL2, IEL31..33 when WMIN > M_p (public source, `github.com/spiesber/DJANGOH` at `e356b84`) — so it unblocks only with a DJANGOH run at IEL31..33 ≠ 0 **and** WMIN ≤ M_p. LiPolGen's W² ≥ 9 cut is on the leptonic W, not on W_h | 16.2–16.5 s | `test_bench_chain_rc.py` (4) | `phase_B3_chain_rc.md` §1 |
| 9 | `t5_est_identity.py` | `est_identity_sources.json` — the EST closed form; COMPASS ⁶LiD (Koivuniemi *et al.*, SPIN 2004) | **PASS** | max \|`spin_temperature_pzz(1, P_z)` − (2 − √(4 − 3P_z²))\| = **1.138e-14** over 19 999 interior points, tolerance **2.31e-14** (bisection bound 2.13e-14; the planned 1e-14 does not hold); COMPASS numbers at the printed digit | 0.14–0.16 s | `test_bench_spin_deuteron.py::test_est_*` (6); doctest *"bench row 9: …"* | `phase_B1_spin_deuteron.md` §2 |
| 10 | `t2_bonus_spectator_shape.py` | `clas_db_deuteron_query_2026-09-23.json` (the CLAS Physics Database catalogue); `clas_deeps_klimenko2006_f2P.json` (Deeps, 115 blocks, vendored, not wired: 60 at the paper's Q² 1.8 / 2.8, 55 with unverified labels, one empty; 60 rows value = stat = 0 are missing, not zeros) | **BLOCKED** | needs Tkachenko *et al.*, PRC 89 (2014) 045206, Supplemental Material (HTTP 401 without an APS login) — REFERENCES.md §1 #59 | 0.03 s | `test_bench_spin_deuteron.py::test_bonus_*` (4; `test_bonus_row_is_blocked_by_name` skips naming it) | `phase_B1_spin_deuteron.md` §4 |

Rows 2 and 8 of the plan have no harness: row 2 was resolved in code on
2026-09-06 (`tests/test_tagged.cpp`'s Cosyn–Weiss gate), and row 8 is a
relabelling (George–Knutson's η is a consistency band on the quadrupole dial,
not a benchmark). Records are under `docs/open_items/run_2026-09-23/`.

### By tier, added 2026-09-26 (BENCHMARK_PLAN.md §4, "Then, by tier")

The "row" is the survey row the harness wires (`02_data_nucleon_deuteron.md`,
`03_data_nuclear.md`, `05_chain.md`). None was tuned; the record is
`docs/open_items/run_2026-09-26/SUMMARY.md` and the harness header.

| row | harness | reference, vendored under `data/` | status | the number, with its tolerance | runtime | tests |
|---|---|---|---|---|---|---|
| T1 F-1 | `t1_nmc_f2d.py` | `nnpdf_NMC_NC_NOTFIXED_D.json` — NMC F₂ᵈ per nucleon, NPB 483 (1997) 3, HEPData ins424154 tables 17–32 (CC0), via the PyPI nnpdf-data 4.1.5 wheel (GPL-3.0-or-later); 158 points | **FAIL, recorded** | p(χ², ndf) ≥ 0.01 on the DIS cut Q² ≥ 1, W² ≥ 4 GeV² (155 points), stat ⊕ sys diagonal: shipped ToyF2 **8307.9/155** (p = 0). Sub-rows: MSTW 416.5/155 (p = 8.2e-26); CT18NLO blocked (environment). Full covariance recorded beside it (21100.7; MSTW 858.5) | 0.78–0.85 s | `test_bench_t1_nucleon.py` (42 for the three T1 rows: `test_f2d_*` 6, `test_dp_*` 5, `test_g1d_*` 8, the rest shared — file guards, blocked LHAPDF configs by name, exit codes, the no-scipy χ² fallback; one needs `$LIPOLGEN_NNPDF_WHEEL`) |
| T1 F-2 | `t1_nmc_f2d_over_f2p.py` | `nnpdf_NMC_NC_NOTFIXED.json` — NMC F₂ᵈ/F₂ᵖ, NPB 487 (1997) 3, HEPData ins426595 tables 2–21, same wheel; 260 points, both the `hepdata` and the five-component `legacy` uncertainties | **FAIL, recorded** | as F-1, 211 points in the cut: ToyF2 **2161.4/211** (p = 0). Sub-rows: MSTW 247.2/211 (p = 0.044, within its own tolerance; 260.7 through the Q²-frozen `f2n_over_f2p` hook); CT18NLO blocked (environment) | 0.78–0.88 s | (above) |
| T1 D-2, D-1, D-6 | `t1_g1d_world.py` | `nnpdf_COMPASS15_NC_NOTFIXED_MUD.json` (15), `nnpdf_HERMES_NC_7GEV_ED.json` (15), `nnpdf_E143_NC_NOTFIXED_ED.json` (28), `nnpdf_SMC_NC_NOTFIXED_MUD.json` (13, recorded only) — g₁ᵈ per nucleon, same wheel; whether the data are D-state corrected is UNVERIFIED on every file | **FAIL, recorded** | p ≥ 0.01 for each of COMPASS, HERMES, E143 on the cut: ToyG1 as g₁_nucleus(d)/2 = 0.9325·(g₁ᵖ + g₁ⁿ)/2 gives **107.8/15, 69.4/15, 54.9/28**. Sub-rows: MSTW-F₁ 53.2 / 47.1 / 50.2 (FAIL); CT18NLO and NNPDFpol11_100 blocked (environment); the no-D-state variant recorded | 0.70–0.81 s | (above) |
| T2 §5.2 | `t2_deuteron_static.py` | `deuteron_static_measured.json` — η_d = 0.0256(4) (Rodning–Knutson, PRC 41 (1990) 898), Q_d = +0.2859(3) fm² (Bishop–Cheung, PRA 20 (1979) 381); a second-hand transcription (02 §5.2, `b1_nuclear.hpp`), the papers not re-read | **FAIL, recorded** | \|η − η_d\| ≤ 2σ for every shipped deuteron: CD-Bonn 0.025571 (−0.07σ), AV18 0.025045 (−1.39σ), **Hulthén** (the tagged channel's default, β = 0.30 GeV, P_D = 0.045) **0.032478, +17.2σ**. Q_d in impulse approximation recorded, not asserted: 0.2705 / 0.2697 / 0.2903 fm² (−5.4 / −5.7 / +1.6 %) | 1.3–1.5 s | `test_bench_static.py::test_deuteron_static_*` (9) |
| 11 | `t2_deeps_spectator_tail.py` | `clas_deeps_klimenko2006_f2P.json` (the row-10 file's Deeps table; the 60 blocks at Q² = 1.8 / 2.8, p_s the paper's bin averages 0.30–0.53 GeV/c; 604 of 643 rows, stat = 0 rows dropped) | **FAIL, recorded** | backward-window (cos θ_pq ≤ −0.3) shape, one free scale per (Q², W*, cos) group: χ²/ndf = **627.85/187** for the shipped Hulthén deuteron (light-cone map coded in the harness); AV18 579.12, CD-Bonn 278.81 recorded; tolerance p ≥ 0.01, declared before any χ²; transverse-window Glauber check recorded, not asserted | ≈ 2 s | `test_bench_deeps.py` (15) | `../../docs/open_items/run_2026-09-27/DECISIONS.md` D3 |
| T3 | `t3_li6_rms_devries.py` | `devries1987_li6_rms.json` — de Vries *et al.*, ADNDT 36 (1987) Table I, 2.54(5) / 2.56(5) / 2.57(10) fm; Angeli–Marinova 2.589(39) fm (closure, recorded); a survey transcription, not re-read | **PASS** (a consistency check of an input, not a prediction: the point radius is a tree input derived from Angeli–Marinova) | r_ch of the shipped `HoSpin1FF` (nucleon-folded, no Darwin–Foldy term) from the Richardson-extrapolated C0 slope: **2.57100 fm** (± 1.4e-9 numerical), within 2σ of each de Vries row ([2.46, 2.64] fm; pulls +0.62 / +0.22 / +0.01σ). The tolerance was chosen after the 2.5710 preview was known; the header says so | 0.10–0.13 s | `test_bench_static.py::test_devries_*` (7) |
| T4 D-1 | `t4_pythia_ep_closure.py` | none vendored: stock PYTHIA `WeakBosonExchange:ff2ff(t:gmZ)` at D-1's settings (`tests/test_pythia.cpp`), same library and default PDF (pSet 13) as the bridge, run in process | **PASS** (opt-in) | \|bridge/stock − 1\| ≤ 0.20 per observable per window (D-1's), at matched (x, Q²), e 10 × p 100 GeV, 20 000 events per side per window: worst **0.023** over 15 sub-rows (⟨n_ch⟩ 1.003 / 0.998) | 9.2–10.2 s | `test_bench_t4_t6.py::test_t4_*` (9, two of them T4 / T6 exit 2 without `lipolgen`; the full row behind `LIPOLGEN_BENCH`) |
| T5 §1.2 | `t5_sum_rules.py` | embedded, sha256-checked: g_A, M_Z, m_b, m_c from PDG 2026 (`pdg.sqlite` in the PyPI `pdg-2026.0` wheel, CC BY 4.0); α_s(M_Z) = 0.1180(9) not in it — from memory, the least-verified input | **FAIL, recorded** (any sub-row fails) | BC / g₂^WW implementation closure **PASS**: 45/45 cells, worst \|resid\|/bound 0.667. Bjorken on ToyG1 **FAIL**: at Q² = 5 the truncated ∫(g₁ᵖ − g₁ⁿ) = 0.6430 against g_A/6·C_Bj = 0.1854 ± 0.0023, not converged (g₁ⁿ ~ x^−1.22 is not integrable); tolerance converged and within 10 %. Bjorken NNPDFpol11 blocked (environment); FSI unitarity blocked (not definable in the tree as built) | 1.4–1.6 s | `test_bench_t5_closures.py::test_sr_*` (13; one skips naming NNPDFpol11_100) |
| T5 §1.1 | `t5_weighted_unweighted.py` | none: the same physics point by the other route (⁶Li inclusive config 1, x ≥ 0.1, y ≥ 0.5, helicity-flip `apar+`) | **PASS** | \|pull\| ≤ 3: unweighted 89.1934 ± 0.0512 pb vs Mode W 89.1674 ± 0.0516 pb, pull **+0.359** (≈ 3.0e6 events per route, seed 20260926; a no-op Mode W would pull +18.3). The window replaced a first design that could not fail; the header records both | 2.2–2.5 s | `test_bench_t5_closures.py::test_wu_*` (5) |
| T6 | `t6_hepmc3_convention.py` | `docs/HEPMC3_CONVENTION.md` (sha256 `56d3a91a…` since D1, 23 clauses, each tied to its line) — a document, not data | **PASS** (FAIL 21/23 until 2026-09-27) | every clause on every event, the doc's 1e-9 — for four-momentum read per component relative to the summed beam energy (200–700× looser on p_x, p_y than `tests/test_hepmc.cpp`'s doctest rule; the verdict is the same under either): **23 of 23** hold on 8 `lipolgen-run` files (800 events: inclusive ⁶Li/⁷Li, tagged ⁶Li-α/⁷Li-α/d-p, coherent at T2; inclusive at T0; one `--rc tensor-band`). C13 (four-momentum) and C14 (charge) failed on the inclusive files, which leave out the (A−1) remnant by design (`docs/USAGE.md` §2); decision D1 (`../../docs/open_items/run_2026-09-27/DECISIONS.md`) scoped the doc: inclusive records balance per nucleon | 3.7–4.2 s | `test_bench_t4_t6.py::test_t6_*` (8) |
