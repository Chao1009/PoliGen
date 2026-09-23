<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# validation/benchmarks — the automated benchmark harnesses

`docs/benchmarking/BENCHMARK_PLAN.md` §3: *a benchmark that is not automated is
not a benchmark.* Every wired row of the plan is one file here. This page states
the conventions once and lists every harness; the evidence for each number is in
the run record named in the last column. (Merged 2026-09-23 from this file and
the nuclear group's `README_nuclear.md`, which it replaces.)

## Conventions

A harness is `<tier>_<name>.py` (`t1_` nucleon, `t2_` deuteron, `t3_` nuclear
inputs, `t4_` generator-vs-generator, `t5_` MC closures, `t6_` chain), and each:

1. carries the three normalisation rules of BENCHMARK_PLAN.md §8 in its header,
   stated **for that row** on both sides (per nucleon / per nucleus / per
   deuteron; the A-scaling assumption; the isoscalar denominator), or says
   explicitly why a rule does not apply to it;
2. reads its reference from a **vendored** table under `data/`, whose file
   carries its provenance (URL, date fetched, licence, the paper, the
   publishing convention; the nuclear-group files also carry a sha256 the
   harness re-checks on every read). No harness and no test needs the network;
3. runs the generator or the library function at the reference's kinematics;
4. prints one `REPORT | name | generator value | reference | tolerance | STATUS`
   line (`PASS` / `FAIL` / `BLOCKED`; the DJANGOH harness prints the same
   fields on a `ROW |` line) and returns the row from `run()`. No `REPORT.md`
   generator exists yet;
5. is exercised by a pytest under `python/tests/test_bench_*.py` (rows 1 and 9
   also by `tests/test_bench_spin.cpp`);
6. **exits 0** whenever it ran and printed its row — `PASS`, recorded `FAIL` and
   `BLOCKED` alike, since the row carries the verdict — and **exits 2** only when
   the harness itself is broken (a missing or altered vendored file, a missing
   dependency such as an LHAPDF set, any exception). A harness is a measurement,
   not a CI gate: the pytests are the gate. (Uniform since 2026-09-23's fix
   stage; before it, five harnesses exited 1 on a `FAIL` and three exited 0.)

**A benchmark that disagrees with the tree is RECORDED, never tuned.** No
harness moves a default, a registry row or a reference JSON; a `FAIL` is a
measurement with its distance, priced for the registry row it bears on.

**Blocked rows.** A row whose input is neither on disk nor open access is
`BLOCKED`: its harness exists, and its test skips with a message naming the
exact document to download (citation and link as in
`docs/references/REFERENCES.md` §1). A row can also be blocked because the
comparison it was planned as does not exist (row 7).

**The opt-in variable.** A row that costs more than 30 s (a full generation, a
chain leg, a sample) is to run only with

    LIPOLGEN_BENCH=1 python -m pytest python/tests -q

and otherwise be reported as skipped with that reason. **No row needs it
today** — the slowest, row 7, takes 17.3–19.6 s — so every row runs in the
ordinary `python -m pytest python/tests -q`, and no test reads the variable yet.

**Environment.**

    cd LiPolGen && source env.sh     # build/python on PYTHONPATH, LHAPDF sets
    python3 validation/benchmarks/t5_epios_source_modes.py
    python -m pytest python/tests -q -k bench -rs

Row 4 needs LHAPDF with `EPPS21nlo_CT18Anlo_Li6` and `CT18NLO` installed
(`env.sh` points `LHAPDF_DATA_PATH` at them); its tests skip, naming the sets,
if they are not.

## The harnesses (status and runtime measured 2026-09-23)

Runtime is the wall time of `python3 validation/benchmarks/<file>` on the
development machine, interpreter start-up included, as the range measured on
2026-09-23 (three runs each by the fix stage on an idle machine, plus the
verification lenses' runs on a loaded one: row 7's 19.5–19.6 s, row 3's 4.4 s).

| row | harness | reference, vendored under `data/` | status | the number, with its tolerance | runtime | tests | record |
|---|---|---|---|---|---|---|---|
| 1 | `t5_epios_source_modes.py` | `epios2026_table2_deuteron_source_modes.json` — EPIOS, PRC 113 (2026) 060501, Table II (COSY/ANKE deuteron source; the P_zz column is ideal, only P_z is measured) | **PASS** | max moment round-trip residual **0** over 8 modes, tolerance 1e-12. Recorded, not asserted: the EST (`ladder`) fill matches the ideal P_zz at 3 of 8 modes | 0.3 s | `test_bench_spin_deuteron.py::test_epios_*` (4); doctest *"bench row 1: …"* | `phase_B1_spin_deuteron.md` §1 |
| 3 | `t2_hermes_b1_table2.py` | `hermes2005_b1d_table2.json` — HERMES, PRL 95 (2005) 242001, Table II (per nucleon); `whitlow1990_r1990.json` — R1990 (the scan's b₁ = 0.635 read as 0.0635) | **FAIL, recorded** | every bin inside stat ⊕ syst: met by no configuration. χ²(b₁ᵈ)/6 = 22.26 (convolution + MSTW, the A = 2 gate's), 5.26 (shipped Miller ×0.5), 5.30 (Miller ×1), 21.80 (b₁ = 0) — HERMES does not decide registry row 2 | 4.0–4.4 s | `test_bench_spin_deuteron.py::test_hermes_*` (6), `::test_r1990_typo_guard` | `phase_B1_spin_deuteron.md` §3 |
| 4 | `t3_nmc_li6_over_d.py` | `hepdata_ins394050_table1.csv` — NMC F₂(⁶Li)/F₂(D), HEPData ins394050 Table 1, CC0 | **PASS** | p(χ², ndf) ≥ 0.01 on both sets: **4.6269/4** on x ≥ 0.30 (p = 0.328), **26.941/15** on the 15 on-grid points (p = 0.029); 9 points below the EPPS21 grid reported, not compared | 1.3–1.4 s | `test_bench_nuclear.py::test_nmc_*` (4) | `phase_B2_nuclear.md` §2 |
| 5 | `t3_li6_charge_ff_fb.py` | `uva_ncd_fb_data_li6.dat` — the ⁶Li row of the UVa NCD archive's `FB_data.dat` (no licence stated; provenance of the row UNVERIFIED) | **FAIL, recorded** | the C0 zero must lie in [2.6944 (FB zero), 2.8284 (Li *et al.* 1971 minimum)] fm⁻¹: `HoSpin1FF`'s is **3.0998**, **+0.2713 fm⁻¹** above; T11's [2.9, 3.3] not moved (open item Q1) | 0.3 s | `test_bench_nuclear.py::test_fb_*` (3) | `phase_B2_nuclear.md` §3 |
| 6 | `t3_li_magnetization_rfy.py` | *none on disk* — awaits `rfy1966_li_magnetization.csv` (format in the header) | **BLOCKED** | needs Rand–Frosch–Yearian, Phys. Rev. 144 (1966) 859, or de Jager *et al.*, ADNDT 14 (1974) 479 Table V — REFERENCES.md §1 #8, #9. Tree side measured, uncompared: r_mag(⁶Li) = 2.9469 fm | 0.3 s | `test_bench_nuclear.py::test_rfy_*` (2; one skips naming both) | `phase_B2_nuclear.md` §4 |
| 7 | `t4_djangoh_rad_noRad.py` | `djangoh_rad_norad_ep18x275.json` — `eic/InclusiveDjangohSamples` at `9869d9a`, four Q² bins, e 18 × p 275 GeV | **BLOCKED — not a benchmark as planned** | none definable: all eight DJANGOH logs have IEL2 = IEL31 = IEL32 = IEL33 = 0 (elastic radiative tail OFF), and `rc_tail` is the elastic (+ nuclear QE) tail only, so the two share no O(α) term. Both sides printed: DJANGOH 1.0725 / 1.0894 / 1.1606 / 1.2280, LiPolGen proton t-peak 1.0162 / 1.0132 / 1.0150 / 1.0207. Unblocks with a DJANGOH run at IEL31..33 ≠ 0 | 17.3–19.6 s | `test_bench_chain_rc.py` (3) | `phase_B3_chain_rc.md` §1 |
| 9 | `t5_est_identity.py` | `est_identity_sources.json` — the EST closed form; COMPASS ⁶LiD (Koivuniemi *et al.*, SPIN 2004) | **PASS** | max \|`spin_temperature_pzz(1, P_z)` − (2 − √(4 − 3P_z²))\| = **1.138e-14** over 19 999 interior points, tolerance **2.31e-14** (bisection bound 2.13e-14; the planned 1e-14 does not hold); COMPASS numbers at the printed digit | 0.35–0.37 s | `test_bench_spin_deuteron.py::test_est_*` (5); doctest *"bench row 9: …"* | `phase_B1_spin_deuteron.md` §2 |
| 10 | `t2_bonus_spectator_shape.py` | `clas_db_deuteron_query_2026-09-23.json` (the CLAS Physics Database catalogue); `clas_deeps_klimenko2006_f2P.json` (Deeps, 115 tables, vendored, not wired) | **BLOCKED** | needs Tkachenko *et al.*, PRC 89 (2014) 045206, Supplemental Material (HTTP 401 without an APS login) — REFERENCES.md §1 #59 | 0.06–0.08 s | `test_bench_spin_deuteron.py::test_bonus_row_is_blocked_by_name` (skips naming it) | `phase_B1_spin_deuteron.md` §4 |

Rows 2 and 8 of the plan have no harness: row 2 was resolved in code on
2026-09-06 (`tests/test_tagged.cpp`'s Cosyn–Weiss gate), and row 8 is a
relabelling (George–Knutson's η is a consistency band on the quadrupole dial,
not a benchmark). Records are under `docs/open_items/run_2026-09-23/`.
