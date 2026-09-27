<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
<!-- GENERATED FILE -- written by validation/benchmarks/make_report.py; DO NOT EDIT. Regenerate it instead. -->
# Benchmark report

> **Generated — do not edit.** Written by `validation/benchmarks/make_report.py` from the harnesses `validation/benchmarks/t[1-6]_*.py`, together with `REPORT.json` beside it (the same rows, machine-readable: the M2 baseline). `BENCHMARK_PLAN.md` §3: *the report is regenerated, never hand-edited* — `make_report.py --check` fails on any difference other than the date, the revision, a runtime or the environment.

- **generated:** 2026-09-27T01:35:01Z
- **git revision:** `f0a8f1e` + uncommitted changes (dirty)
- **environment:** Python 3.11.16; numpy 2.4.6; scipy 1.17.1; PYTHIA 8.312 (pythia8-config --version; env.sh pins 8.317); LHAPDF 6.5.6; HepMC3 3.03.01 (pyhepmc 2.16.1); lipolgen 0.1.0 (LHAPDF tier on, PYTHIA tier on)
- **LHAPDF sets present:** **none** (searched: `/home/user/deps/install/share/LHAPDF`, `/opt/deps/share/LHAPDF`)
- **opt-in (expensive) rows:** run (`--all` / `LIPOLGEN_BENCH=1`)
- **this file was made by:** `cd LiPolGen && source env.sh && python3 validation/benchmarks/make_report.py --all`
- **check it:** `cd LiPolGen && source env.sh && python3 validation/benchmarks/make_report.py --check`; **M2 compare:** `cd LiPolGen && source env.sh && python3 validation/benchmarks/make_report.py --all --compare docs/benchmarking/REPORT.json`

**17 harnesses:** 5 pass, 8 fail, 3 blocked, 1 blocked (environment), 0 error, 0 skipped (opt-in).

## The rows

One row per harness. *plan row* is the `BENCHMARK_PLAN.md` §4 row number when the harness is one of the first ten. The generator column is the harness's own value text; the sub-rows and the configuration are under *Detail*.

| harness | tier | plan row | status | generator (the number) | reference | tolerance | runtime |
|---|---|---|---|---|---|---|---|
| [`t1_g1d_world`](#t1_g1d_world) | T1 | — | **FAIL** | toy (shipped ToyG1, g1_nucleus(d)/2): compass 107.8/15 (p = 4.17e-16), hermes 69.4/15 (p = 5.75e-09), e143 54.9/28 (p = 0.00173) on Q2 &gt;= 1, W2 &gt;= 4 [diagonal]; sub-rows: mstw FAIL compass 53.2/15 (p = 3.52e-06), hermes 47.1/15 (p = 3.5e-05), e143 50.2/28 (p = 0.00607); ct18nlo BLOCKED (environment: LHAPDF set CT18NLO not installed); nnpdfpol BLOCKED (environment: LHAPDF set NNPDFpol11_100 not installed) | g1d per nucleon: COMPASS (D-2), HERMES (D-1), E143 (D-6) via nnpdf-data 4.1.5; SMC recorded only; 15 + 15 + 28 points in the DIS cut | p(chi2, ndf) &gt;= 0.01 on the DIS cut for each of compass, hermes, e143, stat (+) sys diagonal | 0.56 s |
| [`t1_nmc_f2d`](#t1_nmc_f2d) | T1 | — | **FAIL** | toy (shipped ToyF2): chi2/ndf = 8307.9/155, p = 0 on Q2 &gt;= 1, W2 &gt;= 4 [diagonal]; sub-rows: mstw 416.5/155 (p = 8.21e-26); ct18nlo BLOCKED (environment: LHAPDF set CT18NLO not installed) | NMC F2d per nucleon, NPB 483 (1997) 3, HEPData ins424154 tables 17-32 via nnpdf-data 4.1.5 (NMC_NC_NOTFIXED_D), 155 of 158 points in the DIS cut | p(chi2, ndf) &gt;= 0.01 on the DIS cut, stat (+) sys diagonal | 0.71 s |
| [`t1_nmc_f2d_over_f2p`](#t1_nmc_f2d_over_f2p) | T1 | — | **FAIL** | toy (shipped ToyF2): chi2/ndf = 2161.4/211, p = 0 on Q2 &gt;= 1, W2 &gt;= 4 [diagonal]; sub-rows: mstw 247.2/211 (p = 0.0444; frozen hook 260.7); ct18nlo BLOCKED (environment: LHAPDF set CT18NLO not installed) | NMC F2d/F2p (per nucleon), NPB 487 (1997) 3, HEPData ins426595 tables 2-21 via nnpdf-data 4.1.5 (NMC_NC_NOTFIXED), 211 of 260 points in the DIS cut | p(chi2, ndf) &gt;= 0.01 on the DIS cut, stat (+) sys diagonal | 0.72 s |
| [`t2_bonus_spectator_shape`](#t2_bonus_spectator_shape) | T2 | 10 | **BLOCKED** | not evaluated | none machine-readable (see reason) | — | 0.00 s |
| [`t2_deuteron_static`](#t2_deuteron_static) | T2 | — | **FAIL** | eta: cdbonn 0.025571 (-0.07 sigma), av18 0.025045 (-1.39 sigma), hulthen 0.032478 (+17.20 sigma); Q_d (impulse approx., recorded): 0.270493 (-5.39 %) / 0.269670 (-5.68 %) / 0.290349 (+1.56 %) fm^2 | eta_d = 0.0256(4) (Rodning-Knutson PRC 41 (1990) 898), Q_d = +0.2859(3) fm^2 (Bishop-Cheung PRA 20 (1979) 381); second-hand transcription | abs(eta - eta_d) &lt;= 2 sigma for each of cdbonn, av18, hulthen; Q_d recorded, not asserted | 1.26 s |
| [`t2_hermes_b1_table2`](#t2_hermes_b1_table2) | T2 | 3 | **FAIL** | cdks_conv_mstw: chi2/ndf b1 = 22.26/6, A_zz = 22.01/6; 2/6 b1 bins within tol | HERMES PRL 95 (2005) 242001 Table II, 6 bins, per nucleon | per bin sqrt(stat^2+syst^2) | 4.41 s |
| [`t3_li6_charge_ff_fb`](#t3_li6_charge_ff_fb) | T3 | 5 | **FAIL** | HoSpin1FF first C0 zero q0 = 3.0998 fm^-1 (shipped; T11 window [2.9, 3.3] not moved) | band [2.6944 (UVa FB zero), 2.8284 (Li71 minimum)] fm^-1 | q0 inside the band (must not contradict) | 0.06 s |
| [`t3_li6_rms_devries`](#t3_li6_rms_devries) | T3 | — | **PASS** | HoSpin1FF.for_ion(li6()) C0 slope, nucleon-folded (G_E^p + G_E^n): r_ch = 2.57100 fm (+- 1.4e-09 numerical) | de Vries ADNDT 36 (1987) Table I 2.54(5) / 2.56(5) / 2.57(10) fm; Angeli-Marinova 2.589(39) fm (closure, recorded); survey transcription | abs(r - r_i) &lt;= 2 sigma_i for each de Vries row (i.e. r in [2.46, 2.64] fm) | 0.06 s |
| [`t3_li_magnetization_rfy`](#t3_li_magnetization_rfy) | T3 | 6 | **BLOCKED** | r_mag(6Li) = 2.9469 fm (tree F_m slope; UNCOMPARED) | none on disk | n/a | 0.06 s |
| [`t3_nmc_li6_over_d`](#t3_nmc_li6_over_d) | T3 | 4 | **BLOCKED (ENVIRONMENT)** | — (*environment: LHAPDF set EPPS21nlo_CT18Anlo_Li6 not installed*) | — | — | 0.06 s |
| [`t4_djangoh_rad_noRad`](#t4_djangoh_rad_norad) | T4 | 7 | **BLOCKED** | t-peak: 1.01617 / 1.01317 / 1.01497 / 1.02069; polrad-full: 1.0271 / 1.01507 / 1.01514 / 1.02076 | 1.07245 / 1.08941 / 1.16063 / 1.22797 | none definable (the two RC models share no term) | 16.90 s |
| [`t4_pythia_ep_closure`](#t4_pythia_ep_closure) | T4 | — | **PASS** | max \|bridge/stock - 1\| = 0.023 (W2 n_ch_eta[-2,0)) over 15 compared sub-rows, 0 fail; &lt;n_ch&gt; W1 1.003; &lt;n_ch&gt; W2 0.998 | stock PYTHIA 8.312 e p WeakBosonExchange:ff2ff(t:gmZ), D-1 settings, e 10 x p 100, PDF:pSet 13, 20000 + 20000 events | \|bridge/stock - 1\| &lt;= 0.20 per observable per window, at matched (x, Q2) (D-1) | 10.13 s |
| [`t5_epios_source_modes`](#t5_epios_source_modes) | T5 | 1 | **PASS** | max moment round-trip residual 0.000e+00 over 8 modes (min population 0) | EPIOS Table II ideal (P_z, P_zz), 8 modes | 1e-12 | 0.07 s |
| [`t5_est_identity`](#t5_est_identity) | T5 | 9 | **PASS** | max\|spin_temperature_pzz(1,Pz) - closed form\| = 1.138e-14 at Pz = -0.5996 over 19999 interior points | P_zz = 2 - sqrt(4 - 3 P_z^2) (EST, analytic) | 2.30647e-14 | 0.10 s |
| [`t5_sum_rules`](#t5_sum_rules) | T5 | — | **FAIL** | bc_ww closure max\|resid\|/tol = 0.667 over 45 cells (PASS); Bjorken ToyG1 Q2=5: int_{1e-6}^1 (g1p-g1n) = 0.6430, last decade +0.2198 (divergent: g1n ~ x^-1.217) | BC: -x_min int g1/u du (exact truncation of int g2^WW); Bjorken: g_A/6 C_Bj(alpha_s) = 0.1854 +- 0.0023 at Q2=5 | BC: Peano trapezoid bound (h^2/8)L^2 int\|D^2 g1\| + quad errors; Bjorken: converged (last decade &lt;= 1 %) and within 10 % | 1.40 s |
| [`t5_weighted_unweighted`](#t5_weighted_unweighted) | T5 | — | **PASS** | unweighted sigma = 89.1934 +- 0.0512 pb (N = 3032577) vs weighted (Mode W) 89.1674 +- 0.0516 pb (N = 2987706, &lt;w&gt; = 1.01472): pull +0.359 | same physics point, the other route: 6Li inclusive config 1, x &gt;= 0.1, y &gt;= 0.5, helicity-flip apar+ (exact sigma_tot_pb 89.1671 pb; a no-op Mode W would pull +18.3) | \|pull\| &lt;= 3 | 2.25 s |
| [`t6_hepmc3_convention`](#t6_hepmc3_convention) | T6 | — | **FAIL** | 21/23 clauses hold on 8 files, 800 events (lipolgen-run, all CLI channels, T2 + T0 + --rc); failing: C13 (inclusive-6Li-T2, inclusive-7Li-T2, inclusive-6Li-T0, inclusive-6Li-T2-rc); C14 (inclusive-6Li-T2, inclusive-7Li-T2, inclusive-6Li-T2-rc) | docs/HEPMC3_CONVENTION.md sha256 5c6b55dcd9345542; 23 clauses | every clause on every event; numeric identities at the doc's 1e-9 (L137-138), four-momentum per component relative to max(1, summed beam E) | 3.69 s |

## Detail

### t1_g1d_world

- **status:** **FAIL**
- **vendored data (sha256 of every byte):** `data/nnpdf_COMPASS15_NC_NOTFIXED_MUD.json` `35ef52704ee46f4f7cd49803fc452e3a759a10fa91547f4378307e16b45c7892`; `data/nnpdf_E143_NC_NOTFIXED_ED.json` `5b5d7faefef49e3b8c3c17cdac3f277921799c2b3d42c721eb7366eccb623915`; `data/nnpdf_HERMES_NC_7GEV_ED.json` `5b397977b8937866a92bf02693d1e7103b1c871933be0b9df32e3e4f5e0cc0e8`; `data/nnpdf_SMC_NC_NOTFIXED_MUD.json` `8543a2919d8d9b6be0d31b77bf701679498499e054b876149655e01a0911c4fb`

| sub-row | status | values / reason |
|---|---|---|
| `toy` | **FAIL** | label = ToyG1 on ToyF2 (shipped default); q2_floor = 1; kernel_max_rel_diff = 0; world_chi2 = 232.16; world_ndf = 58; world_p = 1.0827e-22 |
| `toy/compass` | **FAIL** | all.n = 15; all.chi2 = 107.849; all.p = 4.16788e-16; all.chi2_cov = n/a; all.p_cov = n/a; all.worst_pull = -6.20059; all.worst_x = 0.00699; all.worst_q2 = 1.39; dis.n = 15; dis.chi2 = 107.849; dis.p = 4.16788e-16; dis.chi2_cov = n/a; dis.p_cov = n/a; dis.worst_pull = -6.20059; dis.worst_x = 0.00699; dis.worst_q2 = 1.39; dis_no_dstate.n = 15; dis_no_dstate.chi2 = 126.37; dis_no_dstate.p = 1.08961e-19; dis_no_dstate.chi2_cov = n/a; dis_no_dstate.p_cov = n/a; dis_no_dstate.worst_pull = -6.66364; dis_no_dstate.worst_x = 0.00699; dis_no_dstate.worst_q2 = 1.39; verdict_differs_under_full_cov = false |
| `toy/hermes` | **FAIL** | all.n = 15; all.chi2 = 69.3839; all.p = 5.74999e-09; all.chi2_cov = 104.46; all.p_cov = 1.85198e-15; all.worst_pull = -3.82896; all.worst_x = 0.3598; all.worst_q2 = 7.24; dis.n = 15; dis.chi2 = 69.3839; dis.p = 5.74999e-09; dis.chi2_cov = 104.46; dis.p_cov = 1.85198e-15; dis.worst_pull = -3.82896; dis.worst_x = 0.3598; dis.worst_q2 = 7.24; dis_no_dstate.n = 15; dis_no_dstate.chi2 = 71.1673; dis_no_dstate.p = 2.7653e-09; dis_no_dstate.chi2_cov = 111.031; dis_no_dstate.p_cov = 1.02198e-16; dis_no_dstate.worst_pull = 3.54256; dis_no_dstate.worst_x = 0.1059; dis_no_dstate.worst_q2 = 2.97; verdict_differs_under_full_cov = false |
| `toy/e143` | **FAIL** | all.n = 28; all.chi2 = 54.927; all.p = 0.00173472; all.chi2_cov = 51.5262; all.p_cov = 0.00434788; all.worst_pull = -2.58094; all.worst_x = 0.182; all.worst_q2 = 4.34; dis.n = 28; dis.chi2 = 54.927; dis.p = 0.00173472; dis.chi2_cov = 51.5262; dis.p_cov = 0.00434788; dis.worst_pull = -2.58094; dis.worst_x = 0.182; dis.worst_q2 = 4.34; dis_no_dstate.n = 28; dis_no_dstate.chi2 = 55.5093; dis_no_dstate.p = 0.00147565; dis_no_dstate.chi2_cov = 54.6984; dis_no_dstate.p_cov = 0.00184783; dis_no_dstate.worst_pull = 2.35372; dis_no_dstate.worst_x = 0.049; dis_no_dstate.worst_q2 = 1.78; verdict_differs_under_full_cov = false |
| `toy/smc` | **PASS** (recorded, not in the verdict) | all.n = 13; all.chi2 = 46.1559; all.p = 1.34037e-05; all.chi2_cov = n/a; all.p_cov = n/a; all.worst_pull = -5.4317; all.worst_x = 0.002; all.worst_q2 = 0.5; dis.n = 12; dis.chi2 = 16.6525; dis.p = 0.163143; dis.chi2_cov = n/a; dis.p_cov = n/a; dis.worst_pull = -1.83408; dis.worst_x = 0.049; dis.worst_q2 = 10; dis_no_dstate.n = 12; dis_no_dstate.chi2 = 16.5038; dis_no_dstate.p = 0.169237; dis_no_dstate.chi2_cov = n/a; dis_no_dstate.p_cov = n/a; dis_no_dstate.worst_pull = 1.90061; dis_no_dstate.worst_x = 0.077; dis_no_dstate.worst_q2 = 14.4; verdict_differs_under_full_cov = false |
| `mstw` | **FAIL** | label = ToyG1 on MstwSF (MSTW2008 LO, PYTHIA 8.312 pdfdata; the --unpol-sf mstw g1); q2_floor = 1; kernel_max_rel_diff = 0; world_chi2 = 150.615; world_ndf = 58; world_p = 3.61869e-10 |
| `mstw/compass` | **FAIL** | all.n = 15; all.chi2 = 53.237; all.p = 3.52132e-06; all.chi2_cov = n/a; all.p_cov = n/a; all.worst_pull = -4.26403; all.worst_x = 0.00699; all.worst_q2 = 1.39; dis.n = 15; dis.chi2 = 53.237; dis.p = 3.52132e-06; dis.chi2_cov = n/a; dis.p_cov = n/a; dis.worst_pull = -4.26403; dis.worst_x = 0.00699; dis.worst_q2 = 1.39; dis_no_dstate.n = 15; dis_no_dstate.chi2 = 66.2986; dis_no_dstate.p = 2.02092e-08; dis_no_dstate.chi2_cov = n/a; dis_no_dstate.p_cov = n/a; dis_no_dstate.worst_pull = -4.5869; dis_no_dstate.worst_x = 0.00699; dis_no_dstate.worst_q2 = 1.39; verdict_differs_under_full_cov = false |
| `mstw/hermes` | **FAIL** | all.n = 15; all.chi2 = 47.1342; all.p = 3.50414e-05; all.chi2_cov = 68.8271; all.p_cov = 7.22068e-09; all.worst_pull = 3.37598; all.worst_x = 0.1059; all.worst_q2 = 2.97; dis.n = 15; dis.chi2 = 47.1342; dis.p = 3.50414e-05; dis.chi2_cov = 68.8271; dis.p_cov = 7.22068e-09; dis.worst_pull = 3.37598; dis.worst_x = 0.1059; dis.worst_q2 = 2.97; dis_no_dstate.n = 15; dis_no_dstate.chi2 = 64.0401; dis_no_dstate.p = 5.02998e-08; dis_no_dstate.chi2_cov = 94.6119; dis_no_dstate.p_cov = 1.35733e-13; dis_no_dstate.worst_pull = 3.92619; dis_no_dstate.worst_x = 0.1059; dis_no_dstate.worst_q2 = 2.97; verdict_differs_under_full_cov = false |
| `mstw/e143` | **FAIL** | all.n = 28; all.chi2 = 50.2437; all.p = 0.0060742; all.chi2_cov = 51.3395; all.p_cov = 0.00456668; all.worst_pull = 2.82365; all.worst_x = 0.259; all.worst_q2 = 5.26; dis.n = 28; dis.chi2 = 50.2437; dis.p = 0.0060742; dis.chi2_cov = 51.3395; dis.p_cov = 0.00456668; dis.worst_pull = 2.82365; dis.worst_x = 0.259; dis.worst_q2 = 5.26; dis_no_dstate.n = 28; dis_no_dstate.chi2 = 61.4966; dis_no_dstate.p = 0.000261336; dis_no_dstate.chi2_cov = 60.9935; dis_no_dstate.p_cov = 0.000303612; dis_no_dstate.worst_pull = 3.22924; dis_no_dstate.worst_x = 0.259; dis_no_dstate.worst_q2 = 5.26; verdict_differs_under_full_cov = false |
| `mstw/smc` | **PASS** (recorded, not in the verdict) | all.n = 13; all.chi2 = 20.6699; all.p = 0.0797187; all.chi2_cov = n/a; all.p_cov = n/a; all.worst_pull = -2.5037; all.worst_x = 0.002; all.worst_q2 = 0.5; dis.n = 12; dis.chi2 = 14.4014; dis.p = 0.275812; dis.chi2_cov = n/a; dis.p_cov = n/a; dis.worst_pull = -1.92466; dis.worst_x = 0.049; dis.worst_q2 = 10; dis_no_dstate.n = 12; dis_no_dstate.chi2 = 14.7566; dis_no_dstate.p = 0.255017; dis_no_dstate.chi2_cov = n/a; dis_no_dstate.p_cov = n/a; dis_no_dstate.worst_pull = 1.84943; dis_no_dstate.worst_x = 0.077; dis_no_dstate.worst_q2 = 14.4; verdict_differs_under_full_cov = false |
| `ct18nlo` | **BLOCKED** | **reason:** environment: LHAPDF set CT18NLO not installed; label = ToyG1(LhapdfSF(CT18NLO)) |
| `nnpdfpol` | **BLOCKED** | **reason:** environment: LHAPDF set NNPDFpol11_100 not installed; label = LhapdfG1(NNPDFpol11_100) |

### t1_nmc_f2d

- **status:** **FAIL**
- **vendored data (sha256 of every byte):** `data/nnpdf_NMC_NC_NOTFIXED_D.json` `d38b3c0c06ca4e123983cffd2f647d82d14f272e8ee04e13f73f9fd098af8e3d`

| sub-row | status | values / reason |
|---|---|---|
| `toy` | **FAIL** | label = ToyF2 (shipped default); q2_floor = 1; all.n = 158; all.chi2 = 9582.57; all.p = 0; all.chi2_cov = 24996; all.p_cov = 0; all.mean_ratio = 1.07418; all.worst_pull = 26.0729; all.worst_x = 0.008; all.worst_q2 = 0.75; dis.n = 155; dis.chi2 = 8307.9; dis.p = 0; dis.chi2_cov = 21100.7; dis.p_cov = 0; dis.mean_ratio = 1.06018; dis.worst_pull = 21.4617; dis.worst_x = 0.0125; dis.worst_q2 = 1.25; kernel_max_rel_diff = 0 |
| `mstw` | **FAIL** | label = MstwSF (MSTW2008 LO, PYTHIA 8.312 pdfdata, i_fit 3); q2_floor = 1; all.n = 158; all.chi2 = 485.451; all.p = 5.46777e-35; all.chi2_cov = 1090.25; all.p_cov = 5.18242e-139; all.mean_ratio = 1.02876; all.worst_pull = 5.93528; all.worst_x = 0.008; all.worst_q2 = 0.75; dis.n = 155; dis.chi2 = 416.459; dis.p = 8.21357e-26; dis.chi2_cov = 858.536; dis.p_cov = 6.94928e-98; dis.mean_ratio = 1.02545; dis.worst_pull = 4.96173; dis.worst_x = 0.0125; dis.worst_q2 = 1.25; kernel_max_rel_diff = 0 |
| `ct18nlo` | **BLOCKED** | **reason:** environment: LHAPDF set CT18NLO not installed; label = LhapdfSF(CT18NLO) |

### t1_nmc_f2d_over_f2p

- **status:** **FAIL**
- **vendored data (sha256 of every byte):** `data/nnpdf_NMC_NC_NOTFIXED.json` `47e833cfb399eb984455e09ed37d8833edd20445607cefe0a103545cd8e43900`

| sub-row | status | values / reason |
|---|---|---|
| `toy` | **FAIL** | label = ToyF2 (shipped default); q2_floor = 1; all.n = 260; all.chi2 = 2251.99; all.p = 0; all.chi2_cov = 1464.97; all.p_cov = 6.80022e-167; all.mean_diff = 0.0351552; all.worst_pull = 7.42661; all.worst_x = 0.35; all.worst_q2 = 14.85; all.chi2_legacy_diag = 2245.92; all.chi2_legacy_cov = 990.572; all.p_legacy_cov = 9.33888e-86; dis.n = 211; dis.chi2 = 2161.44; dis.p = 0; dis.chi2_cov = 1154.29; dis.p_cov = 9.28003e-130; dis.mean_diff = 0.0391701; dis.worst_pull = 7.42661; dis.worst_x = 0.35; dis.worst_q2 = 14.85; dis.chi2_legacy_diag = 2155.29; dis.chi2_legacy_cov = 950.707; dis.p_legacy_cov = 2.45165e-94; all_frozen.n = 260; all_frozen.chi2 = 2251.99; all_frozen.p = 0; all_frozen.chi2_cov = 1464.97; all_frozen.p_cov = 6.80022e-167; all_frozen.mean_diff = 0.0351552; all_frozen.worst_pull = 7.42661; all_frozen.worst_x = 0.35; all_frozen.worst_q2 = 14.85; all_frozen.chi2_legacy_diag = 2245.92; all_frozen.chi2_legacy_cov = 990.572; all_frozen.p_legacy_cov = 9.33888e-86; dis_frozen.n = 211; dis_frozen.chi2 = 2161.44; dis_frozen.p = 0; dis_frozen.chi2_cov = 1154.29; dis_frozen.p_cov = 9.28003e-130; dis_frozen.mean_diff = 0.0391701; dis_frozen.worst_pull = 7.42661; dis_frozen.worst_x = 0.35; dis_frozen.worst_q2 = 14.85; dis_frozen.chi2_legacy_diag = 2155.29; dis_frozen.chi2_legacy_cov = 950.707; dis_frozen.p_legacy_cov = 2.45165e-94; max_abs_frozen_minus_q2 = 1.11022e-16; kernel_max_rel_diff = 0 |
| `mstw` | **PASS** | label = MstwSF (MSTW2008 LO, PYTHIA 8.312 pdfdata, i_fit 3); q2_floor = 1; all.n = 260; all.chi2 = 297.896; all.p = 0.0530191; all.chi2_cov = 305.145; all.p_cov = 0.0284317; all.mean_diff = 0.0016849; all.worst_pull = 3.06959; all.worst_x = 0.14; all.worst_q2 = 46.95; all.chi2_legacy_diag = 297.213; all.chi2_legacy_cov = 269.588; all.p_legacy_cov = 0.328324; dis.n = 211; dis.chi2 = 247.179; dis.p = 0.0444218; dis.chi2_cov = 251.635; dis.p_cov = 0.0290014; dis.mean_diff = 0.000627172; dis.worst_pull = 3.06959; dis.worst_x = 0.14; dis.worst_q2 = 46.95; dis.chi2_legacy_diag = 246.498; dis.chi2_legacy_cov = 224.677; dis.p_legacy_cov = 0.246927; all_frozen.n = 260; all_frozen.chi2 = 310.178; all_frozen.p = 0.0177867; all_frozen.chi2_cov = 317.801; all_frozen.p_cov = 0.00826963; all_frozen.mean_diff = 0.00167267; all_frozen.worst_pull = 3.42309; all_frozen.worst_x = 0.14; all_frozen.worst_q2 = 46.95; all_frozen.chi2_legacy_diag = 309.591; all_frozen.chi2_legacy_cov = 282.396; all_frozen.p_legacy_cov = 0.162479; dis_frozen.n = 211; dis_frozen.chi2 = 260.702; dis_frozen.p = 0.0111915; dis_frozen.chi2_cov = 269.079; dis_frozen.p_cov = 0.00421108; dis_frozen.mean_diff = 0.000400036; dis_frozen.worst_pull = 3.42309; dis_frozen.worst_x = 0.14; dis_frozen.worst_q2 = 46.95; dis_frozen.chi2_legacy_diag = 260.108; dis_frozen.chi2_legacy_cov = 239.138; dis_frozen.p_legacy_cov = 0.0892618; max_abs_frozen_minus_q2 = 0.0158072; kernel_max_rel_diff = 0 |
| `ct18nlo` | **BLOCKED** | **reason:** environment: LHAPDF set CT18NLO not installed; label = LhapdfSF(CT18NLO) |

### t2_bonus_spectator_shape

- **status:** **BLOCKED**
- **reason:** BLOCKED: the CLAS Physics Database (queried 2026-09-23) holds BONuS (E-03-012, eid 135) only as ONE measurement, E135M1 = F2n/F2p vs W\* (Baillie et al., PRL 108 (2012) 142001) -- no spectator-momentum spectrum, no Tkachenko 2014 entry; EG1b (eid 95, 146) is inclusive g1/A1 only. To unblock, download the Supplemental Material (tables of numerical results, their ref. [71]) of S. Tkachenko et al. (CLAS BONuS), Phys. Rev. C 89 (2014) 045206, doi:10.1103/PhysRevC.89.045206, https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.89.045206 (HTTP 401 without an APS login on 2026-09-23; arXiv:1402.2477 has no ancillary files), into validation/benchmarks/data/. Found instead and vendored unwired: Deeps F2N x P(p_s) at 0.30-0.53 GeV/c (eid 90, Klimenko et al., PRC 73 (2006) 035212).
- **vendored data (sha256 of every byte):** `data/clas_db_deuteron_query_2026-09-23.json` `340e04a5ae8ada1eabf4d90c2c76c82b07bda48d28c7f43cdffd881c8ee3fa11`; `data/clas_deeps_klimenko2006_f2P.json` `1d98fd11e1b11b1e070e5504f4919b98b87dbabfe29cdca347522ea4f8737256`

### t2_deuteron_static

- **status:** **FAIL**
- **vendored data (sha256 of every byte):** `data/deuteron_static_measured.json` `6433c31e3c837bd2a88c71417eca5fa4eba51f17495486dd2c21bc99d8f692d5`
- **config:** lipolgen = 0.1.0; numpy = 2.4.6; python = 3.11.16; data_file_sha256 = 6433c31e3c837bd2a88c71417eca5fa4eba51f17495486dd2c21bc99d8f692d5

| sub-row | status | values / reason |
|---|---|---|
| `cdbonn` | **PASS** | label = CD-Bonn closed form (lg.cdbonn_wave(), Machleidt PRC 63 (2001) 024001 App. D); eta = 0.0255714; eta_source = tree: CdBonnWave.eta() = D_1/C_1; eta_numerical_error = 0; eta_check_r60 = 0.0255714; q_d = 0.270493; q_d_quad_error = 7.38298e-15; p_d = 0.0485621; norm = 1; p_d_tree = 0.0485621; a_s = 0.88473; kappa_fm = 0.231538; cross_checks.q_tree_default = 0.270178; cross_checks.q_tree_default_note = lg.alpha_d_quadrupole_fm2 on cdbonn_fdeut_table(), r_max 20 fm, dk 0.1 fm^-1 (tests/test_b1_nuclear.cpp item 5); cross_checks.q_tree_converged = 0.270494; cross_checks.q_tree_converged_note = the same routine on dk = 0.02 fm^-1, r_max = 60 fm, n_r = 6000; config.r_max_fm = 80; config.n_intervals = 100000; eta_pull = -0.0715534; eta_status = pass; q_d_distance = -0.0154066; q_d_rel_distance = -0.0538881 |
| `av18` | **PASS** | label = AV18 table (vmc/deuteron/fdeut.av18, r-space block; the tree reads its k block); eta = 0.0250447; eta_source = this harness: w/(u h(kappa r)) averaged over r in [20, 40] fm of the file's r block; eta_numerical_error = 2.94678e-08; kappa_fm = 0.231607; a_s = 0.885057; q_d = 0.26967; q_d_quad_error = 9.92824e-10; p_d = 0.0575986; norm = 0.999998; p_d_tree = 0.0575999; header.ebind = 2.22457; header.dstate = 0.057599; header.qm = 0.269673; header.mu = 0.84699; header.as = 0.885056; header.eta = 0.025045; header.rd = 1.96736; header.r4 = 55.3761; cross_checks.header_eta = 0.025045; cross_checks.header_qm = 0.269673; cross_checks.header_dstate = 0.057599; cross_checks.header_as = 0.885056; cross_checks.q_tree_kblock = 0.269362; cross_checks.q_tree_kblock_note = lg.alpha_d_quadrupole_fm2 on the k block the generator reads (deuteron_channel(source=VmcAV18)), default numerics (test…; cross_checks.eta_kblock_r = 0.0249839; cross_checks.eta_kblock_note = w/(u h(kappa r)) at r = 20 fm from the k block through ClusterPartialWave's spline and a 20001-point bare-j_L transform; config.file = /home/user/PoliGen/data/vmc/deuteron/fdeut.av18; config.eta_window_fm = 20 / 40; eta_pull = -1.38832; eta_status = pass; q_d_distance = -0.0162297; q_d_rel_distance = -0.0567672 |
| `hulthen` | **FAIL** | label = Hulthen pair (deuteron_channel() default: beta = 0.30 GeV, P_D = 0.045, kappa = 0.045702 GeV), TaggedModel normalisation; eta = 0.0324784; eta_source = this harness: -N_2 A/N_0, the k^2 = -kappa^2 pole residues, N_L read from TaggedModel; eta_numerical_error = 6.32532e-12; eta_check_r80 = 0.0324784; kappa_fm = 0.231604; beta_gev = 0.3; p_d_nominal = 0.045; n0 = 0.298545; n2 = 0.0358768; n0_continuum_over_tree = 0.99997; q_d = 0.290349; q_d_quad_error = 1.22347e-13; p_d = 0.0450804; norm = 1.00014; config.beta_gev = 0.3; config.p_d = 0.045; config.k_grid_gev = 0.0001 / 1.2 / 280; config.r_max_fm = 80; config.n_intervals = 100000; cross_checks.q_tree_tables = 0.289791; cross_checks.q_tree_tables_note = lg.alpha_d_quadrupole_fm2 on the TaggedModel tables (k &lt;= 1.2 GeV), default numerics; eta_pull = 17.1961; eta_status = fail; q_d_distance = 0.00444905; q_d_rel_distance = 0.0155616 |

### t2_hermes_b1_table2

- **status:** **FAIL**
- **vendored data (sha256 of every byte):** `data/hermes2005_b1d_table2.json` `20b9104edf52414822e88305f7144566568ee2979bd531f77dd79f3a92006fc5`; `data/whitlow1990_r1990.json` `ce330eeef6305315a78d66422ff81fe8307bf71fdca4b8ef6e781e71698591d9`

| sub-row | status | values / reason |
|---|---|---|
| `cdks_conv_mstw` | **FAIL** | chi2_b1 = 22.2637; within_b1 = 2; chi2_azz = 22.0087; within_azz = 2; chi2_azz_r1998 = 22.0079; ndf = 6 |
| `cdks_conv_default` | **FAIL** | chi2_b1 = 22.6692; within_b1 = 2; chi2_azz = 22.3919; within_azz = 2; chi2_azz_r1998 = 22.3926; ndf = 6 |
| `miller_shipped` | **FAIL** | chi2_b1 = 5.26186; within_b1 = 4; chi2_azz = 5.54993; within_azz = 4; chi2_azz_r1998 = 5.4104; ndf = 6 |
| `miller_x1` | **FAIL** | chi2_b1 = 5.29707; within_b1 = 5; chi2_azz = 4.5139; within_azz = 5; chi2_azz_r1998 = 4.58653; ndf = 6 |
| `cdks_digitized` | **FAIL** | chi2_b1 = 21.7817; within_b1 = 2; chi2_azz = 21.559; within_azz = 2; chi2_azz_r1998 = 21.5484; ndf = 6 |
| `null_b1_zero` | **FAIL** | chi2_b1 = 21.7994; within_b1 = 2; chi2_azz = 21.5857; within_azz = 2; chi2_azz_r1998 = 21.5857; ndf = 6 |

### t3_li6_charge_ff_fb

- **status:** **FAIL**
- **vendored data (sha256 of every byte):** `data/uva_ncd_fb_data_li6.dat` `6d5174b80f1a0da2197354feb244307ede17d38a1a9169e546abc18e78da7160`

### t3_li6_rms_devries

- **status:** **PASS**
- **vendored data (sha256 of every byte):** `data/devries1987_li6_rms.json` `ffab49a1b67d92a209de3341de994c6b1fc115fe7f69044ef054d3154a53ad3c`
- **config:** lipolgen = 0.1.0; python = 3.11.16; data_file_sha256 = ffab49a1b67d92a209de3341de994c6b1fc115fe7f69044ef054d3154a53ad3c

| sub-row | status | values / reason |
|---|---|---|
| `Su67` | **PASS** | compilation = de Vries, de Jager, de Vries, ADNDT 36 (1987) 495, Table I; model = MI30; rms_fm = 2.54; sigma_fm = 0.05; as_printed = 2.54(5); q_range_fm-1 = 0.69 / 2.52; source_key = Su67; source_key_survey_expansion = Suelzle 1967; transcribed_from = docs/benchmarking/03_data_nuclear.md:186; pull = 0.620016; within = true; within_1sigma = true |
| `Li71a` | **PASS** | compilation = de Vries, de Jager, de Vries, ADNDT 36 (1987) 495, Table I; model = MI30; rms_fm = 2.56; sigma_fm = 0.05; as_printed = 2.56(5); q_range_fm-1 = 0.56 / 3.66; source_key = Li71a; source_key_survey_expansion = Li-Sick-Whitney-Yearian; transcribed_from = docs/benchmarking/03_data_nuclear.md:187; pull = 0.220016; within = true; within_1sigma = true |
| `Bu72` | **PASS** | compilation = de Vries, de Jager, de Vries, ADNDT 36 (1987) 495, Table I; model = MI; rms_fm = 2.57; sigma_fm = 0.1; as_printed = 2.57(10); q_range_fm-1 = 0.09 / 0.9; source_key = Bu72; source_key_survey_expansion = n/a; transcribed_from = docs/benchmarking/03_data_nuclear.md:188; pull = 0.010008; within = true; within_1sigma = true |

### t3_li_magnetization_rfy

- **status:** **BLOCKED**
- **reason:** BLOCKED: reference not on disk and not open access. Download R. E. Rand, R. F. Frosch, M. R. Yearian, Phys. Rev. 144 (1966) 859, https://doi.org/10.1103/PhysRev.144.859 OR C. W. de Jager, H. de Vries, C. de Vries, At. Data Nucl. Data Tables 14 (1974) 479 (Table V), https://doi.org/10.1016/S0092-640X(74)80002-1 (docs/references/REFERENCES.md sec. 1 items 8 and 9), then vendor the lithium rows as validation/benchmarks/data/rfy1966_li_magnetization.csv
- **vendored data (sha256 of every byte):** `data/rfy1966_li_magnetization.csv` **not on disk**

### t3_nmc_li6_over_d

- **status:** **BLOCKED (ENVIRONMENT)**
- **reason:** environment: LHAPDF set EPPS21nlo_CT18Anlo_Li6 not installed
- **vendored data (sha256 of every byte):** `data/hepdata_ins394050_table1.csv` `1f2c36b02972c13482a252c0a13bea734496439db92c60829c337d0a964e254b`

### t4_djangoh_rad_noRad

- **status:** **BLOCKED**
- **reason:** no shared O(alpha) term: DJANGOH Rad=1 ran with IEL2=IEL31=IEL32=IEL33=0 (elastic radiative tail OFF; all eight ePIC logs), LiPolGen's rc_tail is the elastic(+QE) tail ONLY and applies no inelastic/vertex shift -- disjoint additive pieces, so must-not-contradict holds vacuously; the samples' W_h &gt;= 3 GeV cut also excludes the elastic tail (W_h = M_p) whatever IEL is (DJANGOH zeroes IEL2/IEL31..33 when WMIN &gt; M_p), so it unblocks only with a DJANGOH run at IEL31..33 != 0 AND WMIN &lt;= M_p
- **vendored data (sha256 of every byte):** `data/djangoh_rad_norad_ep18x275.json` `2fd4b9b6bb2c0de7c7a9b338ad58cc3cbc221c768d0ce2f257ac6f57139be8fa`
- **config:** beams = e 18 x p 275 (config 2); x_min = 1e-05; y_min = 0.0001; y_max = 0.985; W2_min = 9; rc_ff = zero coherent vertex; qe_kf_gev = 0; clipped_node_fraction.t-peak = 0; clipped_node_fraction.polrad-full = 0.00424328; unpol_sf.t-peak = UnpolSfSource.Toy; unpol_sf.polrad-full = UnpolSfSource.Toy

### t4_pythia_ep_closure

- **status:** **PASS**; expensive (opt-in)
- **config:** expensive = true; n_events = 20000; seed = 20260926; beam_config = 1; pythia_version_used = 8.312; pythia_version_pinned = 8.317; pythia_libraries = /opt/deps/lib/libpythia8.so; pdf_pset = 13; n_sub = 4; eta_edges = -inf / -2 / 0 / 2 / 4 / inf; min_stock_counts = 100; stock_settings_added = SpaceShower:QEDshowerByQ = off

| sub-row | status | values / reason |
|---|---|---|
| `W1/n_ch` | **PASS** | window = W1; observable = n_ch; stock = 8.0267; stock_err = 0.0227026; bridge = 8.05391; bridge_err = 0.0194701; ratio = 1.00339; ratio_err = 0.00373336; distance = 0.00339012; bridge_raw = 8.1546; stock_raw = 8.0267; ratio_raw = 1.01593; n_stock_particles = n/a; compared = true |
| `W1/sum_empz` | **PASS** | window = W1; observable = sum_empz; stock = 2.7229; stock_err = 0.0193708; bridge = 2.74055; bridge_err = 0.00637262; ratio = 1.00648; ratio_err = 0.0075329; distance = 0.00647918; bridge_raw = 2.90458; stock_raw = 2.7229; ratio_raw = 1.06672; n_stock_particles = n/a; compared = true |
| `W1/sum_pt_ch` | **PASS** | window = W1; observable = sum_pt_ch; stock = 3.82709; stock_err = 0.0110632; bridge = 3.84917; bridge_err = 0.00983067; ratio = 1.00577; ratio_err = 0.00387963; distance = 0.00576844; bridge_raw = 3.87739; stock_raw = 3.82709; ratio_raw = 1.01314; n_stock_particles = n/a; compared = true |
| `W1/n_ch_eta[-inf,-2)` | **PASS** | window = W1; observable = n_ch_eta[-inf,-2); stock = 0.0286; stock_err = 0.00136536; bridge = 0.028317; bridge_err = 0.00115488; ratio = 0.990104; ratio_err = 0.0621672; distance = 0.00989624; bridge_raw = 0.0325; stock_raw = 0.0286; ratio_raw = 1.13636; n_stock_particles = 572; compared = true |
| `W1/n_ch_eta[-2,0)` | **PASS** | window = W1; observable = n_ch_eta[-2,0); stock = 0.6916; stock_err = 0.00842306; bridge = 0.680522; bridge_err = 0.00538898; ratio = 0.983983; ratio_err = 0.0142945; distance = 0.0160173; bridge_raw = 0.7344; stock_raw = 0.6916; ratio_raw = 1.06189; n_stock_particles = 13832; compared = true |
| `W1/n_ch_eta[0,2)` | **PASS** | window = W1; observable = n_ch_eta[0,2); stock = 2.4579; stock_err = 0.0127797; bridge = 2.46123; bridge_err = 0.0113746; ratio = 1.00135; ratio_err = 0.00696589; distance = 0.00135466; bridge_raw = 2.49415; stock_raw = 2.4579; ratio_raw = 1.01475; n_stock_particles = 49158; compared = true |
| `W1/n_ch_eta[2,4)` | **PASS** | window = W1; observable = n_ch_eta[2,4); stock = 2.87375; stock_err = 0.0128736; bridge = 2.91045; bridge_err = 0.0128979; ratio = 1.01277; ratio_err = 0.0063818; distance = 0.0127702; bridge_raw = 2.91315; stock_raw = 2.87375; ratio_raw = 1.01371; n_stock_particles = 57475; compared = true |
| `W1/n_ch_eta[4,+inf)` | **PASS** | window = W1; observable = n_ch_eta[4,+inf); stock = 1.97485; stock_err = 0.00830989; bridge = 1.97339; bridge_err = 0.00836494; ratio = 0.999263; ratio_err = 0.00596837; distance = 0.000737182; bridge_raw = 1.9804; stock_raw = 1.97485; ratio_raw = 1.00281; n_stock_particles = 39497; compared = true |
| `W2/n_ch` | **PASS** | window = W2; observable = n_ch; stock = 7.2606; stock_err = 0.0205277; bridge = 7.24281; bridge_err = 0.0182964; ratio = 0.99755; ratio_err = 0.00378214; distance = 0.00244989; bridge_raw = 7.3014; stock_raw = 7.2606; ratio_raw = 1.00562; n_stock_particles = n/a; compared = true |
| `W2/sum_empz` | **PASS** | window = W2; observable = sum_empz; stock = 2.00117; stock_err = 0.0121968; bridge = 2.006; bridge_err = 0.00364611; ratio = 1.00242; ratio_err = 0.00637546; distance = 0.002417; bridge_raw = 2.05983; stock_raw = 2.00117; ratio_raw = 1.02932; n_stock_particles = n/a; compared = true |
| `W2/sum_pt_ch` | **PASS** | window = W2; observable = sum_pt_ch; stock = 6.1278; stock_err = 0.0186538; bridge = 6.13132; bridge_err = 0.0158241; ratio = 1.00057; ratio_err = 0.00399322; distance = 0.00057363; bridge_raw = 6.17419; stock_raw = 6.1278; ratio_raw = 1.00757; n_stock_particles = n/a; compared = true |
| `W2/n_ch_eta[-inf,-2)` | **NOT COMPARED** | window = W2; observable = n_ch_eta[-inf,-2); stock = 0.00015; stock_err = 8.65982e-05; bridge = 0.000238271; bridge_err = 0.000106293; ratio = 1.58848; ratio_err = 1.15894; distance = 0.588476; bridge_raw = 0.00025; stock_raw = 0.00015; ratio_raw = 1.66667; n_stock_particles = 3; compared = false |
| `W2/n_ch_eta[-2,0)` | **PASS** | window = W2; observable = n_ch_eta[-2,0); stock = 0.0381; stock_err = 0.00205247; bridge = 0.0372096; bridge_err = 0.00164641; ratio = 0.97663; ratio_err = 0.0680833; distance = 0.0233704; bridge_raw = 0.03925; stock_raw = 0.0381; ratio_raw = 1.03018; n_stock_particles = 762; compared = true |
| `W2/n_ch_eta[0,2)` | **PASS** | window = W2; observable = n_ch_eta[0,2); stock = 2.5015; stock_err = 0.0151981; bridge = 2.49662; bridge_err = 0.0105711; ratio = 0.998048; ratio_err = 0.00739101; distance = 0.00195212; bridge_raw = 2.56265; stock_raw = 2.5015; ratio_raw = 1.02445; n_stock_particles = 50030; compared = true |
| `W2/n_ch_eta[2,4)` | **PASS** | window = W2; observable = n_ch_eta[2,4); stock = 3.229; stock_err = 0.0135771; bridge = 3.2194; bridge_err = 0.0133563; ratio = 0.997027; ratio_err = 0.00588933; distance = 0.00297341; bridge_raw = 3.2034; stock_raw = 3.229; ratio_raw = 0.992072; n_stock_particles = 64580; compared = true |
| `W2/n_ch_eta[4,+inf)` | **PASS** | window = W2; observable = n_ch_eta[4,+inf); stock = 1.49185; stock_err = 0.00757757; bridge = 1.48935; bridge_err = 0.00745836; ratio = 0.998323; ratio_err = 0.00712089; distance = 0.00167653; bridge_raw = 1.49585; stock_raw = 1.49185; ratio_raw = 1.00268; n_stock_particles = 29837; compared = true |

### t5_epios_source_modes

- **status:** **PASS**
- **vendored data (sha256 of every byte):** `data/epios2026_table2_deuteron_source_modes.json` `8492859b8fe98a61af5a00333da6ef151e98883d934fecaf023c7009dbd1320a`

### t5_est_identity

- **status:** **PASS**
- **vendored data (sha256 of every byte):** `data/est_identity_sources.json` `32d6f70d57a3b5bc816229ded9d9f6dbc3412598e70b4978c28c61c12ef8ac6b`

### t5_sum_rules

- **status:** **FAIL**

| sub-row | status | values / reason |
|---|---|---|
| `bc_ww` | **PASS** | n_cells = 45; n_within = 45; max_ratio = 0.666801; worst_cell.target = p; worst_cell.q2 = 2; worst_cell.x_min = 0.001 |
| `bjorken[toyg1]` | **FAIL** | per_q2.2.0.reference.q2 = 2; per_q2.2.0.reference.nf = 4; per_q2.2.0.reference.alpha_s = 0.358817; per_q2.2.0.reference.c_bj = 0.822753; per_q2.2.0.reference.gamma = 0.174876; per_q2.2.0.reference.err = 0.00462445; per_q2.2.0.reference.err_g_a = 0.000178263; per_q2.2.0.reference.err_alpha_s = 0.00145436; per_q2.2.0.reference.err_truncation = 0.00438618; per_q2.2.0.verdict.last_decade_step = 0.128925; per_q2.2.0.verdict.distance = 0.283232; per_q2.2.0.verdict.distance_rel = 1.61961; per_q2.2.0.verdict.converged = false; per_q2.2.0.verdict.within = false; per_q2.2.0.verdict.ok = false; per_q2.5.0.reference.q2 = 5; per_q2.5.0.reference.nf = 4; per_q2.5.0.reference.alpha_s = 0.285006; per_q2.5.0.reference.c_bj = 0.87219; per_q2.5.0.reference.gamma = 0.185384; per_q2.5.0.reference.err = 0.00232788; per_q2.5.0.reference.err_g_a = 0.000188975; per_q2.5.0.reference.err_alpha_s = 0.00074298; per_q2.5.0.reference.err_truncation = 0.00219802; per_q2.5.0.verdict.last_decade_step = 0.219819; per_q2.5.0.verdict.distance = 0.457638; per_q2.5.0.verdict.distance_rel = 2.46859; per_q2.5.0.verdict.converged = false; per_q2.5.0.verdict.within = false; per_q2.5.0.verdict.ok = false; per_q2.10.0.reference.q2 = 10; per_q2.10.0.reference.nf = 4; per_q2.10.0.reference.alpha_s = 0.247637; per_q2.10.0.reference.c_bj = 0.894197; per_q2.10.0.reference.gamma = 0.190062; per_q2.10.0.reference.err = 0.00153939; per_q2.10.0.reference.err_g_a = 0.000193743; per_q2.10.0.reference.err_alpha_s = 0.000503248; per_q2.10.0.reference.err_truncation = 0.00144184; per_q2.10.0.verdict.last_decade_step = 0.331185; per_q2.10.0.verdict.distance = 0.6651; per_q2.10.0.verdict.distance_rel = 3.49939; per_q2.10.0.verdict.converged = false; per_q2.10.0.verdict.within = false; per_q2.10.0.verdict.ok = false |
| `bjorken[nnpdfpol11]` | **BLOCKED** | **reason:** environment: LHAPDF set NNPDFpol11_100 not installed |
| `fsi_unitarity` | **BLOCKED** | **reason:** not definable in the tree as built (no GlauberFsiWeight option realises int dGamma S[FSI + FSI^2] = 0: the default is absorptive by design, elastic_gain would be a tune, weight_normalised is a division, the table stops at k_max); tree_side.channel = li6_alpha_channel() defaults (Hulthen beta 0.3); tree_side.sigma_xn_mb = 40; tree_side.survival = 0.520239; tree_side.int_w_minus_1 = -0.479761; tree_side.sigma_tot_cluster_mb = 131.045; tree_side.sigma_el_cluster_mb = 35.2237; tree_side.elastic_gain = 1; tree_side.k_max_gev = 1.2; tree_side.w_max = 50; tree_side.clipped_grid_fraction = 0 |

### t5_weighted_unweighted

- **status:** **PASS**

### t6_hepmc3_convention

- **status:** **FAIL**
- **reason:** C13 (doc lines 123-138): inclusive-6Li-T2: event 0: residual 0.82 (status 1) / 0.657 (status 1 + X) of the beam energy
- **config:** n_events = 100; seed = 20260926; tol = 1e-09; lipolgen = 0.1.0; pyhepmc = 2.16.1; pythia_charge_table = 8.312; hepmc3_file_version = 3.03.01

| sub-row | status | values / reason |
|---|---|---|
| `C01` | **PASS** | clause = C01; doc_lines = 16-19; anchor_line = 16; anchor_found = true; what = Asciiv3 only: HepMC::Version / Asciiv3-START header, END trailer, every requested event readable; n_checked = 8; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C02` | **PASS** | clause = C02; doc_lines = 28-29; anchor_line = 28; anchor_found = true; what = every event in GEV and MM; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C03` | **PASS** | clause = C03; doc_lines = 33-37; anchor_line = 33; anchor_found = true; what = exactly two status-4 particles, e (11) and the run's ion code, both incoming at ONE vertex that has no other incoming p…; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C04` | **PASS** | clause = C04; doc_lines = 36-37; anchor_line = 37; anchor_found = true; what = every nuclear PDG code is a well-formed 10LZZZAAAI (L = 0, I = 0, 1 &lt;= Z &lt;= A); n_checked = 1370; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C05` | **PASS** | clause = C05; doc_lines = 39-55; anchor_line = 55; anchor_found = true; what = every non-beam particle has a production vertex in the event; every vertex has an incoming particle; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C06` | **PASS** | clause = C06; doc_lines = 56-60; anchor_line = 58; anchor_found = true; what = e' and gamma\* are outgoing at the primary vertex (and the alpha on an alpha tag); X is produced at its own downstream v…; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C07` | **PASS** | clause = C07; doc_lines = 70-73; anchor_line = 71; anchor_found = true; what = every status is 1, 2, 3 or 4; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C08` | **PASS** | clause = C08; doc_lines = 72-78; anchor_line = 72; anchor_found = true; what = exactly one pdg-92 X per event, always status 3; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C09` | **PASS** | clause = C09; doc_lines = 77-81; anchor_line = 78; anchor_found = true; what = exactly one status-3 gamma\* with negative (spacelike) mass; on coherent files exactly one status-3 990 with mass -sqrt\|…; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C10` | **PASS** | clause = C10; doc_lines = 87-93; anchor_line = 91; anchor_found = true; what = the beam electron and e' carry generated_mass 0.51099895e-3 (rel 1e-9); n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C11` | **PASS** | clause = C11; doc_lines = 95-109; anchor_line = 95; anchor_found = true; what = the role table: status-3 pdg in {22, 990, 92, 2212, 2112, ions}, status-4 in {11, ions}; per channel the rows it writes…; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C12` | **PASS** | clause = C12; doc_lines = 111-121; anchor_line = 111; anchor_found = true; what = alpha tags: the struck cluster's end vertex emits &gt;= 1 status-1 partner (p, n or ion) and exactly one status-3 struck n…; n_checked = 200; files = 2; n_violations = 0; first_violation = n/a; failing_files = |
| `C13` | **FAIL** | clause = C13; doc_lines = 123-138; anchor_line = 130; anchor_found = true; what = four-momentum from the file: sum(status 1) == sum(beams) on a hadronized record, + X when X alone carries the HFS; neve…; n_checked = 800; files = 8; n_violations = 400; first_violation = inclusive-6Li-T2: event 0: residual 0.82 (status 1) / 0.657 (status 1 + X) of the beam energy; failing_files = inclusive-6Li-T2 / inclusive-7Li-T2 / inclusive-6Li-T0 / inclusive-6Li-T2-rc |
| `C14` | **FAIL** | clause = C14; doc_lines = 123-138 (charge); anchor_line = 123; anchor_found = true; what = charge from the file's PDG codes (Z of an ion code, PYTHIA ParticleData otherwise): sum(status 1) == sum(beams) whereve…; n_checked = 700; files = 7; n_violations = 300; first_violation = inclusive-6Li-T2: event 0: charge in +2, out +0; failing_files = inclusive-6Li-T2 / inclusive-7Li-T2 / inclusive-6Li-T2-rc |
| `C15` | **PASS** | clause = C15; doc_lines = 142-170; anchor_line = 142; anchor_found = true; what = all 24 named attributes on every event, each parsing as its listed HepMC3 type; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C16` | **PASS** | clause = C16; doc_lines = 147-159; anchor_line = 147; anchor_found = true; what = spin_J in {1, 1.5}; lam_e in {-1, 0, +1}; struck_cluster_m NaN on inclusive files and finite on tagged ones; channel = …; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C17` | **PASS** | clause = C17; doc_lines = 172-180; anchor_line = 179; anchor_found = true; what = a GenCrossSection on every event, index 0 finite and &gt; 0, its error finite and &gt;= 0; n_checked = 800; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C18` | **PASS** | clause = C18; doc_lines = 182-184; anchor_line = 184; anchor_found = true; what = a particle 'pol' attribute is a double and never the sentinel 9.0; n_checked = 1600; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C19` | **PASS** | clause = C19; doc_lines = 188-190; anchor_line = 188; anchor_found = true; what = ONE GenRunInfo block in the file (one tool line, one weight-name line, both before the first event); every event report…; n_checked = 8; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C20` | **PASS** | clause = C20; doc_lines = 192-193; anchor_line = 192; anchor_found = true; what = exactly one ToolInfo, name LiPolGen, version = lipolgen.__version__; n_checked = 8; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C21` | **PASS** | clause = C21; doc_lines = 194-201; anchor_line = 194; anchor_found = true; what = weight names: 'nominal' at 0, then spin_weight_1..N consecutive; the first event carries exactly that many values; n_checked = 8; files = 8; n_violations = 0; first_violation = n/a; failing_files = |
| `C22` | **PASS** | clause = C22; doc_lines = 202-219; anchor_line = 208; anchor_found = true; what = --rc tensor-band: the rc names APPENDED after the spin block, exactly rc_tensor_lo rc_tensor_hi rc_tail then the _k tri…; n_checked = 1; files = 1; n_violations = 0; first_violation = n/a; failing_files = |
| `C23` | **PASS** | clause = C23; doc_lines = 221-223; anchor_line = 221; anchor_found = true; what = --rc off: no rc_\* name, and 1 + N weight values; n_checked = 7; files = 7; n_violations = 0; first_violation = n/a; failing_files = |

## Planned, not wired

Nothing below is a measurement. The first table is this run's BLOCKED rows and sub-rows, with the reason each one gives; the other two are read verbatim from `BENCHMARK_PLAN.md` (§4's rows that have no harness, and §6, what cannot be benchmarked until data exists).

### Blocked in this run

| row | status | reason |
|---|---|---|
| `t1_g1d_world[ct18nlo]` | **BLOCKED** | environment: LHAPDF set CT18NLO not installed |
| `t1_g1d_world[nnpdfpol]` | **BLOCKED** | environment: LHAPDF set NNPDFpol11_100 not installed |
| `t1_nmc_f2d[ct18nlo]` | **BLOCKED** | environment: LHAPDF set CT18NLO not installed |
| `t1_nmc_f2d_over_f2p[ct18nlo]` | **BLOCKED** | environment: LHAPDF set CT18NLO not installed |
| `t2_bonus_spectator_shape` | **BLOCKED** | BLOCKED: the CLAS Physics Database (queried 2026-09-23) holds BONuS (E-03-012, eid 135) only as ONE measurement, E135M1 = F2n/F2p vs W\* (Baillie et al., PRL 108 (2012) 142001) -- no spectator-momentum spectrum, no Tkachenko 2014 entry; EG1b (eid 95, 146) is inclusive g1/A1 only. To unblock, download the Supplemental Material (tables of numerical results, their ref. [71]) of S. Tkachenko et al. (CLAS BONuS), Phys. Rev. C 89 (2014) 045206, doi:10.1103/PhysRevC.89.045206, https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.89.045206 (HTTP 401 without an APS login on 2026-09-23; arXiv:1402.2477 has no ancillary files), into validation/benchmarks/data/. Found instead and vendored unwired: Deeps F2N x P(p_s) at 0.30-0.53 GeV/c (eid 90, Klimenko et al., PRC 73 (2006) 035212). |
| `t3_li_magnetization_rfy` | **BLOCKED** | BLOCKED: reference not on disk and not open access. Download R. E. Rand, R. F. Frosch, M. R. Yearian, Phys. Rev. 144 (1966) 859, https://doi.org/10.1103/PhysRev.144.859 OR C. W. de Jager, H. de Vries, C. de Vries, At. Data Nucl. Data Tables 14 (1974) 479 (Table V), https://doi.org/10.1016/S0092-640X(74)80002-1 (docs/references/REFERENCES.md sec. 1 items 8 and 9), then vendor the lithium rows as validation/benchmarks/data/rfy1966_li_magnetization.csv |
| `t3_nmc_li6_over_d` | **BLOCKED (ENVIRONMENT)** | environment: LHAPDF set EPPS21nlo_CT18Anlo_Li6 not installed |
| `t4_djangoh_rad_noRad` | **BLOCKED** | no shared O(alpha) term: DJANGOH Rad=1 ran with IEL2=IEL31=IEL32=IEL33=0 (elastic radiative tail OFF; all eight ePIC logs), LiPolGen's rc_tail is the elastic(+QE) tail ONLY and applies no inelastic/vertex shift -- disjoint additive pieces, so must-not-contradict holds vacuously; the samples' W_h &gt;= 3 GeV cut also excludes the elastic tail (W_h = M_p) whatever IEL is (DJANGOH zeroes IEL2/IEL31..33 when WMIN &gt; M_p), so it unblocks only with a DJANGOH run at IEL31..33 != 0 AND WMIN &lt;= M_p |
| `t5_sum_rules[bjorken[nnpdfpol11]]` | **BLOCKED** | environment: LHAPDF set NNPDFpol11_100 not installed |
| `t5_sum_rules[fsi_unitarity]` | **BLOCKED** | not definable in the tree as built (no GlauberFsiWeight option realises int dGamma S[FSI + FSI^2] = 0: the default is absorptive by design, elastic_gain would be a tune, weight_normalised is a division, the table stops at k_max) |

### `BENCHMARK_PLAN.md` §4 rows with no harness (status column, verbatim)

| # | benchmark | tier | status (the plan's words) |
|---|---|---|---|
| 2 | **Cosyn–Weiss II Table II — re-derive the A_T∥ ↔ A_zz^wf mapping** | T2 | **resolved and fixed 2026-09-06** (not a harness; §5.1). The CW TABLE II gate is `tests/test_tagged.cpp` *"the Cosyn-Weiss deuteron tensor gate"* |
| 8 | **George–Knutson η — retitle** | T3 | **recorded, not a harness.** Code comments relabelled 2026-09-23 (`include/lipolgen/cluster_config.hpp`, three `///` blocks); doc sites G-1 … G-16 applied 2026-09-23; the T22 test name now says what it checks; G-17 (the Mäntysaari draft) left to the author, registry row 21 |

### `BENCHMARK_PLAN.md` §6 — what cannot be benchmarked until data exists (verbatim)

| # | observable | closest proxy | distance |
|---|---|---|---|
| U-1 | b₁(⁶Li), b₂(⁶Li) | b₁ᵈ (HERMES; CDKS/Miller) | A = 2 → 6 through `LI6_B1_PER_NUCLEON` = 2/6 — a modelling ASSUMPTION (inert α), not a convention |
| U-2 | tagged tensor asymmetry with an α/t spectator | Cosyn–Weiss A_T∥ with a NUCLEON spectator | different spectator, different wave function; and see §5.1 |
| U-3 | coherent ⁶Li tensor cos 2φ / a₂ | Mäntysaari *et al.* polarized-DEUTERON a₂(\|t\|) | the deuteron's quadrupole is 3.5× ⁶Li's and of the opposite sign; O5 is a closed-form estimate, not a Good–Walker amplitude |
| U-4 | ⁷Li rank-2 sector | Jaffe–Manohar (formalism only); Fu–Sun–Dong GPD forward limits | no number for ⁷Li anywhere; the sector is identically zero in-tree |
| U-5 | b₃, b₄, Δ (double-helicity-flip gluon) | Detmold–Shanahan lattice A₂ for the φ meson | a different hadron |
| U-6 | tensor-sector radiative corrections | Gakh–Shekhovtsova (zero citations, one Q² panel, deuteron) | no calculation for any A > 2 |
| U-7 | cluster-spectator Glauber FSI | Ciofi degli Atti–Kaptari ³He(e,e′d)X | the only cluster-spectator FSI calculation; a d spectator on A = 3 |
| U-8 | the α–d D-wave from data | George–Knutson η | refuted as a benchmark; η is a relabelling of the quadrupole dial |
| U-9 | A_zz(⁶Li) inclusive | the relation A_zz = −(2/3) b₁/F₁ (four sources agree) | a convention gate on the RELATION, not on either side |
| U-10 | ⁶Li elastic FF shape in the RC tail | F_q(0), F_m(0) on measured moments | the C0 SHAPE is a model band (now: the UVa density is a data-derived shape, and it lies OUTSIDE the [ho, vmc-ft] band on every priced tail number — `../open_items/run_2026-09-23/phase_B2_nuclear.md` §3.3) |
| U-11 | the spin bookkeeping | nothing external today | #1 above (wired 2026-09-23): the bookkeeping now has an external anchor, EPIOS Table II, passing at residual 0; the anchor is for an atomic-beam deuteron source, not for lithium |
| U-12 | polarized ⁶Li QE/elastic tails | 133 unpolarized QE points | the tensor QE tail is priced, not computed |
