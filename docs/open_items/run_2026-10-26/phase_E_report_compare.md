<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — phase E: the report regenerated, and its `--compare` against 2026-09-27

The output of the one regeneration of the integration stage, kept verbatim
(PLAN.md phase E1). Written by the integration stage; the board is
`STATUS.md`.

**Command** (from the repository root, after `source env.sh`; 2026-10-27
01:10:19–01:11:36 UTC, 77 s wall, exit 1 — the M2 rule: any change exits 1):

    LIPOLGEN_BENCH=1 python3 validation/benchmarks/make_report.py --all --plots \
        --compare docs/benchmarking/REPORT.json --out-dir docs/benchmarking

`--out-dir docs/benchmarking` is the default output directory spelled out:
`make_report.py` refuses `--plots` with `--compare` unless `--out-dir` is
given (exit 2, "--plots draws beside a written report", phase P's rule), and
`--compare` without `--out-dir` writes nothing. The baseline was the
committed `docs/benchmarking/REPORT.json` of 2026-09-27 (generated
2026-09-27T14:23:30Z at `ed4f9ff` plus uncommitted changes; 18 harnesses: 6 pass, 8 fail, 3 blocked, 1 blocked (environment)), read into
memory before the new report overwrote it. The recorded command line in
`REPORT.md` is the canonical `make_report.py --all --plots`.

**Result:** 23 harnesses — **13 pass, 10 fail, 0 blocked, 0 blocked
(environment), 0 error, 0 skipped (opt-in)**; **41 changes**: 4 STATUS
(`t2_bonus_spectator_shape` BLOCKED → FAIL, `t3_li_magnetization_rfy`
BLOCKED → PASS, `t3_nmc_li6_over_d` BLOCKED (ENVIRONMENT) → PASS,
`t4_djangoh_rad_noRad` BLOCKED → PASS), 8 MOVED, 5 SUBROW status changes
(`t1_g1d_world[ct18nlo]` and `[nnpdfpol]` BLOCKED → FAIL,
`t1_nmc_f2d[ct18nlo]` BLOCKED → FAIL, `t1_nmc_f2d_over_f2p[ct18nlo]`
BLOCKED → PASS, `t5_sum_rules[bjorken[nnpdfpol11]]` BLOCKED → FAIL), 24
SUBMOVED; and, informational (not counted), 5 NEW rows (`t1_hermes_f2d`
FAIL, `t1_r_world` PASS, `t2_boeglin_nk` PASS,
`t3_li6_ff_wiringa_schiavilla` PASS, `t5_delta_moment_bound` PASS) and
the 52 new sub-rows of the re-wired and new rows. Every change traces to a phase record
of this run: the LHAPDF sets installed (phase 0, `phase_A_baseline.md`), the
three re-wired rows (B1, B2, B3), PYTHIA 8.317 (T4; phase 0 / A), the Q_d
reference corrected (D), the new rows (C1, C2, C4, C5). `make_report.py
--check` then passed (2026-10-27 01:12:21–01:13:25 UTC, 64 s, exit 0). The
figures: 19 of 23 rows carry one (`docs/benchmarking/figures/`); none failed
to draw.

## The output, verbatim

The first 23 lines are the harnesses' own `REPORT |` lines as echoed by
`make_report.py` (each whole);
then the summary; then every `--compare` line.

```text
REPORT | t1_g1d_world | toy (shipped ToyG1, g1_nucleus(d)/2): compass 107.8/15 (p = 4.17e-16), hermes 69.4/15 (p = 5.75e-09), e143 54.9/28 (p = 0.00173) on Q2 >= 1, W2 >= 4 [diagonal]; sub-rows: mstw FAIL compass 53.2/15 (p = 3.52e-06), hermes 47.1/15 (p = 3.5e-05), e143 50.2/28 (p = 0.00607); ct18nlo FAIL compass 20.2/11 (p = 0.0429), hermes 64.3/11 (p = 1.46e-09), e143 59.0/24 (p = 8.96e-05); nnpdfpol FAIL compass 45.0/15 (p = 7.67e-05), hermes 52.6/15 (p = 4.53e-06), e143 51.0/28 (p = 0.00505) | g1d per nucleon: COMPASS (D-2), HERMES (D-1), E143 (D-6) via nnpdf-data 4.1.5; SMC recorded only; 15 + 15 + 28 points in the DIS cut | p(chi2, ndf) >= 0.01 on the DIS cut for each of compass, hermes, e143, stat (+) sys diagonal | FAIL
REPORT | t1_hermes_f2d | toy (shipped ToyF2) on F2d: chi2/ndf = 5303.5/68, p = 0 on Q2 >= 1, W2 >= 4 [per-bin diagonal + correlated 7.5 % normalisation]; sub-rows: mstw 318.9/68 (p = 4.03e-34); ct18nlo 1788.2/45 (p = 0); HERMES b1 inversion: F2d_implied/F2d_HERMES = 0.929 / 1.039 / 0.986 / 0.958 / 1.006 / 1.006 -> PER NUCLEON | HERMES F2d (and F2p) per nucleon, JHEP 05 (2011) 126, HEPData ins894309 Tables 2 (d) and 1 (p), 68 of 81 deuteron points in the DIS cut | p(chi2, ndf) >= 0.01 on the DIS cut, per-bin errors diagonal + the overall normalisation fully correlated | FAIL
REPORT | t1_nmc_f2d | toy (shipped ToyF2): chi2/ndf = 8307.9/155, p = 0 on Q2 >= 1, W2 >= 4 [diagonal]; sub-rows: mstw 416.5/155 (p = 8.21e-26); ct18nlo 3401.8/146 (p = 0) | NMC F2d per nucleon, NPB 483 (1997) 3, HEPData ins424154 tables 17-32 via nnpdf-data 4.1.5 (NMC_NC_NOTFIXED_D), 155 of 158 points in the DIS cut | p(chi2, ndf) >= 0.01 on the DIS cut, stat (+) sys diagonal | FAIL
REPORT | t1_nmc_f2d_over_f2p | toy (shipped ToyF2): chi2/ndf = 2161.4/211, p = 0 on Q2 >= 1, W2 >= 4 [diagonal]; sub-rows: mstw 247.2/211 (p = 0.0444; frozen hook 260.7); ct18nlo 215.5/189 (p = 0.0904; frozen hook 239.5) | NMC F2d/F2p (per nucleon), NPB 487 (1997) 3, HEPData ins426595 tables 2-21 via nnpdf-data 4.1.5 (NMC_NC_NOTFIXED), 211 of 260 points in the DIS cut | p(chi2, ndf) >= 0.01 on the DIS cut, stat (+) sys diagonal | FAIL
REPORT | t1_r_world | r_sigma_lt (shipped R) on NMC R: chi2/ndf = 12.54/10, p = 0.25 [diagonal, DIS cut]; r1998 (clip) on NMC 6.20/10 (p = 0.798); h1[r_sigma_lt] 8.82/8 (p = 0.358); h1[r1998] 8.18/8 (p = 0.416); r1990[r_sigma_lt] 209.76/68 (p = 2.28e-16); r1990[r1998] 63.05/68 (p = 0.647); e140[r_sigma_lt] 29.66/9 (p = 0.000502); e140[r1998] 12.71/9 (p = 0.176); bias (1+r1998)/(1+r_sigma_lt) at the sweet spots 1.1848 / 1.1780 / 1.0762 / 0.9563 | NMC R, NPB 483 (1997) 3, HEPData ins424154 Table 33 (PRIMARY, 10 of 10 points in the DIS cut); H1 F_L, PLB 665 (2008) 139, HEPData ins786161 Table 9 (+ Tables 1-8 recorded) (EXTRAPOLATION, 8 of 8 points in the DIS cut); Whitlow et al. R (global SLAC), PLB 250 (1990) 193, HEPData ins296980 Tables 3-13 (CLOSURE, 68 of 89 points in the DIS cut); SLAC E140 R (deuterium), PRD 49 (1994) 5641, HEPData ins360765 Tables 23-25 (CLOSURE, 9 of 11 points in the DIS cut) | p(chi2, ndf) >= 0.01 per data set on Q2 >= 1, W2 >= 4 GeV2, stat (+) sys diagonal; status = r_sigma_lt on the PRIMARY nmc set | PASS
REPORT | t2_boeglin_nk | n(k = p_m) x s, one free scale: hulthen (default) chi2 = 21.98/21 (p = 0.4, s = 0.481); av18 11.54 (p = 0.951); cdbonn 8.89 (p = 0.99) | Boeglin et al. PRL 107 (2011) 262501 Fig. 2, d(e,e'p)n sigma_red at Q^2 = 3.5 GeV^2, theta_nq = 35 deg (10) and 45 deg (12 points), p_m = 0.055-0.538 GeV/c, digitized from the vector PDF | p(chi2, 21) >= 0.01 with one free scale, the bars (+) an allowance A = 0.2 of each point for FSI and the CC1 off-shell prescription (declared in the hashed design block) | PASS
REPORT | t2_bonus_spectator_shape | Hulthen 40 mb GlauberFsiWeight over its backward mean inside the data's 1-sigma flat band in 4/12 (4 GeV) and 0/16 (5 GeV) DIS blocks; worst v/delta 2.62 (4 GeV Q2=1.66,W*=2.44,ps=0.135), 3.44 (5 GeV Q2=1.66,W*=2.03,ps=0.078); backward chi2 flat 739.4, tree shape 226.3 on 84 dof | CLAS BONuS R_D/S vs cos theta_pq, Tkachenko PRC 89 (2014) 045206 Suppl. (100 blocks; 32 in the DIS window W* >= 2, Q2 >= 1, 28 with backward bins) | max over cos <= -0.3 of abs(w/mean_bw(w) - 1) <= delta = S sigma(Rbar)/Rbar (stat, S = max(1, sqrt(chi2/ndf))), every block, both beams | FAIL
REPORT | t2_deeps_spectator_tail | shape chi2/ndf = 627.85/187 (p = 4.93e-49), Hulthen default, LC n(k) E_k/(2-alpha), cos_pq <= -0.3; AV18 579.12/187, CD-Bonn 278.81/187 | CLAS Deeps F2N x P(p_s), Klimenko PRC 73 (2006) 035212 (60 blocks at Q2 = 1.8/2.8, 245 PWIA-window points) | p(chi2, ndf) >= 0.01, default wave, one free scale per (Q2, W*, cos) group | FAIL
REPORT | t2_deuteron_static | eta: cdbonn 0.025571 (-0.07 sigma), av18 0.025045 (-1.39 sigma), hulthen 0.032478 (+17.20 sigma); Q_d (impulse approx., recorded): 0.270493 (-5.42 %) / 0.269670 (-5.71 %) / 0.290349 (+1.52 %) fm^2 | eta_d = 0.0256 +- 0.0004 (Rodning-Knutson PRC 41 (1990) 898, Eq. (10)), Q_d = +0.2860 +- 0.0015 fm^2 (Bishop-Cheung PRA 20 (1979) 381, confirming Reid-Vaida); read from the papers 2026-10-26, Q_d corrected from 0.2859(3) | abs(eta - eta_d) <= 2 sigma for each of cdbonn, av18, hulthen; Q_d recorded, not asserted | FAIL
REPORT | t2_hermes_b1_table2 | cdks_conv_mstw: chi2/ndf b1 = 22.26/6, A_zz = 22.01/6; 2/6 b1 bins within tol | HERMES PRL 95 (2005) 242001 Table II, 6 bins, per nucleon | per bin sqrt(stat^2+syst^2) | FAIL
REPORT | t3_li6_charge_ff_fb | HoSpin1FF first C0 zero q0 = 3.0998 fm^-1 (shipped; T11 window [2.9, 3.3] not moved) | band [2.6944 (UVa FB zero), 2.8284 (Li71 minimum)] fm^-1 | q0 inside the band (must not contradict) | FAIL
REPORT | t3_li6_ff_wiringa_schiavilla | HoSpin1FF C0 zero q0 = 3.0998 fm^-1 (shipped ho); vmc-ft: no zero below 5 fm^-1 | Wiringa-Schiavilla PRL 81 (1998) 4317 Fig. 1, digitized from the vector PDF: the IA+MEC C0 line is drawn non-zero to q_last = 2.9190 fm^-1 (first zero > q_last; not located above it); IA to 3.2000 | q0 >= q_last - (0.0047619 + 0.1) = 2.8143 fm^-1 (the plan's |q0 - q0_WS| <= delta + 0.1 on the zeros the figure allows; the upper side is not tested; declared in the hashed design block) | PASS
REPORT | t3_li6_rms_devries | HoSpin1FF.for_ion(li6()) C0 slope, nucleon-folded (G_E^p + G_E^n): r_ch = 2.57100 fm (+- 1.4e-09 numerical) | de Vries ADNDT 36 (1987) Table I 2.54(5) / 2.56(5) / 2.57(10) fm; Angeli-Marinova 2.589(39) fm (closure, recorded); survey transcription | abs(r - r_i) <= 2 sigma_i for each de Vries row (i.e. r in [2.46, 2.64] fm) | PASS
REPORT | t3_li_magnetization_rfy | HoSpin1FF F_mag^2 (q_z = 1.30 fm^-1, b = 1.85 fm; unfitted starting values) at RFY Table I's 7 points: chi2 = 13.35 / 6 with one free scale s = 1.419 +- 0.206 (p = 0.0378); s = 1: chi2 = 17.49 / 7 (p = 0.0145); 7Li: not implemented in the tree | Rand-Frosch-Yearian, Phys. Rev. 144 (1966) 859, Table I: 6Li F_M^2 at q^2 = 0.72-1.92 fm^-2 (q = 0.85-1.39 fm^-1), 7 points, read from the page image; Table II (7Li, 24 points) recorded | p(chi2, 6) >= 0.01 with one free scale, the paper's errors, diagonal (declared in the hashed design block) | PASS
REPORT | t3_nmc_li6_over_d | chi2/ndf = 4.6269/4 (x >= 0.30, p = 0.328); on-grid 26.941/15 (p = 0.029) | NMC F2(6Li)/F2(D), HEPData ins394050 Table 1 (15 of 24 on grid) | p(chi2, ndf) >= 0.01 on both sets | PASS
REPORT | t4_djangoh_rad_noRad | R = sigma_tail(LiPolGen)/sigma_tail(DJANGOH) per Q2 bin: polrad-full 1.0199 / 1.0261 / 1.0219 / 1.0098; t-peak (shipped) 0.6118 / 0.8968 / 1.0112 / 1.0067 | DJANGOH 4.6.21 elastic radiative tail IEL31+IEL32+IEL33 (fixed alpha), fraction of its own Born in the matched window: 0.04073 / 0.01644 / 0.01710 / 0.02189; control reproduces the 8 ePIC numbers: 8/8 | polrad-full: abs(R - 1) <= 0.03 + 2 delta_b each bin; t-peak (a one-sided ceiling; the tree calls it a lower bound): R <= 1.03 + 2 delta_b each bin; delta_b = DJANGOH's relative integration error (3.7e-04 / 4.8e-04 / 4.4e-04 / 3.9e-04) | PASS
REPORT | t4_pythia_ep_closure | max |bridge/stock - 1| = 0.076 (W2 n_ch_eta[-2,0)) over 15 compared sub-rows, 0 fail; <n_ch> W1 1.007; <n_ch> W2 0.997 | stock PYTHIA 8.317 e p WeakBosonExchange:ff2ff(t:gmZ), D-1 settings, e 10 x p 100, PDF:pSet 13, 20000 + 20000 events | |bridge/stock - 1| <= 0.20 per observable per window, at matched (x, Q2) (D-1) | PASS
REPORT | t5_delta_moment_bound | shipped toy_delta_gluon (scale 0.01) per nucleon, 6Li: M2 = int x Delta dx = 1.2383e-04 / 1.3108e-04 / 1.3805e-04 at Q2 = 2 / 5 / 10 GeV2, |M2|/B = 0.0644 / 0.0803 / 0.0941 (V6, Delta_DS16 = 2 Delta_tree: 0.1289 / 0.1605 / 0.1881); helicity-state rates min W = 0.949400 (A_T in [-0.101, 0.0506]) | gluonic Soffer-type bound, n = 2 moment: Detmold-Shanahan PRD 94 (2016) 014507 Eqs. (13), (32) [(36)]; f1^g CT18NLO, g1^g = P_eff x NNPDFpol11 Delta g, f1LL^g = 0, Tr Q^2 = 2/3, alpha_s CT18NLO; helicity positivity -2 <= A_T <= 1 (Cosyn et al. EPJ A 61 (2025) 83 Eq. (2)) | a bound under assumptions: PASS iff 1 - |M2|/B >= -1e-12 at each Q2 = 2, 5, 10 GeV2 for the shipped toy AND every helicity-state rate >= -1e-12 (-2 <= A_T <= 1) on the vendored grid; the letter's moment_A/B carry their own status, recorded | PASS
REPORT | t5_epios_source_modes | max moment round-trip residual 0.000e+00 over 8 modes (min population 0) | EPIOS Table II ideal (P_z, P_zz), 8 modes | 1e-12 | PASS
REPORT | t5_est_identity | max|spin_temperature_pzz(1,Pz) - closed form| = 1.138e-14 at Pz = -0.5996 over 19999 interior points | P_zz = 2 - sqrt(4 - 3 P_z^2) (EST, analytic) | 2.30647e-14 | PASS
REPORT | t5_sum_rules | bc_ww closure max|resid|/tol = 0.667 over 45 cells (PASS); Bjorken ToyG1 Q2=5: int_{1e-6}^1 (g1p-g1n) = 0.6430, last decade +0.2198 (divergent: g1n ~ x^-1.217) | BC: -x_min int g1/u du (exact truncation of int g2^WW); Bjorken: g_A/6 C_Bj(alpha_s) = 0.1854 +- 0.0023 at Q2=5 | BC: Peano trapezoid bound (h^2/8)L^2 int|D^2 g1| + quad errors; Bjorken: converged (last decade <= 1 %) and within 10 % | FAIL
REPORT | t5_weighted_unweighted | unweighted sigma = 89.1934 +- 0.0512 pb (N = 3032577) vs weighted (Mode W) 89.1674 +- 0.0516 pb (N = 2987706, <w> = 1.01472): pull +0.359 | same physics point, the other route: 6Li inclusive config 1, x >= 0.1, y >= 0.5, helicity-flip apar+ (exact sigma_tot_pb 89.1671 pb; a no-op Mode W would pull +18.3) | |pull| <= 3 | PASS
REPORT | t6_hepmc3_convention | 23/23 clauses hold on 8 files, 800 events (lipolgen-run, all CLI channels, T2 + T0 + --rc) | docs/HEPMC3_CONVENTION.md sha256 56d3a91a9bdf789b; 23 clauses | every clause on every event; numeric identities at the doc's 1e-9 (L137-138), four-momentum per component relative to max(1, summed beam E) | PASS
wrote docs/benchmarking/REPORT.md and docs/benchmarking/REPORT.json
summary: 13 pass, 10 fail, 0 blocked, 0 blocked (environment), 0 error, 0 skipped (opt-in)
MOVED    t1_g1d_world:
    was: toy (shipped ToyG1, g1_nucleus(d)/2): compass 107.8/15 (p = 4.17e-16), hermes 69.4/15 (p = 5.75e-09), e143 54.9/28 (p = 0.00173) on Q2 >= 1, W2 >= 4 [diagonal]; sub-rows: mstw FAIL compass 53.2/15 (p = 3.52e-06), hermes 47.1/15 (p = 3.5e-05), e143 50.2/28 (p = 0.00607); ct18nlo BLOCKED (environment: LHAPDF set CT18NLO not installed); nnpdfpol BLOCKED (environment: LHAPDF set NNPDFpol11_100 not installed)
    now: toy (shipped ToyG1, g1_nucleus(d)/2): compass 107.8/15 (p = 4.17e-16), hermes 69.4/15 (p = 5.75e-09), e143 54.9/28 (p = 0.00173) on Q2 >= 1, W2 >= 4 [diagonal]; sub-rows: mstw FAIL compass 53.2/15 (p = 3.52e-06), hermes 47.1/15 (p = 3.5e-05), e143 50.2/28 (p = 0.00607); ct18nlo FAIL compass 20.2/11 (p = 0.0429), hermes 64.3/11 (p = 1.46e-09), e143 59.0/24 (p = 8.96e-05); nnpdfpol FAIL compass 45.0/15 (p = 7.67e-05), hermes 52.6/15 (p = 4.53e-06), e143 51.0/28 (p = 0.00505)
SUBROW   t1_g1d_world[ct18nlo]: BLOCKED -> FAIL
SUBMOVED t1_g1d_world[ct18nlo]: kernel_max_rel_diff n/a -> 0; q2_floor n/a -> 1.67702; world_chi2 n/a -> 143.433; world_ndf n/a -> 46; world_p n/a -> 6.05607e-12
SUBROW   t1_g1d_world[ct18nlo/compass]: new (informational)
SUBROW   t1_g1d_world[ct18nlo/e143]: new (informational)
SUBROW   t1_g1d_world[ct18nlo/hermes]: new (informational)
SUBROW   t1_g1d_world[ct18nlo/smc]: new (informational)
SUBROW   t1_g1d_world[nnpdfpol]: BLOCKED -> FAIL
SUBMOVED t1_g1d_world[nnpdfpol]: kernel_max_rel_diff n/a -> 0; q2_floor n/a -> 1; world_chi2 n/a -> 148.534; world_ndf n/a -> 58; world_p n/a -> 6.9932e-10
SUBROW   t1_g1d_world[nnpdfpol/compass]: new (informational)
SUBROW   t1_g1d_world[nnpdfpol/e143]: new (informational)
SUBROW   t1_g1d_world[nnpdfpol/hermes]: new (informational)
SUBROW   t1_g1d_world[nnpdfpol/smc]: new (informational)
NEW      t1_hermes_f2d: FAIL (not in the baseline; informational)
MOVED    t1_nmc_f2d:
    was: toy (shipped ToyF2): chi2/ndf = 8307.9/155, p = 0 on Q2 >= 1, W2 >= 4 [diagonal]; sub-rows: mstw 416.5/155 (p = 8.21e-26); ct18nlo BLOCKED (environment: LHAPDF set CT18NLO not installed)
    now: toy (shipped ToyF2): chi2/ndf = 8307.9/155, p = 0 on Q2 >= 1, W2 >= 4 [diagonal]; sub-rows: mstw 416.5/155 (p = 8.21e-26); ct18nlo 3401.8/146 (p = 0)
SUBROW   t1_nmc_f2d[ct18nlo]: BLOCKED -> FAIL
SUBMOVED t1_nmc_f2d[ct18nlo]: all.chi2 n/a -> 3523.99; all.chi2_cov n/a -> 6838.21; all.mean_ratio n/a -> 1.07888; all.n n/a -> 158; all.p n/a -> 0; all.p_cov n/a -> 0; all.worst_pull n/a -> 12.2692; all.worst_q2 n/a -> 2.5; all.worst_x n/a -> 0.05; dis.chi2 n/a -> 3401.76; dis.chi2_cov n/a -> 4553.94; dis.mean_ratio n/a -> 1.08624; dis.n n/a -> 146; dis.p n/a -> 0; dis.p_cov n/a -> 0; dis.worst_pull n/a -> 12.2692; dis.worst_q2 n/a -> 2.5; dis.worst_x n/a -> 0.05; kernel_max_rel_diff n/a -> 0; q2_floor n/a -> 1.67702
MOVED    t1_nmc_f2d_over_f2p:
    was: toy (shipped ToyF2): chi2/ndf = 2161.4/211, p = 0 on Q2 >= 1, W2 >= 4 [diagonal]; sub-rows: mstw 247.2/211 (p = 0.0444; frozen hook 260.7); ct18nlo BLOCKED (environment: LHAPDF set CT18NLO not installed)
    now: toy (shipped ToyF2): chi2/ndf = 2161.4/211, p = 0 on Q2 >= 1, W2 >= 4 [diagonal]; sub-rows: mstw 247.2/211 (p = 0.0444; frozen hook 260.7); ct18nlo 215.5/189 (p = 0.0904; frozen hook 239.5)
SUBROW   t1_nmc_f2d_over_f2p[ct18nlo]: BLOCKED -> PASS
SUBMOVED t1_nmc_f2d_over_f2p[ct18nlo]: all.chi2 n/a -> 295.361; all.chi2_cov n/a -> 317.557; all.chi2_legacy_cov n/a -> 290.171; all.chi2_legacy_diag n/a -> 294.85; all.mean_diff n/a -> 0.000289167; all.n n/a -> 260; all.p n/a -> 0.0649639; all.p_cov n/a -> 0.00848358; all.p_legacy_cov n/a -> 0.0961021; all.worst_pull n/a -> 3.18221; all.worst_q2 n/a -> 0.64; all.worst_x n/a -> 0.008; all_frozen.chi2 n/a -> 321.445; all_frozen.chi2_cov n/a -> 342.406; all_frozen.chi2_legacy_cov n/a -> 300.31; all_frozen.chi2_legacy_diag n/a -> 320.793; all_frozen.mean_diff n/a -> 0.000671017; all_frozen.n n/a -> 260; all_frozen.p n/a -> 0.00560404; all_frozen.p_cov n/a -> 0.000452474; all_frozen.p_legacy_cov n/a -> 0.0433854; all_frozen.worst_pull n/a -> -3.29512; all_frozen.worst_q2 n/a -> 5.49; all_frozen.worst_x n/a -> 0.18; dis.chi2 n/a -> 215.508; dis.chi2_cov n/a -> 226.204; dis.chi2_legacy_cov n/a -> 202.52; dis.chi2_legacy_diag n/a -> 214.689; dis.mean_diff n/a -> -0.00122211; dis.n n/a -> 189; dis.p n/a -> 0.0903941; dis.p_cov n/a -> 0.033293; dis.p_legacy_cov n/a -> 0.237752; dis.worst_pull n/a -> -3.07475; dis.worst_q2 n/a -> 5.49; dis.worst_x n/a -> 0.18; dis_frozen.chi2 n/a -> 239.52; dis_frozen.chi2_cov n/a -> 250.791; dis_frozen.chi2_legacy_cov n/a -> 220.018; dis_frozen.chi2_legacy_diag n/a -> 238.569; dis_frozen.mean_diff n/a -> -0.00127637; dis_frozen.n n/a -> 189; dis_frozen.p n/a -> 0.00752728; dis_frozen.p_cov n/a -> 0.00175431; dis_frozen.p_legacy_cov n/a -> 0.0606314; dis_frozen.worst_pull n/a -> -3.29512; dis_frozen.worst_q2 n/a -> 5.49; dis_frozen.worst_x n/a -> 0.18; kernel_max_rel_diff n/a -> 0; max_abs_frozen_minus_q2 n/a -> 0.0166052; q2_floor n/a -> 1.67702
NEW      t1_r_world: PASS (not in the baseline; informational)
NEW      t2_boeglin_nk: PASS (not in the baseline; informational)
STATUS   t2_bonus_spectator_shape: BLOCKED -> FAIL
MOVED    t2_bonus_spectator_shape:
    was: not evaluated
    now: Hulthen 40 mb GlauberFsiWeight over its backward mean inside the data's 1-sigma flat band in 4/12 (4 GeV) and 0/16 (5 GeV) DIS blocks; worst v/delta 2.62 (4 GeV Q2=1.66,W*=2.44,ps=0.135), 3.44 (5 GeV Q2=1.66,W*=2.03,ps=0.078); backward chi2 flat 739.4, tree shape 226.3 on 84 dof
SUBROW   t2_bonus_spectator_shape[4GeV/av18-40mb]: new (informational)
SUBROW   t2_bonus_spectator_shape[4GeV/hulthen-20mb]: new (informational)
SUBROW   t2_bonus_spectator_shape[4GeV/hulthen-40mb]: new (informational)
SUBROW   t2_bonus_spectator_shape[4GeV/hulthen-40mb stat+syst]: new (informational)
SUBROW   t2_bonus_spectator_shape[5GeV/av18-40mb]: new (informational)
SUBROW   t2_bonus_spectator_shape[5GeV/hulthen-20mb]: new (informational)
SUBROW   t2_bonus_spectator_shape[5GeV/hulthen-40mb]: new (informational)
SUBROW   t2_bonus_spectator_shape[5GeV/hulthen-40mb stat+syst]: new (informational)
SUBROW   t2_bonus_spectator_shape[F2n/F2p ToyF2 vs W*>1.8 table (recorded)]: new (informational)
SUBROW   t2_bonus_spectator_shape[PWIA prescription instant/LC (recorded, not FSI)]: new (informational)
SUBROW   t2_bonus_spectator_shape[forward rise (recorded)]: new (informational)
SUBROW   t2_bonus_spectator_shape[out of model (census)]: new (informational)
MOVED    t2_deuteron_static:
    was: eta: cdbonn 0.025571 (-0.07 sigma), av18 0.025045 (-1.39 sigma), hulthen 0.032478 (+17.20 sigma); Q_d (impulse approx., recorded): 0.270493 (-5.39 %) / 0.269670 (-5.68 %) / 0.290349 (+1.56 %) fm^2
    now: eta: cdbonn 0.025571 (-0.07 sigma), av18 0.025045 (-1.39 sigma), hulthen 0.032478 (+17.20 sigma); Q_d (impulse approx., recorded): 0.270493 (-5.42 %) / 0.269670 (-5.71 %) / 0.290349 (+1.52 %) fm^2
SUBMOVED t2_deuteron_static[av18]: q_d_distance -0.0162297 -> -0.0163297; q_d_rel_distance -0.0567672 -> -0.057097
SUBMOVED t2_deuteron_static[cdbonn]: q_d_distance -0.0154066 -> -0.0155066; q_d_quad_error 7.38298e-15 -> 7.43849e-15; q_d_rel_distance -0.0538881 -> -0.0542189
SUBMOVED t2_deuteron_static[hulthen]: q_d_distance 0.00444905 -> 0.00434905; q_d_quad_error 1.22347e-13 -> 1.22291e-13; q_d_rel_distance 0.0155616 -> 0.0152065
NEW      t3_li6_ff_wiringa_schiavilla: PASS (not in the baseline; informational)
STATUS   t3_li_magnetization_rfy: BLOCKED -> PASS
MOVED    t3_li_magnetization_rfy:
    was: r_mag(6Li) = 2.9469 fm (tree F_m slope; UNCOMPARED)
    now: HoSpin1FF F_mag^2 (q_z = 1.30 fm^-1, b = 1.85 fm; unfitted starting values) at RFY Table I's 7 points: chi2 = 13.35 / 6 with one free scale s = 1.419 +- 0.206 (p = 0.0378); s = 1: chi2 = 17.49 / 7 (p = 0.0145); 7Li: not implemented in the tree
SUBROW   t3_li_magnetization_rfy[closure[6Li Table VII]]: new (informational)
SUBROW   t3_li_magnetization_rfy[closure[7Li Table VIII]]: new (informational)
SUBROW   t3_li_magnetization_rfy[closure_refit[7Li Table VIII, post-hoc]]: new (informational)
SUBROW   t3_li_magnetization_rfy[li6_absolute[s = 1]]: new (informational)
SUBROW   t3_li_magnetization_rfy[li6_shape[free scale]]: new (informational)
SUBROW   t3_li_magnetization_rfy[li6_shape[q2 from E0]]: new (informational)
SUBROW   t3_li_magnetization_rfy[li7[Table II]]: new (informational)
SUBROW   t3_li_magnetization_rfy[q2_column_check]: new (informational)
STATUS   t3_nmc_li6_over_d: BLOCKED (ENVIRONMENT) -> PASS
STATUS   t4_djangoh_rad_noRad: BLOCKED -> PASS
MOVED    t4_djangoh_rad_noRad:
    was: t-peak: 1.01617 / 1.01317 / 1.01497 / 1.02069; polrad-full: 1.0271 / 1.01507 / 1.01514 / 1.02076
    now: R = sigma_tail(LiPolGen)/sigma_tail(DJANGOH) per Q2 bin: polrad-full 1.0199 / 1.0261 / 1.0219 / 1.0098; t-peak (shipped) 0.6118 / 0.8968 / 1.0112 / 1.0067
SUBROW   t4_djangoh_rad_noRad[control[Rad=0, Q2 1-10]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[control[Rad=0, Q2 10-100]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[control[Rad=0, Q2 100-1000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[control[Rad=0, Q2 1000-10000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[control[Rad=1, Q2 1-10]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[control[Rad=1, Q2 10-100]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[control[Rad=1, Q2 100-1000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[control[Rad=1, Q2 1000-10000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[history_epic[Q2 1-10]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[history_epic[Q2 10-100]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[history_epic[Q2 100-1000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[history_epic[Q2 1000-10000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[polrad-full[Q2 1-10]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[polrad-full[Q2 10-100]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[polrad-full[Q2 100-1000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[polrad-full[Q2 1000-10000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[recorded[Q2 1-10]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[recorded[Q2 10-100]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[recorded[Q2 100-1000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[recorded[Q2 1000-10000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[t-peak[Q2 1-10]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[t-peak[Q2 10-100]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[t-peak[Q2 100-1000]]: new (informational)
SUBROW   t4_djangoh_rad_noRad[t-peak[Q2 1000-10000]]: new (informational)
MOVED    t4_pythia_ep_closure:
    was: max |bridge/stock - 1| = 0.023 (W2 n_ch_eta[-2,0)) over 15 compared sub-rows, 0 fail; <n_ch> W1 1.003; <n_ch> W2 0.998
    now: max |bridge/stock - 1| = 0.076 (W2 n_ch_eta[-2,0)) over 15 compared sub-rows, 0 fail; <n_ch> W1 1.007; <n_ch> W2 0.997
SUBMOVED t4_pythia_ep_closure[W1/n_ch]: bridge 8.05391 -> 8.20351; bridge_err 0.0194701 -> 0.019703; bridge_raw 8.1546 -> 8.3075; distance 0.00339012 -> 0.00729505; ratio 1.00339 -> 1.0073; ratio_err 0.00373336 -> 0.0037645; ratio_raw 1.01593 -> 1.02006; stock 8.0267 -> 8.1441; stock_err 0.0227026 -> 0.0233189; stock_raw 8.0267 -> 8.1441
SUBMOVED t4_pythia_ep_closure[W1/n_ch_eta[-2,0)]: bridge 0.680522 -> 0.693477; bridge_err 0.00538898 -> 0.00549713; bridge_raw 0.7344 -> 0.74295; distance 0.0160173 -> 0.0106047; n_stock_particles 13832 -> 13724; ratio 0.983983 -> 1.0106; ratio_err 0.0142945 -> 0.0147439; ratio_raw 1.06189 -> 1.0827; stock 0.6916 -> 0.6862; stock_err 0.00842306 -> 0.00840446; stock_raw 0.6916 -> 0.6862
SUBMOVED t4_pythia_ep_closure[W1/n_ch_eta[-inf,-2)]: bridge 0.028317 -> 0.0274033; bridge_err 0.00115488 -> 0.00114868; bridge_raw 0.0325 -> 0.03175; distance 0.00989624 -> 0.0316846; n_stock_particles 572 -> 566; ratio 0.990104 -> 0.968315; ratio_err 0.0621672 -> 0.0599736; ratio_raw 1.13636 -> 1.12191; stock 0.0286 -> 0.0283; stock_err 0.00136536 -> 0.00129036; stock_raw 0.0286 -> 0.0283
SUBMOVED t4_pythia_ep_closure[W1/n_ch_eta[0,2)]: bridge 2.46123 -> 2.49663; bridge_err 0.0113746 -> 0.0114792; bridge_raw 2.49415 -> 2.53945; distance 0.00135466 -> 0.000490194; n_stock_particles 49158 -> 49957; ratio 1.00135 -> 0.99951; ratio_err 0.00696589 -> 0.00700139; ratio_raw 1.01475 -> 1.01665; stock 2.4579 -> 2.49785; stock_err 0.0127797 -> 0.0132001; stock_raw 2.4579 -> 2.49785
SUBMOVED t4_pythia_ep_closure[W1/n_ch_eta[2,4)]: bridge 2.91045 -> 3.00318; bridge_err 0.0128979 -> 0.0132006; bridge_raw 2.91315 -> 3.0044; distance 0.0127702 -> 0.0190966; n_stock_particles 57475 -> 58938; ratio 1.01277 -> 1.0191; ratio_err 0.0063818 -> 0.00634547; ratio_raw 1.01371 -> 1.01951; stock 2.87375 -> 2.9469; stock_err 0.0128736 -> 0.0129962; stock_raw 2.87375 -> 2.9469
SUBMOVED t4_pythia_ep_closure[W1/n_ch_eta[4,+inf)]: bridge 1.97339 -> 1.98283; bridge_err 0.00836494 -> 0.00839894; bridge_raw 1.9804 -> 1.98895; distance 0.000737182 -> 0.00101778; n_stock_particles 39497 -> 39697; ratio 0.999263 -> 0.998982; ratio_err 0.00596837 -> 0.00597514; ratio_raw 1.00281 -> 1.00207; stock 1.97485 -> 1.98485; stock_err 0.00830989 -> 0.0083818; stock_raw 1.97485 -> 1.98485
SUBMOVED t4_pythia_ep_closure[W1/sum_empz]: bridge 2.74055 -> 2.74493; bridge_err 0.00637262 -> 0.0064404; distance 0.00647918 -> 0.00353297; ratio 1.00648 -> 1.00353; ratio_err 0.0075329 -> 0.00752765; ratio_raw 1.06672 -> 1.0619; stock 2.7229 -> 2.73527; stock_err 0.0193708 -> 0.0194881; stock_raw 2.7229 -> 2.73527
SUBMOVED t4_pythia_ep_closure[W1/sum_pt_ch]: bridge 3.84917 -> 3.88154; bridge_err 0.00983067 -> 0.00983322; bridge_raw 3.87739 -> 3.91528; distance 0.00576844 -> 0.00863876; ratio 1.00577 -> 1.00864; ratio_err 0.00387963 -> 0.00386883; ratio_raw 1.01314 -> 1.0174; stock 3.82709 -> 3.8483; stock_err 0.0110632 -> 0.0110834; stock_raw 3.82709 -> 3.8483
SUBMOVED t4_pythia_ep_closure[W2/n_ch]: bridge 7.24281 -> 7.38264; bridge_err 0.0182964 -> 0.0184182; bridge_raw 7.3014 -> 7.4489; distance 0.00244989 -> 0.00284517; ratio 0.99755 -> 0.997155; ratio_err 0.00378214 -> 0.00373987; ratio_raw 1.00562 -> 1.00611; stock 7.2606 -> 7.4037; stock_err 0.0205277 -> 0.0207336; stock_raw 7.2606 -> 7.4037
SUBMOVED t4_pythia_ep_closure[W2/n_ch_eta[-2,0)]: bridge 0.0372096 -> 0.0358614; bridge_err 0.00164641 -> 0.00155116; bridge_raw 0.03925 -> 0.03875; distance 0.0233704 -> 0.0757371; n_stock_particles 762 -> 776; ratio 0.97663 -> 0.924263; ratio_err 0.0680833 -> 0.0614823; ratio_raw 1.03018 -> 0.998711; stock 0.0381 -> 0.0388; stock_err 0.00205247 -> 0.00196085; stock_raw 0.0381 -> 0.0388
SUBMOVED t4_pythia_ep_closure[W2/n_ch_eta[-inf,-2)]: bridge 0.000238271 -> 0.000182074; bridge_err 0.000106293 -> 9.10119e-05; bridge_raw 0.00025 -> 0.0002; distance 0.588476 -> 0.271706; n_stock_particles 3 -> 5; ratio 1.58848 -> 0.728294; ratio_err 1.15894 -> 0.488459; ratio_raw 1.66667 -> 0.8; stock 0.00015 -> 0.00025; stock_err 8.65982e-05 -> 0.000111792; stock_raw 0.00015 -> 0.00025
SUBMOVED t4_pythia_ep_closure[W2/n_ch_eta[0,2)]: bridge 2.49662 -> 2.53452; bridge_err 0.0105711 -> 0.0105762; bridge_raw 2.56265 -> 2.6071; distance 0.00195212 -> 0.00446787; n_stock_particles 50030 -> 50465; ratio 0.998048 -> 1.00447; ratio_err 0.00739101 -> 0.00737996; ratio_raw 1.02445 -> 1.03323; stock 2.5015 -> 2.52325; stock_err 0.0151981 -> 0.0152584; stock_raw 2.5015 -> 2.52325
SUBMOVED t4_pythia_ep_closure[W2/n_ch_eta[2,4)]: bridge 3.2194 -> 3.31303; bridge_err 0.0133563 -> 0.0135194; bridge_raw 3.2034 -> 3.29615; distance 0.00297341 -> 0.00687143; n_stock_particles 64580 -> 66719; ratio 0.997027 -> 0.993129; ratio_err 0.00588933 -> 0.00575882; ratio_raw 0.992072 -> 0.988069; stock 3.229 -> 3.33595; stock_err 0.0135771 -> 0.0137433; stock_raw 3.229 -> 3.33595
SUBMOVED t4_pythia_ep_closure[W2/n_ch_eta[4,+inf)]: bridge 1.48935 -> 1.49904; bridge_err 0.00745836 -> 0.00745683; bridge_raw 1.49585 -> 1.5067; distance 0.00167653 -> 0.00425726; n_stock_particles 29837 -> 30109; ratio 0.998323 -> 0.995743; ratio_err 0.00712089 -> 0.00703967; ratio_raw 1.00268 -> 1.00083; stock 1.49185 -> 1.50545; stock_err 0.00757757 -> 0.00756283; stock_raw 1.49185 -> 1.50545
SUBMOVED t4_pythia_ep_closure[W2/sum_empz]: bridge 2.006 -> 2.00161; bridge_err 0.00364611 -> 0.00363994; distance 0.002417 -> 0.00313219; ratio 1.00242 -> 1.00313; ratio_err 0.00637546 -> 0.0063713; ratio_raw 1.02932 -> 1.03231; stock 2.00117 -> 1.99536; stock_err 0.0121968 -> 0.0121428; stock_raw 2.00117 -> 1.99536
SUBMOVED t4_pythia_ep_closure[W2/sum_pt_ch]: bridge 6.13132 -> 6.16328; bridge_err 0.0158241 -> 0.0157376; bridge_raw 6.17419 -> 6.20385; distance 0.00057363 -> 0.00496323; ratio 1.00057 -> 0.995037; ratio_err 0.00399322 -> 0.00393079; ratio_raw 1.00757 -> 1.00159; stock 6.1278 -> 6.19403; stock_err 0.0186538 -> 0.0186702; stock_raw 6.1278 -> 6.19403
NEW      t5_delta_moment_bound: PASS (not in the baseline; informational)
SUBROW   t5_sum_rules[bjorken[nnpdfpol11]]: BLOCKED -> FAIL
SUBMOVED t5_sum_rules[bjorken[nnpdfpol11]]: per_q2.10.0.reference.alpha_s n/a -> 0.247637; per_q2.10.0.reference.c_bj n/a -> 0.894197; per_q2.10.0.reference.err n/a -> 0.00153939; per_q2.10.0.reference.err_alpha_s n/a -> 0.000503248; per_q2.10.0.reference.err_g_a n/a -> 0.000193743; per_q2.10.0.reference.err_truncation n/a -> 0.00144184; per_q2.10.0.reference.gamma n/a -> 0.190062; per_q2.10.0.reference.nf n/a -> 4; per_q2.10.0.reference.q2 n/a -> 10; per_q2.10.0.verdict.distance n/a -> 0.018204; per_q2.10.0.verdict.distance_rel n/a -> 0.0957794; per_q2.10.0.verdict.last_decade_step n/a -> 0.00427738; per_q2.2.0.reference.alpha_s n/a -> 0.358817; per_q2.2.0.reference.c_bj n/a -> 0.822753; per_q2.2.0.reference.err n/a -> 0.00462445; per_q2.2.0.reference.err_alpha_s n/a -> 0.00145436; per_q2.2.0.reference.err_g_a n/a -> 0.000178263; per_q2.2.0.reference.err_truncation n/a -> 0.00438618; per_q2.2.0.reference.gamma n/a -> 0.174876; per_q2.2.0.reference.nf n/a -> 4; per_q2.2.0.reference.q2 n/a -> 2; per_q2.2.0.verdict.distance n/a -> 0.0337965; per_q2.2.0.verdict.distance_rel n/a -> 0.19326; per_q2.2.0.verdict.last_decade_step n/a -> 0.00378394; per_q2.5.0.reference.alpha_s n/a -> 0.285006; per_q2.5.0.reference.c_bj n/a -> 0.87219; per_q2.5.0.reference.err n/a -> 0.00232788; per_q2.5.0.reference.err_alpha_s n/a -> 0.00074298; per_q2.5.0.reference.err_g_a n/a -> 0.000188975; per_q2.5.0.reference.err_truncation n/a -> 0.00219802; per_q2.5.0.reference.gamma n/a -> 0.185384; per_q2.5.0.reference.nf n/a -> 4; per_q2.5.0.reference.q2 n/a -> 5; per_q2.5.0.verdict.distance n/a -> 0.0230361; per_q2.5.0.verdict.distance_rel n/a -> 0.124262; per_q2.5.0.verdict.last_decade_step n/a -> 0.00408475
compare: 41 changes (the M2 rule: any change exits 1)
```
