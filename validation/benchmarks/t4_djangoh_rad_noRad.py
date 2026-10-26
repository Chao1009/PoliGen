#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T4 row 7 (BENCHMARK_PLAN.md sec. 4): DJANGOH's ELASTIC RADIATIVE TAIL,
isolated, against LiPolGen's radiative tail, per Q^2 bin, on a PROTON at
e 18 x p 275 GeV -- with the ePIC Rad/noRad four-bin table kept as the
RECORDED HISTORY of why this row was BLOCKED from 2026-09-23 to 2026-10-26.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), FOR THIS ROW
----------------------------------------------------------------------
1. b1.  Neither side carries a b1.  DJANGOH is unpolarised e-p; LiPolGen runs
   an UNPOLARISED proton (spin 1/2, `helicity_flip_plan(0.5, 0, 0)`: P_z = 0,
   P_e = 0), where every rank-2 slot is identically zero.  A = 1 on both
   sides, so no A-scaling assumption enters.
2. Nuclear PDF grids / the isoscalar denominator.  None.  The compared
   quantity -- the elastic radiative tail e p -> e p gamma -- carries NO
   parton distribution on either side: it is the elastic vertex (the proton's
   Sachs form factors, the dipole on both sides) times lepton-line
   bremsstrahlung.  PDFs enter only the DIS Born used to DISPLAY the tail as
   a fraction: DJANGOH's Born uses the free-proton CTEQ6.1M (LHAPDF set
   10150, the ePIC samples' set); LiPolGen's own Born is `InclusiveKernel::
   dsigma_unpol` with the shipped `UnpolSfSource.Toy`.  There is no
   average-nucleon grid to divide by; the harness asserts A == 1.
3. Cross sections vs structure functions; per what.  Both tails are ABSOLUTE
   cross sections in microbarn per e-p collision (A = 1: per nucleus = per
   nucleon), integrated over the same window: DJANGOH's from its integration
   output (nb), LiPolGen's as INT Born_LPG(x, Q^2) r(x, Q^2) over the window
   with r = `RcModel.tail_ratio_at` = sigma_tail/sigma_Born on LiPolGen's OWN
   Born -- i.e. LiPolGen's absolute tail as the generator applies it.  The
   verdict is on R = sigma_tail(LiPolGen)/sigma_tail(DJANGOH), which needs no
   Born at all.  The "elastic-tail fraction" the figure draws divides BOTH
   tails by ONE Born -- DJANGOH's Rad=0 at the same cuts -- so
   f_LiPolGen/f_DJANGOH = R exactly.  LiPolGen's fraction on its OWN Born is
   recorded beside it and never compared: the Toy-vs-CTEQ6.1M Born ratio
   (1.05-1.41 on the ePIC window, 2026-09-23) is a structure-function
   difference, not a radiative correction.

THE DESIGN AND THE TOLERANCE -- fixed 2026-10-26 in this header before ANY
DJANGOH run, the control included, and before any LiPolGen number in the new
window was computed (PLAN.md ground rule 4).  This block, between its two
marker lines, is hashed into the vendored file (`design_block_sha256`).
==== DESIGN BLOCK BEGIN ====
SOURCE.  DJANGOH 4.6.21 = HERACLES 4.6.21 ("VERSION 4.6.21, 05.12.2022", the
banner of every ePIC log) from github.com/spiesber/DJANGOH at e356b84 -- the
file djangoh_h-4.6.21.f (that repository's djangoh_h.f is 4.6.20) with
djangoh_l/u/p/t.f, sophia-dj.f, pythia-6.4.24-dj.f, jetset7409.f,
gmc_random.f and cernlibs.f -- built with gfortran 11.4.0 against LHAPDF
6.5.5, and the manual djangoh_m-4.6.21.pdf.  Nothing of it is vendored; the
numbers are, with every input card verbatim.
(a) IEL2 = 0 in every run.  IEL2 (INT2(3)) integrates HSEL22, elastic
    e p -> e p AT x_l = 1 with its virtual and soft-photon corrections below
    the infrared cutoff XIRDEL = 2e-2 (manual sec. 15.13; HSDELX): the
    NON-radiative elastic peak at W_l = M_p, outside every W_l^2 >= 9 GeV^2
    window.  It is not a tail (and HSPRLG zeroes it anyway whenever
    XMAX < 1 - m_e^2/S - XIRDEL).
(b) The tail is isolated by option (i): INT-OPT-NC 0 0 0 0 0 0 N N N --
    INC2 = INC31..INC34 = 0 and IEL31 = IEL32 = IEL33 = N VEGAS iterations.
    The manual defines each channel's integration separately and the MAIN
    program integrates exactly the channels whose flag is set (its "ELASTIC
    TAIL" blocks), so (i) is allowed and needs no on/off pair: the summary
    SIGTW is a plain sum over channels, so option (ii)'s difference would be
    this sum identically.  IEL31, IEL32, IEL33 are three phase-space mappings
    -- the 1/(2k.l), 1/(2k.l') and 1/t_hadronic peaks: HSELK1, HSELK2,
    HSELCO -- of ONE squared matrix element, lepton-line bremsstrahlung
    e p -> e p gamma with the elastic F1, F2 and the lepton mass kept, so
    only their SUM is physical: the harness compares the sum and prints the
    three only as a record.  The Born behind the displayed fraction is the
    ePIC Rad=0 card (GSW-PARAM 2 0 3 1 0 0 2 1 1 1 1, INT-OPT-NC
    1 0 0 0 0 0 0 0 0) with only KINEM-CUTS changed to the window below.
(c) WMIN = 0.0 GeV on input, which HSPRLG raises to M_p = 0.93828 GeV
    (ICUT = 3: "IF(WMIN.LT.MPRO) WMIN=MPRO"); the elastic channels are zeroed
    only when WMIN > MPRO, so the tail's elastic vertex (W_h = M_p) is inside
    the acceptance.  THE LEPTONIC W.  KINEM-CUTS cuts x_l, y_l, Q^2_l and the
    HADRONIC W_h only, and W_h = M_p on every elastic-tail event, so
    integration alone cannot apply LiPolGen's W_l^2 >= 9 GeV^2.  The window
    is therefore matched BY THE CUT VARIABLES THE MANUAL PROVIDES: per bin
    [Q2_lo, Q2_hi], XMAX_b = Q2_lo/(Q2_lo + 9 GeV^2 - M_p^2) rounded DOWN to
    4 decimals and capped at 0.95 -- 0.1096 / 0.5518 / 0.9249 / 0.95 -- so
    that every point with x_l <= XMAX_b and Q^2_l >= Q2_lo has W_l^2 >= 9
    GeV^2.  The compared window is the rectangle
        Q2_lo <= Q^2_l < Q2_hi,  1e-5 <= x_l <= XMAX_b,  1e-4 <= y_l <= 0.985
    on BOTH sides: y_l <= 0.985 is LiPolGen's tail-table ceiling
    (`Scenario.y_max`), now applied to DJANGOH as well through YMAX, and
    LiPolGen integrates the identical rectangle (its own W^2 >= 9 cut is then
    implied).  HOW IT DIFFERS from the generator's W_l^2 >= 9 window (the
    history rows' window): it drops the strip XMAX_b < x_l <=
    Q^2_l/(Q^2_l + 9 - M_p^2) at Q^2_l > Q2_lo -- high x, low y, where the
    tail fraction is smallest -- and, in bin 4 only, x_l > 0.95: the cap
    keeps clear of DJANGOH's soft boundary x_l ~ 1 - XIRDEL = 0.98 (HSXMAX),
    above which IEL31..33 hand the tail to IEL2.  The harness MEASURES that
    strip on the LiPolGen side (its share of the Born and of the tail):
    recorded, never in the verdict.
    Q2MAX IS NOT APPLIED by the elastic-tail integrands: HSELK1, HSELK2 and
    HSELCO read Q2MIN, XMIN/XMAX and YMIN/YMAX only (HSPRLG's XMAX <=
    Q2MAX/(YMIN S) is a consistency bound, not a cut).  So every DJANGOH
    number per bin is a DIFFERENCE of two runs whose cards differ only in
    Q2MIN = Q2_lo or Q2_hi, both at Q2MAX = 1e8 (inactive):
        sigma_b = sigma(Q^2_l >= Q2_lo) - sigma(Q^2_l >= Q2_hi),
    integration errors added in quadrature.  The Born is built the same way,
    so the two windows are identical by construction; one direct Born run at
    Q2MAX = Q2_hi per bin is a recorded cross-check of the Born's own Q2MAX.
(d) FIXED IN ADVANCE.  PDF: cteq61 = LHAPDF6 set index 10150 (ICODE 10150,
    the ePIC samples' CTEQ6.1M central member), fetched 2026-10-26 from
    lhapdfsets.web.cern.ch/current/cteq61.tar.gz and kept outside the
    repository (its sha256 in the vendored file); the ePIC samples read
    LHAPDF 5.9.1's cteq61.LHgrid -- an interpolation difference, recorded.
    Seeds: RNDM-SEEDS 0 0 in every run (RANLUX with the code's fixed seed
    1213141567, no seed file); the ePIC runs used -1 (date and time), so
    their VEGAS numbers can be reproduced only to their errors.
    Integration: INT-POINTS 100000 (ePIC: 30000).  The control keeps the
    ePIC 18 iterations per radiative channel and INT-ONLY = 1; the isolating
    runs use 20 iterations per IEL channel and INT-ONLY = -1 (integration
    only: HSESTM's maxima search is skipped, which changes no integral).
    START 0 everywhere: no event sampling -- the integrated cross sections
    are what is compared.  RULE, FIXED NOW: a bin whose DJANGOH tail has a
    relative integration error above 1 % is rerun at INT-POINTS 400000 and
    the rerun replaces it; the task's floor (the error of each ratio below
    10 % of its distance from 1) is then met with room to spare.
    Running alpha: the tail runs use GSW-PARAM 2 1 2 1 1 1 2 1 1 1 1 -- the
    ePIC Rad=1 card with LPARIN(3) = 2 instead of 3.  LPARIN(3) >= 3
    multiplies the tail integrand by RUNALP^2 = [1 + Re Sigma_gg(t)/t]^-2 at
    the hadronic vertex (HSELK1, HSELK2, HSELCO), and LiPolGen's tail has a
    fixed alpha = 1/137.036 and no vacuum polarisation, so LPARIN(3) = 2 is
    the like-for-like run; nothing else in the photon-exchange elastic tail
    reads LPARIN(3).  The ePIC value (3) is run too, as a RECORDED sub-row,
    never in the verdict.
    THE CONTROL, FIRST: the eight ePIC cards (runcards/ep.Rad={0,1}.NC.Q2_*
    at eic/InclusiveDjangohSamples 9869d9a) verbatim except RNDM-SEEDS 0 0,
    INT-POINTS 100000 and START 0 -- IEL all 0, WMIN = 3, INT-ONLY = 1.
    "Reproduced" iff each of the eight numbers agrees with the README value
    within 3 x (ePIC's printed relative error (+) ours), the ePIC errors
    being the logs' "Estimated relative error" (Rad=1: 0.43 / 0.47 / 0.45
    e-3 in bins 2-4, and 0.47e-3, the largest, for bin 1, whose log is
    truncated; Rad=0: 1.0e-4, the code's EPSO).  The agreement is recorded
    per number whatever it is; a control that does not reproduce does not
    change the verdict rule -- it is reported with its size and its likely
    cause (LHAPDF 5.9.1 vs 6.5.5 interpolation; bin 1 reaches below the
    set's Q_min = 1.3 GeV).
(e) THE TOLERANCE.  Per bin b, R_b = sigma_tail(LiPolGen model)/
    sigma_tail(DJANGOH IEL31+IEL32+IEL33 at LPARIN(3) = 2), both absolute.
    The accounting, term by term:
    DJANGOH: the exact O(alpha) lepton-line bremsstrahlung e p -> e p gamma
      -- no peaking approximation, lepton mass kept, all three peak regions
      (manual sec. 7 and its Ref. [38], Bucoveanu & Spiesberger, EPJ A 55
      (2019) 57; the matrix element in HSELK1/HSELK2/HSELCO) -- with the
      dipole form factors (IDIPOL = 0: G_E = (1 + Q^2/0.71)^-2,
      G_M = 2.7927 G_E, HSFIE0) and a fixed alpha = 1/137.0359895.  No
      proton-line radiation, no two-photon exchange.
    LiPolGen `polrad-full` (RcTailModel.PolradFull, opt-in): POLRAD Eq. (18),
      the exact tau_A quadrature with all three peaks and the lepton mass,
      applied to the proton through the Eq. (44) quasi-elastic slot at
      Z = 1, N = 0, qe_kf_gev = 0 (`polrad_full_sigma_qe_u`: Im_1 =
      eta G_M^2, Im_2 = (G_E^2 + eta G_M^2)/(1 + eta), `nucleon_ff`'s dipole
      with mu_p = 2.7928473) and a fixed alpha = 1/137.036 -- THE SAME TERM
      as DJANGOH's.  What differs: G_M^2 by 1.1e-4 (mu_p), alpha^3 by 2e-7,
      M_p by 1e-5, (hbar c)^2 by 5.1e-5 (LiPolGen 0.3894e9 pb GeV^2,
      DJANGOH 3.8938e5 nb GeV^2) -- together < 2e-4 -- and the NUMERICS:
      LiPolGen evaluates the tail on an (x, y) table (log-bilinear
      between the sampler's x edges and >= 24
      log-spaced y rows plus 0.9 / 0.95 / 0.97 / 0.98 / 0.985) and CLIPS r
      at `RcOptions.tail_max` = 10, neither of which any test had priced
      against an external exact tail before this row (rc.hpp,
      RcTailModel::PolradFull point (3)); DJANGOH has its VEGAS error.
    So the leading-log vs exact difference (~5-10 %, rc.hpp's accuracy
    statement for `t-peak+ll`, which this row does not wire) is NOT the
    allowance for this pair: it is exact against exact, and the manual's
    sec. 7 and the source say DJANGOH's tail is the exact expression.  The
    allowance is numerical -- 2 % for LiPolGen's table and clip (chosen, not
    measured), plus the constants (< 2e-4), rounded up to 3 % -- plus twice
    DJANGOH's relative integration error delta_b of that bin, carried
    explicitly:
        PASS(polrad-full)  iff  |R_b - 1| <= 0.03 + 2 delta_b  in EACH bin.
    LiPolGen `t-peak` (RcTailModel.TPeak, the SHIPPED default): POLRAD
      Eq. (38)'s spin-1/2 line through Eq. (44) -- the t-peak ALONE,
      ultrarelativistic, no lepton mass -- which the tree itself calls a
      LOWER BOUND, "one peak of a three-peak object" (rc.hpp).  A two-sided
      tolerance on it would test a claim the tree does not make; the claim
      it does make is tested:
        PASS(t-peak)  iff  R_b <= 1 + 0.03 + 2 delta_b  in EACH bin
      (it must not exceed the exact tail beyond the same numerical
      allowance); its shortfall 1 - R_b is RECORDED as the price of the
      shipped default on this window, never asserted.
    The row PASSES iff both hold in all four bins; otherwise it is a
    recorded FAIL with its distance; it is BLOCKED only if no isolating
    configuration integrates.  The DJANGOH LPARIN(3) = 3 (running alpha)
    variant and LiPolGen's own-Born fraction are recorded sub-rows, outside
    the verdict.
==== DESIGN BLOCK END ====

ADDENDA (2026-10-26, after an independent adversarial verification of the
row; the block above is hashed and stays as written):
  * (b)'s "three phase-space mappings of ONE squared matrix element" is
    loose.  In the source IEL31, IEL32 and IEL33 are the 1/A1, 1/A2 and
    1/t_s partial-fraction pieces of that one squared matrix element
    (HSELK1, HSELK2, HSELCO), each integrated with its own mapping, the
    interference split through the 1/(A1 + A2) terms.  "Only their sum is
    physical" stands.
  * HSELK2 (djangoh_h-4.6.21.f, lines 17470-17474) zeroes the FSR piece at
    every point whose hadronic-t range exceeds GS00 = 20 GeV^2 -- an
    undocumented gate, and the cause of IEL32 = 0 at Q2MIN >= 100 GeV^2
    (not "the dipole is dead").  A rebuild with GS00 = 1e5 changes IEL32 by
    +0.0146 / +3.0e-4 / +2.8e-6 / +2.6e-8 nb, i.e. +6.2e-4 / +2.7e-4 / +5e-5
    / +1.7e-5 of the bin tail: below every allowance; not vendored; the
    verdict unchanged.
  * (d)'s fallback "rerun at INT-POINTS 400000" cannot run on this build:
    VEGAS's grid array is dimensioned 50 and the run aborts at 400000 points
    (200000 runs).  It was never triggered (delta_b <= 4.8e-4).
  * TOL_NUMERICAL and N_DELTA (the module constants the verdict uses) sat
    outside the hashed block; since this addendum `design_block_sha256()`
    also READS the tolerance back from the block's two "PASS(...) iff"
    sentences (`tolerance_from_design_block`) and refuses to run (exit 2)
    if the constants differ from it.  The REPORT line's tolerance text is
    built from the constants.

WHAT WAS MEASURED (2026-10-26; docs/open_items/run_2026-10-26/
phase_B3_djangoh.md has every command and runtime)
-----------------------------------------------------------------------------
Vendored in `data/djangoh_eltail_ep18x275.json` (40 runs, every card
verbatim, every log's sha256).  Matched window per bin, nb:
                  Born (Rad=0)   tail IEL31+32+33   fraction   rel. error
    [1, 10]        574.056        23.3809            0.040729   3.7e-4
    [10, 100]       67.6268        1.11147           0.016435   4.8e-4
    [100, 1000]      3.22611       0.0551624         0.017099   4.4e-4
    [1000, 10000]    0.0693644     0.00151826        0.021888   3.9e-4
THE VERDICT, by the rule above: R(polrad-full) = 1.0199 / 1.0261 / 1.0219 /
1.0098 against allowances 0.0307 / 0.0310 / 0.0309 / 0.0308 -- PASS in all
four bins; R(t-peak) = 0.6118 / 0.8968 / 1.0112 / 1.0067 <= 1.0307 /
1.0310 / 1.0309 / 1.0308 -- PASS (shortfall 38.8 % / 10.3 % / -1.1 % /
-0.7 %).  The row PASSES.
RECORDED BESIDE IT, never in the verdict:
  * WHERE THE 2 % COMES FROM (`--diagnostics`, ~110 s).  Straight from the
    quadrature (`RcModel.tail_sigma_at`, no table, no clip) the polrad-full
    proton tail is R = 1.0004 / 1.0010 / 0.9995 / 0.9866 of DJANGOH's: the
    two exact calculations agree to <= 0.1 % in bins 1-3, and bin 4's
    -1.34 % (34 DJANGOH sigma) is NOT numerical on either side as far as
    tested (verification 2026-10-26: LiPolGen's eta quadrature at n_eta 256
    and the direct route with its panels doubled move it < 1e-6; DJANGOH at
    200000 points moves it 2e-4; HSELK2's gate 1.7e-5) -- cause unknown,
    recorded.  The shipped (x, y) TABLE --
    what the generator weights events with -- lifts the window integral by
    +1.9 / +2.5 / +2.2 / +2.3 % (polrad-full) and +2.1 / +2.2 / +2.2 /
    +2.2 % (t-peak), so the verdict's polrad-full margin is mostly the table
    and the t-peak's R > 1 in bins 3-4 is the table too (direct: 0.9895,
    0.9847, below the exact tail, as a lower bound should be).  It is the
    TAIL column's log-bilinear interpolation between the nodes (x1.020 ...
    x1.037), not the Born's (<= 0.1 %); in bin 4 the `tail_max` = 10 clip
    takes 1.3 % of polrad-full's back.  The 2 % allowance for the table was
    chosen before this was known; the table's measured shift is 1.9-2.5 %,
    which EXCEEDS that 2 % term in bins 2-4: those bins pass through the
    rounding up to 3 %, and bin 2 (1.0261 against 1.0310) nearly used all
    of it.  The table route's quadrature moves <= 1.8e-4 when its panels
    are doubled.
  * DJANGOH with running alpha (LPARIN(3) = 3, the ePIC value): the tail
    rises by x1.0152 / 1.0133 / 1.0151 / 1.0206 -- vacuum polarisation at
    the hadronic vertex, which LiPolGen's tail does not carry.
  * The window: the strip of the W_l^2 >= 9 window outside the rectangle
    carries 4.7 % / 0.28 % / 4e-6 / 5e-6 of LiPolGen's Born and <= 5.2e-3
    of its tail (polrad-full, bin 1; <= 1e-5 for t-peak).
  * Q2MAX: a tail run at Q2MAX = Q2_hi returns the Q2MAX = 1e8 integral to
    every digit (ratio 1 exactly, four bins) -- the elastic-tail channels
    ignore it, as the source says; a Born run at Q2MAX = Q2_hi agrees with
    the subtraction to <= 3e-5 (its 5-digit print).
  * The control: all eight ePIC numbers reproduce within the declared
    3-sigma rule (ratio - 1 from -2.6e-4 to +4.0e-4, |pull| <= 1.74; the
    four Born numbers sit +0.2e-4 ... +2.5e-4 high, the LHAPDF 5 -> 6
    interpolation's size).  So this build reproduces the eight ePIC
    integrated cross sections (Born and inelastic channels) to <= 4e-4;
    the elastic-tail channels IEL31..33, which carry the verdict, have no
    external control (every ePIC card runs them off).
  * LiPolGen's OWN Born in the matched window is x1.469 / 1.116 / 1.155 /
    1.065 DJANGOH's (Toy SF vs CTEQ6.1M), so its own-Born fraction differs
    from the compared one by that factor; the 2026-09-23 history numbers
    (below) reproduce on this build to every printed digit.

WHAT THIS ROW DOES NOT VALIDATE.  The COHERENT spin-1 tail (Eqs. (37)-(39)
or (18) with the nuclear form factor: zero on a proton, so the 6Li coherent
tail is untested here); the Pauli-suppressed quasi-elastic tail
(`qe_kf_gev` = 0 here; the shipped 6Li value 0.169 GeV is not exercised);
any tensor or polarised tail; the INELASTIC O(alpha) correction (LiPolGen
applies none -- the ePIC history is its size, 7-23 %); radiation from the
proton line and two-photon exchange (neither side carries them); the
running of alpha (recorded above); the form factor itself (both sides use
the same dipole -- a common-mode input, not a check of it); Q^2 < 1 GeV^2,
y > 0.985, x_l > XMAX_b (the strip, bin 4's x_l > 0.95); the generator's
event sampling and acceptance (this is a window integral of the library's
weight, not a generated sample); `t-peak+ll`; other beams, and 6Li or 7Li.

RECORDED HISTORY: THE ePIC Rad/noRad TABLE, AND WHY IT COULD NOT BE THE
BENCHMARK (BLOCKED 2026-09-23 -> 2026-10-26)
-----------------------------------------------------------------------------
Reference then: sigma(Rad=1)/sigma(Rad=0) of DJANGOH 4.6.21, e- 18 x p 275,
per Q^2 bin [1,10], [10,100], [100,1000], [1000,10000] GeV^2, x_l in
[1e-5, 1], y_l in [1e-4, 1], W_h >= 3 GeV (ICUT = 3) -- vendored in
`data/djangoh_rad_norad_ep18x275.json` (eic/InclusiveDjangohSamples at
9869d9a): 1.072454 / 1.089407 / 1.160626 / 1.227965.  LiPolGen's
1 + INT Born r / INT Born on its W_l^2 >= 9 window (the LEPTONIC W:
the same number as the samples' W_h cut, not the same cut):
t-peak 1.016169 / 1.013172 / 1.014970 / 1.020686, polrad-full 1.027096 /
1.015072 / 1.015138 / 1.020757 (docs/open_items/run_2026-09-23/
phase_B3_chain_rc.md sec. 1; reproduced 2026-10-26).  Why no tolerance
existed: DJANGOH Rad=1 carried the non-radiative O(alpha) piece (INC2 with
virtual + soft), the inelastic lepton-line radiation and the Compton term
(INC31..INC33) and NOT the elastic radiative tail -- all eight ePIC logs
read IEL2 = IEL31 = IEL32 = IEL33 = 0 -- while LiPolGen's `rc_tail` is the
elastic (+ quasi-elastic) tail and nothing else (rc.hpp: the lepton-vertex
piece is priced by the tensor BAND, never applied).  Disjoint additive
pieces of one O(alpha) correction: no value of either could contradict the
other.  And the samples' W_h >= 3 GeV cut excluded the tail whatever IEL
was: W_h = M_p in e p -> e p gamma, and HSPRLG zeroes IEL2 and IEL31..33
whenever WMIN > M_p.  The row therefore needed a run at IEL31..33 != 0 AND
WMIN <= M_p -- the run made 2026-10-26 and designed in the block above.
The ePIC ratio is kept in the report as a recorded sub-row and in the
figure as a grey series: it is the size of the INELASTIC correction
LiPolGen does not apply, never the comparison.

Run:  source env.sh && python3 validation/benchmarks/t4_djangoh_rad_noRad.py
      [--plot]          also writes docs/benchmarking/figures/<NAME>.png
      [--diagnostics]   the table-vs-quadrature record instead (~110 s)
(22-26 s measured 2026-10-26 on a 4-core machine shared with another
worker, the polrad-full tail table 21 s of it there -- rc.hpp says ~10 s
unloaded; not opt-in: no event sample, under 30 s.)

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a design block that no longer hashes to the recorded
digest, a vendored per-bin number that the raw run records do not
reproduce, a missing dependency, any exception).  A harness is a
measurement, not a CI gate: the pytests are the gate.  Both vendored files'
sha256 (every byte) are re-checked on EVERY read, and a missing numpy exits
2 like any other missing dependency.
"""

import hashlib
import json
import math
import pathlib
import re
import sys
import traceback
import time

try:
    import numpy as np
except ImportError:  # a missing dependency is a broken harness: exit 2
    if __name__ == "__main__":
        traceback.print_exc()
        sys.exit(2)
    raise

HERE = pathlib.Path(__file__).resolve().parent
#: the ePIC Rad/noRad table -- the RECORDED HISTORY of this row
DATA = HERE / "data" / "djangoh_rad_norad_ep18x275.json"
#: sha256 of every byte of the ePIC file (recorded 2026-09-26 from the file
#: as committed at f0a8f1e); re-checked on every read.
DATA_SHA256 = "2fd4b9b6bb2c0de7c7a9b338ad58cc3cbc221c768d0ce2f257ac6f57139be8fa"
#: the isolated elastic tail and its Born, from the DJANGOH runs of
#: 2026-10-26 -- the PRIMARY comparison
ELTAIL_DATA = HERE / "data" / "djangoh_eltail_ep18x275.json"
#: sha256 of every byte of that file, as vendored 2026-10-26; re-checked on
#: every read.
ELTAIL_SHA256 = "8c0bd1ee4bb176db9e3d1601c43413ba41bd0d93581d11d20c29478ae05731c7"
#: sha256 of this header's DESIGN BLOCK (its two marker lines included) as
#: written at 2026-10-26T19:11:53Z, before any DJANGOH run; the vendored
#: file records the same digest.  `design_block_sha256()` re-checks both on
#: every run, so the design and the tolerance cannot move after the numbers
#: without the harness refusing to run.
DESIGN_BLOCK_SHA256 = (
    "fbe067892c7d3b8d6de27f32ce1c52138bd80162b24a7b1677dc6d412af2c91b")
DESIGN_BEGIN = "==== DESIGN BLOCK BEGIN ===="
DESIGN_END = "==== DESIGN BLOCK END ===="

NAME = "t4_djangoh_rad_noRad"
#: LiPolGen beam configuration index for e 18 GeV x p 275 GeV/u.
CONFIG = 2
#: The leptonic cuts shared by every window: DJANGOH's XMIN, YMIN, and the
#: history window's W_l^2 floor (the ePIC samples' WMIN = 3 GeV squared --
#: there on the HADRONIC W_h, here on the LEPTONIC W).
X_MIN = 1.0e-5
Y_MIN = 1.0e-4
W2_MIN = 9.0
#: The matched window's x ceilings per Q^2 bin (design block (c)):
#: Q2_lo/(Q2_lo + 9 - M_p^2) rounded DOWN to 4 decimals, capped at 0.95.
XMAX_DECIMALS = 4
XMAX_CAP = 0.95
#: The tolerance (design block (e)): a 3 % numerical allowance plus twice
#: DJANGOH's relative integration error of the bin.
TOL_NUMERICAL = 0.03
N_DELTA = 2.0
#: Gauss-Legendre nodes per panel, and panels, in ln Q^2 and in each of the
#: two y pieces (ln y below Y_SPLIT, -ln(1 - y) above it).  Converged: see
#: `convergence()`.
N_GL = 12
PANELS_Q2 = 4
PANELS_Y = 6
Y_SPLIT = 0.5
TAIL_MODELS = ("t-peak", "polrad-full")
#: The DJANGOH elastic-tail channels, in the MAIN program's order.
EL_CHANNELS = ("IEL31", "IEL32", "IEL33")
#: Not opt-in: the row integrates vendored DJANGOH numbers against a library
#: evaluation (no event sample), 22-26 s measured 2026-10-26 on a shared
#: machine -- under the 30 s line of validation/benchmarks/README.md.
EXPENSIVE = False


# ------------------------------------------------------------ the references

def read_vendored(path, sha256=DATA_SHA256):
    """The bytes of a vendored file -- RuntimeError unless their sha256 is
    the recorded one, so a moved number is refused, never silently read."""
    raw = pathlib.Path(path).read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != sha256:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, sha256))
    return raw


def load_reference(path=DATA):
    """The four ePIC rows (the history) -> (doc, list of dicts) with the
    Rad/noRad ratio computed from the two columns."""
    doc = json.loads(read_vendored(path).decode("utf-8"))
    rows = []
    for lo, hi, rad1, rad0 in doc["rows"]:
        rows.append({"q2_lo": lo, "q2_hi": hi, "sigma_rad1_ub": rad1,
                     "sigma_rad0_ub": rad0, "ratio": rad1 / rad0})
    return doc, rows


_TOL_FULL_RE = re.compile(
    r"PASS\(polrad-full\)\s+iff\s+\|R_b - 1\| <= ([0-9.]+) \+ ([0-9.]+) delta_b")
_TOL_TPK_RE = re.compile(
    r"PASS\(t-peak\)\s+iff\s+R_b <= 1 \+ ([0-9.]+) \+ ([0-9.]+) delta_b")


def tolerance_from_design_block(block):
    """(tol_numerical, n_delta) as the hashed DESIGN BLOCK states them in its
    two "PASS(...) iff" sentences, which must agree -- RuntimeError
    otherwise.  The verdict's constants are checked against this on every
    run (verification of 2026-10-26: they sat outside the hashed text)."""
    m_full, m_tpk = _TOL_FULL_RE.search(block), _TOL_TPK_RE.search(block)
    if m_full is None or m_tpk is None or m_full.groups() != m_tpk.groups():
        raise RuntimeError("the design block's two PASS(...) iff sentences "
                           "do not state one tolerance")
    return float(m_full.group(1)), float(m_full.group(2))


def design_block_sha256(source=None):
    """sha256 of the DESIGN BLOCK of this file's header (both marker lines
    included) -- RuntimeError unless it is `DESIGN_BLOCK_SHA256` AND the
    block's stated tolerance equals (TOL_NUMERICAL, N_DELTA)."""
    text = (pathlib.Path(__file__).read_text(encoding="utf-8")
            if source is None else source)
    i0 = text.index(DESIGN_BEGIN)
    i1 = text.index(DESIGN_END) + len(DESIGN_END)
    got = hashlib.sha256(text[i0:i1].encode("utf-8")).hexdigest()
    if got != DESIGN_BLOCK_SHA256:
        raise RuntimeError("the header's design block sha256 %s != recorded "
                           "%s: the design or the tolerance was edited"
                           % (got, DESIGN_BLOCK_SHA256))
    stated = tolerance_from_design_block(text[i0:i1])
    if stated != (TOL_NUMERICAL, N_DELTA):
        raise RuntimeError("the tolerance constants (TOL_NUMERICAL, N_DELTA) "
                           "= %r differ from the hashed design block's %r: "
                           "the tolerance constants were edited"
                           % ((TOL_NUMERICAL, N_DELTA), stated))
    return got


def xmax_rule(q2_lo, m_p):
    """The matched window's x ceiling for a bin starting at `q2_lo`."""
    x = q2_lo / (q2_lo + W2_MIN - m_p * m_p)
    x = math.floor(x * 10 ** XMAX_DECIMALS) / 10 ** XMAX_DECIMALS
    return min(x, XMAX_CAP)


def _diff(runs, lo_name, hi_name, pick):
    a, b = pick(runs[lo_name]), pick(runs[hi_name])
    return a[0] - b[0], math.hypot(a[1], b[1])


def _born(run):
    return tuple(run["summary"]["sections"]["nc"][
        "VIRTUAL & SOFT CONTRIBUTIONS SIG2"])


def _channel(k):
    return lambda run: (run["vegas"][k]["integral_nb"],
                        run["vegas"][k]["std_dev_nb"])


def load_eltail(path=ELTAIL_DATA, m_p=None):
    """The vendored DJANGOH runs -> (doc, bins, control).

    Every per-bin number is RECOMPUTED here from the raw run records (the
    9-digit VEGAS integrals and the integration summary) as the difference
    of the Q2MIN = Q2_lo and Q2MIN = Q2_hi runs (design block (c)) and must
    equal the file's own `bins` entry; the window's x ceilings must equal
    the design rule; the design block digest must match.  Any mismatch is a
    RuntimeError (the harness exits 2)."""
    doc = json.loads(read_vendored(path, ELTAIL_SHA256).decode("utf-8"))
    own = design_block_sha256()
    if doc["provenance"]["design"]["design_block_sha256"] != own:
        raise RuntimeError("%s records design block %s, the harness header "
                           "hashes to %s" % (path, doc["provenance"]["design"]
                                             ["design_block_sha256"], own))
    runs = doc["runs"]
    bins = []
    for stored in doc["bins"]:
        lo, hi = stored["q2_lo"], stored["q2_hi"]
        tag = "%d_%d" % (lo, hi)
        if m_p is not None and stored["xmax"] != xmax_rule(lo, m_p):
            raise RuntimeError("bin %s: vendored xmax %r != the design rule %r"
                               % (tag, stored["xmax"], xmax_rule(lo, m_p)))
        b = {"q2_lo": lo, "q2_hi": hi, "xmax": stored["xmax"]}
        b["born_nb"], b["born_err_nb"] = _diff(
            runs, "born_%s_q2min_lo" % tag, "born_%s_q2min_hi" % tag, _born)
        b["born_direct_q2max_nb"] = _born(runs["born_%s_direct_q2max" % tag])[0]
        for var in ("novp", "vp"):
            chans = {}
            for k, ch in enumerate(EL_CHANNELS):
                chans[ch] = _diff(runs, "tail_%s_%s_q2min_lo" % (var, tag),
                                  "tail_%s_%s_q2min_hi" % (var, tag),
                                  _channel(k))
            b["tail_%s_nb" % var] = sum(v[0] for v in chans.values())
            b["tail_%s_err_nb" % var] = math.sqrt(
                sum(v[1] ** 2 for v in chans.values()))
            b["tail_%s_channels_nb" % var] = {c: v[0] for c, v in chans.items()}
        b["tail_novp_q2min_lo_nb"] = sum(
            _channel(k)(runs["tail_novp_%s_q2min_lo" % tag])[0]
            for k in range(3))
        b["tail_novp_direct_q2max_nb"] = sum(
            _channel(k)(runs["tail_novp_%s_direct_q2max" % tag])[0]
            for k in range(3))
        for key in ("born_nb", "born_err_nb", "tail_novp_nb", "tail_novp_err_nb",
                    "tail_vp_nb", "tail_vp_err_nb", "born_direct_q2max_nb",
                    "tail_novp_direct_q2max_nb"):
            if not math.isclose(b[key], stored[key], rel_tol=1e-12):
                raise RuntimeError("bin %s: %s recomputed %r != vendored %r"
                                   % (tag, key, b[key], stored[key]))
        b["delta"] = b["tail_novp_err_nb"] / b["tail_novp_nb"]
        bins.append(b)
    return doc, bins, doc["control"]


# ------------------------------------------------------------ LiPolGen side

def _gl(n):
    t, w = np.polynomial.legendre.leggauss(n)
    return t, w


def _panel_nodes(a, b, panels, n):
    """Composite Gauss-Legendre nodes/weights on [a, b]."""
    t, w = _gl(n)
    edges = np.linspace(a, b, panels + 1)
    nodes, weights = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        h = 0.5 * (hi - lo)
        nodes.append(lo + h * (t + 1.0))
        weights.append(h * w)
    return np.concatenate(nodes), np.concatenate(weights)


def build_pipeline(tail_model, qe_kf_gev=0.0):
    """LiPolGen's proton, RC on, at 18 x 275 with DJANGOH's x/y floors.

    The ONLY non-default RC settings: a zero coherent form factor (the proton
    has no spin-1 coherent vertex; `HoSpin1FF::for_ion` throws for it) and,
    unless `qe_kf_gev` is passed, no Pauli factor (a free proton).  The
    `tail_model` is the knob under test.  Nothing here is a physics default
    of the library; this is a harness configuration."""
    import lipolgen as lg
    from lipolgen import _lipolgen as _l

    class _NoCoherentVertex(_l.Spin1ElasticFF):
        def fc(self, t):
            return 0.0

        def fm(self, t):
            return 0.0

        def fq(self, t):
            return 0.0

        def provenance(self):
            return ("t4_djangoh_rad_noRad: zero coherent spin-1 vertex (a "
                    "proton has none); the proton elastic tail is the "
                    "Eq. (44) Z = 1, N = 0 slot")

    cfg = lg.make_config(isotope="p", channel="inclusive", config=CONFIG,
                         events=10, seed=1, rc="tensor-band",
                         rc_tail_model=tail_model)
    g = cfg.grid
    g.x_min = X_MIN
    g.q2_min = 1.0
    g.q2_max = 1.0e4
    cfg.grid = g
    sc = cfg.scenario
    sc.w2_min = W2_MIN
    sc.y_min = Y_MIN
    sc.q2_min = 1.0
    # DJANGOH has no scattered-electron acceptance: open LiPolGen's.
    sc.eta_min, sc.eta_max, sc.e_prime_min = -50.0, 50.0, 0.0
    cfg.scenario = sc
    ff = _NoCoherentVertex()
    cfg.rc_options.ff = ff
    cfg.rc_options.qe_kf_gev = qe_kf_gev
    p = lg.Pipeline(cfg, lg.helicity_flip_plan(0.5, 0.0, 0.0))
    assert p.dis_sampler.kernel.ion.A == 1, "rule 2/3: A must be 1"
    # `ff` is returned so the caller holds the Python trampoline for as long
    # as the model lives.
    return p, cfg, ff


def _y_nodes(y_lo, y_hi):
    """Nodes/weights in y on [y_lo, y_hi]: ln y below Y_SPLIT (the Born's
    low-x peak) and -ln(1 - y) above it (the tail's 1/(1 - y) growth)."""
    ys, ws = [], []
    a, b = y_lo, min(y_hi, Y_SPLIT)
    if b > a:
        v, w = _panel_nodes(math.log(a), math.log(b), PANELS_Y, N_GL)
        ys.append(np.exp(v))
        ws.append(w * np.exp(v))
    a, b = max(y_lo, Y_SPLIT), y_hi
    if b > a:
        t, w = _panel_nodes(-math.log(1.0 - a), -math.log(1.0 - b),
                            PANELS_Y, N_GL)
        ys.append(1.0 - np.exp(-t))
        ws.append(w * np.exp(-t))
    if not ys:
        return np.zeros(0), np.zeros(0)
    return np.concatenate(ys), np.concatenate(ws)


def integrate_window(p, q2_lo, q2_hi, y_max, x_max=1.0, tail_max=None):
    """INT Born and INT Born * r over one window, in (ln Q^2, y) [pb].

    Born = d2sigma/(dx dQ2) [pb/GeV^2], so dx dQ2 = (Q^2/(s y^2)) dy Q^2
    dlnQ^2.  The cuts, all on the LEPTONIC variables: Q2_lo <= Q^2 <= Q2_hi,
    X_MIN <= x <= x_max, Y_MIN <= y <= y_max, W^2 >= W2_MIN.  With
    x_max = 1 this is the history window (W_l^2 >= 9); with the bin's
    XMAX_b it is the matched rectangle, the W cut then implied.  Also
    returned: the Born that falls outside the tail table (`uncovered`) and
    the Born and tail weight at nodes where r hit `tail_max` (clipped)."""
    import lipolgen as lg
    ds = p.dis_sampler
    k, m, s = ds.kernel, p.rc_model, ds.s
    m_p2 = lg.PROTON_MASS ** 2
    uq, wq = _panel_nodes(math.log(q2_lo), math.log(q2_hi), PANELS_Q2, N_GL)
    born = tail = uncovered = clip_born = clip_tail = 0.0
    for u, wu in zip(uq, wq):
        q2 = math.exp(u)
        y_lo = max(Y_MIN, q2 / (s * x_max), (W2_MIN - m_p2 + q2) / s)
        y_hi = min(y_max, q2 / (s * X_MIN))
        for y, wy in zip(*_y_nodes(y_lo, y_hi)):
            x = q2 / (s * y)
            b = (k.dsigma_unpol(x, q2, s) * q2 * q2 / (s * y * y)) * wu * wy
            born += b
            try:
                r = m.tail_ratio_at(x, q2, 0.0)
            except Exception:                 # outside the tail table
                uncovered += b
                continue
            tail += b * r
            if tail_max is not None and r >= tail_max:
                clip_born += b
                clip_tail += b * r
    return {"born_pb": float(born), "tail_pb": float(tail),
            "uncovered_pb": float(uncovered), "clipped_born_pb": float(clip_born),
            "clipped_tail_pb": float(clip_tail)}


def convergence(p, bins, y_max, factor=2):
    """The matched-window tails at the shipped node counts and at `factor` x
    the panels in both variables -> max |relative change| over the bins."""
    global PANELS_Q2, PANELS_Y
    base = (PANELS_Q2, PANELS_Y)
    vals = []
    try:
        for pq, py in (base, (base[0] * factor, base[1] * factor)):
            PANELS_Q2, PANELS_Y = pq, py
            vals.append([integrate_window(p, b["q2_lo"], b["q2_hi"], y_max,
                                          b["xmax"])["tail_pb"] for b in bins])
    finally:
        PANELS_Q2, PANELS_Y = base
    return max(abs(b / a - 1.0) for a, b in zip(*vals))


def lipolgen_side(tail_models=TAIL_MODELS, bins=None, with_convergence=False):
    """Per tail model: the matched-window tail and Born per bin, and the
    history window's (W_l^2 >= 9) tail and Born, both [microbarn]."""
    out = {}
    for tm in tail_models:
        t0 = time.time()
        p, cfg, _ff = build_pipeline(tm)
        build_s = time.time() - t0
        # The shipped Scenario.y_max, read and never retyped: it is the top of
        # the tail table's y support (`RcModel` refuses to extrapolate).
        ym = cfg.scenario.y_max
        tmax = cfg.rc_options.tail_max
        per_bin = []
        for b in bins:
            rect = integrate_window(p, b["q2_lo"], b["q2_hi"], ym, b["xmax"],
                                    tmax)
            hist = integrate_window(p, b["q2_lo"], b["q2_hi"], ym, 1.0, tmax)
            # pb -> microbarn
            per_bin.append({
                "tail_ub": rect["tail_pb"] * 1e-6,
                "born_ub": rect["born_pb"] * 1e-6,
                "uncovered_born_fraction": rect["uncovered_pb"] / rect["born_pb"],
                "clipped_born_fraction": rect["clipped_born_pb"] / rect["born_pb"],
                "clipped_tail_fraction": (rect["clipped_tail_pb"] / rect["tail_pb"]
                                          if rect["tail_pb"] > 0 else 0.0),
                "hist_tail_ub": hist["tail_pb"] * 1e-6,
                "hist_born_ub": hist["born_pb"] * 1e-6,
                # the 2026-09-23 definition, kept: 1 + INT Born r / INT Born
                # over the covered Born
                "hist_ratio": 1.0 + hist["tail_pb"] / (hist["born_pb"]
                                                       - hist["uncovered_pb"]),
                "hist_uncovered_born_fraction":
                    hist["uncovered_pb"] / hist["born_pb"],
            })
        conv = convergence(p, bins, ym) if with_convergence else None
        out[tm] = {"bins": per_bin, "build_s": build_s, "y_max": ym,
                   "tail_max": tmax,
                   "clipped_node_fraction": p.rc_model.clipped_cell_fraction,
                   "ff_provenance": p.rc_model.ff_provenance,
                   "unpol_sf": str(cfg.unpol_sf),
                   "quadrature_max_rel_change": conv}
    return out


def _table_sigma_qe(m, x, y):
    """The shipped table's quasi-elastic column at (x, y), interpolated
    exactly as `RcModel::tail_ratio_at` does (log-bilinear in ln x, ln y
    between the bracketing nodes) -- for the attribution in `diagnostics`
    only.  On this row's proton the coherent columns are zero, so the tail
    numerator is this column (times `qe_suppression`)."""
    def locate(nodes, v):
        k = int(np.searchsorted(nodes, v, side="right"))
        k = min(max(k, 1), len(nodes) - 1)
        lo, hi = math.log(nodes[k - 1]), math.log(nodes[k])
        f = (math.log(v) - lo) / (hi - lo) if hi > lo else 0.0
        return k - 1, min(max(f, 0.0), 1.0)
    tx_nodes, ty_nodes = m.table_x, m.table_y
    ix, tx = locate(tx_nodes, x)
    iy, ty = locate(ty_nodes, y)
    ny = len(ty_nodes)
    q = m.sigma_tail_qe
    corners = ((ix, iy, (1 - tx) * (1 - ty)), (ix, iy + 1, (1 - tx) * ty),
               (ix + 1, iy, tx * (1 - ty)), (ix + 1, iy + 1, tx * ty))
    return math.exp(sum(w * math.log(max(q[i * ny + j], 1e-300))
                        for i, j, w in corners))


def direct_window_tail(p, q2_lo, q2_hi, y_max, x_max, qe_suppression):
    """The same window integral as `integrate_window`'s tail, two ways
    [pb]: `direct`, the tail STRAIGHT FROM THE QUADRATURE at every node
    (`RcModel.tail_sigma_at`: no (x, y) table, no `tail_max` clip), and
    `sigma_interp`, the table's interpolated tail over the EXACT Born and
    without the clip (the table route with its Born interpolation and its
    `tail_max` clip removed).  A diagnostic of the table, never in the
    verdict: Born * r = (Q^2/(s y^2)) sigma_tail GEV2_TO_PB/(x s) per unit
    dy Q^2 dlnQ^2."""
    import lipolgen as lg
    ds = p.dis_sampler
    m, s = p.rc_model, ds.s
    m_p2 = lg.PROTON_MASS ** 2
    uq, wq = _panel_nodes(math.log(q2_lo), math.log(q2_hi), PANELS_Q2, N_GL)
    direct = interp = 0.0
    for u, wu in zip(uq, wq):
        q2 = math.exp(u)
        y_lo = max(Y_MIN, q2 / (s * x_max), (W2_MIN - m_p2 + q2) / s)
        y_hi = min(y_max, q2 / (s * X_MIN))
        for y, wy in zip(*_y_nodes(y_lo, y_hi)):
            x = q2 / (s * y)
            el_u, _el_t, qe, u_sp, qe_sp = m.tail_sigma_at(x, q2)
            sig = el_u + u_sp + qe_suppression * (qe + qe_sp)
            jac = (q2 * q2 / (s * y * y)) * wu * wy * lg.GEV2_TO_PB / (x * s)
            direct += jac * sig
            interp += jac * qe_suppression * _table_sigma_qe(m, x, y)
    return {"direct_pb": float(direct), "sigma_interp_pb": float(interp)}


def diagnostics(tail_models=TAIL_MODELS):
    """`python3 t4_djangoh_rad_noRad.py --diagnostics`: per model and bin,
    the matched-window tail through the shipped table (what the generator
    applies, the verdict's number) over the same integral straight from the
    quadrature, split into (the table's Born interpolation and its
    `tail_max` clip) x (its tail interpolation); and the table route's own
    convergence (panels doubled).
    ~110 s (the polrad-full quadrature is ~2 ms a node); a record of where
    the table moves the number, never in the verdict."""
    import lipolgen as lg
    _doc, bins, _control = load_eltail(m_p=lg.PROTON_MASS)
    out = {}
    for tm in tail_models:
        p, cfg, _ff = build_pipeline(tm)
        ym = cfg.scenario.y_max
        conv = convergence(p, bins, ym)
        print("DIAG | %s | matched-window tail, panels doubled in ln Q^2 and "
              "y: max relative change %.2e" % (tm, conv), flush=True)
        rows = []
        for b in bins:
            tab = integrate_window(p, b["q2_lo"], b["q2_hi"], ym, b["xmax"],
                                   cfg.rc_options.tail_max)["tail_pb"]
            dw = direct_window_tail(p, b["q2_lo"], b["q2_hi"], ym, b["xmax"],
                                    cfg.rc_options.qe_suppression)
            direct = dw["direct_pb"]
            dj = b["tail_novp_nb"] * 1e3
            rows.append({"bin": (b["q2_lo"], b["q2_hi"]),
                         "table_over_direct": tab / direct,
                         "born_interp_and_clip_factor":
                             tab / dw["sigma_interp_pb"],
                         "tail_interp_factor": dw["sigma_interp_pb"] / direct,
                         "R_table": tab / dj, "R_direct": direct / dj})
            print("DIAG | %s[%s] | table/direct %.5f = (Born interpolation "
                  "and clip) %.5f x tail interpolation %.5f | R(table) %.5f | "
                  "R(direct) %.5f" % (
                      tm, _bin_label(b), tab / direct,
                      tab / dw["sigma_interp_pb"],
                      dw["sigma_interp_pb"] / direct, tab / dj, direct / dj),
                  flush=True)
        out[tm] = {"bins": rows, "quadrature_max_rel_change": conv}
    return out


def log_slope(q2_mid, excess):
    """d ln(excess) / d ln Q^2 by least squares -- the Q^2 TREND, as a sign
    and a size, of (ratio - 1).  Used on the history rows only."""
    xs = np.log(np.asarray(q2_mid))
    ys = np.log(np.asarray(excess))
    return float(np.polyfit(xs, ys, 1)[0])


# ------------------------------------------------------------ the verdict

def verdict(dj_bins, lpg):
    """Per bin and model: R = sigma_tail(LiPolGen)/sigma_tail(DJANGOH), the
    allowance 0.03 + 2 delta_b, and the declared test (design block (e)):
    polrad-full two-sided, t-peak one-sided (the lower-bound claim)."""
    out = []
    for i, b in enumerate(dj_bins):
        tail_dj_ub = b["tail_novp_nb"] * 1e-3
        allow = TOL_NUMERICAL + N_DELTA * b["delta"]
        for tm in TAIL_MODELS:
            lb = lpg[tm]["bins"][i]
            r = lb["tail_ub"] / tail_dj_ub
            if tm == "polrad-full":
                ok = abs(r - 1.0) <= allow
                test = "abs(R - 1) <= %.4f" % allow
            else:
                ok = r <= 1.0 + allow
                test = "R <= %.4f (lower bound)" % (1.0 + allow)
            sane = math.isfinite(r) and r >= 0.0
            out.append({"bin": i, "model": tm, "R": r, "allowance": allow,
                        "test": test, "pass": bool(ok and sane), "sane": sane})
    return out


def _bin_label(b):
    return "Q2 %g-%g" % (b["q2_lo"], b["q2_hi"])


def run(tail_models=TAIL_MODELS, verbose=True):
    import lipolgen as lg
    if tuple(tail_models) != TAIL_MODELS:
        raise ValueError("the verdict is declared for %s exactly" % (TAIL_MODELS,))
    epic_doc, epic_rows = load_reference()
    doc, bins, control = load_eltail(m_p=lg.PROTON_MASS)
    gen = lipolgen_side(tail_models, bins)
    ver = verdict(bins, gen)
    status = "pass" if all(v["pass"] for v in ver) else "fail"

    detail, subrows = [], []
    for i, b in enumerate(bins):
        lab = _bin_label(b)
        born_dj_ub = b["born_nb"] * 1e-3
        f_dj = b["tail_novp_nb"] / b["born_nb"]
        d = {"bin": (b["q2_lo"], b["q2_hi"]), "xmax": b["xmax"],
             "djangoh_born_ub": born_dj_ub,
             "djangoh_tail_ub": b["tail_novp_nb"] * 1e-3,
             "djangoh_tail_rel_err": b["delta"],
             "djangoh_fraction": f_dj,
             "djangoh_fraction_err": f_dj * math.hypot(
                 b["delta"], b["born_err_nb"] / b["born_nb"]),
             "djangoh_channels_ub": {c: v * 1e-3 for c, v in
                                     b["tail_novp_channels_nb"].items()},
             "djangoh_fraction_running_alpha": b["tail_vp_nb"] / b["born_nb"],
             "running_alpha_over_fixed": b["tail_vp_nb"] / b["tail_novp_nb"],
             "epic_rad_over_norad": epic_rows[i]["ratio"],
             "allowance": TOL_NUMERICAL + N_DELTA * b["delta"]}
        for tm in tail_models:
            lb = gen[tm]["bins"][i]
            v = [x for x in ver if x["bin"] == i and x["model"] == tm][0]
            d[tm] = {"R": v["R"], "pass": v["pass"],
                     "tail_ub": lb["tail_ub"],
                     "fraction_on_djangoh_born": lb["tail_ub"] / born_dj_ub,
                     "fraction_on_own_born": lb["tail_ub"] / lb["born_ub"],
                     "own_born_ub": lb["born_ub"],
                     "strip_born_share": 1.0 - lb["born_ub"] / lb["hist_born_ub"],
                     "strip_tail_share": 1.0 - lb["tail_ub"] / lb["hist_tail_ub"],
                     "hist_ratio": lb["hist_ratio"],
                     "uncovered_born_fraction": lb["uncovered_born_fraction"],
                     "clipped_born_fraction": lb["clipped_born_fraction"],
                     "clipped_tail_fraction": lb["clipped_tail_fraction"]}
            sub = {"name": "%s[%s]" % (tm, lab),
                   "status": "pass" if v["pass"] else "fail",
                   "in_verdict": True,
                   "R": v["R"], "allowance": v["allowance"],
                   "test": v["test"],
                   "tail_lipolgen_ub": lb["tail_ub"],
                   "tail_djangoh_ub": b["tail_novp_nb"] * 1e-3,
                   "fraction_lipolgen": lb["tail_ub"] / born_dj_ub,
                   "fraction_djangoh": f_dj,
                   "fraction_lipolgen_own_born": lb["tail_ub"] / lb["born_ub"],
                   "clipped_tail_fraction": lb["clipped_tail_fraction"]}
            if tm == "t-peak":
                sub["shortfall_1_minus_R"] = 1.0 - v["R"]
            subrows.append(sub)
        subrows.append({
            "name": "recorded[%s]" % lab, "status": "recorded",
            "in_verdict": False,
            "djangoh_running_alpha_over_fixed": b["tail_vp_nb"] / b["tail_novp_nb"],
            "djangoh_fraction_running_alpha": b["tail_vp_nb"] / b["born_nb"],
            "strip_born_share": d["t-peak"]["strip_born_share"],
            "strip_tail_share_t-peak": d["t-peak"]["strip_tail_share"],
            "strip_tail_share_polrad-full": d["polrad-full"]["strip_tail_share"],
            "born_q2max_direct_over_subtraction":
                b["born_direct_q2max_nb"] / b["born_nb"],
            "tail_q2max_direct_over_q2min_only":
                b["tail_novp_direct_q2max_nb"] / b["tail_novp_q2min_lo_nb"],
            "lipolgen_own_born_over_djangoh_born":
                d["t-peak"]["own_born_ub"] / born_dj_ub})
        subrows.append({
            "name": "history_epic[%s]" % lab, "status": "recorded",
            "in_verdict": False,
            "reason": ("ePIC Rad/noRad: inelastic O(alpha) with the elastic "
                       "tail OFF and W_h >= 3 GeV; LiPolGen ratios on the "
                       "W_l^2 >= 9 window -- disjoint pieces, no tolerance"),
            "epic_rad_over_norad": epic_rows[i]["ratio"],
            "lipolgen_t-peak_ratio": d["t-peak"]["hist_ratio"],
            "lipolgen_polrad-full_ratio": d["polrad-full"]["hist_ratio"]})
        detail.append(d)
    for c in control:
        subrows.append({
            "name": "control[Rad=%d, Q2 %g-%g]" % (c["rad"], c["q2_lo"], c["q2_hi"]),
            "status": "pass" if c["reproduced"] else "fail",
            "in_verdict": False,
            "sigma_nb": c["sigma_nb"], "epic_readme_nb": c["epic_readme_nb"],
            "ratio_minus_1": c["ratio_minus_1"], "pull": c["pull"]})

    q2_mid = [math.sqrt(b["q2_lo"] * b["q2_hi"]) for b in bins]
    r_full = [d["polrad-full"]["R"] for d in detail]
    r_tpk = [d["t-peak"]["R"] for d in detail]
    generator = ("R = sigma_tail(LiPolGen)/sigma_tail(DJANGOH) per Q2 bin: "
                 "polrad-full %s; t-peak (shipped) %s" % (
                     " / ".join("%.4f" % r for r in r_full),
                     " / ".join("%.4f" % r for r in r_tpk)))
    reference = ("DJANGOH 4.6.21 elastic radiative tail IEL31+IEL32+IEL33 "
                 "(fixed alpha), fraction of its own Born in the matched "
                 "window: %s; control reproduces the 8 ePIC numbers: %d/8" % (
                     " / ".join("%.5f" % d["djangoh_fraction"] for d in detail),
                     sum(1 for c in control if c["reproduced"])))
    tolerance = ("polrad-full: abs(R - 1) <= %g + %g delta_b each bin; "
                 "t-peak (a one-sided ceiling; the tree calls it a lower "
                 "bound): R <= %g + %g delta_b each bin; delta_b = DJANGOH's "
                 "relative integration error (%s)" % (
                     TOL_NUMERICAL, N_DELTA, 1.0 + TOL_NUMERICAL, N_DELTA,
                     " / ".join("%.1e" % b["delta"] for b in bins)))
    fails = [v for v in ver if not v["pass"]]
    if status == "pass":
        reason = ("polrad-full (two-sided) within %g %% + %g delta of "
                  "DJANGOH's exact O(alpha) proton elastic tail in all four "
                  "bins (R %s); the shipped t-peak passes its one-sided "
                  "CEILING only (R %s): it is %.0f %% / %.0f %% below the "
                  "exact tail at Q2 = 1-10 / 10-100 GeV2, and its R > 1 in "
                  "bins 3-4 is the tail table's +2 %% (--diagnostics)" % (
                      100.0 * TOL_NUMERICAL, N_DELTA,
                      " / ".join("%.4f" % r for r in r_full),
                      " / ".join("%.4f" % r for r in r_tpk),
                      100.0 * (1.0 - r_tpk[0]), 100.0 * (1.0 - r_tpk[1])))
    else:
        reason = "; ".join(
            "%s %s: R = %.4f outside %s" % (v["model"], _bin_label(bins[v["bin"]]),
                                            v["R"], v["test"]) for v in fails)
    row = {
        "name": NAME,
        "generator": generator,
        "reference": reference,
        "tolerance": tolerance,
        "status": status,
        "reason": reason,
        "generator_value": {"polrad-full": [round(r, 8) for r in r_full],
                            "t-peak": [round(r, 8) for r in r_tpk]},
        "reference_value": [round(d["djangoh_fraction"], 8) for d in detail],
        "subrows": subrows,
        "detail": detail,
        "control": control,
        # The 2026-09-23 row, kept as RECORDED HISTORY (never the verdict):
        # the ePIC Rad/noRad ratios and LiPolGen's 1 + INT Born r / INT Born
        # on the W_l^2 >= 9 window, with their Q^2 trends.
        "history": {
            "epic_ratio": [round(float(r["ratio"]), 8) for r in epic_rows],
            "lipolgen_ratio": {tm: [round(float(d[tm]["hist_ratio"]), 8)
                                    for d in detail] for tm in tail_models},
            "trend_dlnexcess_dlnQ2": dict(
                [("djangoh", log_slope(q2_mid, [r["ratio"] - 1.0
                                                for r in epic_rows]))] +
                [(tm, log_slope(q2_mid, [d[tm]["hist_ratio"] - 1.0
                                         for d in detail]))
                 for tm in tail_models])},
        "config": {
            "beams": "e- 18 x p 275 GeV (LiPolGen config 2; DJANGOH EL-BEAM 18, "
                     "PR-BEAM 275)",
            "window": "Q2_lo <= Q2_l < Q2_hi, 1e-5 <= x_l <= xmax_b, "
                      "1e-4 <= y_l <= 0.985 (W_l^2 >= 9 implied)",
            "xmax_b": [b["xmax"] for b in bins],
            "x_min": X_MIN, "y_min": Y_MIN, "y_max": gen["t-peak"]["y_max"],
            "W2_min": W2_MIN,
            "tol_numerical": TOL_NUMERICAL, "n_delta": N_DELTA,
            "design_block_sha256": DESIGN_BLOCK_SHA256,
            "djangoh": "4.6.21 (spiesber/DJANGOH e356b84), cteq61 = LHAPDF "
                       "10150, RNDM-SEEDS 0 0, INT-POINTS 100000, LPARIN(3) = 2",
            "rc_ff": "zero coherent vertex", "qe_kf_gev": 0.0,
            "tail_max": {tm: gen[tm]["tail_max"] for tm in tail_models},
            "clipped_node_fraction": {tm: gen[tm]["clipped_node_fraction"]
                                      for tm in tail_models},
            "tail_build_s": {tm: round(gen[tm]["build_s"], 1)
                             for tm in tail_models},
            "unpol_sf": {tm: gen[tm]["unpol_sf"] for tm in tail_models}},
    }
    if verbose:
        print_row(row)
    return row


def print_row(row):
    w = sys.stdout.write
    w("%-16s %-7s %-21s %-21s %-12s %-12s\n" % (
        "Q2 bin [GeV2]", "xmax", "DJANGOH f (+- err)", "polrad-full f (R)",
        "t-peak f (R)", "allowance"))
    for d in row["detail"]:
        w("%-16s %-7g %.5f (+- %.1e)  %.5f (%.4f)     %.5f (%.4f) %.4f\n" % (
            "[%g, %g]" % d["bin"], d["xmax"], d["djangoh_fraction"],
            d["djangoh_fraction_err"], d["polrad-full"]["fraction_on_djangoh_born"],
            d["polrad-full"]["R"], d["t-peak"]["fraction_on_djangoh_born"],
            d["t-peak"]["R"], d["allowance"]))
    for s in row["subrows"]:
        vals = {k: v for k, v in s.items()
                if k not in ("name", "status", "in_verdict", "reason")}
        w("SUBROW | %s[%s] | %s | %s\n" % (
            row["name"], s["name"], "; ".join(
                "%s %s" % (k, ("%.6g" % v) if isinstance(v, float) else v)
                for k, v in vals.items()),
            s["status"].upper() + ("" if s["in_verdict"] else " (not in verdict)")))
    w("REPORT | %s | %s | %s | %s | %s -- %s\n" % (
        row["name"], row["generator"], row["reference"], row["tolerance"],
        row["status"].upper(), row["reason"]))


# ------------------------------------------------------------ the figure

def plot(report, out_dir=None):
    """The elastic-tail fraction per Q^2 bin (DJANGOH points with their
    integration errors, LiPolGen t-peak and polrad-full as markers, the
    tolerance band around DJANGOH, the ePIC Rad/noRad - 1 greyed as the
    recorded history), over the ratio R = LiPolGen/DJANGOH with the bands."""
    sys.path.insert(0, str(HERE))
    from _plots import (new_figure, save, style_series, status_title, MUTED,
                        BAND)
    det = report["detail"]
    q = np.array([math.sqrt(d["bin"][0] * d["bin"][1]) for d in det])
    f_dj = np.array([d["djangoh_fraction"] for d in det])
    e_dj = np.array([d["djangoh_fraction_err"] for d in det])
    allow = np.array([d["allowance"] for d in det])
    fig, (ax, axr) = new_figure(two_panel=True)
    # the tolerance band, drawn per bin as a box half a decade wide
    for qi, fi, ai in zip(q, f_dj, allow):
        lo, hi = qi / 10 ** 0.25, qi * 10 ** 0.25
        ax.fill_between([lo, hi], [fi * (1 - ai)] * 2, [fi * (1 + ai)] * 2,
                        color=BAND, alpha=0.8, linewidth=0, zorder=1)
        axr.fill_between([lo, hi], [1 - ai] * 2, [1 + ai] * 2, color=BAND,
                         alpha=0.8, linewidth=0, zorder=1)
    from matplotlib.patches import Patch
    ax.errorbar(q, f_dj, yerr=e_dj, **style_series(
        "reference", label="DJANGOH 4.6.21 IEL31+32+33 (fixed alpha)"))
    for tm, slot, mk, off in (("polrad-full", "alt2", "s", 10 ** 0.06),
                              ("t-peak", "default", "D", 10 ** -0.06)):
        fr = np.array([d[tm]["fraction_on_djangoh_born"] for d in det])
        rr = np.array([d[tm]["R"] for d in det])
        lab = ("LiPolGen polrad-full (opt-in)" if tm == "polrad-full"
               else "LiPolGen t-peak (shipped default)")
        ax.plot(q * off, fr, **style_series(slot, label=lab, marker=mk))
        axr.plot(q * off, rr, **style_series(slot, label=None, marker=mk))
    ep = np.array([d["epic_rad_over_norad"] - 1.0 for d in det])
    ax.plot(q, ep, color=MUTED, marker="o", markersize=5, linestyle=":",
            linewidth=1.0, zorder=2, alpha=0.8,
            label="ePIC Rad/noRad - 1 (inelastic, IEL off; recorded)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-2, 1.0)
    ax.set_ylabel("tail / Born (DJANGOH)")
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Patch(color=BAND, alpha=0.8, linewidth=0))
    labels.append("tolerance 3 % + 2 delta (polrad-full)")
    ax.legend(handles, labels, loc="upper left", ncol=1)
    axr.axhline(1.0, color=MUTED, linewidth=1.0, zorder=0)
    axr.set_xscale("log")
    axr.set_xlabel("Q^2 bin (geometric centre) [GeV^2]")
    axr.set_ylim(0.5, 1.15)
    axr.set_ylabel("R = LPG / DJ")
    status_title(fig, NAME, report["status"],
                 "proton, e 18 x p 275, W_l^2 >= 9, y <= 0.985;  polrad-full "
                 "|R-1| <= 3% + 2 delta, t-peak R <= 1.03 + 2 delta")
    return save(fig, NAME, out_dir)


if __name__ == "__main__":
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    try:
        if "--diagnostics" in sys.argv[1:]:
            diagnostics()
        else:
            rep = run()
            if "--plot" in sys.argv[1:]:
                print("figure: %s" % plot(rep))
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        sys.exit(2)
    sys.exit(0)
