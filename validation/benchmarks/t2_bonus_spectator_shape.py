#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 10 -- CLAS BONuS: the cos theta_pq SHAPE of
R_D/S (data over BONuS's own spectator simulation) in the DIS blocks, held
against the tree's Glauber FSI weight on the deuteron control (tier T2).
BLOCKED from 2026-09-23 to 2026-10-26 (the history is kept below, dated);
a measurement since 2026-10-26 (run 2026-10-26, phase B1).

======================================================================
WHAT IS COMPARED (2026-10-26) -- A NEW DESIGN, PRESENTED AS SUCH
======================================================================
The 2026-09-23 row wanted the spectator-momentum SPECTRUM at A = 2.  The
document it waited for does not hold one (below: R_D/S is normalised per
p_s bin), so the row is redesigned as a check of the DIRECTION and SIZE of
the tree's spectator FSI weight at low p_s, where the paper finds the ratio
flattest (Sec. V A 1).  Nothing of the 2026-09-23 plan is reused but the
frame caveat.

THE REFERENCE.  The Supplemental Material (the paper's ref. [71]) of
S. Tkachenko et al. (CLAS), Phys. Rev. C 89 (2014) 045206,
doi:10.1103/PhysRevC.89.045206, BONuSprcSupMat2.txt (88 845 bytes, 5452
lines, sha256 18c7eb38...b3f65f0e), fetched 2026-10-26 with curl's own
identity through the institution's APS licence and kept OUTSIDE the
repository; its NUMBERS are vendored in
`data/bonus2014_tkachenko_supplemental.json` (provenance block inside:
URL, fetch date, 'facts quoted with attribution; APS copyright on the
document; the licence question is listed for the author', the convention
of every table, the source line of every block).  Sections 1 and 2 of the
supplement: R_D/S against cos theta_pq, 100 blocks of 10 rows (cos, R_D/S,
stat, syst):
  * 4 GeV beam (4.23 GeV, the paper's Fig. 12): 3 Q2 (0.93, 1.66, 3.38
    GeV^2) x 5 W* (1.17, 1.48, 1.73, 2.03, 2.44 GeV) x 4 p_s (0.078, 0.093,
    0.11, 0.135 GeV/c) = 60 blocks;
  * 5 GeV beam (5.25 GeV, Figs. 10-11): 2 Q2 (1.66, 3.38) x 5 W* x 4 p_s
    = 40 blocks.
CORRECTED 2026-10-26: this header said '6 W* x 5 Q2 x 4 p_s' until today --
the paper's text (Sec. V A 1), which counts the quasi-elastic W* region and
two Q2 bins (0.34, 0.61) that are not in the supplement.  The labels are
bin centres: W* = 2.03 is the paper's region 1.85-2.20 GeV and 2.44 is
2.2-2.68 GeV (Sec. V A); Q2 = 1.66 is 1.10-2.23 GeV^2 (Figs. 10-13
captions); p_s = 70-85, 85-100, 100-120, 120-150 MeV/c; cos theta_pq ten
even bins over [-1, 1].  A row whose four columns all read 0 is a MISSING
bin (22 rows, all in the four 4 GeV blocks at Q2 = 3.38, W* = 2.44), never
a measured zero.  Section 5 (DIS-region F2n/F2d and F2n/F2p for three W*
cuts, 13 / 12 / 10 rows) is vendored for the F2n/F2p sub-row below;
sections 3-4 (F2n/F2d per bin) and the two F2n_eff(p_s) ratio tables
(Figs. 15-16) are not vendored (listed in the file).

WHAT THE DENOMINATOR IS (read in arXiv:1402.2477v2 Sec. IV C before any
design).  BONuS's Monte Carlo: a PLANE-WAVE spectator model -- the Paris
deuteron wave function rescaled with the LIGHT-CONE formalism of Frankfurt
and Strikman, the struck neutron off shell with E = M_d - sqrt(M^2 + p_s^2),
p = -p_s, electron kinematics generated in the struck nucleon's rest frame
and corrected for its off-shell mass, Bosted-Christy F2n, Mo-Tsai
radiative effects, GEANT4 RTPC + GSIM CLAS acceptance; and in the paper's
words 'Our generator did not contain additional processes like coherent
pion production, final-state interactions and other two-nucleon effects;
therefore, the ratio of measured to simulated tagged data can be
interpreted as a direct test of the spectator picture' (Sec. IV C 1).  So
the denominator is PWIA: NO FSI.  The normalisation (Sec. IV C 3) is ONE
factor N(p_s, E_b) per p_s bin and beam energy, fitted on -1 < cos theta_pq
< -0.2 at W* = 2.0-2.2 GeV in the lowest fully accepted Q2 bin and common
to every cos theta_pq and (W*, Q2) bin of that (p_s, E_b); its +-10 % scale
is not in the syst column.  If the spectator picture holds, R_D/S is
F2n_eff(W*, Q2, p_s, cos theta_pq) / F2n_model(W*, Q2) (Eq. (15)), CONSTANT
in cos theta_pq at fixed (W*, Q2, p_s) -- 'Any overall trend, such as a
monotonic increase or decrease with cos theta_pq, would indicate a
shortcoming of the deuteron wave function model, while FSI effects are
expected to give rise to more complicated structures' (Sec. V A 1).

WHAT R_D/S CAN AND CANNOT SAY.  It CAN bound a cos theta_pq dependence at
fixed (W*, Q2, p_s) beyond the light-cone Paris PWIA: an FSI, off-shell or
fragmentation effect that varies with angle.  It CANNOT give (i) the p_s
spectrum (N(p_s, E_b) is fitted per p_s bin -- the RTPC efficiency falls
from ~0.6 to ~0.23 across the range and is normalised away), so the
spectator-momentum shape this row was planned for is NOT in the document;
(ii) an angle-independent FSI suppression (absorbed by the same
normalisation), so the weight's absolute level is untested; (iii) a
separation of FSI from the PWIA prescription: the denominator is
light-cone, the tree's spectator density is instant-form (recorded below,
as a number).

----------------------------------------------------------------------
THE RULE AND ITS TOLERANCE -- written 2026-10-26 BEFORE any fit, band, chi2
or tree number of this row was computed.  The supervising session previewed
the R_D/S tables before this design (run_2026-10-26/PLAN.md, ground rules
4), and this worker read every value while transcribing them into the
vendored file; the window below is defined on the block LABELS alone.
----------------------------------------------------------------------
(a) THE DIS WINDOW, declared first: W* >= 2 GeV AND Q2 >= 1 GeV^2 (block
    centres).  It keeps 32 blocks -- at EACH beam the W* = 2.03 and 2.44
    blocks at Q2 = 1.66 and 3.38, x 4 p_s.  Of these, the four 4 GeV blocks
    at (Q2 = 3.38, W* = 2.44) have NO measured backward bin (their rows at
    cos theta_pq <= -0.1 are MISSING: the 4 GeV beam does not reach W* =
    2.44 at Q2 = 3.38 with a backward spectator); they carry no verdict and
    are counted.  So 28 VERDICT BLOCKS: 12 at 4 GeV, 16 at 5 GeV.  OUT OF
    MODEL, recorded by count and never compared: the 68 other blocks -- the
    20 at Q2 = 0.93 (4 GeV; the bin lies mostly below 1 GeV^2) and the 48
    resonance-region blocks (W* = 1.17, 1.48, 1.73 at both Q2 and both
    beams).  Why: `GlauberFsiWeight` models the rescattering of DIS DEBRIS
    on the spectator at sigma_XN = 40 mb (a free-hadron cross section); a
    resonance or quasi-elastic final state is not that, and the generator
    never produces such an event -- its own scenario floor is W^2 >= 8
    GeV^2 (`generator_scenario`, include/lipolgen/sampler.hpp).  NO BONuS
    block reaches that floor: the highest W* bin is 2.2-2.68 GeV (W*^2 <=
    7.2 GeV^2), and the 2.03 bin is 1.85-2.20 GeV, its lower half in the
    third-resonance tail.  The window is as close to the generator's floor
    as the data allow, and the row tests the weight BELOW the generator's
    own W range -- an extrapolation of its 40 mb input to W* ~ 2-2.7 GeV.
(b) THE DATA'S FLATNESS BAND.  Backward window cos theta_pq <= -0.3: bin
    centres -0.9, -0.7, -0.5, -0.3, i.e. the bins [-1, -0.2] -- the paper's
    'backward' coarse bin and its normalisation region.  Per verdict block
    (>= 2 measured backward bins): the weighted constant fit with STAT
    errors, Rbar = sum(r/s^2)/sum(1/s^2), sigma = sum(1/s^2)^(-1/2),
    chi2_flat on ndf = n - 1; the scale factor S = max(1, sqrt(chi2_flat /
    ndf)); the band half-width relative to the level, delta = S sigma /
    Rbar.  STAT because the paper draws R_D/S with statistical error bars
    and its systematic as a band ('a relatively smooth function of the
    kinematic variables', Sec. IV C 4), and the +-10 % scale is common to
    every cos bin of a (p_s, beam): stat is the error of a within-block
    cos theta_pq shape.  S widens the band to the data's own departure from
    flat when that exceeds the errors: the departure is what bounds an
    effect.
(c) THE TREE -- a LIBRARY EVALUATION, no generation: `lipolgen.
    GlauberFsiWeight(TaggedModel(deuteron_channel()), GlauberFsiOptions())
    .weight(FsiKinematics(k = p_s, cos_theta_k = cos theta_pq, phi_k = 0,
    w = W*, q2 = Q2, x = x*, spectator_z = 1, spectator_a = 1, channel =
    Channel.TaggedDeuteronP))` (python/bindings.cpp `bind_fsi` exposes the
    weight itself, so no Pipeline is built and the row is not EXPENSIVE), at
    the block's p_s and at the cos bin centres; x* = Q2 / (W*^2 - M^2 + Q2).
    Default options: sigma_XN = 40 mb, eps = -0.5, B_XN = 6 GeV^-2,
    formation ramp OFF -- so the weight depends on (k, |cos theta_k|) alone
    (it reads its (k_z, k_T) table at |k_z|, symmetric under cos -> -cos,
    include/lipolgen/fsi.hpp), and W*, Q2, x reach it only through the ramp,
    which is off.  k = p_s: BONuS is a fixed-target experiment, so the lab
    IS the deuteron rest frame, the frame of the tree's k.  Normalised to its
    backward-window mean: u(c) = w(c) / mean_bw(w), the plain mean over the
    block's measured backward bin centres; the tree's variation v =
    max_bw |u(c) - 1|.
(d) TOLERANCE (must not contradict): a block PASSES iff v <= delta -- the
    tree's normalised weight stays inside the data's 1-sigma flat band at
    every backward bin centre.  A beam sub-row passes iff every one of its
    verdict blocks passes.  THE ROW is PASS iff the DEFAULT configuration
    passes at both beams: the Hulthen deuteron (deuteron_channel() default,
    P_D = 0.045, the tagged channel's shipped wave) at sigma_XN = 40 mb.
(e) RECORDED beside the verdict, never in it, each by the same rule where a
    rule applies:
    * AV18 (deuteron_channel(source = VmcAV18)) at 40 mb;
    * Hulthen at sigma_XN = 20 mb, the other end of the documented 20-40 mb
      band (fsi.hpp: 'band it; never quote one row alone');
    * the default with a stat (+) syst band (syst read as uncorrelated point
      to point -- the loose reading);
    * per block, chi2 of the backward data against s w(c), one free scale
      s, beside chi2_flat (same ndf): which shape the data prefer;
    * THE FORWARD RISE (cos theta_pq > 0: bins 0.1 ... 0.9), one free scale
      per block (BENCHMARK_PLAN.md sec. 8 rule 3 -- the ratio is normalised
      to a simulation): rho_data = Rbar_fwd / Rbar_bw (stat-weighted means)
      against rho_tree = mean_fwd(w) / mean_bw(w), their SIGN and SIZE per
      block and per (beam, p_s).  Said before the number: a weight that is
      symmetric under cos -> -cos cannot rise forward beyond the mirror of
      its own backward shape (the dropped principal-value piece, fsi.hpp);
    * THE PWIA PRESCRIPTION, not FSI: the tree's Hulthen density as an
      instant-form P = n(k = p_s) over the light-cone one n(k_LC) E_k /
      (2 - alpha_s) (the map t2_deeps_spectator_tail.py codes, coded here
      the same), normalised to its backward mean -- its max deviation per
      p_s is the size of what the prescription alone moves in this ratio's
      shape.  Not the tree's prediction of R_D/S (the tree is not run
      through BONuS's binning or acceptance);
    * F2n/F2p: the shipped `ToyF2.f2n_over_f2p(x)` at x*_aver against the
      section-5 W* > 1.8 GeV table (the paper's DIS-most cut, Fig. 21; 10
      rows, ratio method, cross-normalised to the CJ fit at x = 0.32),
      stat (+) syst: chi2 absolute and with one free scale, tree/data per
      row, and the 5 rows with W*_aver >= 2 GeV separately.  No tolerance.
(f) THE FIGURE (`--plot`; a view, never the record): one panel, R_D/S over
    its block's backward Rbar against cos theta_pq for the two highest-W*
    blocks per beam at Q2 = 1.66 GeV^2, p_s = 0.093 GeV/c (85-100 MeV/c,
    the upper half of the paper's VIP region) -- chosen here, before any
    number -- with the three configurations' u(c) at that p_s as lines.

MEASURED (2026-10-26, this tree, AFTER the rule above was fixed): FAIL.
  Default (Hulthen, 40 mb), stat band: the weight stays inside the data's
  1-sigma flat band in 4 of 12 blocks at 4 GeV (the four at Q2 = 3.38,
  W* = 2.03, where sigma(Rbar) is 6-12 %) and 0 of 16 at 5 GeV; worst
  v/delta 2.62 (4 GeV) and 3.44 (5 GeV).  The weight's variation across the
  backward window is v = 0.034 / 0.055 / 0.085 / 0.141 at p_s = 0.078 /
  0.093 / 0.11 / 0.135 GeV/c -- w(-0.9)/w(-0.3) = 1.071 / 1.114 / 1.179 /
  1.310, more absorption toward 90 deg, growing with p_s -- against data
  bands whose narrowest per p_s are delta = 0.0100 / 0.0255 / 0.0505 /
  0.0536.
  WHAT THE FAIL MEASURES (read before quoting it).  The rule's premise --
  that the backward window is flat within its errors, so that its flatness
  bounds any FSI effect -- does NOT hold at the stat level: chi2_flat =
  739.4 on 84 dof over the 28 blocks (205.7/36 at 4 GeV, 533.7/48 at
  5 GeV).  The data FALL toward 90 deg in the backward window, as the
  weight does: chi2 against the weight's shape, one free scale per block,
  is 226.3/84 for Hulthen 40 mb, 221.6 for AV18 40 mb, 403.3 for Hulthen
  20 mb (recorded, (e)) -- the 40 mb shape is preferred to flat by
  Delta chi2 = 513 and to 20 mb by 177, and none fits at stat.  Over the
  whole cos range the data dip near cos theta_pq = 0 and rise on both
  sides, the dip deepening with p_s -- what the paper reports ('a
  depression at angles around 90 deg develops at higher p_s', Sec. V A 1)
  and the shape of the symmetric Glauber weight.  So the FAIL says the
  weight varies MORE than the 1-sigma level of a flat fit; it does not say
  the data contradict the weight's direction, which the data share.
  RECORDED: stat (+) syst band (syst ~12 % per point, chi2_flat 57.4/84):
  8/12 and 8/16 blocks inside; Hulthen 20 mb: 8/12 and 13/16; AV18 40 mb:
  3/12 and 0/16.  Forward rise: rho_data = 1.036 / 0.973 / 1.020 / 1.026
  (4 GeV) and 1.031 / 0.994 / 0.974 / 1.030 (5 GeV) per p_s against
  rho_tree = 0.990 / 0.984 / 0.974 / 0.955 (the symmetric weight's forward
  mean includes the cos = 0.1 bin of its dip); same sign in 12 of 28
  blocks; the data's forward excess beyond the mirror (e.g. 1.355 at
  cos = 0.9 against the weight's 1.126; 5 GeV, Q2 = 1.66, W* = 2.03,
  p_s = 0.135) is what a cos -> -cos symmetric weight cannot carry.  The
  PWIA prescription (NOT FSI): the Hulthen instant-form / light-cone factor
  varies by 0.023 / 0.031 / 0.041 / 0.057 across the backward window, in
  the SAME direction as the weight (f(-0.9)/f(-0.3) = 1.046 / 1.063 /
  1.084 / 1.118): in the backward window alone R_D/S cannot separate the
  weight's dip from a light-cone vs instant-form difference of 40-66 % of
  its size; the full-range shape, a dip with a rise on BOTH sides, is not
  the monotonic prescription's.  F2n/F2p: ToyF2 / BONuS = 1.190 at
  x* = 0.177 ... 1.103 at 0.32 ... 0.990 at 0.52 ... 1.286 at 0.62;
  chi2 40.5/10 absolute (p = 1.4e-5), 8.30/9 (p = 0.50) with one free
  scale 0.902; on the 5 rows with W*_aver >= 2 GeV, 35.5/5 and 1.70/4 at
  scale 0.878 (the table is cross-normalised to CJ at x = 0.32).
  CHECKS (asserted in python/tests/test_bench_bonus.py): the weight's v is
  converged in its (k_z, k_T) table (the 161 x 161 default against 321 x
  321 on [0, 0.6] GeV: |dv| < 5e-4); it is evaluated at the bin labels, and
  averaging it uniformly over each p_s bin and each cos bin moves v by
  -0.4 ... -1.9 % of itself (< 2 %); it is exactly symmetric,
  max |w(c) - w(-c)| = 0.
D6 NOTE (2026-10-29; decision D6 of 2026-10-27,
docs/open_items/run_2026-10-27/DECISIONS.md): the tagged channels' shipped
wave is AV18 since (CLUSTER_WAVE_DEFAULT).  Rule (d) above DEFINES the
verdict configuration as the Hulthen deuteron at 40 mb, so DEFAULT_CONFIG
stays "hulthen-40mb"; AV18 at 40 mb, recorded under (e), fails too (3/12
and 0/16; tree shape 221.6/84 against Hulthen's 226.3/84).  The row is FAIL
either way; re-pointing it to the shipped wave belongs to the author's
pending re-declaration of this row (run_2026-10-26/SUMMARY.md, "For the
author", item 6).

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE (unpolarized; the weight is spin
     independent by construction and never enters A_zz).
  2. isoscalar denominator -- NOT APPLICABLE; its analogue, F2n of the
     struck neutron at (W*, Q2), is FIXED within a block (x* depends on
     |p_s|, W*, Q2, not on the angle) and cancels in the flatness test.
     The F2n/F2p sub-row compares a ratio to a ratio: F2d in the vendored
     F2n/F2d is PER DEUTERON (F2n/F2d ~ r/(1 + r)), not read here.
  3. per nucleus / per nucleon -- R_D/S is a ratio of two tagged
     D(e,e'p_s)X yields per deuteron, dimensionless; its level is a
     normalisation to a SIMULATION, so only the shape across cos theta_pq
     is compared (one free scale per block), never an absolute factor.

THE FRAME CAVEAT (02_data_nucleon_deuteron.md F-9 / sec. 5.5), kept from
2026-09-23 and riding with every number here: BONuS tags BACKWARD in the
target rest frame at p_s < 150 MeV/c; the EIC tags a boosted spectator in
the far-forward acceptance.  Same physics, different frames; acceptance-
driven systematics do not transfer.  LiPolGen's k IS the rest-frame
momentum, and BONuS's spectrum is light-cone-corrected and acceptance-
folded (in the denominator, which the ratio divides out).

WHAT THIS ROW DOES NOT VALIDATE: the spin dependence of the rescattering
(the weight is spin-blind by construction -- it never moves A_zz^tag, only
the spectator spectrum); the weight's absolute level (the survival
probability -- normalised away); the weight at p_s > 0.15 GeV/c (Deeps,
row 11, transverse window, recorded) or at EIC W (formation lengths of many
fm, the sigma_XN(W) ramp is off and untested); the p_s spectrum of n(k);
the tagged lithium channels (A = 6, 7: alpha and triton spectators, other
profiles); the generator's event sampling, binning or acceptance (no event
is generated).

======================================================================
HISTORY (dated record, 2026-09-23 .. 2026-10-26; kept verbatim in substance)
======================================================================
WHAT THE ROW WANTED.  A machine-readable BONuS (Baillie et al., PRL 108
(2012) 142001; Tkachenko et al., PRC 89 (2014) 045206) spectator-momentum
spectrum -- backward protons, p_s < 100 MeV/c, theta_pq > 100 deg, in the
deuteron REST frame -- to compare against the SHAPE of LiPolGen's deuteron-
control tagged spectator spectrum (`deuteron_channel()` -> `TaggedModel`,
n(k, cos theta_k) in the ion rest frame, averaged over the ion's M).

WHAT THE DATABASE RETURNED (queried 2026-09-23, read-only;
`data/clas_db_deuteron_query_2026-09-23.json`).  The CLAS Physics Database
(https://clas.sinp.msu.ru/cgi-bin/jlab/db.cgi) holds, for a deuteron target,
918 measurements in 19 experiments.  BONuS (E-03-012, eid 135) is ONE
measurement, E135M1: the ratio F2n/F2p against W* (Baillie 2012).  There is
no BONuS spectator-momentum spectrum and no Tkachenko 2014 entry.  EG1b
(eid 95 and 146) is inclusive (g1, A1, g1/F1, Gamma_1, ...).  So the table
this row needs is NOT machine-readable in the database.  Both papers are on
disk (refs/1110.2770.pdf, refs/1402.2477.pdf) and neither prints a p_s table
in its body.  Tkachenko et al. ref. [71] puts "tables of numerical results"
in the journal's Supplemental Material, which answered HTTP 401 (APS login)
on 2026-09-23; arXiv:1402.2477 has no ancillary files.  (This paragraph
said the tables were 'for 6 W* x 5 Q2 x 4 p_s bins' -- corrected above.)
Those tables are R_D/S, data over BONuS's OWN spectator simulation, so they
give a p_s SHAPE only together with that simulation.  Unblocked 2026-10-26:
the Supplemental Material was fetched through the institution's licence.

WHAT IT DID FIND.  The only deuteron spectator-momentum distribution in the
database is Deeps (eid 90; Klimenko et al., PRC 73 (2006) 035212,
nucl-ex/0510032): F2N x P(p_s, cos theta_pq) at p_s = 0.30-0.53 GeV/c,
VENDORED in `data/clas_deeps_klimenko2006_f2P.json` and wired since
2026-09-27 as its own row (t2_deeps_spectator_tail.py, row 11).  WHAT THE
FILE HOLDS (counted by `findings()` on every run; tabulated 2026-09-26):
115 blocks (database measurement pages, mid 1-115), of which
  * 60 (mid 1-60) are the paper's bins: Q2 = 1.8 and 2.8 GeV^2 x six W*
    (0.94, 1.25, 1.5, 1.73, 2.02, 2.4 GeV) x five p_s -- one block each;
  * 55 (mid 61-115) carry database Q2 labels 0.18 (5 blocks) and 0.25,
    0.35, 0.45, 0.55, 0.65 (10 each), all at W* = 2.4, that are NOT the
    paper's Q2 bins (kept verbatim, NOT interpreted -- the file's own
    provenance).  Their (q2, W*, p_s) labels are not even unique: 25 label
    triples occur TWICE with different data (e.g. mid 62 and 67 are both
    (0.25, 2.4, 300 MeV), 11 and 8 rows), so nothing may key on them;
  * one block, mid 78 (labelled 0.25, 2.4, 340 MeV), has NO rows;
  * 60 rows (in 59 blocks, 55 of them a block's first cos theta_pq bin)
    read value = 0 with stat = 0 and a nonzero syst: vendored as the
    database lists them, but to be read as MISSING, not measured zeros.

INTEGRITY.  The sha256 of every byte of each vendored file is recorded
below and re-checked on EVERY read: an altered file is refused with a
RuntimeError, so `__main__` exits 2.  QUERY_SHA256 and DEEPS_SHA256 were
taken from the files as committed at f0a8f1e (2026-09-26); BONUS_SHA256
from the file as written on 2026-10-26, which also records the original
txt's sha256 (checked against ORIGINAL_TXT_SHA256 on every read).  The
catalogue's quantity names carry Greek letters and are decoded as UTF-8
whatever the locale; every printout escapes what the terminal cannot
encode, so a non-UTF-8 stdout still gets the REPORT row and exit 0.

Run:  source env.sh && python3 validation/benchmarks/t2_bonus_spectator_shape.py [--plot]

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import collections
import hashlib
import json
import math
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
BONUS = os.path.join(HERE, "data", "bonus2014_tkachenko_supplemental.json")
QUERY = os.path.join(HERE, "data", "clas_db_deuteron_query_2026-09-23.json")
DEEPS = os.path.join(HERE, "data", "clas_deeps_klimenko2006_f2P.json")
#: sha256 of every byte of each vendored file; re-checked on every read.
#: QUERY / DEEPS recorded 2026-09-26 from the files as committed at f0a8f1e;
#: BONUS recorded 2026-10-26 from the file as written by the vendoring.
BONUS_SHA256 = "3636ff08adbd1d7b526b30076e03512164e036c31392e291fbd0fcf6f88d3201"
QUERY_SHA256 = "340e04a5ae8ada1eabf4d90c2c76c82b07bda48d28c7f43cdffd881c8ee3fa11"
DEEPS_SHA256 = "1d98fd11e1b11b1e070e5504f4919b98b87dbabfe29cdca347522ea4f8737256"
#: sha256 of the ORIGINAL BONuSprcSupMat2.txt (88 845 bytes, outside the
#: repository), as the vendored file's provenance records it
ORIGINAL_TXT_SHA256 = "18c7eb384c59c854e9e7926a1f96ec6840a9cc05674142bdfefee837b3f65f0e"
#: the Deeps paper's two Q2 bins (Klimenko 2006 Sec. IV averages), GeV^2
DEEPS_PAPER_Q2 = (1.8, 2.8)

NAME = "t2_bonus_spectator_shape"

# ---- the design (module docstring, 'THE RULE AND ITS TOLERANCE') ----------
#: (a) the DIS window on block centres
W_STAR_MIN = 2.0          # GeV
Q2_MIN = 1.0              # GeV^2
#: (b) backward window: cos theta_pq bin centres <= this
BACKWARD_COS_MAX = -0.3
#: a verdict block needs at least this many measured backward bins
MIN_BACKWARD_POINTS = 2
#: (e) forward window: cos theta_pq bin centres > this
FORWARD_COS_GT = 0.0
#: (e) the F2n/F2p sub-row reads this cut of supplement section 5
F2N_F2P_CUT = "W > 1.8 GeV"
#: (f) the figure's blocks (chosen before any number)
FIGURE_Q2 = 1.66
FIGURE_PS = 0.093
#: the supplement's two ratio sections -> beam labels
BEAMS = {1: "4GeV", 2: "5GeV"}
#: configurations: (wave, sigma_XN [mb], in the verdict)
CONFIGS = collections.OrderedDict([
    ("hulthen-40mb", ("hulthen", 40.0, True)),
    ("av18-40mb", ("av18", 40.0, False)),
    ("hulthen-20mb", ("hulthen", 20.0, False)),
])
DEFAULT_CONFIG = "hulthen-40mb"
WAVE_LABELS = {
    "hulthen": "Hulthen pair, P_D = 0.045 (deuteron_channel() default)",
    "av18": "AV18 (deuteron_channel(source=VmcAV18), fdeut.av18)",
}

#: the row's reason from 2026-09-23 to 2026-10-26 -- the dated record, no
#: longer the row's status
BLOCKED_REASON = (
    "BLOCKED: the CLAS Physics Database (queried 2026-09-23) holds BONuS "
    "(E-03-012, eid 135) only as ONE measurement, E135M1 = F2n/F2p vs W* "
    "(Baillie et al., PRL 108 (2012) 142001) -- no spectator-momentum "
    "spectrum, no Tkachenko 2014 entry; EG1b (eid 95, 146) is inclusive "
    "g1/A1 only.  To unblock, download the Supplemental Material (tables "
    "of numerical results, their ref. [71]) of S. Tkachenko et al. (CLAS "
    "BONuS), Phys. Rev. C 89 (2014) 045206, doi:10.1103/PhysRevC.89.045206, "
    "https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.89.045206 "
    "(HTTP 401 without an APS login on 2026-09-23; arXiv:1402.2477 has no "
    "ancillary files), into validation/benchmarks/data/.  Found instead and "
    "vendored unwired: Deeps F2N x P(p_s) at 0.30-0.53 GeV/c (eid 90, "
    "Klimenko et al., PRC 73 (2006) 035212).")
#: when that reason stopped being the row's status
UNBLOCKED = ("2026-10-26: the Supplemental Material fetched through the "
             "institution's APS licence; its numbers vendored in "
             "data/bonus2014_tkachenko_supplemental.json")


# ------------------------------------------------------------------ files

def read_vendored(path, sha256):
    """The bytes of a vendored file -- RuntimeError unless their sha256 is
    the recorded one, so a moved number is refused, never silently read."""
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != sha256:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, sha256))
    return raw


def load_query(path=QUERY):
    return json.loads(read_vendored(path, QUERY_SHA256).decode("utf-8"))


def load_deeps(path=DEEPS):
    return json.loads(read_vendored(path, DEEPS_SHA256).decode("utf-8"))


def load_bonus(path=BONUS):
    """The vendored supplement, sha256-checked, and the original txt's
    digest it records checked against ORIGINAL_TXT_SHA256."""
    doc = json.loads(read_vendored(path, BONUS_SHA256).decode("utf-8"))
    got = doc["_provenance"]["original_file"]["sha256"]
    if got != ORIGINAL_TXT_SHA256:
        raise RuntimeError("%s: records original sha256 %s != %s"
                           % (path, got, ORIGINAL_TXT_SHA256))
    return doc


def deeps_census(d):
    """What the Deeps file holds, counted (the HISTORY census)."""
    blocks = d["blocks"]
    paper = [b for b in blocks if b["q2"] in DEEPS_PAPER_Q2]
    labels = collections.Counter((b["q2"], b["w_star"], b["ps_MeV"])
                                 for b in blocks if b["q2"] not in DEEPS_PAPER_Q2)
    zero = [(b["mid"], i) for b in blocks for i, r in enumerate(b["rows"])
            if r[1] == 0.0 and r[2] == 0.0]
    return dict(
        deeps_blocks=len(blocks),
        deeps_blocks_paper_q2=len(paper),
        deeps_blocks_unverified_q2=len(blocks) - len(paper),
        deeps_duplicate_label_triples=sum(1 for n in labels.values() if n > 1),
        deeps_empty_blocks=[b["mid"] for b in blocks if not b["rows"]],
        deeps_zero_value_zero_stat_rows=len(zero),
    )


def _printable(text):
    """`text` as the current stdout can encode it (backslash escapes for
    what it cannot), so a non-UTF-8 terminal still gets the REPORT row."""
    enc = getattr(sys.stdout, "encoding", None) or "ascii"
    return text.encode(enc, "backslashreplace").decode(enc)


def findings(query_path=QUERY, deeps_path=DEEPS):
    """The 2026-09-23 database query and the Deeps census (HISTORY)."""
    q = load_query(query_path)
    exps = {e["eid"]: e for e in q["experiments"]}
    d = load_deeps(deeps_path)
    return dict(
        n_experiments=len(exps),
        n_measurements=sum(e["n_measurements"] for e in exps.values()),
        bonus=exps.get(135),
        eg1b=[exps[k] for k in (95, 146) if k in exps],
        deeps=exps.get(90),
        **deeps_census(d),
    )


# ------------------------------------------------------------------ blocks

def is_missing(row):
    """A row whose four columns all read 0 is a MISSING bin."""
    return all(v == 0.0 for v in row)


def blocks(doc):
    """Every R_D/S block as a dict: beam, (q2, w_star, ps), src_line, the
    measured rows as (cos, r, stat, syst) and the count of missing rows.
    The (beam, q2, w_star, ps) key is asserted unique."""
    out, seen = [], set()
    for sec in doc["ratio_blocks"]["sections"]:
        beam = BEAMS[sec["section"]]
        for b in sec["blocks"]:
            key = (beam, b["q2"], b["w_star"], b["ps"])
            if key in seen:
                raise RuntimeError("duplicate R_D/S block %r" % (key,))
            seen.add(key)
            rows = [tuple(r) for r in b["rows"] if not is_missing(r)]
            out.append(dict(beam=beam, q2=b["q2"], w_star=b["w_star"],
                            ps=b["ps"], src_line=b["src_line"], rows=rows,
                            n_missing=len(b["rows"]) - len(rows)))
    return out


def in_dis_window(b):
    return b["w_star"] >= W_STAR_MIN and b["q2"] >= Q2_MIN


def backward(b):
    return [r for r in b["rows"] if r[0] <= BACKWARD_COS_MAX]


def forward(b):
    return [r for r in b["rows"] if r[0] > FORWARD_COS_GT]


def census(doc):
    """The block counts the header states, recomputed on every run."""
    bl = blocks(doc)
    dis = [b for b in bl if in_dis_window(b)]
    verdict = [b for b in dis if len(backward(b)) >= MIN_BACKWARD_POINTS]
    no_bw = [b for b in dis if len(backward(b)) < MIN_BACKWARD_POINTS]
    out = [b for b in bl if not in_dis_window(b)]
    return dict(
        blocks=len(bl),
        blocks_by_beam={k: sum(1 for b in bl if b["beam"] == k)
                        for k in BEAMS.values()},
        rows=sum(len(b["rows"]) + b["n_missing"] for b in bl),
        missing_rows=sum(b["n_missing"] for b in bl),
        blocks_with_missing=[(b["beam"], b["q2"], b["w_star"], b["ps"])
                             for b in bl if b["n_missing"]],
        dis_window=len(dis),
        verdict_blocks=len(verdict),
        verdict_blocks_by_beam={k: sum(1 for b in verdict if b["beam"] == k)
                                for k in BEAMS.values()},
        dis_no_backward=[(b["beam"], b["q2"], b["w_star"], b["ps"])
                         for b in no_bw],
        out_of_model=len(out),
        out_of_model_q2_below=sum(1 for b in out if b["q2"] < Q2_MIN),
        out_of_model_resonance=sum(1 for b in out if b["q2"] >= Q2_MIN
                                   and b["w_star"] < W_STAR_MIN),
        f2n_ratio_tables={t["cut_verbatim"]: len(t["rows"])
                          for t in doc["dis_f2n_ratios"]["tables"]},
    )


# ------------------------------------------------------------------ stats

def chi2_sf_fallback(chi2, ndf):
    """Q(ndf/2, chi2/2) without scipy -- t1_nmc_f2d.py's `chi2_sf_fallback`
    (Numerical Recipes' gammq: the series for P when chi2/2 < ndf/2 + 1,
    else the continued fraction for Q itself), copied unchanged."""
    chi2, ndf = float(chi2), int(ndf)
    if ndf <= 0 or not chi2 >= 0.0:
        raise ValueError("chi2_sf: need ndf > 0 and chi2 >= 0 (got %r, %r)"
                         % (chi2, ndf))
    if chi2 == 0.0:
        return 1.0
    if math.isinf(chi2):
        return 0.0
    a, xx = ndf / 2.0, chi2 / 2.0
    log_pref = -xx + a * math.log(xx) - math.lgamma(a)
    if xx < a + 1.0:
        term = s = 1.0 / a
        k = 1
        while abs(term) > 1e-17 * abs(s):
            term *= xx / (a + k)
            s += term
            k += 1
        return min(1.0, max(0.0, 1.0 - math.exp(log_pref) * s))
    tiny = 1e-300
    b = xx + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 100000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    q = math.exp(log_pref) * h
    return 0.0 if q < 2.2250738585072014e-308 else min(1.0, q)


def chi2_sf(chi2, ndf):
    """Upper-tail probability of chi2 with ndf dof (scipy if present, else
    `chi2_sf_fallback`)."""
    chi2, ndf = float(chi2), int(ndf)
    try:
        from scipy.stats import chi2 as _c
    except ImportError:                             # pragma: no cover
        return chi2_sf_fallback(chi2, ndf)
    return float(_c.sf(chi2, ndf))


def const_fit(vals, errs):
    """Weighted constant: (mean, sigma, chi2, ndf = n - 1)."""
    w = [1.0 / e ** 2 for e in errs]
    sw = sum(w)
    mean = sum(wi * v for wi, v in zip(w, vals)) / sw
    chi2 = sum(((v - mean) / e) ** 2 for v, e in zip(vals, errs))
    return mean, 1.0 / math.sqrt(sw), chi2, len(vals) - 1


def scale_fit(vals, errs, model):
    """One free scale s on a model shape: (s, chi2, ndf = n - 1).
    Invariant under model -> b model (s absorbs b)."""
    swm = sum(m * m / e ** 2 for m, e in zip(model, errs))
    s = sum(v * m / e ** 2 for v, m, e in zip(vals, model, errs)) / swm
    chi2 = sum(((v - s * m) / e) ** 2 for v, m, e in zip(vals, model, errs))
    return s, chi2, len(vals) - 1


def flat_band(vals, errs):
    """(b): the data's 1-sigma flat band of one block's backward window --
    Rbar, sigma, chi2_flat, ndf, S = max(1, sqrt(chi2/ndf)), delta =
    S sigma / Rbar."""
    mean, sig, chi2, ndf = const_fit(vals, errs)
    s_fac = max(1.0, math.sqrt(chi2 / ndf)) if ndf > 0 else 1.0
    return dict(rbar=mean, sigma=sig, chi2=chi2, ndf=ndf, scale_factor=s_fac,
                delta=s_fac * sig / mean)


# ------------------------------------------------------------------ tree

class Tree:
    """The tree side: GlauberFsiWeight on the deuteron control (a library
    evaluation), the waves' unpolarized densities, the light-cone map, and
    the shipped ToyF2."""

    def __init__(self):
        import lipolgen as lg                   # env.sh puts it on the path
        self.lg = lg
        self.m = lg.M_NUCLEON
        self.m_p = lg.PROTON_MASS
        self.m_d = lg.nuclear_mass(1, 2)
        self.models = {
            "hulthen": lg.TaggedModel(lg.deuteron_channel()),
            "av18": lg.TaggedModel(lg.deuteron_channel(
                source=lg.ClusterWaveSource.VmcAV18)),
        }
        self.fsi = {}
        for name, (wave, sigma, _) in CONFIGS.items():
            opt = lg.GlauberFsiOptions()
            opt.sigma_xn_mb = sigma
            self.fsi[name] = lg.GlauberFsiWeight(self.models[wave], opt)
        self.toy = lg.ToyF2()
        self._dens = {}

    def x_star(self, w, q2):
        return q2 / (w * w - self.m ** 2 + q2)

    def weight(self, config, ps, c, w, q2):
        """GlauberFsiWeight.weight at the spectator (k = p_s, cos theta_k =
        cos theta_pq) of a BONuS bin; (w, q2, x) reach it only through the
        formation ramp, which is off."""
        lg = self.lg
        return self.fsi[config].weight(lg.FsiKinematics(
            k=ps, cos_theta_k=c, phi_k=0.0, w=w, q2=q2,
            x=self.x_star(w, q2), spectator_z=1, spectator_a=1,
            channel=lg.Channel.TaggedDeuteronP))

    def density(self, wave, k):
        """The M-averaged (unpolarized, isotropic) momentum density per
        deuteron, sum_L (a_L psi_L(k))^2 / 4 pi with the model's OWN
        normalisation a_L = radial_table(L) / Wave.radial on its grid: it
        equals the M-averaged `TaggedModel.n_of_kc` table on every grid point
        (tested) and is exact between them, where n_of_kc is a nearest-cell
        lookup (the cell below), too coarse for a 1-6 MeV shift of k."""
        if wave not in self._dens:
            mod = self.models[wave]
            ch = mod.channel
            kappa = ch.base.kappa()
            kg = list(mod.k)
            i0 = len(kg) // 4
            comps = []
            for wv in ch.waves:
                tab = list(mod.radial_table(wv.l))
                comps.append((wv, tab[i0] / wv.radial(kg[i0], kappa)))
            self._dens[wave] = (comps, kappa)
        comps, kappa = self._dens[wave]
        return sum((a * wv.radial(k, kappa)) ** 2
                   for wv, a in comps) / (4.0 * math.pi)

    def lc(self, ps, c):
        """(alpha_s, p_T, k_LC, flux E_k / (2 - alpha_s)) for a spectator
        proton of lab momentum ps at cos theta_pq = c -- the map
        t2_deeps_spectator_tail.Tree.lc codes, coded the same."""
        e_s = math.sqrt(self.m_p ** 2 + ps ** 2)
        alpha = (e_s - ps * c) / (0.5 * self.m_d)
        pt2 = ps * ps * max(0.0, 1.0 - c * c)
        k2 = (self.m ** 2 + pt2) / (alpha * (2.0 - alpha)) - self.m ** 2
        k = math.sqrt(max(k2, 0.0))
        return alpha, math.sqrt(pt2), k, math.sqrt(self.m ** 2 + k * k) / (2.0 - alpha)

    def instant_over_lc(self, wave, ps, c):
        """n(k = p_s) / [n(k_LC) E_k / (2 - alpha_s)]: the PWIA prescription
        factor (instant form over light cone), up to a constant."""
        _, _, k, flux = self.lc(ps, c)
        return self.density(wave, ps) / (self.density(wave, k) * flux)


# ------------------------------------------------------------------ the row

def judge_block(b, wfun, errs="stat"):
    """One verdict block under one weight: the data's band (b), the tree's
    normalised variation (c), the verdict (d), and what (e) records."""
    bw = backward(b)
    cs = [r[0] for r in bw]
    vals = [r[1] for r in bw]
    if errs == "stat":
        es = [r[2] for r in bw]
    else:
        es = [math.hypot(r[2], r[3]) for r in bw]
    band = flat_band(vals, es)
    w_bw = [wfun(c) for c in cs]
    mean_w = sum(w_bw) / len(w_bw)
    u = [x / mean_w for x in w_bw]
    v = max(abs(x - 1.0) for x in u)
    _, chi2_tree, _ = scale_fit(vals, es, w_bw)
    out = dict(beam=b["beam"], q2=b["q2"], w_star=b["w_star"], ps=b["ps"],
               n_bw=len(bw), rbar=band["rbar"], sigma=band["sigma"],
               chi2_flat=band["chi2"], ndf=band["ndf"],
               scale_factor=band["scale_factor"], delta=band["delta"],
               v=v, v_over_delta=v / band["delta"],
               u_range=max(u) - min(u), chi2_tree=chi2_tree,
               w_bw_mean=mean_w, passed=v <= band["delta"])
    fw = forward(b)
    if errs == "stat" and fw:
        f_vals, f_errs = [r[1] for r in fw], [r[2] for r in fw]
        fmean, fsig, _, _ = const_fit(f_vals, f_errs)
        rho = fmean / band["rbar"]
        rho_err = rho * math.hypot(fsig / fmean, band["sigma"] / band["rbar"])
        w_fw = [wfun(r[0]) for r in fw]
        rho_tree = (sum(w_fw) / len(w_fw)) / mean_w
        out.update(n_fw=len(fw), rho_data=rho, rho_data_err=rho_err,
                   rho_tree=rho_tree,
                   same_sign=(rho - 1.0) * (rho_tree - 1.0) > 0.0,
                   rise_pull=(rho - rho_tree) / rho_err)
    return out


def summarise(judged):
    """A beam sub-row from its judged blocks."""
    n_pass = sum(1 for j in judged if j["passed"])
    worst = max(judged, key=lambda j: j["v_over_delta"])
    by_ps = collections.OrderedDict()
    for ps in sorted({j["ps"] for j in judged}):
        sel = [j for j in judged if j["ps"] == ps]
        by_ps[ps] = dict(v=max(j["v"] for j in sel),
                         min_delta=min(j["delta"] for j in sel),
                         n=len(sel), n_pass=sum(1 for j in sel if j["passed"]))
    return dict(
        status="pass" if n_pass == len(judged) else "fail",
        n_blocks=len(judged), n_pass=n_pass,
        worst_block="Q2=%g,W*=%g,ps=%g" % (worst["q2"], worst["w_star"],
                                            worst["ps"]),
        worst_v_over_delta=worst["v_over_delta"],
        chi2_flat_sum=sum(j["chi2_flat"] for j in judged),
        chi2_tree_sum=sum(j["chi2_tree"] for j in judged),
        ndf_sum=sum(j["ndf"] for j in judged),
        v_by_ps=[d["v"] for d in by_ps.values()],
        min_delta_by_ps=[d["min_delta"] for d in by_ps.values()],
        n_pass_by_ps=[d["n_pass"] for d in by_ps.values()],
        ps=list(by_ps),
    )


def forward_record(judged):
    """(e) the forward rise per (beam, p_s): stat-weighted mean of rho_data
    over the blocks, the tree's rho (one per p_s with the ramp off)."""
    out = collections.OrderedDict()
    for beam in BEAMS.values():
        for ps in sorted({j["ps"] for j in judged if j["beam"] == beam}):
            sel = [j for j in judged if j["beam"] == beam and j["ps"] == ps
                   and "rho_data" in j]
            if not sel:
                continue
            wts = [1.0 / j["rho_data_err"] ** 2 for j in sel]
            rd = sum(w * j["rho_data"] for w, j in zip(wts, sel)) / sum(wts)
            out["%s/ps=%g" % (beam, ps)] = dict(
                n=len(sel), rho_data=rd, rho_data_err=1.0 / math.sqrt(sum(wts)),
                rho_tree_min=min(j["rho_tree"] for j in sel),
                rho_tree_max=max(j["rho_tree"] for j in sel),
                same_sign_blocks=sum(1 for j in sel if j["same_sign"]))
    return out


def f2n_f2p_record(doc, tree):
    """(e) ToyF2.f2n_over_f2p against the section-5 W* > 1.8 GeV table."""
    tab = [t for t in doc["dis_f2n_ratios"]["tables"]
           if t["cut_verbatim"] == F2N_F2P_CUT][0]
    rows = tab["rows"]
    out = {}
    for label, sel in (("all", rows),
                       ("w_ge_2", [r for r in rows if r[0] >= W_STAR_MIN])):
        x = [r[2] for r in sel]
        d = [r[6] for r in sel]
        e = [math.hypot(r[7], r[8]) for r in sel]
        m = [tree.toy.f2n_over_f2p(xi) for xi in x]
        chi2_abs = sum(((di - mi) / ei) ** 2 for di, mi, ei in zip(d, m, e))
        s, chi2_s, ndf_s = scale_fit(d, e, m)
        out[label] = dict(n=len(sel), chi2_abs=chi2_abs,
                          p_abs=chi2_sf(chi2_abs, len(sel)),
                          scale=s, chi2_scaled=chi2_s, ndf_scaled=ndf_s,
                          p_scaled=chi2_sf(chi2_s, ndf_s),
                          tree_over_data=[mi / di for mi, di in zip(m, d)],
                          x=x)
    return out


def prescription_record(tree, ps_values, wave="hulthen"):
    """(e) the instant-form / light-cone factor per p_s, normalised to its
    backward mean at the backward bin centres: its max deviation, and its
    end-to-end ratio f(-0.9)/f(-0.3) beside the default weight's
    w(-0.9)/w(-0.3) (same direction when both exceed 1)."""
    cs = [-0.9, -0.7, -0.5, -0.3]
    out = collections.OrderedDict()
    for ps in ps_values:
        f = [tree.instant_over_lc(wave, ps, c) for c in cs]
        mean = sum(f) / len(f)
        w = [tree.weight(DEFAULT_CONFIG, ps, c, 2.44, 1.66) for c in cs]
        out["ps=%g" % ps] = dict(max_dev=max(abs(x / mean - 1.0) for x in f),
                                 ends_ratio=f[0] / f[-1],
                                 weight_ends_ratio=w[0] / w[-1])
    return out


def measure(path=BONUS):
    doc = load_bonus(path)
    cen = census(doc)
    bl = [b for b in blocks(doc) if in_dis_window(b)
          and len(backward(b)) >= MIN_BACKWARD_POINTS]
    tree = Tree()
    judged = {}
    for name in CONFIGS:
        judged[name] = [judge_block(b, lambda c, b=b, name=name: tree.weight(
            name, b["ps"], c, b["w_star"], b["q2"])) for b in bl]
    judged_syst = [judge_block(b, lambda c, b=b: tree.weight(
        DEFAULT_CONFIG, b["ps"], c, b["w_star"], b["q2"]), errs="stat+syst")
        for b in bl]
    ps_values = sorted({b["ps"] for b in bl})
    return dict(
        census=cen,
        judged=judged,
        judged_stat_syst=judged_syst,
        forward=forward_record(judged[DEFAULT_CONFIG]),
        prescription=prescription_record(tree, ps_values),
        f2n_f2p=f2n_f2p_record(doc, tree),
        weight_symmetric=max(abs(tree.weight(DEFAULT_CONFIG, ps, c, 2.44, 1.66)
                                 - tree.weight(DEFAULT_CONFIG, ps, -c, 2.44, 1.66))
                             for ps in ps_values for c in (0.1, 0.5, 0.9)),
        fsi_sigma_xn_mb={k: tree.fsi[k].sigma_eff_mb(2.44) for k in CONFIGS},
        fsi_clipped_grid_fraction=tree.fsi[DEFAULT_CONFIG].clipped_grid_fraction(),
        figure=figure_series(doc, tree),
    )


def figure_series(doc, tree):
    """(f) what the figure draws, computed with the row (so REPORT.json
    carries every number the figure shows)."""
    cs = [-0.95 + 0.05 * i for i in range(39)]
    lines = collections.OrderedDict()
    for name in CONFIGS:
        w = [tree.weight(name, FIGURE_PS, c, 2.44, FIGURE_Q2) for c in cs]
        w_bw = [tree.weight(name, FIGURE_PS, c, 2.44, FIGURE_Q2)
                for c in (-0.9, -0.7, -0.5, -0.3)]
        mean = sum(w_bw) / len(w_bw)
        lines[name] = [x / mean for x in w]
    pts = collections.OrderedDict()
    for b in blocks(doc):
        if b["q2"] != FIGURE_Q2 or b["ps"] != FIGURE_PS or not in_dis_window(b):
            continue
        bw = backward(b)
        if len(bw) < MIN_BACKWARD_POINTS:
            continue
        band = flat_band([r[1] for r in bw], [r[2] for r in bw])
        pts["%s/W*=%g" % (b["beam"], b["w_star"])] = dict(
            cos=[r[0] for r in b["rows"]],
            r_over_rbar=[r[1] / band["rbar"] for r in b["rows"]],
            err=[r[2] / band["rbar"] for r in b["rows"]],
            delta=band["delta"])
    return dict(cos=cs, lines=lines, blocks=pts)


def _beam_subrow(s, config, in_verdict, errs="stat"):
    wave, sigma, _ = CONFIGS[config]
    return dict(status=s["status"], in_verdict=in_verdict,
                wave=WAVE_LABELS[wave], sigma_xn_mb=sigma, errors=errs,
                n_blocks=s["n_blocks"], n_pass=s["n_pass"],
                worst_block=s["worst_block"],
                worst_v_over_delta=s["worst_v_over_delta"],
                chi2_flat_sum=s["chi2_flat_sum"],
                chi2_tree_sum=s["chi2_tree_sum"], ndf_sum=s["ndf_sum"],
                ps=s["ps"], v_by_ps=s["v_by_ps"],
                min_delta_by_ps=s["min_delta_by_ps"],
                n_pass_by_ps=s["n_pass_by_ps"])


def run(path=BONUS, verbose=True):
    m = measure(path)
    sub = collections.OrderedDict()
    summ = {}
    for name, (_, _, in_v) in CONFIGS.items():
        for beam in BEAMS.values():
            s = summarise([j for j in m["judged"][name] if j["beam"] == beam])
            summ[(name, beam)] = s
            sub["%s/%s" % (beam, name)] = _beam_subrow(s, name, in_v)
    for beam in BEAMS.values():
        s = summarise([j for j in m["judged_stat_syst"] if j["beam"] == beam])
        summ[("stat+syst", beam)] = s
        sub["%s/%s stat+syst" % (beam, DEFAULT_CONFIG)] = _beam_subrow(
            s, DEFAULT_CONFIG, False, errs="stat (+) syst")
    fwd = m["forward"]
    sub["forward rise (recorded)"] = dict(
        status="recorded", in_verdict=False,
        **{k: "rho_data %.3f +- %.3f, rho_tree %.3f..%.3f, same sign in %d/%d"
           % (v["rho_data"], v["rho_data_err"], v["rho_tree_min"],
              v["rho_tree_max"], v["same_sign_blocks"], v["n"])
           for k, v in fwd.items()})
    pr = m["prescription"]
    sub["PWIA prescription instant/LC (recorded, not FSI)"] = dict(
        status="recorded", in_verdict=False, wave=WAVE_LABELS["hulthen"],
        max_dev_by_ps=[v["max_dev"] for v in pr.values()],
        ends_ratio_by_ps=[v["ends_ratio"] for v in pr.values()],
        weight_ends_ratio_by_ps=[v["weight_ends_ratio"] for v in pr.values()])
    f2 = m["f2n_f2p"]
    sub["F2n/F2p ToyF2 vs W*>1.8 table (recorded)"] = dict(
        status="recorded", in_verdict=False,
        n=f2["all"]["n"], chi2_abs=f2["all"]["chi2_abs"],
        p_abs=f2["all"]["p_abs"], scale=f2["all"]["scale"],
        chi2_scaled=f2["all"]["chi2_scaled"],
        ndf_scaled=f2["all"]["ndf_scaled"],
        tree_over_data_first=f2["all"]["tree_over_data"][0],
        tree_over_data_last=f2["all"]["tree_over_data"][-1],
        n_w_ge_2=f2["w_ge_2"]["n"], chi2_abs_w_ge_2=f2["w_ge_2"]["chi2_abs"],
        chi2_scaled_w_ge_2=f2["w_ge_2"]["chi2_scaled"])
    c = m["census"]
    sub["out of model (census)"] = dict(
        status="recorded", in_verdict=False, blocks=c["blocks"],
        dis_window=c["dis_window"], verdict_blocks=c["verdict_blocks"],
        dis_no_backward=len(c["dis_no_backward"]),
        out_of_model=c["out_of_model"],
        out_of_model_q2_below_1=c["out_of_model_q2_below"],
        out_of_model_resonance=c["out_of_model_resonance"])

    heads = [summ[(DEFAULT_CONFIG, beam)] for beam in BEAMS.values()]
    status = "pass" if all(h["status"] == "pass" for h in heads) else "fail"
    h4, h5 = heads
    report = dict(
        name=NAME,
        generator=(
            "Hulthen (rule (d)'s configuration; AV18 ships since D6 and fails too) "
            "40 mb GlauberFsiWeight over its backward mean inside the "
            "data's 1-sigma flat band in %d/%d (4 GeV) and %d/%d (5 GeV) DIS "
            "blocks; worst v/delta %.2f (4 GeV %s), %.2f (5 GeV %s); backward "
            "chi2 flat %.1f, tree shape %.1f on %d dof"
            % (h4["n_pass"], h4["n_blocks"], h5["n_pass"], h5["n_blocks"],
               h4["worst_v_over_delta"], h4["worst_block"],
               h5["worst_v_over_delta"], h5["worst_block"],
               h4["chi2_flat_sum"] + h5["chi2_flat_sum"],
               h4["chi2_tree_sum"] + h5["chi2_tree_sum"],
               h4["ndf_sum"] + h5["ndf_sum"])),
        reference=(
            "CLAS BONuS R_D/S vs cos theta_pq, Tkachenko PRC 89 (2014) 045206 "
            "Suppl. (%d blocks; %d in the DIS window W* >= 2, Q2 >= 1, %d with "
            "backward bins)" % (c["blocks"], c["dis_window"],
                                c["verdict_blocks"])),
        tolerance=("max over cos <= -0.3 of abs(w/mean_bw(w) - 1) <= delta = "
                   "S sigma(Rbar)/Rbar (stat, S = max(1, sqrt(chi2/ndf))), "
                   "every block, both beams"),
        status=status,
        config=dict(wave=WAVE_LABELS["hulthen"], sigma_xn_mb=40.0,
                    w_star_min=W_STAR_MIN, q2_min=Q2_MIN,
                    backward_cos_max=BACKWARD_COS_MAX,
                    forward_cos_gt=FORWARD_COS_GT, errors="stat",
                    weight_entry="lipolgen.GlauberFsiWeight.weight (library "
                                 "evaluation, no generation)"),
        subrows=sub,
        measured=m,
        history=dict(blocked_2026_09_23=BLOCKED_REASON, unblocked=UNBLOCKED),
        findings=findings(),
    )
    if verbose:
        _print(report)
    return report


def _print(report):
    say = lambda line: print(_printable(line))    # noqa: E731
    m, f = report["measured"], report["findings"]
    c = m["census"]
    say("HISTORY (2026-09-23): CLAS DB, deuteron target: %d experiments, %d "
        "measurements" % (f["n_experiments"], f["n_measurements"]))
    b = f["bonus"]
    say("  BONuS eid 135: %d measurement(s), quantities %s" % (
        b["n_measurements"], b["quantities"]))
    for e in f["eg1b"]:
        say("  EG1b eid %d: %d measurements, quantities %s" % (
            e["eid"], e["n_measurements"], e["quantities"]))
    d = f["deeps"]
    say("  Deeps eid 90: %d measurements %s; %d F2xP(p_s) blocks vendored "
        "(row 11): %d at the paper's Q2 = 1.8/2.8, %d with unverified Q2 "
        "labels (%d label triples duplicated), empty block(s) mid %s, %d "
        "value = stat = 0 rows (missing, not zeros)"
        % (d["n_measurements"], d["quantities"], f["deeps_blocks"],
           f["deeps_blocks_paper_q2"], f["deeps_blocks_unverified_q2"],
           f["deeps_duplicate_label_triples"], f["deeps_empty_blocks"],
           f["deeps_zero_value_zero_stat_rows"]))
    say("BONuS supplement: %d R_D/S blocks (%s), %d rows, %d missing in %d "
        "blocks; DIS window W* >= %g, Q2 >= %g: %d blocks, %d with backward "
        "bins (%s), %d without %s; out of model %d (%d at Q2 < 1, %d "
        "resonance region); F2n ratio tables %s"
        % (c["blocks"], c["blocks_by_beam"], c["rows"], c["missing_rows"],
           len(c["blocks_with_missing"]), W_STAR_MIN, Q2_MIN, c["dis_window"],
           c["verdict_blocks"], c["verdict_blocks_by_beam"],
           len(c["dis_no_backward"]), c["dis_no_backward"], c["out_of_model"],
           c["out_of_model_q2_below"], c["out_of_model_resonance"],
           c["f2n_ratio_tables"]))
    say("GlauberFsiWeight: sigma_XN %s mb, clipped grid fraction %.4g, "
        "max |w(c) - w(-c)| = %.3g (symmetric by construction)"
        % (m["fsi_sigma_xn_mb"], m["fsi_clipped_grid_fraction"],
           m["weight_symmetric"]))
    say("%-30s %-6s %7s %18s %22s %16s" % (
        "sub-row", "", "blocks", "worst v/delta", "chi2 flat / tree / ndf",
        "verdict"))
    for name, s in report["subrows"].items():
        if "n_blocks" not in s:
            continue
        say("%-30s %-6s %3d/%-3d %8.3f %-9s %7.2f %7.2f %4d %16s" % (
            name, "", s["n_pass"], s["n_blocks"], s["worst_v_over_delta"], "",
            s["chi2_flat_sum"], s["chi2_tree_sum"], s["ndf_sum"],
            s["status"].upper() if s["in_verdict"] else
            "(%s)" % s["status"]))
        say("    p_s %s: tree v %s, data min delta %s, pass %s" % (
            s["ps"],
            " ".join("%.4f" % x for x in s["v_by_ps"]),
            " ".join("%.4f" % x for x in s["min_delta_by_ps"]),
            s["n_pass_by_ps"]))
    say("default, per block (beam Q2 W* p_s: Rbar +- sigma, chi2_flat/ndf, "
        "S, delta | tree v, chi2_tree | pass):")
    for j in m["judged"][DEFAULT_CONFIG]:
        say("  %s %4.2f %4.2f %5.3f: %.3f +- %.3f  %5.2f/%d  S %.2f  delta "
            "%.4f | v %.4f  %5.2f | %s" % (
                j["beam"], j["q2"], j["w_star"], j["ps"], j["rbar"],
                j["sigma"], j["chi2_flat"], j["ndf"], j["scale_factor"],
                j["delta"], j["v"], j["chi2_tree"],
                "pass" if j["passed"] else "FAIL"))
    say("forward rise (recorded; cos > 0 over cos <= -0.3, one scale per "
        "block):")
    for k, v in m["forward"].items():
        say("  %-14s rho_data %.3f +- %.3f   rho_tree %.4f..%.4f   same sign "
            "%d/%d" % (k, v["rho_data"], v["rho_data_err"], v["rho_tree_min"],
                       v["rho_tree_max"], v["same_sign_blocks"], v["n"]))
    say("PWIA prescription (recorded, NOT FSI): Hulthen instant-form / "
        "light-cone over the backward window, max deviation | f(-0.9)/f(-0.3)"
        " | the weight's w(-0.9)/w(-0.3):")
    for k, v in m["prescription"].items():
        say("  %-9s %.4f | %.4f | %.4f" % (k, v["max_dev"], v["ends_ratio"],
                                          v["weight_ends_ratio"]))
    for label, r in m["f2n_f2p"].items():
        say("F2n/F2p (recorded) ToyF2 vs BONuS %s [%s, %d rows]: chi2 %.1f "
            "(p = %.3g); one free scale %.4f -> chi2 %.2f/%d (p = %.3g); "
            "tree/data %s" % (
                F2N_F2P_CUT, label, r["n"], r["chi2_abs"], r["p_abs"],
                r["scale"], r["chi2_scaled"], r["ndf_scaled"], r["p_scaled"],
                " ".join("%.3f" % x for x in r["tree_over_data"])))
    say("REPORT | %s | %s | %s | %s | %s" % (
        NAME, report["generator"], report["reference"], report["tolerance"],
        report["status"].upper()))


# ------------------------------------------------------------------ figure

def plot(report, out_dir=None):
    """(f) R_D/S over its block's backward Rbar against cos theta_pq for the
    two highest-W* blocks per beam at Q2 = 1.66, p_s = 0.093 (ink markers,
    stat errors), the three configurations' normalised weight u(c) as
    lines, the backward window shaded."""
    sys.path.insert(0, HERE)
    from _plots import BAND, new_figure, save, status_title, style_series
    fig, ax = new_figure(two_panel=False, height=4.0)
    fs = report["measured"]["figure"]
    ax.axvspan(-1.0, -0.2, color=BAND, alpha=0.6, linewidth=0, zorder=0)
    markers = {"4GeV/W*=2.03": "o", "4GeV/W*=2.44": "s",
               "5GeV/W*=2.03": "^", "5GeV/W*=2.44": "D"}
    offs = {"4GeV/W*=2.03": -0.024, "4GeV/W*=2.44": -0.008,
            "5GeV/W*=2.03": 0.008, "5GeV/W*=2.44": 0.024}
    for key, blk in fs["blocks"].items():
        beam, w = key.split("/")
        kw = style_series("reference", marker=markers.get(key, "o"),
                          label="BONuS %s, %s" % (beam.replace("GeV", " GeV"),
                                                  w.replace("W*=", "W* = ")))
        if markers.get(key) in ("s", "D"):
            kw.update(markerfacecolor="white")
        ax.errorbar([c + offs.get(key, 0.0) for c in blk["cos"]],
                    blk["r_over_rbar"], yerr=blk["err"], **kw)
    slots = {"hulthen-40mb": ("default", "Glauber, Hulthen, 40 mb (verdict)"),
             "av18-40mb": ("alt", "Glauber, AV18, 40 mb"),
             "hulthen-20mb": ("alt2", "Glauber, Hulthen, 20 mb")}
    for name in sorted(fs["lines"], key=lambda n: n == DEFAULT_CONFIG):
        slot, lab = slots[name]                 # the verdict's line on top
        kw = style_series(slot, label=lab)
        if name == DEFAULT_CONFIG:
            kw.update(zorder=4)
        ax.plot(fs["cos"], fs["lines"][name], **kw)
    ax.axhline(1.0, color="#c3c2b7", linewidth=1.0, zorder=0)
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(0.66, 1.22)
    ax.set_xlabel("cos theta_pq  (spectator proton vs q, deuteron rest frame)")
    ax.set_ylabel("R_D/S, and the weight w, / backward mean")
    ax.legend(loc="lower center", ncol=2, fontsize=7)
    ax.text(-0.98, 1.205, "backward window", fontsize=7, va="top",
            color="#52514e")
    ax.text(0.98, 1.205, "shown: Q2 = 1.66 GeV2, p_s = 85-100 MeV/c",
            fontsize=7, va="top", ha="right", color="#52514e")
    h = report["subrows"]
    a, b = h["4GeV/hulthen-40mb"], h["5GeV/hulthen-40mb"]
    status_title(fig, NAME, report["status"],
                 "W*>=2, Q2>=1, cos<=-0.3: in flat 1-sigma band %d/%d (4 GeV), "
                 "%d/%d (5 GeV); chi2 flat %.0f, tree %.0f /%d"
                 % (a["n_pass"], a["n_blocks"], b["n_pass"], b["n_blocks"],
                    a["chi2_flat_sum"] + b["chi2_flat_sum"],
                    a["chi2_tree_sum"] + b["chi2_tree_sum"],
                    a["ndf_sum"] + b["ndf_sum"]))
    return save(fig, NAME, out_dir)


def main(argv=None):
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    argv = sys.argv[1:] if argv is None else argv
    try:
        rep = run()
        if "--plot" in argv:
            print(_printable("figure: %s" % plot(rep)))
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
