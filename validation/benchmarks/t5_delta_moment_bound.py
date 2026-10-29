#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""docs/benchmarking/08_channel_priority.md sec. 4 ingredient C (rank 1, the
letter's claim 1) / run 2026-10-26 PLAN.md phase C2 / 04_theory.md sec. 4,
T-18 -- the one bound that exists on the Delta sector: the n = 2 moment of
the gluonic Soffer-type bound, applied to LiPolGen's double-helicity-flip
structure function Delta(x, Q2) on the 6Li target, AS A BOUND UNDER STATED
ASSUMPTIONS (tier T5).

A BOUND SAYS A SCENARIO IS ALLOWED, NEVER THAT IT IS RIGHT (04_theory.md
sec. 11 item 2).  A PASS here says the shipped toy Delta does not ask for
more linearly polarized glue than positivity permits under the assumptions
below; it says nothing about whether exotic glue in 6Li has that size, that
shape or that sign.

TWO PASSES.  The first pass (2026-10-26) wrote this header down to TOLERANCE
before computing any saturation, ran the row and recorded it.  A second pass
(the same day, AFTER every first-pass number had been computed and seen)
re-read the sources and found two things, written under SECOND PASS below:
the variation V4 is a misreading, and the first pass's normalisation A7
disagrees with DS16's own definition of Delta by a factor 2.  The second
pass changes no central assumption, no central number, not the verdict rule
and not its tolerance; it adds the source KS20, the variation V6 and the
fields that carry it, moves V4 out of the LOOSE edge and V6 into the TIGHT
one.  D4 (2026-10-27, the author's decision) then RELABELS, moving no rule, tolerance, edge or number of rule (a): the D4 note above main().

WHAT IS COMPARED.  For each Delta model, at Q2 = 2, 5, 10 GeV2, the n = 2
Mellin moment (Detmold-Shanahan [DS16] Eq. (12) at n = 2)

    M2(Q2) = int_0^1 x Delta(x, Q2) dx            (= the plain integral of
                                                     the x Delta the letter
                                                     extracts)

against the bound's right-hand side

    B(Q2) = (alpha_s(Q2) / (8 pi)) Tr Q^2 * (1/2) int_0^1 x [ f1^g(x, Q2)
              + (1/2) f1LL^g(x, Q2) + g1^g(x, Q2) ] dx ,

with the saturation fraction S = |M2|/B (how much of the allowed room the
scenario uses -- the number the letter can quote) and the margin 1 - S.
B follows from two equations of DS16, read off the page images:

  * Eq. (13), DS16's restatement of the Jaffe-Manohar parton-model relation:
        Delta(x, Q2) = -(alpha_s/(2 pi)) Tr Q^2 x^2 int_x^1 (dy/y^3) deltaG(y)
    whose n = 2 moment (exchange the integrals: int_0^y x^3 dx = y^4/4) is
        M2 = -(alpha_s/(8 pi)) Tr Q^2 int_0^1 y deltaG(y) dy ,
    the same as DS16 Eqs. (9) + (11) + (12) at n = 2 (M2 = C_2 A_2/2,
    C_2 = -(alpha_s/(2 pi)) Tr Q^2 * 2/4, A_2 = int x deltaG dx);
  * Eq. (32), the x-space bound on the gluon transversity of a spin-1 target,
        |deltaG(x)| <= (1/2) [ f1(x) + (1/2) f1LL(x) + g1(x) ] ,
    published as a positivity bound by Cotogno, van Daal, Mulders [CvDM17]
    Eq. (36) (deltaG = h1TT, DS16 p. 7).  Multiplying by x >= 0 and
    integrating, |A_2| <= int x |deltaG| <= (1/2) int x [...] dx: that is
    the "n = 2 moment form" this row applies.  DS16's own Eq. (36),
    |A_2| <= (5 B_{2,1} - 6 B_{2,2})/24, is the same statement written for
    the phi meson in lattice reduced matrix elements (the g1 term vanishes
    for the phi by operator symmetry; the map from B_{2,1}, B_{2,2} to the
    moments of f1 and f1LL is not printed), so it is quoted, not evaluated.

THE MODELS (which Delta, in which normalisation -- read from the kernel).
The kernel builds the cos 2phi term as a2 += -(1/2) ((1-y)/y^2) c_eff sin^2(theta_S)
Delta/D_phi since D4, 2026-10-27 (no 1/2 before: the D4 note) (src/core/xsec.cpp `tensor_amplitudes`, c_eff = 3 Q_NN), with
Delta = delta_func(x, Q2, f1) evaluated on f1 = F1A/A, the PER-NUCLEON F1
(`InclusiveKernel::tables`), and D_phi per nucleon.  This row therefore takes
Delta PER NUCLEON exactly as the kernel carries it, in the per-nucleon x.
  1. SHIPPED (the verdict): `default_inclusive_kernel(li6())`'s Delta --
     toy_delta_gluon(x, Q2, F1, 1e-2) = 1e-2 F1 x^0.3 (1-x)^4, the pipeline's
     hard-coded 1e-2 "discovery scale" (src/core/pipeline.cpp, the Delta slot
     after the b1 branches; `rank2_input_report`), NOT the function's default
     1e-3; F1 the shipped ToyF2 per nucleon.  The tree applies NO dilution to
     it: P_zz is the whole-nucleus fill and nothing multiplies Delta.  The
     harness reads Delta from the kernel's own tables and records the scale
     it finds (Delta / (F1 x^0.3 (1-x)^4) at every node).
  2. THE LETTER'S moment_A and moment_B (PolarizedLithiumSim plans/07 sec.
     7.1: "moment_A, P_zz = 0.6, dilution 1/3"; the sibling's
     polli_fastsim/delta_models.py, read, never run), through the tree's
     port `make_moment_a` / `make_moment_b` (include/lipolgen/sf.hpp):
     variant mid_x, c = C_BAG = -0.012 (Sather-Schmidt [SS90] Eq. (3.5), the
     bag value for the Delta++), alpha_s the CT18NLO table (the letter's
     production convention), dilution 1/3 -- the counting fraction the
     sibling folds into P_zz = 0.8/3 and the tree expresses on Delta with
     P_zz whole-nucleus (sf.hpp, the `dilution` comment), so with dilution
     1/3 the model IS the letter's physical per-nucleon Delta in the tree's
     normalisation; moment_A solved on the shipped kernel's F1 at q2_ref =
     4.474 GeV2, the letter's toy-column <Q2> (plans/07 addendum: A =
     -0.2924; the harness reproduces it and records the check).  Each is
     RECORDED beside the verdict with its own outcome of rule (a) under the
     central assumptions and under V6 (no run path of the generator uses
     them: the pipeline hard-codes the toy and no CLI flag selects a Delta
     model).
  3. RECORDED: moment_A and moment_B at dilution 1 (the library default:
     the bag value per nucleon, undiluted).

THE ASSUMPTIONS, DECLARED HERE BEFORE ANY NUMBER (first pass; central values
first; the second pass's findings on A1 and A7 are under SECOND PASS).
  A1  the relation: DS16 Eq. (13) as above, with the deltaG that DS16 Eq.
      (32) bounds (internal consistency of one paper).  [First pass: "a
      second printed normalisation exists: SS90 Eq. (2.11) writes Delta =
      (alpha_s/(2 pi)) Tr Q^2 x^2 int (dy/y^3) a(y) with a = -2 (g_xhat -
      g_yhat), TWICE DS16 Eq. (13) for the same g_xhat - g_yhat ... recorded
      as variation V4, never the central reading."  The second pass shows
      the two are ONE relation: S1.]
  A2  Tr Q^2 = 2/3: three flavours, Q = diag(2/3, -1/3, -1/3) -- SS90 above
      Eq. (3.5) (the same paper C_BAG comes from) and LOI18 Eq. (7).  DS16
      gives Tr Q^2 no value.  Four flavours (10/9) is variation V3.
  A3  alpha_s(Q2) = CT18NLO member 0 (LHAPDF's alphasQ2, the set's ipol
      table), the same set as f1^g; the letter's moment models carry the
      same function.
  A4  f1^g = the CT18NLO gluon (pid 21, member 0) per nucleon, FREE nucleon:
      no nuclear modification centrally (EPPS21's 6Li gluon ratio is
      variation V5).  Read from LHAPDF 6.5.5 itself through PYTHIA 8.317's
      LHAPDF6 plugin (`_report.import_pythia8`, then PDF::xfUpdate +
      PDF::xfRaw, so PYTHIA's max(0, .) clip in PDF::xf never acts): the
      tree's LhapdfSF / LhapdfG1 expose F2 and g1, never a parton density,
      and the lhapdf Python module is not built.  Guard: the plugin's
      sum e_q^2 x(q + qbar) reproduces LhapdfSF('CT18NLO').f2p, and its
      3-flavour g1 LhapdfG1.g1p, at six points each (|rel| <= 1e-12, else
      the harness is broken: exit 2); the set files' sha256 must equal the
      vendored record (else BLOCKED, environment: a different set version).
  A5  g1^g = P_eff * Delta g, Delta g = NNPDFpol11_100 member 0 (the replica
      average) per nucleon, P_eff = (Z eff_pol_p + N eff_pol_n)/A of the
      tree's li6() = 0.2704 (the cluster picture's 0.81123/3, the convention
      `PolSF::g1_nucleus` uses), i.e. the 6Li m = +1 state; integrated on
      the set's own grid, x >= 1e-5.  Delta g = 0 is the other band edge (V1).
  A6  f1LL^g = 0.  Band (V2): f1LL^g = +- f1^g, the largest band symmetric
      about zero that CvDM17 Eqs. (34)-(35) allow (it needs |g1^g| <= f1^g/2:
      checked on the fine rule's nodes; the first x where it fails, the first
      where |g1^g| > f1^g -- the mixed central inputs leaving CvDM17 (34)
      themselves -- and the share of <x>_g above, recorded).  Positivity
      alone confines f1LL^g to [-2 (f1^g - |g1^g|), f1^g], at whose lower
      edge the bound collapses to (|g1| + g1)/2: positivity does not bound
      f1LL^g usefully, so the band is a statement of ignorance, not a
      calculation.
  A7  normalisation (first pass): the tree's Delta is the Delta of DS16
      Eqs. (4)-(5) (taken by DS16 and LOI18 from Jaffe-Manohar): LOI18 Eq.
      (9) is dsigma(lambda_m) ~ F1 + ... - ((1-y)/y^2) b_m sin^2(theta_m)
      Delta cos(2 phi), b_m = 3|lambda_m| - 2 -- the kernel's a2 term with
      c_eff = 3m^2 - 2 (docs/PHYSICS_CHANNELS.md sec. 2).  [The second pass
      finds DS16's Eqs. (3)-(5) give Delta_DS16 = 2 Delta_tree instead: S2,
      variation V6.]
  A8  per nucleon on both sides: for A nucleons, x = A x_A; per-nucleon
      F1-like densities are D(x) = D_A(x_A)/A^2 (the kernel's F1 = F1A/A
      with F1A built on the per-nucleon x is F1_A(x_A)/A^2), the same for
      Delta, deltaG, f1^g, g1^g.  Both Eq. (13) and Eq. (32) are homogeneous
      of degree one in the densities and x^2 int dy/y^3 is invariant under
      x -> A x, so they hold unchanged per nucleon in the per-nucleon x
      (upper limit A for 1: the support is x < 1 here).
  The integrals: the composite Gauss-Legendre rule of the vendored grid
  (144 ln-x panels on [1e-9, 0.1], 72 linear panels on [0.1, 1], 8 nodes
  each; the set's own x_min for NNPDFpol11), repeated on twice the panels;
  the finer value is used, the difference is its numerical error, and a
  relative error above QUAD_REL_TOL = 1e-7 makes the harness broken (exit 2).
  (The first rule, 72 + 36 panels, failed that check on the harness's first
  run -- 1.7e-7 on <x>_g at Q2 = 2, LHAPDF's interpolant being only C1 at
  its knots -- before any saturation was computed; the panels were doubled.
  Numerics only: the plateau is ~1e-9 relative.)  Below x_lo the tail is
  not integrated (x_lo |f(x_lo)| is recorded).  Closure: moment_B's M2 is
  dilution x C_BAG x alpha_s(Q2) in closed form, and the integrated value
  must meet it to 1e-9 relative (else broken).

SECOND PASS (written after the first pass's numbers were seen; it moves no
central number and not the verdict rule):
  S1  V4 is a misreading, not a second normalisation.  SS90 write Eq. (2.11)
      for a spin-1 target circularly polarized ("aligned") along x in the
      infinite-momentum frame (SS90 p. 1425); DS16's deltaG is read by KS20
      as the linear-polarization difference in the LINEARLY polarized hadron
      -- KS20 Eq. (68) lists DS16's deltaG and SS90's function on the line
      of h1TT = -g_xhat/xhat + g_yhat/xhat, hadron polarization xhat being
      the linear E_x of KS20 Eq. (22), the one state in which DS16's own
      identification deltaG = h1TT holds in magnitude.  For spin 1 the
      transverse tensor polarization
      T_xx - T_yy is +1/2 in the circular (m_x = +-1) state and -1 in the
      linear (E_x) one, so (g_x - g_y)_circular = -(1/2)(g_x - g_y)_linear,
      SS90's a = -2 (g_x - g_y)_circular = (g_x - g_y)_linear, and SS90
      (2.11) is DS16 (13).  V4 stays a recorded number (what the misreading
      gives) and leaves the LOOSE edge.
  S2  A7 disagrees with DS16's own definition of Delta by a factor 2.  DS16
      Eq. (3) (projection operators P(hH,h'H')_{mu nu,alpha beta} =
      eps*_mu(h') E*_alpha(H') E_beta(H) eps_nu(h)) and Eq. (4) (Delta =
      A_{+-,-+} = A_{-+,+-}) give Eq. (5)'s tensor (W_xx = -W_yy = Delta/2
      for E = xhat) and, for an unpolarized target, W_ij = delta_ij (1/3)
      sum_H A_{+H,+H}: the standard -F1 g_{mu nu}.  Contracted with the
      unpolarized lepton tensor, Eq. (5) gives for a pure m = +-1 state
      along an axis perpendicular to q the cos 2phi term -(1/2) ((1-y)/y^2)
      Delta relative to F1 + (1-y) F2/(x y^2) -- exactly LOI18 Eq. (6), the
      Bjorken limit the LOI prints for a target with its spin along xhat --
      while the pre-D4 kernel and LOI18 Eq. (9) (which A7 read) carry
      -((1-y)/y^2) Delta there (c_eff = 1).  Per unit cos 2phi amplitude,
      Delta_DS16 = DS16_OVER_TREE x Delta_tree with DS16_OVER_TREE = 2 (the pre-D4 kernel's; 1 since D4).  The
      pytest re-derives the factor numerically from Eqs. (3)-(5).  Recorded
      as V6: every model's M2 in DS16's normalisation is 2 M2, so every
      saturation doubles.  The shipped toy passes under V6 too; which
      reading the tree's (and the letter's) Delta follows is a property of
      the kernel's cos 2phi coefficient and the author's to settle (the
      phase record says what each costs) -- this row records both.
  S3  the edges follow: TIGHT = V1 + V2(-) + V6 on the smaller of the free
      and EPPS21 gluon moments; LOOSE = V2(+) + V3 on the larger, with the
      NNPDF Delta g term only if it widens B.  (First pass: TIGHT without
      V6, LOOSE with V4; its edge values are in the phase record.)

VARIATIONS, RECORDED (never the verdict): V1 Delta g = 0; V2 f1LL^g = -f1^g
and +f1^g; V3 Tr Q^2 = 10/9; V4 SS90's coefficient applied to DS16's deltaG
(x2 room; a misreading, S1 -- in no edge); V5 f1^g = R_g^EPPS21(6Li) x
CT18NLO gluon (lipolgen's Epps21Ratio, proton baseline CT18NLO); V6
Delta_DS16 = 2 Delta_tree (S2: every saturation x2; since D4 the pre-D4 kernel's / the letter's printed reading).  For the letter's models
the saturation with the library-default LO alpha_s (`alpha_s_lo`, Lambda
0.22, n_f 4) on the model side is recorded too; and for every model the
constant fraction f1LL^g/f1^g at which |M2| = B (all else central), with and
without V6: the scenario is allowed iff the 6Li gluon's tensor polarization
lies above it.  Registry row 23 (`LI6_CLUSTER_POLARIZATION` against the naive
1/3) enters only through P_eff in g1^g: B with P_eff = 1/3 is recorded.

TOLERANCE, DECLARED BEFORE THE VERDICT (rule (a) fixed by the task brief,
PLAN.md C2; rule (b) the first pass's, with the second, sourced sub-row the
brief asked for; both before any number was computed; unchanged by the
second pass):

    the row PASSES iff, under the central assumptions A1-A8,
      (a) the SHIPPED Delta satisfies |M2| <= B with margin 1 - |M2|/B
          >= -1e-12 at EACH of Q2 = 2, 5, 10 GeV2, and
      (b) every helicity-state rate of the shipped 6Li kernel is >= -1e-12
          relative to the unpolarized rate (i.e. -2 <= A_T <= 1) at every
          point of the vendored positivity grid -- the SECOND, sourced
          sub-row: Cosyn, Roldan Tomei, Sosa, Zec [CRSZ25] Eq. (2), "A^T =
          [dsigma(+1) + dsigma(-1) - 2 dsigma(0)]/[dsigma(+1) + dsigma(0) +
          dsigma(-1)], [-2 <= A^T <= 1]" -- the helicity-state positivity of
          the spin-1 cross sections, here on the generator's own per-state
          W_m(phi') = 1 + w_avg + a1 cos phi' + a2 cos 2phi' (unpolarized
          electrons, pure m = +1, 0, -1 along an axis at theta_S = 0, 45, 90
          deg; x, Q2, y, phi' on the vendored grid), which carries b1, b2 AND
          Delta.
    The letter's moment_A / moment_B get the outcome of rule (a), recorded,
    under the central assumptions and under V6.

PREVIEW, STATED AS THE RULES REQUIRE.  Before this header was written the
first pass made a back-of-envelope estimate (gluon momentum fraction ~0.42 at
these Q2): toy saturation of order 0.1; for the bag-sum-rule models of order
0.7 at dilution 1/3 and 2 at dilution 1, with Tr Q^2 = 2/3 and DS16's
normalisation.  The machinery probes (the plugin against LhapdfSF /
LhapdfG1, the kernel's tables) were run before the header.  The tolerance,
the central assumptions and the band edges were fixed before the first pass
computed any saturation; Tr Q^2 = 2/3 follows SS90 and LOI18, and no central
assumption was changed after a number was seen.  The second pass read every
first-pass number before writing SECOND PASS; S1-S3 change no central
number, no verdict and no tolerance.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- the moment bound has no b1.  The sub-row (b) reads
     the shipped b1(6Li) only through the kernel's own per-state rates,
     whatever its normalisation (Li6B1(MillerB1): 0.921947 x (2/6) x Miller
     per nucleon, per nucleon like F1), so no b1 is compared to another.
  2. isoscalar denominator -- NOT APPLICABLE as a structure-function
     denominator; the gluon is isospin symmetric (f1^g, Delta g of the
     proton = of the neutron) and 6Li has Z = N, so per nucleon = the
     proton's gluon (A4, A5).
  3. per nucleus / per nucleon -- every structure function and density on
     BOTH sides is PER NUCLEON in the per-nucleon x (A8); the per-state
     rates of (b) are ratios, normalisation-free.  (The factor of V6 is not
     a per-nucleon question: it is the cross-section normalisation of Delta
     itself, S2.)

WHAT THIS DOES NOT VALIDATE.  It does not say Delta(6Li) has the toy's size,
shape or sign (a bound, not a calculation; the lattice A_2 of DS16 is a phi
meson's and is not compared: 04_theory.md sec. 4); it does not decide which
of A7 and V6 the tree's Delta follows (the author did, D4: A7; S2 records both); it applies the
MOMENT form only -- the x-space bound Eq. (32) itself is stronger, and a
model can pass the moment and fail pointwise (any Delta falling slower than
the gluon as x -> 1 does); it propagates no PDF uncertainty (CT18NLO's 58
eigenvectors, NNPDFpol1.1's 100 replicas -- the declared Delta g = 0 edge and
the f1LL band are the ones the result is read against); the n = 4, 6, ...
moments; any nuclear modification beyond V5; the RC of the Delta term (none
exists, `RcScope::TensorAll`); the finite-gamma leakage of b1 into the
cos 2phi amplitude; 7Li (its rank-2 sector is identically zero, and T-13's
spin-3/2 bound is not this row); b3/b4.

Run:  source env.sh && python3 validation/benchmarks/t5_delta_moment_bound.py [--plot]

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike (BLOCKED,
with an `environment:` reason, when PYTHIA's python module or its LHAPDF6
plugin is missing or cannot load a set, or a set's member file is missing
or a set file is not the vendored version); 2 when the harness itself is
broken (a missing or altered vendored file, lipolgen missing, a failed
same-grid guard, an unconverged integral, any exception) -- and when one of
the two sets is not installed where LHAPDF searches: the row's only
configuration needs both, so it raises LHAPDF's own "Info file not found for
PDF set" message and the report records it as blocked (environment), as row
4 (README convention 6; a BLOCKED row with exit 0 until phase F of run
2026-10-26).  A harness is a measurement, not a CI gate: the pytests are the
gate.
"""

import hashlib
import json
import math
import os
import subprocess
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "delta_moment_bound_sources.json")
#: sha256 of every byte of the vendored file, and of its canonical body
#: (json.dumps(body, sort_keys=True, separators=(",", ":"))), which the file
#: also records.  Both are re-checked on EVERY read.
FILE_SHA256 = "794d1429662a41bf717a7269071d9e3bdafa6c8f4933a92557e38807912218e7"
BODY_SHA256 = "d355287d1cbb0fd06269ed4ee5af58e5b1d8790b8aa314bce112cdf109364594"
#: the sources, in the vendored order (KS20 added by the second pass)
SOURCE_KEYS = ["DS16", "SS90", "LOI18", "CvDM17", "CRSZ25", "KS20"]

NAME = "t5_delta_moment_bound"
#: a few seconds, no event generation, no generator-vs-generator sample
EXPENSIVE = False

#: the tolerance (declared in the header before any number): margin >= this
MARGIN_MIN = -1e-12
#: (b): every helicity-state rate W_m >= this (relative to the unpolarized)
RATE_MIN = -1e-12
#: A2 central / V3
TR_Q2 = 2.0 / 3.0
TR_Q2_NF4 = 10.0 / 9.0
#: V4: SS90 Eq. (2.11)'s coefficient over DS16 Eq. (13)'s for the same
#: symbols g_xhat - g_yhat (a misreading, SECOND PASS S1: recorded only)
SS90_FACTOR = 2.0
#: Delta in DS16's normalisation (Eqs. (3)-(5)) per unit Delta of the kernel, at equal cos 2phi
#: amplitude -- a DERIVED ratio, re-derived numerically by the pytest; never set from a result.
#: 2 until D4 (SECOND PASS S2); 1 since D4, 2026-10-27 (V6's factor is PRE_D4_OVER_JM: the D4 note)
DS16_OVER_TREE = 1.0
#: V2: f1LL^g = +- F1LL_BAND * f1^g
F1LL_BAND = 1.0
#: the pipeline's toy scale, CHECKED against the kernel, never set by it
SHIPPED_SCALE = 1e-2
#: the quadrature's convergence requirement (fine vs coarse, relative)
QUAD_REL_TOL = 1e-7
#: the same-grid guard (plugin vs LhapdfSF / LhapdfG1), relative
SAME_GRID_TOL = 1e-12
GUARD_POINTS = ((1e-4, 2.0), (1e-2, 2.0), (0.1, 5.0), (0.5, 5.0), (1e-3, 10.0),
                (0.3, 10.0))
#: the letter's printed A (plans/07 addendum) is quoted to 4 decimals
A_PRINTED_HALF_ULP = 5e-5
#: the F2 / g1 flavour weights of src/lhapdf/lhapdf_sf.cpp (F2: d u s c b;
#: g1: d u s only)
E2_PROTON = {1: 1.0 / 9.0, 2: 4.0 / 9.0, 3: 1.0 / 9.0, 4: 4.0 / 9.0, 5: 1.0 / 9.0}
#: registry row 23's alternative effective polarization (beams.hpp
#: LI6_NAIVE_ONE_THIRD), for the recorded B with P_eff = 1/3
P_EFF_NAIVE = 1.0 / 3.0


class EnvironmentBlocked(Exception):
    """The environment cannot run this row (a missing module or plugin, a
    set's missing member file, or a set that is not the vendored version) --
    a BLOCKED row, exit 0.  A set that is not installed at all is not this:
    check_set raises LHAPDF's own message (exit 2, blocked (environment))."""


# ---------------------------------------------------------------- reference

def read_vendored(path=DATA):
    """The vendored document -- RuntimeError unless the file's bytes AND its
    canonical body hash to the recorded digests."""
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != FILE_SHA256:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, FILE_SHA256))
    doc = json.loads(raw.decode("ascii"))
    canon = json.dumps(doc["body"], sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False)
    body = hashlib.sha256(canon.encode("ascii")).hexdigest()
    if not body == doc["_body_sha256"] == BODY_SHA256:
        raise RuntimeError("%s: body sha256 %s; file records %s, harness %s"
                           % (path, body, doc["_body_sha256"], BODY_SHA256))
    return doc


def load_reference(path=DATA):
    """The body, with its shape checked (the keys this harness reads)."""
    body = read_vendored(path)["body"]
    keys = [s["key"] for s in body["sources"]]
    if keys != SOURCE_KEYS:
        raise RuntimeError("%s: sources %s, expected %s" % (path, keys,
                                                            " ".join(SOURCE_KEYS)))
    sets = [s["name"] for s in body["pdf_sets"]]
    if sets != ["CT18NLO", "NNPDFpol11_100"]:
        raise RuntimeError("%s: pdf_sets %s" % (path, sets))
    q2 = body["grid"]["q2_gev2"]
    if len(q2) != 3 or q2 != sorted(q2):
        raise RuntimeError("%s: the Q2 grid %s is not 3 increasing values"
                           % (path, q2))
    return body


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- LHAPDF sets

def lhapdf_dirs():
    """LHAPDF's own search order (LHAPDF::paths()): $LHAPDF_DATA_PATH -- or
    $LHAPATH when that is unset -- split on ':', then the install's datadir
    UNLESS the variable ends in '::'.  (Until phase F of run 2026-10-26 the
    datadir was appended always, so check_set could certify a copy LHAPDF
    never loads.)"""
    var = os.environ.get("LHAPDF_DATA_PATH")
    if var is None:
        var = os.environ.get("LHAPATH", "")
    dirs = [d for d in var.split(":") if d]
    if not var.endswith("::"):
        try:
            p = subprocess.run(["lhapdf-config", "--datadir"], capture_output=True,
                               text=True, timeout=10)
            if p.returncode == 0 and p.stdout.strip():
                dirs.append(p.stdout.strip())
        except (OSError, subprocess.SubprocessError):
            pass
    return dirs


def check_set(entry):
    """The set's .info and member file as LHAPDF would find them, each sha256
    against the vendored record -- EnvironmentBlocked (a BLOCKED row) if the
    member file is absent or either file differs; RuntimeError with LHAPDF's
    own message if the set is not installed (exit 2, blocked (environment))."""
    dirs = lhapdf_dirs()
    for d in dirs:
        info = os.path.join(d, entry["name"], entry["name"] + ".info")
        dat = os.path.join(d, entry["name"], entry["member_file"])
        if os.path.isfile(info):
            if not os.path.isfile(dat):
                raise EnvironmentBlocked("environment: LHAPDF set %s member file %s "
                                         "not installed" % (entry["name"], dat))
            got = (sha256_file(info), sha256_file(dat))
            want = (entry["info_sha256"], entry["member_file_sha256"])
            if got != want:
                raise EnvironmentBlocked(
                    "environment: LHAPDF set %s at %s is not the version this row "
                    "was measured with (sha256 info %s / member %s, vendored %s / %s)"
                    % (entry["name"], d, got[0][:12], got[1][:12], want[0][:12],
                       want[1][:12]))
            return dict(dir=d, info_sha256=got[0], member_sha256=got[1])
    # A set that is not installed is the environment, not a measurement: the
    # row's only configuration needs it, so LHAPDF's own message, uncaught --
    # main() exits 2 and make_report records 'blocked (environment)' (README
    # convention 6, as row 4; phase F of run 2026-10-26).
    raise RuntimeError("Info file not found for PDF set '%s' (searched %s)"
                       % (entry["name"], ", ".join(dirs) or "nothing"))


class PluginPdf:
    """One LHAPDF set read by LHAPDF itself through PYTHIA 8's LHAPDF6 plugin.
    `xf` is the RAW x f(x, Q2) (xfUpdate + xfRaw: no max(0, .) clip)."""

    def __init__(self, p8, setname):
        self.setname = setname
        self._pythia = p8.Pythia("", False)
        self._pythia.readString("Print:quiet = on")
        if not self._pythia.readString("PDF:pSet = LHAPDF6:%s" % setname):
            raise EnvironmentBlocked("environment: PYTHIA refused PDF:pSet = "
                                     "LHAPDF6:%s" % setname)
        pdf = self._pythia.getPDFPtr(2212)
        if pdf is None or not pdf.isSetup():
            raise EnvironmentBlocked("environment: PYTHIA's LHAPDF6 plugin could not "
                                     "load %s" % setname)
        self.pdf = pdf
        self.version = self._pythia.settings.parm("Pythia:versionNumber")

    def xf(self, pid, x, q2):
        self.pdf.xfUpdate(pid, x, q2)
        return self.pdf.xfRaw(pid)

    def alphas(self, q2):
        return self.pdf.alphaS(q2)


def load_sets(body):
    """{'ct18': PluginPdf, 'nnpdf': PluginPdf}, the set records, the PYTHIA
    version -- or EnvironmentBlocked (a BLOCKED row), or check_set's
    RuntimeError for a set that is not installed."""
    import lipolgen as lg  # noqa: WPS433 (env.sh puts it on the path)
    if not getattr(lg, "HAVE_LHAPDF", False):
        raise EnvironmentBlocked("environment: lipolgen built without the LHAPDF tier "
                                 "(HAVE_LHAPDF false): no LhapdfSF / LhapdfG1 / "
                                 "Epps21Ratio for the same-grid guard and V5")
    lg.lhapdf_quiet()      # one libLHAPDF in the process: silences the plugin's loads too
    records = [check_set(e) for e in body["pdf_sets"]]
    sys.path.insert(0, HERE)
    try:   # a missing _report.py is the harness broken (exit 2), not BLOCKED
        from _report import import_pythia8  # noqa: WPS433 (PYTHIA 8.317 cwd defect)
    finally:
        if sys.path and sys.path[0] == HERE:
            sys.path.pop(0)
    try:
        p8 = import_pythia8()
    except ImportError as exc:
        raise EnvironmentBlocked("environment: PYTHIA 8's python module (pythia8) "
                                 "is not importable: %s" % str(exc).splitlines()[0])
    ct = PluginPdf(p8, "CT18NLO")
    nn = PluginPdf(p8, "NNPDFpol11_100")
    return dict(ct18=ct, nnpdf=nn), records, ct.version


def same_grid_guard(sets):
    """The plugin reads the grids the tree's backends read: sum e_q^2 x(q+qbar)
    (5 flavours) = LhapdfSF.f2p and (1/2x) sum e_q^2 x(dq+dqbar) (d u s) =
    LhapdfG1.g1p.  RuntimeError (broken) beyond SAME_GRID_TOL."""
    import lipolgen as lg  # noqa: WPS433
    sf = lg.LhapdfSF("CT18NLO", 0)
    g1 = lg.LhapdfG1("NNPDFpol11_100", 0)
    ct, nn = sets["ct18"], sets["nnpdf"]
    worst = dict(f2p=0.0, g1p=0.0)
    for x, q2 in GUARD_POINTS:
        f2 = sum(e * (ct.xf(q, x, q2) + ct.xf(-q, x, q2)) for q, e in E2_PROTON.items())
        g = 0.5 * sum(e * (nn.xf(q, x, q2) + nn.xf(-q, x, q2))
                      for q, e in E2_PROTON.items() if q <= 3) / x
        for key, mine, tree in (("f2p", f2, sf.f2p(x, q2)), ("g1p", g, g1.g1p(x, q2))):
            rel = abs(mine - tree) / max(abs(tree), 1e-300)
            worst[key] = max(worst[key], rel)
    if max(worst.values()) > SAME_GRID_TOL:
        raise RuntimeError("same-grid guard: the PYTHIA plugin's grids differ from "
                           "LhapdfSF/LhapdfG1 (worst relative %s)" % worst)
    return dict(worst_rel=worst, points=len(GUARD_POINTS))


# ---------------------------------------------------------------- quadrature

def rule(x_lo, spec, refine=1):
    """(x, w) of the composite Gauss-Legendre rule of the vendored grid on
    [x_lo, 1]: ln-x panels on [x_lo, x_mid], linear panels on [x_mid, 1]."""
    import numpy as np
    t, w = np.polynomial.legendre.leggauss(int(spec["nodes_per_panel"]))
    xs, ws = [], []
    ue = np.linspace(math.log(x_lo), math.log(spec["x_mid"]),
                     refine * int(spec["log_panels"]) + 1)
    for a, b in zip(ue[:-1], ue[1:]):
        x = np.exp(0.5 * (a + b) + 0.5 * (b - a) * t)
        xs.append(x)
        ws.append(0.5 * (b - a) * w * x)
    xe = np.linspace(spec["x_mid"], 1.0, refine * int(spec["lin_panels"]) + 1)
    for a, b in zip(xe[:-1], xe[1:]):
        xs.append(0.5 * (a + b) + 0.5 * (b - a) * t)
        ws.append(0.5 * (b - a) * w)
    return np.concatenate(xs), np.concatenate(ws)


def integrate(f, x_lo, spec, label):
    """int_{x_lo}^1 f dx: (value on 2x the panels, |fine - coarse|, x_lo |f(x_lo)|).
    RuntimeError (broken) if the relative difference exceeds QUAD_REL_TOL."""
    import numpy as np
    vals = []
    for refine in (1, 2):
        x, w = rule(x_lo, spec, refine)
        fx = np.array([f(float(v)) for v in x], dtype=float)
        if not np.all(np.isfinite(fx)):
            raise RuntimeError("%s: a non-finite integrand value" % label)
        vals.append(float(np.dot(w, fx)))
    fine, err = vals[1], abs(vals[1] - vals[0])
    if err > QUAD_REL_TOL * max(abs(fine), 1e-300):
        raise RuntimeError("%s: quadrature not converged (%r vs %r)" % (
            label, vals[1], vals[0]))
    return fine, err, x_lo * abs(f(x_lo))


# ---------------------------------------------------------------- the bound

def bound_side(body, sets, ion):
    """Per Q2: alpha_s, the moments <x>_g (free and EPPS21), <x>_dg, P_eff,
    and B under the central assumptions and every recorded variation."""
    import lipolgen as lg  # noqa: WPS433
    spec = body["grid"]["x_integration"]
    ct, nn = sets["ct18"], sets["nnpdf"]
    nn_xmin = [s for s in body["pdf_sets"] if s["name"] == "NNPDFpol11_100"][0]["x_min"]
    p_eff = (ion.Z * ion.eff_pol_p + ion.N * ion.eff_pol_n) / float(ion.A)
    epps = lg.Epps21Ratio()
    out = []
    for q2 in body["grid"]["q2_gev2"]:
        a_s = ct.alphas(q2)
        xg, xg_err, xg_tail = integrate(lambda x: ct.xf(21, x, q2), spec["x_lo"], spec,
                                        "<x>_g(Q2=%g)" % q2)
        xg_a, xg_a_err, _ = integrate(lambda x: epps.ratio(21, x, q2) * ct.xf(21, x, q2),
                                      spec["x_lo"], spec, "<x>_g^EPPS21(Q2=%g)" % q2)
        xdg, xdg_err, xdg_tail = integrate(lambda x: nn.xf(21, x, q2), nn_xmin, spec,
                                           "<x>_dg(Q2=%g)" % q2)
        # CvDM17 (34) at the f1LL = -f1 edge needs |g1| <= f1/2 pointwise:
        # the first fine-rule node where it fails (and where |g1| > f1, i.e.
        # where the mixed central inputs leave positivity themselves), and
        # the share of <x>_g above it
        xs, ws = rule(nn_xmin, spec, 2)
        x_half = x_full = None
        for x in xs:
            g, dg = ct.xf(21, float(x), q2), p_eff * nn.xf(21, float(x), q2)
            if x_half is None and abs(dg) > 0.5 * g:
                x_half = float(x)
            if x_full is None and abs(dg) > g:
                x_full = float(x)
        share = (sum(float(w) * ct.xf(21, float(x), q2) for x, w in zip(xs, ws)
                     if x_half is not None and x >= x_half) / xg)
        pref = a_s / (8.0 * math.pi)
        g1m = p_eff * xdg

        def b(tr=TR_Q2, xg_=xg, ll=0.0, g1=g1m, norm=1.0):
            return pref * tr * norm * 0.5 * (xg_ * (1.0 + 0.5 * ll) + g1)

        xg_lo, xg_hi = min(xg, xg_a), max(xg, xg_a)
        out.append(dict(
            q2=q2, alpha_s=a_s, alpha_s_lo=lg.alpha_s_lo(q2), p_eff=p_eff,
            xg=xg, xg_err=xg_err, xg_tail=xg_tail, xg_epps21=xg_a, xg_epps21_err=xg_a_err,
            xdg=xdg, xdg_err=xdg_err, xdg_tail=xdg_tail, g1_moment=g1m,
            x_g1_above_half_f1=x_half, x_g1_above_f1=x_full,
            xg_share_above_x_half=share,
            B=b(), B_dg0=b(g1=0.0), B_ll_minus=b(ll=-F1LL_BAND), B_ll_plus=b(ll=F1LL_BAND),
            B_nf4=b(tr=TR_Q2_NF4), B_ss90=b(norm=SS90_FACTOR), B_epps21=b(xg_=xg_a),
            B_peff_naive=b(g1=P_EFF_NAIVE * xdg),
            # the edges (SECOND PASS S3): V6 acts on the model side (x2 on M2, PRE_D4_OVER_JM since D4)
            B_tight=b(xg_=xg_lo, ll=-F1LL_BAND, g1=0.0),
            B_loose=b(tr=TR_Q2_NF4, xg_=xg_hi, ll=F1LL_BAND, g1=max(g1m, 0.0))))
    return out


# ---------------------------------------------------------------- the models

def shipped_kernel():
    import lipolgen as lg  # noqa: WPS433
    return lg.default_inclusive_kernel(lg.li6())


def model_moments(body, kernel, bound, ct_alphas):
    """M2 = int x Delta dx per Q2 for every model, and the shipped scale check.
    `ct_alphas` is CT18NLO's alpha_s(Q2) (A3), handed to the letter's models."""
    import lipolgen as lg  # noqa: WPS433
    spec = body["grid"]["x_integration"]
    letter = body["letter_models"]

    def f1(x, q2):
        return kernel.tables(x, q2).f1

    la, lb = letter["moment_A"], letter["moment_B"]
    dil = la["dilution_num"] / float(la["dilution_den"])
    if lb["dilution_num"] / float(lb["dilution_den"]) != dil:
        raise RuntimeError("the vendored letter models disagree on the dilution")
    built = dict(
        moment_A_letter=lg.make_moment_a(f1, la["q2_ref_gev2"], la["variant"], lg.C_BAG,
                                         dil, ct_alphas),
        moment_B_letter=lg.make_moment_b(lb["variant"], lg.C_BAG, dil, ct_alphas),
        moment_A_dilution1=lg.make_moment_a(f1, la["q2_ref_gev2"], la["variant"], lg.C_BAG,
                                            1.0, ct_alphas),
        moment_B_dilution1=lg.make_moment_b(lb["variant"], lg.C_BAG, 1.0, ct_alphas),
        moment_A_letter_lo=lg.make_moment_a(f1, la["q2_ref_gev2"], la["variant"], lg.C_BAG,
                                            dil, None),
        moment_B_letter_lo=lg.make_moment_b(lb["variant"], lg.C_BAG, dil, None))
    a_solved = lg.solve_A_interp_a(f1, la["q2_ref_gev2"], la["variant"], lg.C_BAG)
    closure = 0.0          # moment_B's moment is dilution C_BAG alpha_s exactly
    out = {k: [] for k in ["toy_shipped"] + list(built)}
    scale_dev = 0.0
    scale_ref = None
    for b in bound:
        q2 = b["q2"]

        def toy(x, q2=q2):
            t = kernel.tables(x, q2)
            return x * t.delta

        m, e, tail = integrate(toy, spec["x_lo"], spec, "M2[toy](Q2=%g)" % q2)
        out["toy_shipped"].append(dict(q2=q2, m2=m, err=e, tail=tail))
        for k, mod in built.items():
            def g(x, q2=q2, mod=mod):
                return x * mod(x, q2, kernel.tables(x, q2).f1)
            m, e, tail = integrate(g, spec["x_lo"], spec, "M2[%s](Q2=%g)" % (k, q2))
            out[k].append(dict(q2=q2, m2=m, err=e, tail=tail))
        for k, d in (("moment_B_letter", dil), ("moment_B_dilution1", 1.0)):
            exact = d * lg.C_BAG * ct_alphas(q2)
            closure = max(closure, abs(out[k][-1]["m2"] / exact - 1.0))
        # the shipped scale, measured at every node of the fine rule
        xs, _ = rule(spec["x_lo"], spec, 2)
        for x in xs:
            t = kernel.tables(float(x), q2)
            den = t.f1 * float(x) ** 0.3 * (1.0 - float(x)) ** 4
            if den > 0.0:
                s = t.delta / den
                scale_dev = max(scale_dev, abs(s / SHIPPED_SCALE - 1.0))
        t = kernel.tables(0.1, q2)
        scale_ref = t.delta / (t.f1 * 0.1 ** 0.3 * 0.9 ** 4)
    if closure > 1e-9:
        raise RuntimeError("moment_B's integrated M2 misses dilution C_BAG alpha_s by "
                           "%.2e relative: the quadrature is broken" % closure)
    return out, dict(scale_at_x0p1=scale_ref, max_rel_dev_from_shipped_scale=scale_dev,
                     moment_b_closure_rel=closure,
                     a_solved=a_solved, a_printed=la["A_as_printed"],
                     a_reproduced=abs(a_solved - la["A_as_printed"]) <= A_PRINTED_HALF_ULP)


# ---------------------------------------------------------------- (b) positivity

def helicity_positivity(body, kernel):
    """min over the vendored grid of the helicity-state rates W_m (m = +1, 0,
    -1; relative to the unpolarized rate) and the range of A_T."""
    import lipolgen as lg  # noqa: WPS433
    g = body["grid"]["positivity"]
    w_min, a_lo, a_hi, at_min = float("inf"), float("inf"), -float("inf"), None
    n = 0
    for x in g["x"]:
        for q2 in g["q2_gev2"]:
            t = kernel.tables(x, q2)
            for y in g["y"]:
                s = q2 / (x * y)
                for th in g["theta_s_deg"]:
                    amps = {m: kernel.amplitudes(
                        t, x, q2, s, lg.EventSpinState(lam_e=0, pe=0.0, j=1.0, m=float(m),
                                                       theta_s=math.radians(th),
                                                       phi_s=0.0))
                            for m in (1, 0, -1)}
                    for ph in g["phi_prime_deg"]:
                        w = {m: lg.InclusiveKernel.density(a, math.radians(ph))
                             for m, a in amps.items()}
                        a_t = (w[1] + w[-1] - 2.0 * w[0]) / (w[1] + w[0] + w[-1])
                        n += 1
                        a_lo, a_hi = min(a_lo, a_t), max(a_hi, a_t)
                        for m in (1, 0, -1):
                            if w[m] < w_min:
                                w_min = w[m]
                                at_min = dict(x=x, q2=q2, y=y, theta_s_deg=th,
                                              phi_prime_deg=ph, m=m)
    ok = w_min >= RATE_MIN
    return dict(name="helicity_positivity[shipped kernel]",
                status="pass" if ok else "fail", in_verdict=True,
                source="Cosyn, Roldan Tomei, Sosa, Zec, EPJ A 61 (2025) 83, Eq. (2): "
                       "-2 <= A_T <= 1",
                n_points=n, n_rates=3 * n, w_min=w_min, a_t_min=a_lo, a_t_max=a_hi,
                w_min_at="x=%g Q2=%g y=%g theta_S=%g phi'=%g m=%+d" % (
                    at_min["x"], at_min["q2"], at_min["y"], at_min["theta_s_deg"],
                    at_min["phi_prime_deg"], at_min["m"]))


# ---------------------------------------------------------------- verdict

def f1ll_critical(m2, b, with_g1=True, factor=1.0):
    """The constant fraction f1LL^g/f1^g at which factor x |M2| = B (central
    Tr Q^2, alpha_s, f1^g; g1^g central or 0; factor = PRE_D4_OVER_JM for V6):
    the scenario is allowed iff f1LL^g >= this fraction x f1^g.  Below
    -2 (1 - |g1|/f1) it is outside positivity; above +1 (CvDM17 (35)) no
    positive tensor polarization can reach it."""
    k = b["alpha_s"] / (8.0 * math.pi) * TR_Q2 * 0.5
    g1 = b["g1_moment"] if with_g1 else 0.0
    return 2.0 * ((factor * abs(m2) / k - g1) / b["xg"] - 1.0)


def rule_a(sat):
    """Rule (a) of the tolerance on a list of saturations."""
    return "pass" if all(1.0 - s >= MARGIN_MIN for s in sat) else "fail"


def saturation_row(name, rows, bound, in_verdict, extra=None):
    """One model's sub-row: M2, B, S, margin per Q2 (central), V6, the TIGHT
    and LOOSE edges, every variation; the outcome of rule (a) under the
    central assumptions (`status` for the verdict model, `status_central`
    beside a recorded one) and under V6 (`status_v6`; x PRE_D4_OVER_JM since D4)."""
    m2 = [r["m2"] for r in rows]
    sat = [abs(m) / b["B"] for m, b in zip(m2, bound)]
    sat_v6 = [PRE_D4_OVER_JM * s for s in sat]
    central, v6 = rule_a(sat), rule_a(sat_v6)
    d = dict(name=name, status=central if in_verdict else "recorded",
             in_verdict=in_verdict, status_central=central, status_v6=v6,
             q2=[b["q2"] for b in bound], m2=m2, bound=[b["B"] for b in bound],
             saturation=sat, margin=[1.0 - s for s in sat], saturation_v6=sat_v6,
             saturation_tight=[PRE_D4_OVER_JM * abs(m) / b["B_tight"]
                               for m, b in zip(m2, bound)],
             saturation_loose=[abs(m) / b["B_loose"] for m, b in zip(m2, bound)],
             f1ll_critical=[f1ll_critical(m, b) for m, b in zip(m2, bound)],
             f1ll_critical_dg0=[f1ll_critical(m, b, False) for m, b in zip(m2, bound)],
             f1ll_critical_v6=[f1ll_critical(m, b, True, PRE_D4_OVER_JM)
                               for m, b in zip(m2, bound)],
             m2_numerical_error=[r["err"] for r in rows],
             m2_tail_x_lo=[r["tail"] for r in rows])
    for key in ("B_dg0", "B_ll_minus", "B_ll_plus", "B_nf4", "B_ss90", "B_epps21",
                "B_peff_naive"):
        d["saturation_" + key[2:]] = [abs(m) / b[key] for m, b in zip(m2, bound)]
    if extra:
        d.update(extra)
    return d


def measure(path=DATA):
    import lipolgen as lg  # noqa: WPS433
    body = load_reference(path)
    sets, set_records, pythia_version = load_sets(body)
    guard = same_grid_guard(sets)
    ion = lg.li6()
    kernel = shipped_kernel()
    bound = bound_side(body, sets, ion)
    moments, checks = model_moments(body, kernel, bound, sets["ct18"].alphas)
    letter = body["letter_models"]["moment_A"]
    subrows = [
        saturation_row("toy_1e-2[shipped]", moments["toy_shipped"], bound, True,
                       dict(scale_measured=checks["scale_at_x0p1"])),
        helicity_positivity(body, kernel),
        saturation_row("moment_A[letter]", moments["moment_A_letter"], bound, False, dict(
            variant=letter["variant"], q2_ref=letter["q2_ref_gev2"],
            dilution="%d/%d" % (letter["dilution_num"], letter["dilution_den"]),
            a_solved=checks["a_solved"], a_printed=checks["a_printed"],
            a_reproduced=checks["a_reproduced"],
            saturation_model_lo_alpha_s=[abs(r["m2"]) / b["B"] for r, b in
                                         zip(moments["moment_A_letter_lo"], bound)])),
        saturation_row("moment_B[letter]", moments["moment_B_letter"], bound, False, dict(
            variant=body["letter_models"]["moment_B"]["variant"], dilution="1/3",
            saturation_model_lo_alpha_s=[abs(r["m2"]) / b["B"] for r, b in
                                         zip(moments["moment_B_letter_lo"], bound)])),
        saturation_row("moment_A[dilution 1]", moments["moment_A_dilution1"], bound, False),
        saturation_row("moment_B[dilution 1]", moments["moment_B_dilution1"], bound, False),
    ]
    keys = ("alpha_s", "xg", "xg_epps21", "xdg", "g1_moment", "B", "B_dg0", "B_ll_minus",
            "B_ll_plus", "B_nf4", "B_ss90", "B_epps21", "B_peff_naive", "B_tight",
            "B_loose", "x_g1_above_half_f1", "x_g1_above_f1", "xg_share_above_x_half")
    side = dict(name="bound_side[central]", status="recorded", in_verdict=False,
                q2=[b["q2"] for b in bound], p_eff=bound[0]["p_eff"], tr_q2=TR_Q2,
                ds16_over_tree=DS16_OVER_TREE, pre_d4_over_jm=PRE_D4_OVER_JM)
    side.update({k: [b[k] for b in bound] for k in keys})
    subrows.append(side)
    verdict = [s for s in subrows if s.get("in_verdict")]
    status = "pass" if all(s["status"] == "pass" for s in verdict) else "fail"
    return dict(body=body, bound=bound, subrows=subrows, status=status,
                checks=dict(checks, same_grid=guard), set_records=set_records,
                pythia_version=pythia_version)


def _fmt3(v, f="%.4g"):
    return " / ".join(f % x for x in v)


def run(path=DATA, verbose=True):
    t0 = time.perf_counter()
    import lipolgen as lg  # noqa: WPS433 (missing lipolgen = broken, exit 2)
    load_reference(path)   # an altered file is broken (exit 2), never BLOCKED
    tolerance = ("a bound under assumptions: PASS iff 1 - |M2|/B >= %g at each Q2 = "
                 "2, 5, 10 GeV2 for the shipped toy AND every helicity-state rate >= "
                 "%g (-2 <= A_T <= 1) on the vendored grid; the letter's moment_A/B "
                 "carry their own status, recorded" % (MARGIN_MIN, RATE_MIN))
    reference = ("gluonic Soffer-type bound, n = 2 moment: Detmold-Shanahan PRD 94 "
                 "(2016) 014507 Eqs. (13), (32) [(36)]; f1^g CT18NLO, g1^g = P_eff x "
                 "NNPDFpol11 Delta g, f1LL^g = 0, Tr Q^2 = 2/3, alpha_s CT18NLO; "
                 "helicity positivity -2 <= A_T <= 1 (Cosyn et al. EPJ A 61 (2025) 83 "
                 "Eq. (2))")
    config = dict(lipolgen=getattr(lg, "__version__", None),
                  python=sys.version.split()[0], data_file_sha256=FILE_SHA256,
                  tr_q2=TR_Q2, f1ll_band=F1LL_BAND, ds16_over_tree=DS16_OVER_TREE,
                  pre_d4_over_jm=PRE_D4_OVER_JM, quad_rel_tol=QUAD_REL_TOL)
    try:
        m = measure(path)
    except EnvironmentBlocked as exc:
        report = dict(name=NAME, status="blocked", reason=str(exc),
                      generator="not run (%s)" % exc, reference=reference,
                      tolerance=tolerance, config=config, subrows=[])
        report["runtime_s"] = time.perf_counter() - t0
        if verbose:
            print("REPORT | %s | %s | %s | %s | BLOCKED" % (
                NAME, report["generator"], reference, tolerance))
        return report
    toy = m["subrows"][0]
    pos = m["subrows"][1]
    generator = ("shipped toy_delta_gluon (scale %.3g) per nucleon, 6Li: M2 = int x "
                 "Delta dx = %s at Q2 = 2 / 5 / 10 GeV2, |M2|/B = %s (Jaffe-Manohar Delta, the "
                 "kernel's since D4; V6, %s: %s); helicity-state rates min W = %.6f (A_T in [%.3g, "
                 "%.3g])" % (
                     toy["scale_measured"], _fmt3(toy["m2"], "%.4e"),
                     _fmt3(toy["saturation"], "%.4f"), V6_LABEL, _fmt3(toy["saturation_v6"], "%.4f"),
                     pos["w_min"], pos["a_t_min"], pos["a_t_max"]))
    config.update(pythia=m["pythia_version"], lhapdf_sets=[
        "%s member 0 (sha256 %s)" % (e["name"], e["member_file_sha256"][:12])
        for e in m["body"]["pdf_sets"]])
    report = dict(name=NAME, status=m["status"], generator=generator,
                  reference=reference, tolerance=tolerance, config=config,
                  subrows=m["subrows"], bound_detail=m["bound"], checks=m["checks"])
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        _print(report, m)
    return report


def _print(report, m):
    b = m["bound"]
    print("bound side (per nucleon; Tr Q^2 = %.6g; alpha_s, <x>_g from CT18NLO m0; "
          "<x>_dg from NNPDFpol11_100 m0; P_eff = %.10f):" % (TR_Q2, b[0]["p_eff"]))
    for r in b:
        print("  Q2 = %4g: alpha_s %.6f (LO fallback %.6f); <x>_g %.8f (+- %.1e; "
              "EPPS21 6Li %.8f); <x>_dg %.8f (+- %.1e) -> P_eff <x>_dg %.8f" % (
                  r["q2"], r["alpha_s"], r["alpha_s_lo"], r["xg"], r["xg_err"],
                  r["xg_epps21"], r["xdg"], r["xdg_err"], r["g1_moment"]))
        print("            B = %.6e  [dg0 %.6e, f1LL -/+ %.6e / %.6e, nf4 %.6e, SS90 "
              "(V4, a misreading) %.6e, EPPS21 %.6e, P_eff = 1/3 %.6e; tight %.6e (V6 on "
              "M2: %s), loose %.6e]" % (
                  r["B"], r["B_dg0"], r["B_ll_minus"], r["B_ll_plus"], r["B_nf4"],
                  r["B_ss90"], r["B_epps21"], r["B_peff_naive"], r["B_tight"],
                  V6_LABEL, r["B_loose"]))
        print("            |g1^g| > f1^g/2 from x = %s, > f1^g from x = %s (share of "
              "<x>_g above: %.2e)" % (
                  "%.4f" % r["x_g1_above_half_f1"] if r["x_g1_above_half_f1"] else "never",
                  "%.4f" % r["x_g1_above_f1"] if r["x_g1_above_f1"] else "never",
                  r["xg_share_above_x_half"]))
    c = m["checks"]
    print("checks: shipped scale %.12g (max rel dev from 1e-2 over the nodes %.1e); "
          "moment_A A = %.6f vs the letter's printed %.4f (reproduced: %s); moment_B "
          "closure %.1e; same-grid guard worst %s" % (
              c["scale_at_x0p1"], c["max_rel_dev_from_shipped_scale"], c["a_solved"],
              c["a_printed"], c["a_reproduced"], c["moment_b_closure_rel"],
              ", ".join("%s %.1e" % kv for kv in sorted(c["same_grid"]["worst_rel"].items()))))
    for s in report["subrows"]:
        if "saturation" in s:
            word = (s["status"].upper() if s["in_verdict"] else
                    "RECORDED: %s centrally (A7, the kernel's own Delta since D4), %s under V6 (%s)" % (
                        s["status_central"].upper(), s["status_v6"].upper(), V6_LABEL))
            if s["in_verdict"]:
                word += " (%s under V6, %s)" % (s["status_v6"].upper(), V6_LABEL)
            print("SUBROW | %s[%s] | M2 = %s | B = %s | 1 - |M2|/B >= %g | %s "
                  "(|M2|/B = %s; V6 %s; tight %s; loose %s; saturates at f1LL^g/f1^g = "
                  "%s, V6 %s)" % (
                      NAME, s["name"], _fmt3(s["m2"], "%.4e"), _fmt3(s["bound"], "%.4e"),
                      MARGIN_MIN, word, _fmt3(s["saturation"]), _fmt3(s["saturation_v6"]),
                      _fmt3(s["saturation_tight"]), _fmt3(s["saturation_loose"]),
                      _fmt3(s["f1ll_critical"], "%.3f"),
                      _fmt3(s["f1ll_critical_v6"], "%.3f")))
        elif "w_min" in s:
            print("SUBROW | %s[%s] | min W_m = %.6f over %d rates (A_T in [%.4g, %.4g]; "
                  "min at %s) | -2 <= A_T <= 1 (%s) | %s" % (
                      NAME, s["name"], s["w_min"], s["n_rates"], s["a_t_min"],
                      s["a_t_max"], s["w_min_at"], s["source"], s["status"].upper()))
        else:
            print("SUBROW | %s[%s] | B = %s at Q2 = 2 / 5 / 10 (alpha_s %s; <x>_g %s; "
                  "P_eff <x>_dg %s) | not a verdict | RECORDED" % (
                      NAME, s["name"], _fmt3(s["B"], "%.6e"), _fmt3(s["alpha_s"], "%.6f"),
                      _fmt3(s["xg"], "%.8f"), _fmt3(s["g1_moment"], "%.8f")))
    print("REPORT | %s | %s | %s | %s | %s" % (
        NAME, report["generator"], report["reference"], report["tolerance"],
        report["status"].upper()))


# ------------------------------------------------------------------ figure

def plot(report, out_dir=None):
    """|M2|/B per Q2 and model as grouped markers: filled = the central
    assumptions (A7, the kernel's own Delta since D4), open = V6 (the pre-D4 kernel / the letter as
    printed, x2); the whisker spans the LOOSE to the TIGHT edge (SECOND PASS S3, V6 kept); the bound at 1."""
    sys.path.insert(0, HERE)
    from _plots import INK, MUTED, new_figure, save, status_title, style_series
    fig, ax = new_figure(height=3.9)
    if report.get("status") == "blocked":
        ax.text(0.5, 0.5, "BLOCKED: %s" % report.get("reason"), ha="center",
                va="center", transform=ax.transAxes, fontsize=8, color=MUTED, wrap=True)
        status_title(fig, NAME, report["status"], "not run")
        return save(fig, NAME, out_dir)
    subs = {s["name"]: s for s in report["subrows"] if "saturation" in s}
    series = (("toy_1e-2[shipped]", "default", "o", "toy_delta_gluon @ 1e-2 (shipped; verdict)"),
              ("moment_A[letter]", "alt", "s", "moment_A, letter (mid_x, dilution 1/3)"),
              ("moment_B[letter]", "alt2", "D", "moment_B, letter (mid_x, dilution 1/3)"),
              ("moment_A[dilution 1]", "alt3", "^", "moment_A at dilution 1 (recorded)"))
    q2 = subs["toy_1e-2[shipped]"]["q2"]
    for i, (key, slot, mk, lab) in enumerate(series):
        s = subs[key]
        xs = [j + (i - 1.5) * 0.14 for j in range(len(q2))]
        kw = style_series(slot, label=lab, marker=mk)
        ax.vlines(xs, s["saturation_loose"], s["saturation_tight"], color=kw["color"],
                  linewidth=1.1, zorder=2)
        ax.plot(xs, s["saturation"], **kw)
        open_kw = dict(kw, label=None, markerfacecolor="none", linestyle="none")
        ax.plot(xs, s["saturation_v6"], **open_kw)
    ax.plot([], [], linestyle="none", marker="o", markerfacecolor="none", color=MUTED,
            label="open: V6, %s" % V6_LABEL)
    ax.axhline(1.0, color=INK, linewidth=1.4, zorder=4, label="the bound (|M2| = B)")
    ax.set_yscale("log")
    ax.set_ylim(1.5e-3, 60.0)
    ax.set_xticks(range(len(q2)))
    ax.set_xticklabels(["%g" % v for v in q2])
    ax.set_xlabel("Q2 [GeV2]")
    ax.set_ylabel("saturation |M2| / B  (n = 2 moment)")
    ax.set_xlim(-0.6, len(q2) - 0.4)
    ax.legend(loc="lower right", fontsize=6.2, ncol=2)
    ax.text(0.01, 0.985, "filled: central (f1^g CT18NLO; g1^g = P_eff NNPDFpol11 Delta g; "
            "f1LL^g = 0; Tr Q^2 = 2/3; Jaffe-Manohar Delta, D4)\nwhiskers: loose edge (f1LL^g = +f1^g, "
            "Tr Q^2 = 10/9) to tight edge (Delta g = 0, f1LL^g = -f1^g, V6: pre-D4 x2)",
            transform=ax.transAxes, fontsize=6.3, color=MUTED, va="top")
    status_title(fig, NAME, report["status"],
                 "a bound under assumptions: PASS iff |M2|/B <= 1 at each Q2 (shipped "
                 "toy) and all helicity-state rates >= 0")
    return save(fig, NAME, out_dir)


# ------------------------------------------------------------------ D4
#
# D4 (2026-10-27, after both passes and after the author's decision D4,
# docs/open_items/run_2026-10-27/DECISIONS.md; plan_D4_jm_convention.md and
# its adversarial check).  Kept here, after every line the run's records cite,
# so that none of them moves.  The kernel's cos 2phi coefficient is Jaffe-
# Manohar's since this date, a2 = -(1/2) ((1-y)/y^2) c_eff sin^2(theta_S)
# Delta/D_phi (include/lipolgen/constants.hpp, DELTA_COS2PHI_NORM = 1/2), which
# is what DS16 Eqs. (3)-(5) contract to.  So the first-pass identification A7
# is now the CONSISTENT one: DS16_OVER_TREE, re-derived by the pytest from
# Eqs. (3)-(5) against the kernel, is 1, and the central column is the shipped
# toy's GENERATED amplitude in DS16's units.  V6 is kept as declared,
# numerically unchanged, and only RELABELLED: it is the reading of the PRE-D4
# kernel and of the letter's printed amplitudes (the LOI's Eq. (9), twice
# Jaffe-Manohar's coefficient: PRE_D4_OVER_JM = 2).  The TIGHT edge keeps V6
# as S3 declared it: whether it should lose V6 now is the author's question
# (plan_D4 sec. 12 item 3), and no declared edge is dropped before the answer.
# Nothing else moves: not M2 (the toy keeps its 1e-2 scale and enters the
# kernel unconverted), not B, not the rule or its tolerance, not the verdict;
# sub-row (b)'s rates carry the halved Delta term (min W 0.9493997515 ->
# 0.9493997990 at the same point; A_T max 0.0506002485 -> 0.0506002010; A_T
# min, at theta_S = 0, unmoved).  The vendored file's LOI18 note ("(9) is the
# tree's kernel term") is the second pass's reading of the pre-D4 kernel,
# dated 2026-10-26, and is not re-hashed.
#
#: V6, relabelled by D4: the PRE-D4 kernel's coefficient (the LOI's Eq. (9),
#: still the letter's printed Eqs. (1) and (3)) over Jaffe-Manohar's -- the
#: reading V6 records, numerically unchanged, and the factor the TIGHT edge
#: carries as S3 declared it (module-level, read by the functions above at
#: call time)
PRE_D4_OVER_JM = 2.0
#: how the reports name V6 since D4
V6_LABEL = "the pre-D4 kernel / the letter as printed, x2"


def main(argv=None):
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    argv = sys.argv[1:] if argv is None else argv
    try:
        rep = run()
        if "--plot" in argv:
            print("figure: %s" % plot(rep))
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
