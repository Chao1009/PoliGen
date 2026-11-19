<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Phase B1 — row 10, CLAS BONuS: R_D/S against cos θ_pq in the DIS blocks, held against the tree's Glauber FSI weight

Written 2026-10-26 by the B1 worker against base commit `87787dd` (HEAD
moved to `372b915`, the supervising session's B3 commit, during the work; the
shared tree dirty with other workers' files), after the runs, from output
quoted with the command that produced it. **Nothing in the tree was tuned**:
no physics default, constant, registry status, tolerance or reference JSON
moved; nothing under `include/`, `src/` or `python/lipolgen/` was edited.
The plan item is `PLAN.md` phase B1; the ranking it serves is
`../../benchmarking/08_channel_priority.md` §4 ingredient **I** (spectator
FSI, "the Deeps transverse window and BONuS (DIS blocks only) as direction
checks") under ranks 6 (tagged ⁶Li α) and 7 (the tagged d–p control).

Environment for every command below:

```bash
cd /home/cpeng/Projects/polli/PoliGen && source env.sh
```

## 0. Summary

| | |
|---|---|
| row | `BENCHMARK_PLAN.md` §4 row 10, `validation/benchmarks/t2_bonus_spectator_shape.py` — **BLOCKED → FAIL, recorded** |
| reference | the Supplemental Material of Tkachenko *et al.* (CLAS BONuS), PRC 89 (2014) 045206 — 100 blocks of R_D/S (data over BONuS's own **PWIA, no-FSI** spectator simulation) against cos θ_pq, vendored as numbers in `validation/benchmarks/data/bonus2014_tkachenko_supplemental.json` |
| window (declared first) | W* ≥ 2 GeV and Q² ≥ 1 GeV² on the block labels: 32 blocks, **28 verdict blocks** (12 at 4 GeV, 16 at 5 GeV) with measured backward bins; 68 out of model (20 at Q² = 0.93, 48 resonance region) |
| tolerance (declared before any number) | per block, the tree's `GlauberFsiWeight` over its backward-window mean must stay inside the data's 1σ flat band (stat, PDG scale factor) at every cos θ_pq ≤ −0.3 bin; every block, both beams; default = Hulthén, σ_XN = 40 mb |
| measured | inside the band in **4/12 blocks (4 GeV), 0/16 (5 GeV)**; worst v/δ = 2.62 / 3.44. The weight's variation across the backward window v = **0.034 / 0.055 / 0.085 / 0.141** at p_s = 0.078 / 0.093 / 0.11 / 0.135 GeV/c |
| what the FAIL means | the rule's premise (a flat backward window) fails, not the weight's direction: the backward data are **not flat** at stat (χ²_flat **739.4/84**) and **fall toward 90° as the weight does** — χ² against the weight's shape (one free scale per block) **226.3/84** at 40 mb, 221.6 for AV18, 403.3 at 20 mb (all recorded, all declared before the numbers) |
| recorded beside it | stat ⊕ syst band: 8/12, 8/16; Hulthén 20 mb: 8/12, 13/16; AV18 40 mb: 3/12, 0/16; the forward rise the symmetric weight cannot carry; the instant-form / light-cone PWIA factor (same direction as the weight in the backward window, 40–66 % of its size); ToyF2's F₂ⁿ/F₂ᵖ against BONuS's DIS table: 40.5/10 absolute, 8.30/9 with one free scale (0.902) |
| tests | `python/tests/test_bench_bonus.py`, **26 passed** (13.0 s); the harness itself runs in 2.3–2.4 s, not behind `LIPOLGEN_BENCH` |
| figure | `docs/benchmarking/figures/t2_bonus_spectator_shape.png` |

## 1. Files

New: `validation/benchmarks/data/bonus2014_tkachenko_supplemental.json`
(59 194 bytes, sha256 `3636ff08adbd1d7b526b30076e03512164e036c31392e291fbd0fcf6f88d3201`),
`python/tests/test_bench_bonus.py`,
`docs/benchmarking/figures/t2_bonus_spectator_shape.png` (written by
`--plot` runs of the harness: one of an intermediate version, then one of the
final code, which is the file on disk), and this record. Rewritten:
`validation/benchmarks/t2_bonus_spectator_shape.py` (the 2026-09-23 header
kept as a dated HISTORY block, its "6 W* × 5 Q² × 4 p_s" corrected in
place with a note). No other repository file was edited. Kept outside the
repository, in the session scratch directory: the original
`BONuSprcSupMat2.txt` (`…/scratchpad/refs/`, as fetched by the supervising
session), and this worker's dev-time scripts (`…/scratchpad/b1w/scripts/`:
the vendoring parse, an independent cross-check, a grid-convergence probe).

## 2. The document and what was vendored

The supervising session fetched `BONuSprcSupMat2.txt` on 2026-10-26 from
<https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.89.045206/BONuSprcSupMat2.txt>
with curl's own identity through the institution's APS licence (nothing
spoofed; `phase_0_environment.md`): 88 845 bytes, 5452 LF-terminated lines,
ASCII, tab-separated, sha256
`18c7eb384c59c854e9e7926a1f96ec6840a9cc05674142bdfefee837b3f65f0e`
(re-measured here). Its sections, by the line their header starts on:

| line | content | vendored |
|---|---|---|
| 7 | R_D/S vs cos θ_pq, **4 GeV** beam (4.23 GeV in the paper's Fig. 12): 3 Q² (0.93, 1.66, 3.38) × 5 W* (1.17, 1.48, 1.73, 2.03, 2.44) × 4 p_s (0.078, 0.093, 0.11, 0.135) = **60 blocks** × 10 rows | yes |
| 1695 | the same, **5 GeV** beam (5.25 GeV, Figs. 10–11): 2 Q² (1.66, 3.38) × 5 W* × 4 p_s = **40 blocks** × 10 rows | yes |
| 2823, 3009 | F₂ⁿ_eff(p_s)/F₂ⁿ_eff(lowest p_s) vs x* at Q² = 1.66 and 3.38 (Figs. 15–16) — two sections the task brief did not list | no (listed in the file) |
| 3173, 4047 | F₂ⁿ/F₂ᵈ per (W*², Q², x*) bin, 4 and 5 GeV (Fig. 17) | no (listed) |
| 5363 | F₂ⁿ/F₂ᵈ and F₂ⁿ/F₂ᵖ in the DIS region for W* > 1.4 / 1.6 / 1.8 GeV (13 / 12 / 10 rows; the header says "see Fig. 20", the three-cut F₂ⁿ/F₂ᵖ is the paper's Fig. 21) | yes |

**The binning, corrected.** The old harness header (and `REFERENCES.md`
item 59) said "6 W* × 5 Q² × 4 p_s": that is the paper's text (Sec. V A 1),
which counts the quasi-elastic W* region and the Q² = 0.34, 0.61 bins the
supplement does not hold. The block labels are bin centres: W* = 2.03 and
2.44 are the paper's regions 1.85–2.20 and 2.2–2.68 GeV (the mid-points of
all five check: (1.0 + 1.35)/2 = 1.175 … (2.2 + 2.68)/2 = 2.44); Q² = 1.66
is 1.10–2.23 GeV² (Figs. 10–13 captions); p_s = 70–85, 85–100, 100–120,
120–150 MeV/c; cos θ_pq ten even bins over [−1, 1].

**Missing rows.** 22 rows read "0 0 0 0" — all in the four 4 GeV blocks at
(Q² = 3.38, W* = 2.44), whose backward bins the 4 GeV beam cannot reach;
kept verbatim, read as MISSING, never as measured zeros. One measured value
is negative (−0.447 ± 0.974 at cos θ_pq = 0.3 in the 4 GeV (3.38, 2.44,
0.135) block); it has a stat error and is kept (that block carries no
verdict).

**The parse.** A dev-time script read the original and wrote the JSON (every
number the printed decimal, written as the shortest round-trip double;
the only tokens whose text differs are the integers "0" ×88, "1" ×3 — value
identical); the JSON carries the `src_line` of every section, block and
table, so the parse is re-checkable against the original. An independent
regex pass over the raw text, in file order, against the JSON flattened in
the same order: **2580 / 1720 / 318** numbers in sections 1 / 2 / 5, all
equal. The provenance block names the DOI, URL, fetch date, "facts quoted
with attribution; APS copyright on the document; the licence question …
listed for the author", the original's size and sha256, the paper and its
corpus file's sha256, the convention of every table and what was not
vendored.

## 3. What the denominator is — read before any design

arXiv:1402.2477v2 (the author's corpus copy, `pdftotext` of it), Sec. IV C:
BONuS's simulation is a **plane-wave spectator model** — the Paris deuteron
wave function rescaled with the **light-cone** formalism of Frankfurt and
Strikman; the struck neutron off shell, E = M_d − √(M² + p_s²), p = −p_s;
electrons generated in the struck nucleon's rest frame and corrected for its
off-shell mass; Bosted–Christy F₂ⁿ; Mo–Tsai radiative effects (RCSLACPOL);
GEANT4 RTPC and GSIM CLAS. In the paper's words: "Our generator did not
contain additional processes like coherent pion production, final-state
interactions and other two-nucleon effects; therefore, the ratio of measured
to simulated tagged data can be interpreted as a direct test of the
spectator picture." **No FSI in the denominator.** The normalisation
N(p_s, E_b) is ONE factor per p_s bin and beam energy, fitted on
−1 < cos θ_pq < −0.2 at W* = 2.0–2.2 GeV in the lowest fully accepted Q²
bin and common to every cos θ_pq and (W*, Q²) bin; its ±10 % is not in the
syst column (Sec. IV C 3–4). Table I's point-to-point systematic is 12.5 %
(the supplement's syst column reads about 11–15 % of each value), "a relatively
smooth function of the kinematic variables"; the paper draws R_D/S with
statistical error bars and the systematic as a band.

**So R_D/S can** bound a cos θ_pq dependence at fixed (W*, Q², p_s) beyond
the light-cone Paris PWIA (x* is fixed in a block: W*² = p² + Q²(1 − x*)/x*
with p² depending on |p_s| only). **It cannot** give the p_s spectrum (the
normalisation is per p_s bin — the row's 2026-09-23 goal is not in the
document), an angle-independent FSI suppression (absorbed by the same
factor), or a separation of FSI from the PWIA prescription (light-cone in
the denominator; the tree's density is instant-form, §5.4).

## 4. The design, declared before any number

The harness header block "THE RULE AND ITS TOLERANCE" was written before any
fit, band, χ² or tree evaluation of this row was computed; the window is
defined on the block labels. What had been seen before it, and is said in the
header: the supervising session previewed the tables (`PLAN.md` ground rules,
rule 4); this worker read every value while transcribing them, and the block
census (counts, the 22 missing rows).

- **(a) DIS window** W* ≥ 2 GeV and Q² ≥ 1 GeV²; why: `GlauberFsiWeight`
  models DIS debris rescattering at σ_XN = 40 mb, and the generator never
  produces a resonance event — its own scenario floor is W² ≥ 8 GeV²
  (`include/lipolgen/sampler.hpp:130`, `generator_scenario`, "Q2 >= 0.7, y in
  [0.004, 0.985], W2 >= 8"). **No BONuS block reaches that floor** (the top
  W* bin ends at 2.68 GeV, W*² ≤ 7.2 GeV²; the 2.03 bin is 1.85–2.20 GeV),
  so the window is as close as the data allow and the row tests the weight
  below the generator's own W range.
- **(b)** per verdict block, the weighted constant fit of the four backward
  bins (cos θ_pq = −0.9 … −0.3, i.e. [−1, −0.2], the paper's backward bin and
  normalisation region) with **stat** errors: R̄, σ, χ²_flat (ndf 3),
  S = max(1, √(χ²/ndf)), band half-width δ = S σ / R̄.
- **(c)** the tree: `lipolgen.GlauberFsiWeight(TaggedModel(deuteron_channel()),
  GlauberFsiOptions()).weight(FsiKinematics(k = p_s, cos_theta_k = cos θ_pq,
  …, spectator (1, 1), Channel.TaggedDeuteronP))` — a **library evaluation**
  through `python/bindings.cpp` `bind_fsi` (no Pipeline, no generation, so
  not EXPENSIVE). With the formation ramp off (`GlauberFsiOptions` defaults:
  σ_XN = 40 mb, ε = −0.5, B = 6 GeV⁻²) the weight reads its (k_z, k_T)
  table at |k_z|, so it is symmetric under cos → −cos
  (`include/lipolgen/fsi.hpp:40`, "symmetric under theta_k -> pi - theta_k")
  and W*, Q², x do not reach it. u(c) = w(c)/mean_bw(w); v = max_bw |u − 1|.
- **(d)** a block passes iff v ≤ δ; a beam sub-row iff every block passes;
  the row iff the default (Hulthén, P_D = 0.045, 40 mb) passes at both beams.
- **(e)** recorded, never in the verdict: AV18; σ_XN = 20 mb (the
  documented band's other end); the stat ⊕ syst band; χ² of the backward
  data against s·w(c) beside χ²_flat; the forward rise; the instant-form /
  light-cone factor; F₂ⁿ/F₂ᵖ.
- **(f)** the figure's blocks (Q² = 1.66, p_s = 0.093, the two highest W*
  per beam) fixed in the header before any number.

The four 4 GeV (3.38, 2.44) blocks have no backward bin and carry no verdict
(counted), leaving the 28 verdict blocks.

## 5. The measurement

`python3 validation/benchmarks/t2_bonus_spectator_shape.py` (exit 0,
2.3–2.4 s), the sub-row table as printed:

| sub-row | blocks inside | worst v/δ | Σχ²_flat | Σχ² vs the weight's shape | ndf | |
|---|---|---|---|---|---|---|
| `4GeV/hulthen-40mb` | 4/12 | 2.622 | 205.71 | 109.74 | 36 | **FAIL** (verdict) |
| `5GeV/hulthen-40mb` | 0/16 | 3.441 | 533.66 | 116.59 | 48 | **FAIL** (verdict) |
| `4GeV/av18-40mb` | 3/12 | 2.807 | 205.71 | 112.43 | 36 | (fail) |
| `5GeV/av18-40mb` | 0/16 | 3.560 | 533.66 | 109.19 | 48 | (fail) |
| `4GeV/hulthen-20mb` | 8/12 | 1.263 | 205.71 | 129.18 | 36 | (fail) |
| `5GeV/hulthen-20mb` | 13/16 | 1.663 | 533.66 | 274.12 | 48 | (fail) |
| `4GeV/hulthen-40mb stat+syst` | 8/12 | 2.324 | 23.49 | 19.02 | 36 | (fail) |
| `5GeV/hulthen-40mb stat+syst` | 8/16 | 1.807 | 33.89 | 9.13 | 48 | (fail) |

Per p_s (0.078 / 0.093 / 0.11 / 0.135 GeV/c): the default weight's v =
0.0344 / 0.0548 / 0.0848 / 0.1406 (AV18 0.0355 / 0.0555 / 0.0871 / 0.1505;
20 mb 0.0166 / 0.0264 / 0.0407 / 0.0678); the narrowest data band δ =
0.0138 / 0.0255 / 0.0506 / 0.0536 at 4 GeV and 0.0100 / 0.0257 / 0.0505 /
0.0672 at 5 GeV. The four 4 GeV blocks inside the band are the four at
(Q² = 3.38, W* = 2.03), where σ/R̄ is 6–12 %.

### 5.1 What the FAIL measures

The rule assumed the backward window flat within its errors, so that its
flatness would bound any FSI effect. At the stat level it is not
(χ²_flat = 739.4 on 84 dof). Read with the shape, the data fall from
cos θ_pq = −0.9 toward 90° in most blocks and the fall deepens with p_s —
e.g. (5 GeV, Q² = 1.66, W* = 2.03, p_s = 0.135) R/R̄ = 1.201, 1.091,
0.968, 0.840 against the weight's u = 1.126, 1.053, 0.962, 0.859 — and the
χ² against the weight's shape with one free scale per block drops to
226.3/84 (40 mb), against 403.3 for 20 mb. Over the whole range the data dip
near cos θ_pq = 0 and rise on **both** sides (same block: 0.811 and 0.799
at ∓0.1, 1.355 at +0.9), which is the paper's own reading ("a depression at
angles around 90° develops at higher p_s … consistent with expectations from
some FSI models", Sec. V A 1) and the shape of the symmetric eikonal weight.
**So the FAIL says that the shipped weight varies across the backward window
by more than the 1σ level of a flat fit; it does not say that the data
contradict the weight's direction — they share it — and none of the shapes
fits at stat (p ≪ 0.01 for 226.3/84).** Recorded, never tuned: the
tolerance stays as declared.

### 5.2 The forward rise (recorded)

ρ = R̄(cos > 0)/R̄(cos ≤ −0.3), stat-weighted over the blocks of each (beam,
p_s):

| p_s | 4 GeV ρ_data | 5 GeV ρ_data | ρ_tree |
|---|---|---|---|
| 0.078 | 1.036 ± 0.012 | 1.031 ± 0.008 | 0.9902 |
| 0.093 | 0.973 ± 0.010 | 0.994 ± 0.007 | 0.9838 |
| 0.11 | 1.020 ± 0.013 | 0.974 ± 0.008 | 0.9741 |
| 0.135 | 1.026 ± 0.017 | 1.030 ± 0.012 | 0.9552 |

The symmetric weight's forward mean sits below its backward one (the forward
set holds the cos = 0.1 bin of its dip); same sign in 12 of 28 blocks. The
data's excess beyond the mirror at the most forward bins is what a cos → −cos
symmetric weight cannot carry by construction (the principal-value piece the
eikonal pole approximation drops, `include/lipolgen/fsi.hpp`); the paper
reads it as target fragmentation.

### 5.3 F₂ⁿ/F₂ᵖ (recorded)

The shipped `ToyF2.f2n_over_f2p(x)` = clip(1 − 0.75x, 0.25, 1)
(`src/core/sf.cpp:110`, `ToyF2::f2n_over_f2p`) at x*_aver against the
W* > 1.8 GeV table, stat ⊕ syst: tree/data = 1.190 at x* = 0.177, 1.103 at
0.32, 0.990 at 0.52, 1.286 at 0.62; χ² = 40.5/10 (p = 1.4e−5) absolute,
8.30/9 (p = 0.50) with one free scale 0.902; on the 5 rows with W*_aver ≥
2 GeV, 35.5/5 and 1.70/4 at scale 0.878. The table is cross-normalised to
the CJ fit at x = 0.32 (Fig. 21 caption), so the absolute miss is ToyF2's
level against CJ's; the shape agrees. A second anchor beside
`t1_nmc_f2d_over_f2p`, no tolerance.

### 5.4 The PWIA prescription (recorded, not FSI)

The tree's Hulthén density, instant-form n(k = p_s) over light-cone
n(k_LC) E_k/(2 − α_s) (the map `t2_deeps_spectator_tail.py` codes; a test
checks the two agree to 1e−14), normalised to its backward mean: max
deviation 0.0227 / 0.0309 / 0.0409 / 0.0571, and f(−0.9)/f(−0.3) =
1.046 / 1.063 / 1.084 / 1.118 against the weight's 1.071 / 1.114 / 1.179 /
1.310 — **the same direction**, 40–66 % of the weight's variation. In the
backward window alone R_D/S cannot tell the weight's dip from a
light-cone-vs-instant-form difference; only the full-range shape (a dip with
a rise on both sides, where the prescription factor is monotonic) separates
them. The density is the tree's own: `TaggedModel.n_of_kc` is a nearest-cell
lookup (it reads the cell below the query), so the harness rebuilds
Σ_L (a_L ψ_L(k))²/4π from `Wave.radial` with the model's own a_L, and a
test checks it equals the M-averaged `n_of_kc` table on the grid to 1e−10.

### 5.5 Checks (asserted in the tests)

The weight's v is converged in its (k_z, k_T) table: the 161 × 161 default
against 321 × 321 on [0, 0.6] GeV moves v by < 5e−4 (measured < 2e−4 in a
scratch probe, 481 × 481 on [0, 0.4] GeV likewise). Evaluated at the bin
labels; averaged uniformly over each p_s and cos bin, v moves by −1.27 /
−1.88 / −0.42 / −0.86 % of itself (asserted < 2 %). Exactly symmetric:
max |w(c) − w(−c)| = 0. Clipped grid fraction of the default weight 0.01451
(the high-k_T cells, far from k ≤ 0.135 GeV).

## 6. What the verdict means for the channels it sits under

Ingredient I of `08_channel_priority.md` §4 — the spectator FSI weight —
is **spin-blind by construction** (`include/lipolgen/fsi.hpp`: the
unpolarized-shape distortion, "Quote it as an unpolarized-shape SYSTEMATIC,
never as a correction to A_zz"): nothing here moves A_zz^tag or any tensor
observable of rank 6. What the weight moves is the **tagged spectator
spectrum** — its k and angle distribution, hence the acceptance-folded
tagged rates. For the d–p control (rank 7) this row says, at p_s =
70–150 MeV/c and W* ≈ 2–2.7 GeV: the shipped 40 mb weight has the data's
**direction** (more absorption toward 90°, deepening with p_s, symmetric),
its backward-window size is closer to the data than the 20 mb end
(Δχ² = 177), and the backward window that BONuS uses to extract the free
neutron is **not** FSI-free at the stat level in the DIS blocks — the
must-not-contradict rule written on that premise fails. Beside row 11
(Deeps, p_s = 0.30–0.53 GeV/c: the weight points the right way at transverse
angles in 4 of 5 p_s bins but overshoots 2–6×), the two rows read the same
weight at low and high k; they are recorded side by side, not combined. For
tagged ⁶Li α (rank 6) the row is indirect: a nucleon spectator, not the α
(whose cluster profile, σ_Xα = 131 mb at 40 mb, is a different input), and
no α-spectator datum exists (`BENCHMARK_PLAN.md` §6 U-2, U-7).

## 7. What this row does NOT validate

The spin dependence of the rescattering (none is modelled; Cosyn–Weiss list
it as open); the weight's absolute level (the survival probability is
normalised away with N(p_s, E_b)); the p_s spectrum of n(k); the weight at
p_s > 0.15 GeV/c (row 11) or at EIC W, where formation lengths argue for the
lower end of the 20–40 mb band (the σ_XN(W) ramp is off and untested here,
and BONuS's W* sits below the generator's W² ≥ 8 GeV² floor); the α and
triton spectators of the lithium channels; the generator's sampling, binning
or acceptance (no event is generated); and the frame: BONuS tags backward in
the target rest frame, the EIC a boosted spectator in the far-forward
acceptance (the 2026-09-23 frame caveat, kept in the header).

## 8. Tests

`python/tests/test_bench_bonus.py` (new, 26): file guards —
`test_bonus_file_is_the_record_and_refuses_a_flipped_byte`,
`test_the_recorded_original_digest_is_checked_on_its_own`,
`test_history_files_are_the_vendored_ones`; the blocks —
`test_census_is_the_header_s`,
`test_window_is_on_the_labels_and_missing_rows_are_dropped`,
`test_a_duplicate_block_is_refused`; statistics — `test_flat_band_by_hand`,
`test_scale_fit_absorbs_a_common_factor`, `test_chi2_sf_fallback_is_t1s`;
the tree — `test_weight_is_live_symmetric_and_reads_only_k_and_abs_cos`,
`test_weight_variation_is_converged_in_its_table`,
`test_bin_averaging_moves_the_variation_by_under_two_percent`,
`test_density_is_the_trees_own_n_of_k`, `test_light_cone_map_is_the_deeps_rows`,
`test_toyf2_hook_is_the_shipped_formula`; the row — `test_report_row_shape`,
`test_verdict_follows_from_the_numbers`,
`test_measured_values_pinned_a_recorded_fail`,
`test_forward_rise_recorded_not_asserted`,
`test_prescription_and_f2n_f2p_recorded`,
`test_header_declares_the_rule_the_code_applies`,
`test_report_normalises_for_make_report`; exit codes —
`test_harness_exits_0_and_prints_one_report_row`,
`test_harness_prints_its_row_on_a_non_utf8_stdout`,
`test_harness_exits_2_on_an_altered_or_missing_file`; the figure —
`test_plot_writes_a_png`.

`python/tests/test_bench_spin_deuteron.py` (not this worker's file) carries
four `test_bonus_*` tests: `test_bonus_files_are_the_vendored_ones`,
`test_bonus_deeps_census_matches_the_docstring` and
`test_bonus_harness_prints_its_row_on_a_non_utf8_stdout` still pass (the
harness keeps `load_query`, `load_deeps`, `QUERY`, `DEEPS`, `findings()`, the
history printout and its escaping); **`test_bonus_row_is_blocked_by_name`
now FAILS** (`assert rep["status"] == "blocked"` meets `'fail'`) and must be
retired in the same commit — the integrator's edit. `test_bench_deeps.py`
(15) and `test_spdx_headers.py` pass unchanged.

## 9. Commands

```bash
cd /home/cpeng/Projects/polli/PoliGen && source env.sh
python3 validation/benchmarks/t2_bonus_spectator_shape.py            # exit 0, FAIL row
python3 validation/benchmarks/t2_bonus_spectator_shape.py --plot     # + the figure
python3 -m pytest -p no:cacheprovider python/tests/test_bench_bonus.py -q          # 26 passed
python3 -m pytest -p no:cacheprovider python/tests/test_bench_spin_deuteron.py -k bonus \
        python/tests/test_bench_deeps.py -q   # -k over both files: 4 passed (3 test_bonus_*,
                                              # the Deeps digest test), 1 failed (test_bonus_row_is_blocked_by_name)
python3 -m pytest -p no:cacheprovider python/tests/test_bench_deeps.py \
        python/tests/test_spdx_headers.py -q  # 21 passed
python3 validation/benchmarks/make_report.py --only t2_bonus_spectator_shape \
        --out-dir <scratch>/report_B1                       # 0 pass, 1 fail; the row renders
python3 validation/check_spdx_headers.py                    # 130 files, all carry the header
sha256sum validation/benchmarks/data/bonus2014_tkachenko_supplemental.json
```

## 10. Open items

- The licence question for the supplemental tables (APS copyright on the
  document; numbers vendored as facts with attribution) is the author's,
  with the filing of the original into `PolarizedLithiumSim/refs/`.
- `BENCHMARK_PLAN.md` §7 maps `fsi.*` to T3-cluster and T4-tagged only, so a
  change to `GlauberFsiWeight` re-runs neither T2 row that evaluates it (this
  one and row 11) — an integrator's edit, with `tiers.py` mirrored.
- Whether the tolerance's premise should be re-declared for a future run
  (e.g. a shape χ² of the weight with a declared p-value) is a design choice
  for the author; this run records the declared rule's FAIL and the
  pre-declared χ² beside it, and changes nothing.
