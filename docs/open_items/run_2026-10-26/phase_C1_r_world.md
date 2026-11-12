<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Phase C1 — R = σ_L/σ_T against the world's R, with NMC's R at the sweet spots primary; HERMES F₂ᵈ

Written 2026-10-26 by the C1 worker against base commit `87787dd` (the
shared tree, dirty with other workers' files), after the runs, from output
quoted with the command that produced it. **Nothing in the tree was tuned**:
no physics default, constant, registry status, tolerance or reference JSON
moved; nothing under `include/`, `src/` or `python/lipolgen/` was edited.
The plan item is `PLAN.md` phase C1; the ranking it serves is
`../../benchmarking/08_channel_priority.md` §4 ingredient A (rank 1, the
cos 2φ amplitude) and, for HERMES, rank 4 and `../../benchmarking/02_data_nucleon_deuteron.md`
§10 rank 2.

Environment for every command below:

```bash
cd /home/cpeng/Projects/polli/PoliGen && source env.sh
```

## 0. Summary

| harness (`validation/benchmarks/`) | row status | the number |
|---|---|---|
| `t1_r_world.py` (new) | **PASS** | the shipped `r_sigma_lt` on NMC's R (PRIMARY, 10 points): χ²/ndf = **12.54/10, p = 0.250** (stat ⊕ sys diagonal, the declared tolerance p ≥ 0.01). Recorded beside it: `r1998` 6.20/10 (p = 0.798); with NMC's "largely correlated" systematic taken 100 % correlated, `r_sigma_lt` **60.22/10 (p = 3.3e−9)** and `r1998` 21.28/10 (p = 0.019). H1 F_L (EXTRAPOLATION): both pass (8.82/8, 8.18/8). Whitlow's SLAC R (CLOSURE): `r_sigma_lt` **FAILS 209.76/68 (p = 2.3e−16)**, `r1998` passes 63.05/68. E140 (CLOSURE): `r_sigma_lt` **FAILS 29.66/9 (p = 5.0e−4)**, `r1998` passes 12.71/9 |
| ″ — the bias table | recorded | (1 + r1998)/(1 + r_sigma_lt) at the letter's four sweet-spot cells = **1.1848 / 1.1780 / 1.0762 / 0.9563** (x = 0.0112 / 0.0282 / 0.0708 / 0.141); the generator's own cos 2φ amplitude ratio a₂(R1998)/a₂(shipped), read from the kernel = **0.8441 / 0.8489 / 0.9292 / 1.0457**, i.e. **−15.6 / −15.1 / −7.1 / +4.6 %** |
| `t1_hermes_f2d.py` (new) | **FAIL, recorded** | the shipped ToyF2 on HERMES F₂ᵈ per nucleon: **5303.5/68 (p = 0)** with the per-bin errors diagonal and HERMES's 7.5 % overall normalisation as one correlated component. Sub-rows: MSTW 318.9/68 (p = 4.0e−34; ⟨th/data⟩ 1.013), CT18NLO at LO 1788.2/45 (p = 0; 1.059). Proton recorded: 5678.5 / 247.8 / 1760.9 |
| ″ — the HERMES b₁ inversion | recorded | F₂ᵈ implied by HERMES's b₁ Table II over HERMES's own measured F₂ᵈ = **0.929 / 1.039 / 0.986 / 0.958 / 1.006 / 1.006** in the six bins → **per nucleon** (rule declared first: all six in [0.8, 1.25]) |

Tests: `python/tests/test_bench_t1_r_world.py`, **51 passed** in 8.9 s
(§9). Both harnesses run in < 2 s; neither is behind `LIPOLGEN_BENCH`.

## 1. Files

New, all this worker's: `validation/benchmarks/t1_r_world.py`,
`validation/benchmarks/t1_hermes_f2d.py`, five vendored files under
`validation/benchmarks/data/` — `hepdata_ins424154_nmc_r.json`,
`hepdata_ins786161_h1_fl.json`, `hepdata_ins296980_whitlow_r.json`,
`hepdata_ins360765_e140_r.json`, `hepdata_ins894309_hermes_f2.json` —
`python/tests/test_bench_t1_r_world.py`, the two figures
`docs/benchmarking/figures/t1_r_world.png` and `t1_hermes_f2d.png` (each
written by one `--plot` run of the final harness; the `t1_r_world` one twice,
the first with a superseded legend), and this file. **No existing file was
edited.** `t1_hermes_f2d.py` also reads two older files,
`hermes2005_b1d_table2.json` and `whitlow1990_r1990.json`, under the sha256
`t2_hermes_b1_table2.py` records for them (asserted by a test).

## 2. The data

**Fetched** 2026-10-26 18:54–19:29 UTC with the standard Python `requests`
2.25.1 client, default identity (curl's default identity gets 403 from
HEPData and was not used; nothing was spoofed). Raw downloads stay in the
session scratch directory (`…/scratchpad/hepdata/<ins>/`); only numbers are
vendored. 28 of the E140 per-table downloads (yaml and json of Tables 8–19,
json of 20, 22, 24, 27) first came back as HEPData's **HTTP 429 "Too many
requests"** page; the yaml of those twelve tables and the json of the three
tables used were re-fetched one at a time, three seconds apart (no other
change), and every file was checked to be the table, not an HTML page,
before it was read.

Each vendored file carries a `_provenance` block — the record (INSPIRE id,
HEPData DOI, version, catalogue URL and sha256), every table used with its
yaml and json download URL, sha256, size and fetch time, the licence (**CC0**,
read from each table's `table_license` field), the paper, the publishing
convention, the errors, the paper cross-check and whether the data are an
input of the R1998 fit — and a `_body_sha256`; the harness pins both digests.

| file | record (HEPData DOI) | tables | points | checked against the paper |
|---|---|---|---|---|
| `hepdata_ins424154_nmc_r.json` | ins424154 (`10.17182/hepdata.32752.v1`) | 33 | 10 | all ten rows = NPB 483 (1997) 3 Table 7 (hep-ph/9610231 p. 14) |
| `hepdata_ins786161_h1_fl.json` | ins786161 (`…45340.v1`) | 9 (averaged), 1–8 (per bin) | 8 + 25 | the eight averaged rows = PLB 665 (2008) 139 Table 2 (arXiv:0805.2809 p. 14); their stat ⊕ uncorr ⊕ corr = the paper's total to the printed digit |
| `hepdata_ins296980_whitlow_r.json` | ins296980 (`…29555.v1`) | 3–13 (global R), 1 (R_d − R_p) | 89 | not possible: the letter prints no table ("given in ref. 11", a thesis) |
| `hepdata_ins360765_e140_r.json` | ins360765 (`…22468.v1`) | 23–25 (deuterium R) | 11 | all eleven = PRD 49 (1994) 5641 Table VII, D rows (SLAC-PUB-5814 p. 30) |
| `hepdata_ins894309_hermes_f2.json` | ins894309 (`…66147.v1`) | 1 (p), 2 (d) | 81 + 81 | every row of both (bins, ⟨x⟩, ⟨Q²⟩, F₂ and the five per-bin uncertainty columns) = JHEP 05 (2011) 126 Tables 5, 6 (arXiv:1103.5704 pp. 30–34), checked by machine |

Facts found while vendoring, recorded in the files:

- **HEPData record defect (ins894309):** Table 2's d²σ/dxdQ² column (the
  deuteron's) is identical, in all 81 rows, to Table 1's (the proton's), and
  equals the paper's proton Table 5 (81/81), not its deuteron Table 6 (1/81).
  The F₂ᵈ column and every uncertainty column agree with Table 6. Only F₂ is
  vendored. Worth reporting to HEPData.
- **The R NMC extracted its F₂ with** (an open item of `t1_nmc_f2d.py`'s
  header, owned by task D): in the F₂ᵈ Tables 17–26 (x = 0.0045 … 0.11) the
  independent variable SIG(L)/SIG(T) is one value per table, equal to NMC's
  own measured R (Table 33) at that x in all ten; in Tables 27–32 (x ≥ 0.14)
  it is Q²-dependent — the paper's sec. 4: "At x > 0.12 … we used the R
  parametrisation given in ref. [2]" (R1990). Checked against the 16 yaml
  downloads, whose sha256 are in the NMC file's `f2d_tables_checked`.
- **NMC's R is an R1998 fit input.** Abe *et al.*, PLB 452 (1999) 194
  (hep-ex/9808028 p. 3): "We have performed a new fit using the present data
  [3,11–14]", with [12] = "NMC, M. Arneodo et al., Nucl. Phys. B483, 3
  (1997)" — the table of the PRIMARY sub-row — and [3] = Whitlow *et al.*
  1990, whose record carries E140's R as its Table 2. So for `r1998` the
  nmc, r1990 and e140 sub-rows are all closures; only h1 is independent of
  the fit, and it lies outside its support. For `r_sigma_lt`, which is not a
  fit, all four are tests. Stated in the harness header before any number.
- **Whitlow's global extraction excludes E140's deuterium cross sections**
  (Whitlow p. 6), so its Tables 3–13 and E140 are independent data; its
  Table 2 (E140 in Whitlow's analysis, R_d/R_Fe averaged) is not used.
- **E140's R carries an additional normalisation error of 0.03** "primarily
  due to radiative corrections" (Table VII caption) that is not in the
  record's numbers; kept in the file and used only in a recorded χ².
- **COMPASS (PLAN item C1 (3), optional):** the final COMPASS g₁ᵈ, HEPData
  ins1501480 Table 1 (`10.17182/hepdata.78374.v1/t1`, yaml sha256
  `23cc1f83b41a24f60e1945b7ca4b8c027676622219d613674758d6305f075844`), equals
  the vendored `nnpdf_COMPASS15_NC_NOTFIXED_MUD.json` that `t1_g1d_world`'s
  `compass` sub-row reads: **15 points each, and all 15 (x, Q², g₁, stat,
  sys) identical** — first (0.0046, 1.1, −0.13 ± 0.17 ± 0.11), last (0.567,
  60.8, 0.0203 ± 0.0035 ± 0.0035). No COMPASS harness was added.

## 3. `t1_r_world` — R against data

### 3.1 Design (the header, written before any comparison)

The tree's R backends at each data point's own (x, Q²): `r_sigma_lt` =
0.18/(1 + Q²/50) (declared at `include/lipolgen/sf.hpp:46`, defined at
`src/core/sf.cpp:52`; the null `RFunc` hook resolves to it — the default R of
the whole program), `r1998` with `clip = True` as shipped (its support
0.005 ≤ x ≤ 0.86, 0.5 ≤ Q² ≤ 130 GeV², from `R1998_X_MIN` at
`include/lipolgen/sf.hpp:53` and the three constants after it; the harness
copies the four numbers and checks them against the clip at run time), and
`r1998` unclipped (recorded). Sub-rows in this order: **nmc PRIMARY**, **h1
EXTRAPOLATION**, **r1990 CLOSURE**, **e140 CLOSURE**.

**Tolerance, declared before any comparison was computed:** per (sub-row,
backend) cell p(χ², ndf) ≥ 0.01, ndf = the number of points, χ² diagonal
with every published component in quadrature, on Q² ≥ 1 GeV² and W² ≥ 4 GeV²
— the T1 rows' rule applied unchanged. **The row's status = (nmc,
r_sigma_lt).** What had been seen first is stated in the header: every data
table (so NMC's 0.537 at x = 0.0045 against r_sigma_lt's 0.175 was visible),
one evaluation of each backend at one point and the kernel's amplitude
ratio at one cell while probing the API, and `docs/USAGE.md`'s r1998 prices
at x = 0.05 / 0.10 / 0.30. No χ², pull or p had been computed.

**H1 is compared as F_L, not R.** The record carries F_L only — no F₂, no
covariance — so R = F_L/(F₂ − F_L) cannot be formed from it. The tree side
is the tree's own proton `flp` = F₂ᵖ R/(1 + R) with R from the backend under
test (through `set_r_func`; the null hook is r_sigma_lt). F₂ is the nuisance;
the verdict cells take F₂ᵖ from MstwSF (MSTW2008 LO, a fit that includes
HERA F₂), and record ToyF2's (the shipped generator's own proton F_L) and
CT18NLO's (at LO) beside it. The verdict set is the eight averaged points
(Table 9); the 25 per-bin points are recorded. The figure shows H1 in R as
F_L/(F₂ᵖ_MSTW − F_L), the verdict's own F₂.

### 3.2 The measurement

```bash
python3 validation/benchmarks/t1_r_world.py --plot
```

Points in the DIS cut: nmc 10 of 10, h1 8 of 8, r1990 68 of 89, e140 9 of 11.
Outside R1998's support (in the cut): NMC's x = 0.0045 bin and all eight H1
points.

| sub-row | backend | χ²/ndf | p | cell | recorded beside it |
|---|---|---|---|---|---|
| **nmc** (PRIMARY) | **r_sigma_lt** | **12.54/10** | **0.250** | **PASS** (the row) | sys 100 % correlated: **60.22/10, p = 3.3e−9**; worst pull −2.81 at x = 0.0045 |
| nmc | r1998 | 6.20/10 | 0.798 | pass | sys correlated: 21.28/10, p = 0.019 |
| nmc | r1998 unclipped | 6.19/10 | 0.799 | pass | 21.22/10, p = 0.020 |
| h1 (EXTRAPOLATION; F₂ᵖ MSTW) | r_sigma_lt | 8.82/8 | 0.358 | pass | corr. sys correlated 13.37/8 (p = 0.100); 25 per-bin points 38.70/25; ToyF2-F₂ 6.96/8; CT18NLO-F₂ 8.48/8 |
| h1 | r1998 (frozen at x = 0.005) | 8.18/8 | 0.416 | pass | 13.43/8; per bin 37.65/25; ToyF2 6.48/8; CT18NLO 7.85/8 |
| h1 | r1998 unclipped | 8.11/8 | 0.423 | pass | 13.53/8; per bin 37.48/25 |
| r1990 (CLOSURE) | r_sigma_lt | **209.76/68** | **2.3e−16** | **fail** | worst pull −4.04 at x = 0.175, Q² = 1.25 |
| r1990 | r1998 | 63.05/68 | 0.647 | pass | (= unclipped: every point inside the support) |
| e140 (CLOSURE) | r_sigma_lt | **29.66/9** | **5.0e−4** | **fail** | + E140's 0.03 normalisation correlated: 29.46/9 |
| e140 | r1998 | 12.71/9 | 0.176 | pass | 12.07/9 |

Recorded, not compared: Whitlow's R_d − R_p = −0.001 ± 0.009 ± 0.009 and
E140's average R_A − R_D = 0.001 ± 0.018 (stat) ± 0.016 (syst) (its sec.
IV.B) — the tree's R has no target dependence, ⁶Li included. Without scipy
(the no-scipy χ² fallback) every p agrees to 1e−14 relative and every status
is the same (checked).

## 4. The bias table — the price of R on the cos 2φ amplitude

**The four cells** are the letter's sweet spots as
`PolarizedLithiumSim/plans/07` §7.4 WP1 addendum lists them (the toy → grid
drift table, e 10 × ⁶Li 99.5 GeV/u): spot 1 (x 0.0282, Q² 1.14), spot 2 (x
0.0112, Q² 1.14), spot 3 (x 0.0708, Q² 3.13), spot 4 (x 0.141, Q² 14.3);
§7.1's YR addendum quotes the same cells rounded, (0.028, 1.14), (0.011,
1.14), (0.071, 3.13), (0.141, 14.3).

**The generator's own factor, read from the kernel.** The cos 2φ amplitude is
a₂ = −(1−y)/y² c_eff sin²θ_S Δ/D_φ with D_φ = F₁ + (1−y)/(x y²) F₂
(include/lipolgen/xsec.hpp:16 and include/lipolgen/xsec.hpp:21; the term
itself at src/core/xsec.cpp:281), and the shipped Δ is `toy_delta_gluon` at
scale 10⁻² (`src/core/pipeline.cpp:1265`), = scale · F₁ · x^0.3 (1−x)^4
(`src/core/sf.cpp:377`). With F₂ = 2x(1 + R)F₁ this is a₂ ∝ F₁/D_φ =
1/[1 + 2(1−y)(1+R)/y²]: **F₂ cancels exactly, and at the sweet spots' y
(0.010–0.026) the amplitude moves as 1/(1 + R)** — because the toy Δ is
pinned to F₁, not through D_φ's own R dependence, which is suppressed there
by y²/[2(1−y)(1+R)] ≈ 4×10⁻⁵–3×10⁻⁴. The harness evaluates a₂ through
`default_inclusive_kernel(li6())` (shipped R) and
`default_inclusive_kernel(li6(), Li6Convolution, r_source=R1998)` (the only
path that accepts an R source; the Δ slot is bit-identical across b₁ models,
asserted), at y = Q²/(s x) with s = 3980.0 GeV² (`default_configs('6Li')[1]`).

| cell | x | Q² | y | R_lt | R_98 (spread; fit err) | (1+R_98)/(1+R_lt) | kernel a₂(R1998)/a₂(shipped) | nearest NMC R (x, ⟨Q²⟩) |
|---|---|---|---|---|---|---|---|---|
| spot 2 | 0.0112 | 1.14 | 0.02557 | 0.1760 | 0.3933 (0.082; 0.031) | **1.1848** | **0.84409** (−15.6 %) | 0.246 ± 0.116 (0.0125, 2.2) |
| spot 1 | 0.0282 | 1.14 | 0.01016 | 0.1760 | 0.3853 (0.062; 0.028) | **1.1780** | **0.84893** (−15.1 %) | 0.099 ± 0.096 (0.025, 4.5) |
| spot 3 | 0.0708 | 3.13 | 0.01111 | 0.1694 | 0.2585 (0.008; 0.016) | **1.0762** | **0.92917** (−7.1 %) | 0.113 ± 0.090 (0.07, 14.1) |
| spot 4 | 0.141 | 14.3 | 0.02548 | 0.1400 | 0.0902 (0.005; 0.008) | **0.9563** | **1.04566** (+4.6 %) | 0.043 ± 0.079 (0.11, 20.6) |

The kernel ratio equals the closed form [1 + 2(1−y)(1+R_lt)/y²]/[1 +
2(1−y)(1+R_98)/y²] to 1e−15 and 1/bias to ≤ 3.7×10⁻⁵ (the O(y²) term). The
sibling's own fast-sim measurement of the same swap ("R alone",
`plans/07` §7.4 WP1 addendum: 0.849 / 0.842 / 0.930 / 1.039 at spots 1–4,
or 0.855 / 0.848 / 0.936 / 1.046 after its bag re-solve factor 0.993 is
divided out) and this kernel's 0.849 / 0.844 / 0.929 / 1.046 agree within
0.75 % with either set — a different code and Δ model, so a consistency,
not a closure. **Over the grid**
x = 0.001 … 0.7 × Q² = 1 … 100 GeV², the ratio is 0.921 … 1.192 inside the
support (W² ≥ 4) and 1.004 … 1.192 in the extrapolated region x < 0.005,
where it is frozen in x at its x = 0.005 value (asserted); the largest
ratios sit at Q² ≈ 1–2 and the smallest at x ≥ 0.2, Q² ≈ 10–20.

## 5. What it means for rank 1 (the cos 2φ amplitude, Δ)

- **The R lever is about twice what `08_channel_priority.md` quotes for the
  two low-x sweet spots.** Its rank-1 row
  (docs/benchmarking/08_channel_priority.md:91 as of 2d89a64; phase E corrected its numbers) and ingredient A
  (docs/benchmarking/08_channel_priority.md:308 as of 2d89a64, line 327 since phase E corrected its numbers) give "−8.3 / −6.7 / −3.1 %
  at the sweet spots generator-side"; those three numbers are
  `docs/USAGE.md:345`'s, at x =
  0.05 / 0.10 / 0.30, y = 0.5, Q² = 2.5 — not the sweet spots. At the four
  cells the generator's own factor is **−15.6 / −15.1 / −7.1 / +4.6 %**
  (r1998 against the shipped r_sigma_lt), in line with the sibling's ≤ 18 %.
  The text belongs to the integrator (notes below).
- **The data do not pin R at the cells' Q².** NMC's R covers the cells' x,
  not their Q²: at x ≈ 0.011–0.028 its points sit at ⟨Q²⟩ = 2.2–4.5 GeV²
  against the cells' 1.14. Its errors there are ±0.08–0.12 in R, i.e.
  ±7.5–9.3 % on (1 + R) — larger than the swap's 7.6 % and 4.4 % on (1 + R)
  at spots 3 and 4, and about half of its 17.8 % and 18.5 % at spots 1 and 2. A letter that quotes Δ extracted at
  those two cells carries an R systematic of at least that size unless the
  analysis measures its own R.
- **The shipped R is not excluded at the sweet-spot x on the declared
  reading, and it is the wrong shape.** On NMC's R with the diagonal χ² both
  backends pass (r_sigma_lt p = 0.25, r1998 p = 0.80). Read with NMC's
  systematic fully correlated — the extreme the paper's "largely correlated"
  allows; the paper gives no covariance — r_sigma_lt is excluded (p = 3.3e−9)
  and r1998 is marginal (p = 0.019): NMC's R falls from 0.54 to 0.04 across
  the bins while 0.18/(1 + Q²/50) has no x dependence and falls only from
  0.175 to 0.128 across the bins' ⟨Q²⟩. At x ≥ 0.1 the SLAC data
  exclude r_sigma_lt outright (Whitlow p = 2.3e−16, E140 p = 5.0e−4), where
  r1998 passes — as a closure on its own fit input, which tests the
  transcription, not the physics. The PASS of the row is therefore the
  declared diagonal reading of the PRIMARY set; it is not evidence that
  r_sigma_lt is right.
- **Below x = 0.005 nothing distinguishes the two.** H1's F_L (errors
  0.10–0.39, 26 % to 140 % of the value) is consistent with r_sigma_lt and with r1998 frozen at x = 0.005
  (p = 0.36, 0.42); the extrapolation the clip makes is not contradicted, and
  not tested at better than those errors.
- Registry row 3 (`Li6ConvolutionOptions` keeps `r_sigma_lt`) is priced by
  these numbers, not decided; the default R is unchanged.

## 6. `t1_hermes_f2d` — HERMES F₂ᵈ and the b₁ inversion

### 6.1 Design

The `t1_nmc_f2d` structure re-used (nothing imported from it): ToyF2
(shipped, the verdict on the deuteron), MstwSF and CT18NLO (sub-rows,
blocked (environment) by name when absent), each against F₂ᵈ = (F₂ᵖ + F₂ⁿ)/2
(rule 2, asserted against `InclusiveKernel(deuteron).tables().f2`: max
relative difference 0 for all three) and F₂ᵖ (recorded), at each point's own
(⟨x⟩, ⟨Q²⟩), on the same DIS cut (Q² ≥ 1, W² ≥ 4, and ≥ CT18NLO's grid floor
1.677). **The error model differs from `t1_nmc_f2d`'s, declared first and
argued:** HERMES states its 7.5 % (d) / 7.6 % (p) normalisation is an overall
one, "dominated by the uncertainty of the year-dependent luminosity
constant" (sec. 5.1), so the verdict χ² takes the five per-bin components
(stat, PID, model, misalignment, Bethe–Heitler) diagonal and the
normalisation as one fully correlated component; recorded beside it: the
per-bin diagonal alone (strictest) and the normalisation in quadrature per
point (`t1_nmc_f2d`'s literal convention, most lenient), and the
normalisation pull. **Per nucleon** on both sides: HERMES normalises to "the
integrated luminosity L per nucleon" (its Eq. (4.3), factor A/Z); its F₂ᵈ is
0.036–0.354.

### 6.2 The measurement

```bash
python3 validation/benchmarks/t1_hermes_f2d.py --plot
```

| config | target | n | χ²/ndf (verdict) | p | normalisation pull | per-bin only | norm. in quadrature | ⟨th/data⟩ |
|---|---|---|---|---|---|---|---|---|
| **toy** | **deuteron** | 68 | **5303.5** | 0 | +0.69 σ | 5718.4 | 455.2 (p = 1.2e−58) | 0.965 |
| toy | proton | 68 | 5678.5 | 0 | +1.08 σ | 6692.9 | 515.1 | 0.944 |
| mstw | deuteron | 68 | 318.9 | 4.0e−34 | −0.28 σ | 385.9 | 31.4 (p = 1.000) | 1.013 |
| mstw | proton | 68 | 247.8 | 2.9e−22 | −0.45 σ | 421.4 | 31.7 | 1.025 |
| ct18nlo | deuteron | 45 | 1788.2 | 0 | −1.01 σ | 2535.0 | 154.5 | 1.059 |
| ct18nlo | proton | 45 | 1760.9 | 0 | −1.27 σ | 2963.0 | 171.9 | 1.083 |

Worst per-bin pulls: toy −26.6 at x = 0.49, Q² = 6.0 (ToyF2 is 35 % low at
high x); MSTW −9.3 at x = 0.66, Q² = 12.8; CT18NLO +16.5 at x = 0.135, Q² =
4.6. Every configuration fails; the order is that of phase 0's NMC F₂ᵈ row
(toy worst, MSTW best; CT18NLO at LO 5.9 % high here, 8.6 % on NMC), and on
HERMES too no shipped or opt-in backend passes absolute F₂ᵈ. With the normalisation thrown
into every point's error MSTW would "pass" at p = 1.000 — the reason that
reading is not the verdict.

### 6.3 The b₁ inversion (02 §10 rank 2)

HERMES's b₁ paper (PRL 95 (2005) 242001, Eq. (5)) extracted b₁ᵈ =
−(3/2) A_zzᵈ F₁ᵈ with a per-nucleon F₂ᵈ = F₂ᵖ(1 + F₂ⁿ/F₂ᵖ)/2 (ALLM97 × NMC)
and R1990. Inverting its Table II, F₁_implied = −(2/3) b₁/A_zz and F₂ᵈ_implied
= 2x(1 + R1990) F₁_implied/(1 + γ²), against this record's measured
per-nucleon F₂ᵈ at the nearest HERMES bin:

| b₁ bin (x, Q²) | F₁_implied | R1990 | F₂ᵈ_implied | HERMES F₂ᵈ (⟨x⟩, ⟨Q²⟩) | ratio |
|---|---|---|---|---|---|
| 0.012, 0.51 | 7.044 | 0.243 | 0.2099 | 0.226 (0.0117, 0.477) | 0.929 |
| 0.032, 1.06 | 3.427 | 0.383 | 0.3022 | 0.291 (0.0327, 1.144) | 1.039 |
| 0.063, 1.65 | 1.929 | 0.341 | 0.3234 | 0.328 (0.0596, 1.675) | 0.986 |
| 0.128, 2.33 | 1.018 | 0.262 | 0.3209 | 0.335 (0.1334, 2.158) | 0.958 |
| 0.248, 3.11 | 0.496 | 0.190 | 0.2735 | 0.272 (0.2724, 3.479) | 1.006 |
| 0.452, 4.69 | 0.161 | 0.129 | 0.1428 | 0.142 (0.4605, 4.224) | 1.006 |

**Per nucleon** by the rule declared first (all six in [0.8, 1.25]; per
deuteron would need [1.6, 2.5]): geometric mean 0.986. The premise of
`constants.hpp`'s B1_PER_DEUTERON_TO_PER_NUCLEON block — that HERMES's
published b₁ᵈ is per nucleon — is now checked against a published HERMES
table, not only argued from an inversion against ALLM97. It says nothing
about Miller's own normalisation (registry row 2, an author decision) and
nothing about b₁. One bin had been estimated by hand (≈ 1.05) before the
rule was written; the header says so.

**For rank 4 (A_zz / b₁(⁶Li)):** the per-nucleon convention HERMES's b₁ and
the tree's F₁ share is confirmed; the F₂ that A_zz divides by fails HERMES's
own F₂ᵈ for every backend at the per-bin level (the toy by ×5300/68), so the
default-backend decision of ingredient A stays the M2 question
`08_channel_priority.md` §4 states.

## 7. What these rows do NOT validate

R on a nucleus (no Li data; the tree's R has no A dependence); the tensor
L/T separation; R's Q² dependence inside an NMC x bin; R at the sweet-spot
cells' own Q²; the Δ model (the bias table prices R at fixed Δ/F₁); F₂ below
Q² = 1 GeV²; HERMES's σ_d/σ_p (its Table 3); b₁ itself. For h1 the F₂ is
MSTW2008 LO's, so that sub-row carries MSTW's F₂ error (the ToyF2 and
CT18NLO variants bracket it: 6.96 … 8.82 for r_sigma_lt).

## 8. Findings for files this worker does not own (exact text in the return)

1. docs/benchmarking/08_channel_priority.md:91 and
   docs/benchmarking/08_channel_priority.md:308 (both as of 2d89a64: phase E corrected the numbers and row A is at line 327 since) quoted "−8.3 / −6.7 / −3.1 %
   at the sweet spots generator-side"; those are docs/USAGE.md:345's x =
   0.05 / 0.10 / 0.30 at y = 0.5, Q² = 2.5. At the sweet spots: −15.6 /
   −15.1 / −7.1 / +4.6 %. docs/benchmarking/08_channel_priority.md:117's (as of 2d89a64; reworded in phase F) "R
   through D_φ" is right only with Δ/F₁ held fixed (the shipped toy); §4.
2. docs/references/REFERENCES.md:136 (as of 2d89a64; line 137 since phase E corrected the wording) said R1998 "is LiPolGen's default R
   … *Serves:* `f1_from_f2` everywhere"; the default is `r_sigma_lt`
   (`include/lipolgen/sf.hpp:46`), and r1998 is reached only through
   `--r-source r1998` on the `li6-convolution` path and the CDKS gate.
3. docs/references/REFERENCES.md:140 (as of 2d89a64; line 141 since phase E corrected the wording) and `02_data_nucleon_deuteron.md` §4
   R-2 call E140
   "the dominant input" / "the largest single R data set" of R1998. Abe *et
   al.* fit 237 points from [3, 11–14] and their own; E140 enters through
   Whitlow's [3] with 11 deuterium (and 11 Fe/Au) R values, against Whitlow's
   89 global points — not shown to be dominant. And NMC's R is itself an
   R1998 input (§2).
4. `t1_nmc_f2d.py`'s header (task D's): "the R NMC extracted F2 with is not
   stated in the file" — now known (§2).

## 9. Tests

```bash
python3 -m pytest -p no:cacheprovider python/tests/test_bench_t1_r_world.py -q
```

**51 passed** (8.9 s): the five vendored files (point counts, first and last
values, a flipped byte in the body and one in the provenance refused, the
body digest checked on its own, the two older files guarded by
`t2_hermes_b1_table2.py`'s digests); the R1998 support constants against
the clip; the DIS-cut counts; the theory side is the tree's backends and
`flp`; the h1 cells block by name without MSTW; the verdict cell's χ² and
the correlated χ² recomputed from the file; every cell pinned (recorded
FAILs pinned as FAILs); verdicts follow p ≥ 0.01; the bias cells, the
kernel = closed form, the frozen extrapolated region; HERMES's theory is the
isoscalar nucleon and the kernel, its covariance χ² the profiled
normalisation, its pins, the b₁ inversion and the rule both ways; one
REPORT line and every SUBROW line; `_report.normalise_row`; exit 0 on a run,
exit 2 on four altered files; `plot()` writes a PNG; the χ² fallback is
`t1_nmc_f2d`'s and scipy's. Also run: `python3 validation/check_spdx_headers.py`
(130 files, all carry the header) and `python/tests/test_spdx_headers.py`
(6 passed).

## 10. Commands

```bash
python3 validation/benchmarks/t1_r_world.py --plot
python3 validation/benchmarks/t1_hermes_f2d.py --plot
python3 validation/benchmarks/make_report.py --only t1_r_world,t1_hermes_f2d \
    --out-dir /tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/report_C1
python3 -m pytest -p no:cacheprovider python/tests/test_bench_t1_r_world.py -q
```

The partial report: **2 harnesses: 1 pass, 1 fail, 0 blocked, 0 blocked
(environment), 0 error**. Runtimes measured here: `t1_r_world` 1.6–1.8 s,
`t1_hermes_f2d` 0.5–0.6 s (wall, interpreter start-up included for the
latter's 0.6 s).
