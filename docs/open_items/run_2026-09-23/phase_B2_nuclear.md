<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Phase B2 — the NUCLEAR-INPUTS benchmark harnesses (BENCHMARK_PLAN.md §4 rows 4, 5, 6, 8)

Written 2026-09-23 against base commit `d735abc`, after the runs, from output
quoted with the command that produced it. Nothing here is an estimate; where a
number is *not* measured the sentence says so. **Nothing in the tree was
tuned**: no physics default, no registry row of
`run_2026-09-03/AUTHOR_DECISIONS.md`, no rtol-1e-12 reference JSON moved. The
only edit under `include/` is comment-only (§5.2) and changes no line count.

Environment for every command below:

```bash
cd /home/cpeng/Projects/polli/LiPolGen && source env.sh
```

## 0. Summary

| row | harness (`validation/benchmarks/`) | verdict | the number |
|---|---|---|---|
| 4 NMC F₂(⁶Li)/F₂(D) | `t3_nmc_li6_over_d.py` | **pass** | χ² = **4.6269/4** on x ≥ 0.30 (p = 0.328) — the survey's 4.63/4 **reproduced**; **26.941/15** on the whole on-grid set (p = 0.029) |
| 5 UVa FB ⁶Li charge density | `t3_li6_charge_ff_fb.py` | **fail, recorded** | tree C0 zero **3.0998** fm⁻¹ vs band **[2.6944, 2.8284]** fm⁻¹: **+0.2713 fm⁻¹ (+9.6 %)** above the band's upper edge; T11's [2.9, 3.3] untouched |
| 6 Rand–Frosch–Yearian magnetization | `t3_li_magnetization_rfy.py` | **blocked** | neither Phys. Rev. 144 (1966) 859 nor ADNDT 14 (1974) 479 is on disk; tree side measured and uncompared: r_mag(⁶Li) = 2.9469 fm |
| 8 George–Knutson relabel | — (no harness: not a benchmark) | **3 comment sites fixed** in `include/lipolgen/cluster_config.hpp`; **17 doc edits (G-1 … G-17), 3 dated-record notes and 1 test rename listed** for the Documents stage (§8) | — |

Test file: `python/tests/test_bench_nuclear.py`, 9 tests, all cheap (1.9 s):
**8 passed, 1 skipped** (the skip is row 6's, naming both documents).

## 1. Files

New, all mine: `validation/benchmarks/t3_nmc_li6_over_d.py`,
`validation/benchmarks/t3_li6_charge_ff_fb.py`,
`validation/benchmarks/t3_li_magnetization_rfy.py`,
`validation/benchmarks/README_nuclear.md` (merged into `README.md` there and removed by the Documents stage, 2026-09-23; another agent owns `README.md`
there), `validation/benchmarks/data/hepdata_ins394050_table1.csv`,
`validation/benchmarks/data/uva_ncd_fb_data_li6.dat`,
`python/tests/test_bench_nuclear.py`, this file.
Edited: `include/lipolgen/cluster_config.hpp`, comments only (§5.2).

Each harness carries BENCHMARK_PLAN.md §8's three rules in its header, stated
for its own row on both sides; reads its reference from the vendored file
(which carries URL, fetch date, licence, the paper, and a sha256 the harness
re-checks on every read — a one-byte change is refused, and a test asserts
that); prints one `REPORT | name | generator | reference | tolerance | STATUS`
line and returns it as a dict from `run()`.

## 2. Row 4 — NMC F₂(⁶Li)/F₂(D) against the EPPS21 baseline

### 2.1 The vendored table

Fetched once, 2026-09-23: `https://hepdata.net/download/table/ins394050/Table%201/1/csv`
(HTTP 200, 1667 bytes; the `www.` host is behind the Cloudflare challenge
`06_critic.md` §4.4 recorded, the bare host is not).
sha256 = `3d8b92785ac82f3fccc3b440e66e8eb0852375e6bdc7fc31dfef1b14e4c1b376`.
Table DOI `10.17182/hepdata.47955.v1/t1`; `table_license` **CC0**; the table's
own description: "Additional normalization error of 0.004 not included."
The file is the download byte for byte under a provenance header. Read in the
held preprint (`PolarizedLithiumSim/refs/hep-ex_9504002.pdf`): the ⁶Li targets
carried 4.5 % ⁷Li, the D targets 1.5 % H, and "a correction was also applied
to take into account the slight non-isoscalarity of the targets" — so the
published ratio is already isoscalar-corrected ⁶Li/D. All 24 values agree with
the survey's table (`03_data_nuclear.md` §1.2) to every printed digit.

### 2.2 The measurement

```bash
python3 validation/benchmarks/t3_nmc_li6_over_d.py
```

Tree side: `lipolgen.Epps21Ratio()` (grid `EPPS21nlo_CT18Anlo_Li6`, proton
`CT18NLO`) `.f2_per_nucleon(x, Q²)` ÷ ½(F₂ᵖ + F₂ⁿ) of `LhapdfSF("CT18NLO")`,
at each point's own (x, Q²). Errors stat ⊕ sys; the 0.004 normalisation is
not included. Points with Q² < 1.69 GeV² (EPPS21's `q2_min`; CT18NLO's is
1.677) are reported and **not** compared: **9 of 24**.

| subset | χ² / ndf | p | verdict (p ≥ 0.01) |
|---|---|---|---|
| x ≥ 0.30 (the valence window) | **4.6269 / 4** | 0.328 | pass |
| all on-grid, x = 0.0125 … 0.65 | **26.941 / 15** | 0.029 | pass |
| x ≥ 0.30, **WRONG** free-proton denominator | 85.17 / 4 | — | never used |

Pulls (data − tree)/σ on the 15 on-grid points, x ascending: +2.55, +2.34,
+1.59, −2.23, −0.24, −1.50, −0.43, −0.44, −0.05, −0.23, −0.31, +0.21, −1.45,
+0.28, −1.55. The whole-set χ² is carried by the four points x = 0.0125 …
0.035 — the shadowing / anti-shadowing crossover — not by the valence window.

**The survey's 4.63/4 is reproduced** (4.626864), at NMC's own Q², isoscalar
denominator, CT18NLO proton. (At Q² = 5 instead it would be 4.4399; the
harness uses the data's Q².) **One correction to the survey's table**: it
printed x = 0.0085, Q² = 1.4 GeV² as on-grid (0.9469); 1.4 < 1.69, so that
value is an LHAPDF extrapolation below the EPPS21 grid, and the harness
excludes it — 9 below-grid points, not 8.

The tolerance p ≥ 0.01 on **both** sets is this harness's choice, made after
the survey's number was known; it is written in the harness header so it can
be argued with. It was not tuned: both p-values are printed.

### 2.3 The depletion, like for like and not

| quantity | value | command |
|---|---|---|
| NMC weighted mean, the 4 points x ≥ 0.30 | 0.9618 ± 0.0219 | harness |
| tree R_iso, same 4 points, same Q², same weights | 0.9795 | harness |
| ⇒ like-for-like distance | −0.81 σ (data deeper) | (0.9618 − 0.9795)/0.0219 |
| ⟨1 − R_iso⟩, 301 x in [0.35, 0.65], Q² = 5, **CT18NLO** proton | 0.029803 | harness |
| `EMC_VALENCE_DEPLETION_EPPS21` (same definition, **CT18ANLO** proton, not installed here) | 0.031052 | `lipolgen.EMC_VALENCE_DEPLETION_EPPS21` |
| ⟨1 − R⟩ with the **WRONG** free-proton denominator | 0.2676 (**9.0×**) | harness, printed once |
| R_u, R_d of the EPPS21 grid at x = 0.65, Q² = 5 | 0.578, 2.829 | `Epps21Ratio.ratio(pid, 0.65, 5)`, pid 2 and 1 |
| WRONG ratio at x = 0.65 | 0.6880 at Q² = 5; 0.6924 at Q² = 39 | harness / §2.2 run |

The survey's "0.038 ± 0.022 against 0.031052, 0.3 σ" compares a four-point
weighted mean at Q² = 23–39 with a 301-point window mean at Q² = 5 on a
different proton baseline; the like-for-like distance at the data's own points
is **−0.81 σ**, still well inside one standard deviation. The survey's "0.692
at x = 0.65" for the wrong denominator is the Q² = 39 value.

**What this does not validate** (stated in the harness header, from
`03_data_nuclear.md` §1.4): the deuteron denominator — NMC divides by D, the
tree by the free isoscalar nucleon; no deuteron correction is applied (the
tree has none), and the direction of the omission is to make the true ⁶Li
depletion **larger**, i.e. up on 0.031. Nor the polarized EMC effect, ⁷Li, or
shadowing point by point.

### 2.4 Price for a registry row

No registry row names `EMC_VALENCE_DEPLETION_EPPS21` or the EPPS21 baseline
(measured: `grep -n -i 'epps\|emc' docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md
docs/open_items/run_2026-09-03/STATUS.md` → no hit). Nothing to price: the row
passes, and the data's ±0.022 on ⟨1 − R⟩ (±70 % of it) bounds a gross error
without refining the constant.

## 3. Row 5 — the UVa ⁶Li Fourier–Bessel charge density

### 3.1 The vendored row

Fetched once, 2026-09-23:
`https://discovery.phys.virginia.edu/research/groups/ncd/dldata/FB_data.dat`,
HTTP 200, 18 306 bytes (the size `06_critic.md` §4.5 recorded),
sha256(FB_data.dat) = `e977046004c97f03cbe77c908cad501e5d454d6eece0107351aef3e7b3577c5d`.
Vendored: the column-header line and the ⁶Li row only, byte for byte; the
row's own sha256 = `b9655bdc97151a58614cd949c5d5a03c54abfee8ae4609c2053c65e482b24a00`.
The seven coefficients and R = 6.0 fm are exactly those `06_critic.md` §4.5
quotes. **Licence: none stated** by the archive (its index page, read the same
day, says only that it collects "data from Atomic and Nuclear Data Tables,
Volumes 14, 36 and 60"); what is vendored is eight published fit parameters
with attribution. **Provenance of the row stays UNVERIFIED** (⁶Li is in neither
the 1974 nor the 1987 Table IV, `06_critic.md` §4.5). One trap, recorded in
the file: the archive's header labels the columns "Z A" but the values are
A, Z (⁶Li reads 6, 3; ³He reads 3, 2).

### 3.2 The measurement

```bash
python3 validation/benchmarks/t3_li6_charge_ff_fb.py
```

FB side, closed form (each term's transform analytic; checked in the test
against a 20 000-point midpoint quadrature of ρ(r), rel. 1e−6):

| quantity | value | survey |
|---|---|---|
| 4π∫ρr²dr | 2.991170 | 2.991170 |
| ⟨r²⟩^½ | **2.5206 fm** | 2.5206 |
| first zero of F_C(q) | **2.694413 fm⁻¹** | 2.6950 |

The survey's fourth decimal of the zero was off by 6e−4 (0.02 %); 2.6944 is
the closed-form value, bisected to 1e−12.

Tree side, `HoSpin1FF.for_ion(li6())`, shipped defaults:

| quantity | value |
|---|---|
| first zero of `fc` (= of F_point; G_E^p + G_E^n has none) | **3.099771 fm⁻¹** |
| minimum of the tree's own \|F_L\|² = F_C0² + F_C2² | 3.099771 fm⁻¹ (\|q_min − q₀\| = 1.7e−11, \|F_L\|² = 6.7e−30 there) |
| charge rms of the FOLDED `fc`, from its q → 0 slope | **2.5710 fm** |

The band and the distance:

| | fm⁻¹ |
|---|---|
| band lower edge, UVa FB zero | 2.6944 |
| band upper edge, Li *et al.* 1971 \|F_L\|² minimum, √8 | 2.8284 |
| tree q₀ | 3.0998 |
| q₀ − upper edge | **+0.2713 (+9.6 %)** |
| q₀ − FB zero | +0.4054 (+15.0 %) |
| gap between the band and T11's window [2.9, 3.3] | 0.072 — **they do not meet** |

REPORT line: `FAIL (recorded, +0.2713 fm^-1 above the band; open item Q1;
nothing moved)`. The test pins the FB numbers and asserts only that the
verdict follows from the numbers, so a future refit that lands in the band
turns the row to pass without a test edit.

Three findings beyond the survey's:

1. **The C2 escape does not exist inside the model.** `03_data_nuclear.md`
   §2.4 called the gap "not a contradiction yet" because Li *et al.*'s
   minimum is in F_C0² + F_C2². In `HoSpin1FF` F_q shares F_point with F_c,
   so C2 vanishes at the same q₀ and the model's own \|F_L\|² minimum is q₀
   (to 1.7e−11 fm⁻¹). What still separates the model zero from the measured
   minimum is Coulomb distortion (not computed here) and nothing in the model.
2. **The folded charge radius is 2.5710 fm, not the 2.589 fm it is built
   from.** `(a, α)` reproduce ⟨r²⟩_point = 6.0788 fm² exactly, which rc.hpp
   derives from r_ch = 2.589 with ⟨r²⟩_p = 0.7071, −⟨r²⟩_n = 0.1155 and a
   Darwin–Foldy 0.033 fm²; but `fc` folds with the dipole G_E^p (⟨r²⟩ =
   0.658 fm²) and Galster G_E^n (−0.127 fm²) and no Darwin–Foldy term, which
   gives 6.610 fm² → 2.5710 fm (measured from the slope, 1e−4 stable in the
   step). It sits between the FB 2.5206 and Angeli's 2.589, and just above
   de Vries's 2.54–2.57 spread. No document claims otherwise; recorded.
3. **The "uncited 0.31 GeV²" the tree withdrew has a source.**
   `run_2026-09-03/phase_B_numbers.md` §B1.9 withdrew `PHYSICS_CHANNELS.md`'s
   "dips near \|t\| ≈ 0.31 GeV²" as having "no source anywhere in the
   repository". 0.31 GeV² is q = 2.82 fm⁻¹ — Li *et al.* 1971's q² = 8 fm⁻²
   minimum, which `03_data_nuclear.md` §2 (written three days later) cites.
   The withdrawal stands as a record of its date; the sentence
   `OPEN_ITEMS_SOLUTIONS.md` carries about it is now wrong (§8, edit D-5).

### 3.3 The price — what adopting the FB shape would move (scratch only)

No registry row names `LI6_FF_HO_*` (measured: `grep -n 'LI6_FF_HO\|Q1'
docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md` → no hit). The HO refit is
open item **Q1** (`OPEN_ITEMS_SOLUTIONS.md`, "Q1 and Q10 stay open"), and its
nearest registry relative is **§B15** of `AUTHOR_DECISIONS.md` (promotion of
`c0_shape`), which has deliberately no row. This prices Q1 for the Documents
stage.

Scratch script, never in the tree:
`/tmp/claude-1000/-home-cpeng-Projects-polli/3ac6e55b-3bf7-4951-8170-ad408d0b635b/scratchpad/fb_swap.py`
— a Python `Spin1ElasticFF` subclass with F_c = Z·F_FB(q)/F_FB(0) (the FB
charge FF; no extra nucleon fold, it is already a charge density),
F_q = F_c·F_q^HO(0)/F_c^HO(0) (the tree's convention: C0 and C2 share the
monopole), F_m = the shipped `HoSpin1FF.fm`; set as
`make_config(isotope="6Li", channel="inclusive", config=1, events=2000,
seed=99, rc="tensor-band").rc_options.ff`, read through
`Pipeline(cfg, tensor_thirds_plan(0.0, 0.6)).rc_model` exactly as
`phase_B_numbers.md` §B1.10 does. **Every `ho` row reproduces §B1.5 to every
printed digit**, which licenses the `fb` rows beside it.

```bash
python /tmp/claude-1000/-home-cpeng-Projects-polli/3ac6e55b-3bf7-4951-8170-ad408d0b635b/scratchpad/fb_swap.py        # raw FB, all q
python /tmp/claude-1000/-home-cpeng-Projects-polli/3ac6e55b-3bf7-4951-8170-ad408d0b635b/scratchpad/fb_swap.py 3.66   # FB zeroed above Li71's 3.66 fm^-1
```

| x | Q² | `(1/6)σ^el_T/σ^el_U` ho | vmc-ft (§B1.5) | **fb** | `σ^el_U` fb/ho | `ΔA_zz` ho (shipped b₁) | **`ΔA_zz` fb** |
|---|---|---|---|---|---|---|---|
| 0.01 | 5 | −5.094304e−04 | −5.470470e−04 | **−5.612501e−04** | 1.0194 | −1.769796e−07 | **−2.240798e−07 (+26.6 %)** |
| 0.03 | 5 | −7.358312e−04 | −8.042604e−04 | **−8.262816e−04** | 1.0384 | — | — |
| 0.10 | 5 | **+1.559536e−04** | −4.070157e−05 | **−8.166613e−05** | 1.1521 | +2.127451e−09 | **+2.033221e−09 (−4.4 %)** |
| 0.30 | 5 | +1.328679e−02 | +1.555578e−02 | **+1.573541e−02** | **3.5322** | +5.130156e−12 | **+1.093263e−11 (+113 %)** |

(Q² = 2 and 10 rows, and σ^el_T ratios fb/ho = +1.1231 / +1.1661 / −0.6033 /
+4.1832 at x = 0.01 / 0.03 / 0.10 / 0.30, in the scratch output
`fb_swap.out`. ΔA_zz uses §8.1c-corr's shipped A_zz(Born) column of
`run_2026-09-02/phase_C_numbers.md`.) Zeroing the FB form factor above
3.66 fm⁻¹ moves nothing beyond the fourth significant digit (x = 0.30:
σ^el_U ratio 3.5314 against 3.5322), so its truncation tail (F_FB = 1.5e−4 at
q = 5 fm⁻¹ against the HO's 4.6e−11) is not what moves these.

What that says, and only that:

* **The data-derived shape lies OUTSIDE the tree's [ho, vmc-ft] band** on
  every row above: the tensor fraction is more negative than `vmc-ft` at
  x = 0.01, 0.03 and 0.10, more positive at x = 0.30, and σ^el_U at x = 0.30
  is ×3.53 against `vmc-ft`'s ×3.03. The band does not bracket the one
  data-derived C0 shape in hand.
* **The x = 0.10 sign question** (`phase_B_numbers.md` §B1.6: "INDETERMINATE
  across the band") gets a data-derived entry, and it lands on the **negative**
  side, with `vmc-ft`: −8.17e−05. It does not settle the sign (the FB row's
  provenance is unverified; `fq_scale` still spans zero), it adds a third
  point off one edge.
* ΔA_zz from the whole tail moves by ≤ 2.2e−07 in absolute value at every
  point priced — three orders below the band half-width 4.4e−04 §B1.6 quotes.
  "The tail is not the leading RC systematic" is untouched.

## 4. Row 6 — Rand–Frosch–Yearian magnetization: BLOCKED

Checked 2026-09-23: `ls PolarizedLithiumSim/refs` (204 entries) — no PDF
for Phys. Rev. 144 (1966) 859 or for ADNDT 14 (1974) 479 (both have `refs_dict.json` entries, `rand-frosch-yearian1966-magnetization`
and `dejager1974-adndt14-charge-magnetization`, each with `file: null`; corrected by the fix stage, this line read "no … `refs_dict.json` key"); the six `SLAC-PUB-*` files there were opened (page 1)
and are Mo–Tsai, Tsai, Tsai, Stein *et al.*, E143/E140x — none is RFY.
`docs/references/REFERENCES.md` §1 items 8 and 9 list both as PLEASE
DOWNLOAD, confirmed. The survey read ADNDT 14 Table V from a third-party
mirror (`06_critic.md` §4.6); that copy is not a licensed open-access source
and was **not** fetched or vendored, so its transcribed rows are not used.

```bash
python3 validation/benchmarks/t3_li_magnetization_rfy.py     # prints BLOCKED + the two documents
python -m pytest python/tests/test_bench_nuclear.py -q -rs   # the skip reason names both, with DOIs
```

The harness header writes out the comparison it will make: the tree's
q → 0 slope radius ⟨r²⟩_mag = −6 d[F_m/F_m(0)]/dq² = 6/q_z² + (3/2)b² of
rc.hpp's F_mag, against the envelope of all vendored ⁶Li rows ± their errors
(both model rows — the band, not a number); ⁷Li recorded only. It runs
unchanged once `validation/benchmarks/data/rfy1966_li_magnetization.csv`
lands (format in the header).

Measured now and **not compared**: r_mag(⁶Li) = **2.9469 fm** (numerical slope
and closed form agree to 1e−5; q_z = 1.30 fm⁻¹, b = 1.85 fm, both unfitted
placeholders). ⁷Li: `HoSpin1FF.for_ion(li7())` throws — the tree has no ⁷Li
F_m to compare. (For orientation only, and UNVERIFIED because the document is
not held: the survey's transcription puts ⁶Li at 2.81–3.42 fm.)

## 5. Row 8 — George–Knutson: relabel, not a benchmark

### 5.1 What is verified about the technique

The title, verified via Crossref by the reference series
(`docs/references/REFERENCES.md` §1 #31, `01_li-nuclear-data.md` §7.2):
*"Determination of the ⁶Li → α + d asymptotic D- to S-state ratio by a
restricted phase shift analysis"*, PRC 59 (1999) 598. So: a **restricted
phase-shift analysis**, not a direct measurement of d + α tensor analysing
powers. **Not verified in this tree**: which scattering system and data set
the analysis fits (BENCHMARK_PLAN.md §4 row 8 says "a ⁶Li + ⁴He phase-shift
analysis"; the paper is not held, and no abstract was fetched — the network
rule of this run allows fetching only for vendoring). Every replacement
wording below therefore names the technique by the title and nothing more.
Why it is not a D-wave benchmark: `asymptotic_ds_ratio()` is exactly linear in
the quadrupole dial (T22b), so η relabels the dial; it is a **consistency band
on a dial** (`06_critic.md` §3.1 C-1).

A second over-claim found by the sweep: `REFERENCES.md` #30,
`02_data_nucleon_deuteron.md` §5.2, `01_nucleon-deuteron-data.md` F-4 and
`00_corpus.md` call Rodning–Knutson's deuteron η_d "the same technique and
group" as George–Knutson. Rodning–Knutson is sub-Coulomb d + ⁴He tensor
analysing powers (as those entries say); George–Knutson is, by its title, a
restricted phase-shift analysis. Same group, not the same technique.

### 5.2 Fixed in `include/` — comments only, line count unchanged

`include/lipolgen/cluster_config.hpp`, three comment blocks (the file stays
733 lines, so none of the ~30 `cluster_config.hpp` line citations in the tree
move; `diff` shows only `///` lines):

* the first paragraph of the comment above `LI6_ETA_DS_GK` (lines 145–151
  of the file): "The MEASURED … from d + alpha elastic tensor analysing powers" →
  "The EMPIRICAL … '… by a restricted phase shift analysis' (the title;
  REFERENCES.md sec. 1 #31) -- NOT a direct d + alpha tensor-analysing-power
  measurement, as this line said until 2026-09-23 … it is a CONSISTENCY BAND
  ON THE QUADRUPOLE DIAL, not a D-wave benchmark (06_critic.md C-1)". The
  earlier sentence "NOT physics_literature.md, which predates this citation"
  was dropped to hold the line count; `[GK99]` of `PHYSICS_CHANNELS.md`
  stays named as the record.
* `asymptotic_ds_ratio()`'s doc comment: "against the MEASURED
  `LI6_ETA_DS_GK` …: a real but MODERATE D-wave excess, ~2x" → "against GK's
  phase-shift-analysis `LI6_ETA_DS_GK` …: the undialled tables sit ~2x off,
  not 5-15x -- a consistency band on the quadrupole dial, not a D-wave
  measurement."
* the same comment: "a budget leg anchored on a MEASUREMENT rather than on
  GFMC" → "anchored on an EMPIRICAL number, not on GFMC".

No `src/` site names George–Knutson (`grep -rn -i 'george\|knutson\|ETA_DS_GK' src`
→ none). `python/bindings.cpp` and `python/lipolgen/configs.py` carry the
constant and a band sentence that is already correct; nothing owed there.

## 6. Suite tallies

In §9, measured on the final tree.

## 7. Problems

* HEPData's `www.` host still serves a Cloudflare challenge; the bare
  `hepdata.net` host did not. Recorded so the next fetch does not repeat
  `06_critic.md` §4.4's conclusion that HEPData is unreachable from here.
* The UVa FB row's source paper remains unidentified; the harness says so.
* Li *et al.* 1971's minimum (q² = 8 fm⁻²) is used at the precision `03` §2
  quotes it; the paper is not held and its own precision is unknown.
* No abstract of George–Knutson was fetched (network rule); the scattering
  system of their phase-shift analysis is not verified by this record.

## 8. Doc edits requested for the Documents stage (file: site: exact replacement text)

Nothing below was edited by this phase. "Dated record" = a file under
`docs/open_items/run_*`: the request there is an **appended dated note**,
never a rewrite.

### Row 8 — George–Knutson

* **G-1** `docs/PHYSICS_CHANNELS.md`: reference [GK99] (item 55):
  "the measured asymptotic D/S ratio of the ⁶Li → α + d overlap, η = −0.025 ± 0.006 ± 0.010, from d + α elastic tensor analysing powers. The §10 wave functions give ≈ −0.05, about 2× that." →
  "the asymptotic D/S ratio of the ⁶Li → α + d overlap, η = −0.025 ± 0.006 ± 0.010, determined by a restricted phase-shift analysis (the paper's title) — not a direct measurement of d + α tensor analysing powers. The §10 wave functions give ≈ −0.05, about 2× that; because η is exactly linear in the quadrupole dial, this is a consistency band on that dial, not a benchmark of the α–d D wave (`benchmarking/06_critic.md` §3.1, C-1)."
* **G-2** `docs/benchmarking/00_in_tree_checks.md`: row C-1, first cell's "from d+α elastic tensor analysing powers" → "by a restricted phase-shift analysis (the paper's title; not a direct d+α tensor-analysing-power measurement)"; and the row's class: move C-1 out of Table C (DATA) into Table E with the verdict cell "a **consistency band on the quadrupole dial**, not a D/S benchmark: η is exactly linear in the dial (T22b), and T22 asserts `1.2 < |η/η_GK| < 3.0`, pinning a documented 1.93× distance of the UNDIALLED tables".
* **G-3** `docs/benchmarking/00_in_tree_checks.md`: the spot-check sentence "C-1 (T22 η)" → "C-1 (T22 η, a consistency band on the quadrupole dial)"; and the closing paragraph "several of its "gates" (C-1, C-2, E-13/T6, E-16)" → "several of its "gates" (C-1 — a consistency band on a dial, not a gate on data —, C-2, E-13/T6, E-16)".
* **G-4** `docs/benchmarking/00_in_tree_checks.md`: row E-18 title "George–Knutson band ↔ quadrupole dial ↔ a₂ map" → "the η consistency band on the quadrupole dial ↔ a₂ map (η from George–Knutson's restricted phase-shift analysis)". Its "does not validate" cell already says "bookkeeping validation of a model band, not a physics check" — keep.
* **G-5** `docs/benchmarking/BENCHMARK_PLAN.md` §4 row 8: "(it is a ⁶Li+⁴He phase-shift analysis, not d+α tensor analysing powers, and η is exactly linear in the quadrupole dial)" → "(it is a restricted phase-shift analysis — PRC 59 (1999) 598's own title; the scattering system it analyses is not verified in this tree, the paper is not held — not a direct d+α tensor-analysing-power measurement, and η is exactly linear in the quadrupole dial)"; and append "— **done 2026-09-23 for the code comments** (`run_2026-09-23/phase_B2_nuclear.md` §5.2); the doc sites are G-1 … G-17 there."
* **G-6** `docs/benchmarking/03_data_nuclear.md` §3 "What they do NOT validate" and §8 item 6: "The 1.93× η discrepancy of C-1" / "The 1.93× η discrepancy (C-1) has no better external anchor than George–Knutson" → "The 1.93× distance of the undialled tables from George–Knutson's phase-shift-analysis η (C-1, a consistency band on the quadrupole dial)" / "… has no better external anchor than that phase-shift analysis".
* **G-7** `docs/benchmarking/03_data_nuclear.md` §5 table, η row, last cell "Already C-1." → "Already C-1 — a consistency band on the quadrupole dial, not a D-wave benchmark."
* **G-8** `docs/benchmarking/04_theory.md` T-23 row: "that C-1's η gate currently reads only from the George–Knutson *measurement*" → "that C-1's η consistency band currently reads only from the George–Knutson *phase-shift analysis*"; and the "An α–d D/S ratio from an ab-initio calculation that could replace C-1's "gate a documented 1.93× discrepancy"" bullet → "… that could turn C-1's consistency band on a dial into an independent comparison — pending T-23."
* **G-9** `docs/benchmarking/02_data_nucleon_deuteron.md` §5.2 table, value row: "— **the same technique and the same group** as George & Knutson's η(⁶Li → α+d) = −0.025 ± 0.006 ± 0.010, which `00_in_tree_checks.md` row C-1 already gates" → "— the same group as George & Knutson's η(⁶Li → α+d) = −0.025 ± 0.006 ± 0.010, which was determined by a restricted phase-shift analysis (a different technique) and which `00_in_tree_checks.md` row C-1 carries as a consistency band on the quadrupole dial"; "does NOT validate" cell: "which the tree already knows is off by 1.93× (row C-1)" → "whose undialled model sits 1.93× from the phase-shift-analysis value (row C-1, a band on a dial)".
* **G-10** `docs/benchmarking/01_generators.md` (the TAGWF paragraph): "the D/S ratio η that `C-1` and `E-1` live on" → "the D/S ratio η (C-1's consistency band) and the tensor observable `E-1` lives on".
* **G-11** `docs/references/REFERENCES.md` §1 #30: "— same technique and group as the George–Knutson η(⁶Li → α+d) the tree gates on" → "— same group as the George–Knutson η(⁶Li → α+d), which is a restricted phase-shift analysis and which the tree carries as a consistency band on the quadrupole dial"; #31: "`LI6_ETA_DS_GK` and gate **C-1**" → "`LI6_ETA_DS_GK` and **C-1** (a consistency band on the quadrupole dial; retitled in the code comments 2026-09-23)"; §3 note "**#31** George–Knutson 1999 (`LI6_ETA_DS_GK`, gate **C-1**)" → "(`LI6_ETA_DS_GK`, consistency band **C-1**)".
* **G-12** `docs/references/00_corpus.md` row G9: "the `C-1` gate in `benchmarking/00_in_tree_checks.md`" → "the `C-1` consistency band in `benchmarking/00_in_tree_checks.md`"; and the Rodning & Knutson row "the deuteron analogue of the George–Knutson technique" → "by the same group as George–Knutson, by a different technique (sub-Coulomb tensor analysing powers vs a restricted phase-shift analysis)".
* **G-13** `docs/references/01_nucleon-deuteron-data.md` F-4: "the same technique and group as the George–Knutson η(⁶Li → α+d) the tree already gates on" → "the same group as the George–Knutson η(⁶Li → α+d) (a restricted phase-shift analysis, not this technique), which the tree carries as a consistency band on the quadrupole dial".
* **G-14** `docs/references/01_li-nuclear-data.md` §9 item 1 ("What this domain does *not* constrain") and §10.4: "The 1.93× η discrepancy (C-1) is untouched by all of it." → "The 1.93× distance of the undialled tables from the phase-shift-analysis η (C-1, a consistency band on a dial) is untouched by all of it."; "C-1's band is as good as the external world gets" — keep (already "band").
* **G-15** `docs/OPEN_ITEMS_SOLUTIONS.md` §11.2 (the O1 paragraph, "gives **η = −0.048** against the measured **η = −0.025 ± 0.006 ± 0.010** (George & Knutson, PRC 59, 598 (1999)): a real but *moderate* ≈ 2× D-wave excess") → "gives **η = −0.048** against **η = −0.025 ± 0.006 ± 0.010** from George & Knutson's restricted phase-shift analysis (PRC 59, 598 (1999)): the undialled tables sit ≈ 2× from it, not 5–15× — a consistency band on the quadrupole dial rather than a measurement of the D wave, since η is exactly linear in that dial"; and "the D-wave excess is moderate *and measured*" → "the D-wave excess is moderate *on the phase-shift-analysis η*"; "which lets a *measured* observable, not the GFMC number, supply the first leg" → "which lets an *empirical* number, not the GFMC one, supply the first leg".
* **G-16** `docs/USAGE.md` §9 table row "asymptotic η": "(measured `LI6_ETA_DS_GK` = …" → "(George–Knutson phase-shift analysis `LI6_ETA_DS_GK` = …"; "The asymptotic D/S ratio η says the excess is a real but *moderate* ≈ 2× effect" → "The asymptotic D/S ratio η puts the undialled tables ≈ 2× from George–Knutson's phase-shift-analysis value — a consistency band on the dial, not a D-wave measurement"; "the dial setting that matches the *measured* η" → "the dial setting that matches the *empirical* η"; "the measured η carries ±0.011662" → "the empirical η carries ±0.011662"; "are the single home of the measurement" → "are the single home of that number".
* **G-17** `docs/open_items/run_2026-09-03/mantysaari_collaboration_draft.md` — **OUTWARD-FACING DRAFT; the author decides** (registry row 21). Site: "anchored on the *measured* asymptotic α–d D/S ratio η = −0.025 ± 0.006 ± 0.010, George & Knutson" → "anchored on the asymptotic α–d D/S ratio η = −0.025 ± 0.006 ± 0.010 from George & Knutson's restricted phase-shift analysis". Nothing else in the letter moves.
* **Dated records — appended note only**, text for each: "*Note 2026-09-23 (`run_2026-09-23/phase_B2_nuclear.md` §5): George–Knutson's η is from a restricted phase-shift analysis (the paper's title), not from d + α elastic tensor analysing powers, and it is a consistency band on the quadrupole dial, not a D-wave benchmark. This record's wording is left as written on its date.*" — at `docs/open_items/run_2026-09-02/design_G_cluster_config.md` §2.7 ("it is measured … from d + α elastic tensor analysing powers") and §10 ("from d + α elastic tensor analysing powers"); `docs/open_items/run_2026-09-03/phase_C_numbers.md` §C1 ("George & Knutson's *measured* η"). `run_2026-09-03/PLAN.md`, `AUTHOR_DECISIONS.md` §B5/§B8 and `phase_C_numbers.md`'s σ-pull tables state distances on the band, not a benchmark claim: **no note owed**.
* **Tests (not in this phase's partition)** — `tests/test_cluster_config.cpp`: `TEST_CASE("T22 the asymptotic D/S ratio, Whittaker-divided")` → `TEST_CASE("T22 the asymptotic D/S ratio, Whittaker-divided -- a consistency band on the dial, not a D-wave benchmark")`; its `MESSAGE(... << " (FitRescaled), measured " << ...)` → `" (FitRescaled), George-Knutson phase-shift analysis "`; its comment "~2x the measured value, not 5-15x: a real but MODERATE D-wave excess." → "~2x the phase-shift-analysis value, not 5-15x: a consistency band on the quadrupole dial (eta is linear in it, T22b), not a D-wave measurement." The assertions are unchanged; the name is the only thing that claims more than it checks. E-18's cases (T22b, T23, "the design's sec. 8 table") are named for what they check — **no rename owed**. `python/tests/test_cluster_config.py`'s "THE BAND IS THE PHYSICS" docstring is already correct.

### Row 4 — NMC

* **N-1** `docs/benchmarking/03_data_nuclear.md` §1.2 table, row x = 0.0085: "0.9469" in the "at NMC's Q²" column → "*below grid* (1.4 < 1.69 GeV²; the LHAPDF extrapolation reads 0.9469)".
* **N-2** same section, first bullet: after "EPPS21's ⁶Li valence depletion is consistent with the only measurement of it." add "Over all 15 on-grid points χ² = 26.94/15 (p = 0.029), carried by the shadowing/anti-shadowing crossover (pulls +2.55, +2.34, +1.59, −2.23 at x = 0.0125 … 0.035); wired as `validation/benchmarks/t3_nmc_li6_over_d.py` 2026-09-23 (pass at p ≥ 0.01 on both sets)."
* **N-3** same section, second bullet: after "The agreement is 0.3 σ" add "(against the library's window mean at Q² = 5 on CT18ANLO; like for like — the tree at the same four points, Q² and weights, 0.9795 — the distance is −0.81 σ)".
* **N-4** §1.3: "gives 0.692 at x = 0.65" → "gives 0.6880 at x = 0.65, Q² = 5 (0.6924 at NMC's Q² = 39)".
* **N-5** `docs/benchmarking/BENCHMARK_PLAN.md` §4 row 4: "(χ² = 4.63/4 measured in the survey)" → "(χ² = 4.6269/4 on x ≥ 0.30, p = 0.33, and 26.94/15 on all 15 on-grid points, p = 0.029 — wired 2026-09-23, `validation/benchmarks/t3_nmc_li6_over_d.py`, pass)".
* **N-6** `docs/benchmarking/06_critic.md` §4.4: after "remains the only route from this environment" add "— superseded 2026-09-23: the bare host `hepdata.net` (no `www.`) served ins394050 Table 1 with HTTP 200 to a plain user agent (`run_2026-09-23/phase_B2_nuclear.md` §2.1)."
* **N-7** `docs/references/REFERENCES.md` §8 (the HEPData line): "`free, not fetched`, CC0" → "`free`, CC0; ins394050 Table 1 **vendored 2026-09-23** as `validation/benchmarks/data/hepdata_ins394050_table1.csv`".
* **N-8** (code comment, not in this phase's named sites) `include/lipolgen/lhapdf_sf.hpp`, the `Epps21Ratio` class comment: "the two proton baselines differ from each other by a few tenths of a percent at DIS kinematics" → "the two proton baselines differ by a few tenths of a percent in R, which is 4.0 % relative on the valence depletion ⟨1 − R⟩ (0.029803 on CT18NLO vs 0.031052 on CT18ANLO, `03_data_nuclear.md` §1.2)". Comment-only; must keep the line count (the file is cited).

### Row 5 — FB charge density

* **F-1** `docs/benchmarking/06_critic.md` §4.5 code block: "first zero of F_C(q)= 2.6950 fm^-1" → "first zero of F_C(q)= 2.6944 fm^-1   (closed form; 2.6950 in this file's first printing)"; and "**Provenance is the one open item**" paragraph: add "FB_data.dat re-fetched 2026-09-23, 18 306 bytes, sha256 `e9770460…3577c5d`; the ⁶Li row is vendored in `validation/benchmarks/data/uva_ncd_fb_data_li6.dat`. The archive states no licence."
* **F-2** `docs/benchmarking/BENCHMARK_PLAN.md` §4 row 5: "its zero at 2.695 fm⁻¹ and the diffraction minimum at 2.828 both sit OUTSIDE the tree's asserted [2.9, 3.3] window" → "its zero at 2.694 fm⁻¹ and the diffraction minimum at 2.828 both sit below the tree's asserted [2.9, 3.3] window — wired 2026-09-23, `validation/benchmarks/t3_li6_charge_ff_fb.py`: the shipped q₀ = 3.0998 is +0.2713 fm⁻¹ above the band [2.694, 2.828] (recorded FAIL, nothing moved; Q1 priced in `run_2026-09-23/phase_B2_nuclear.md` §3.3)"; §6 U-10: "(now: UVa density is a data-derived shape — #5 above)" → "(now: the UVa density is a data-derived shape, and it lies OUTSIDE the [ho, vmc-ft] band on every priced tail number — `run_2026-09-23/phase_B2_nuclear.md` §3.3)".
* **F-3** `docs/benchmarking/03_data_nuclear.md` §2 "What that means for Q1" item 4: after "the resolution is arithmetic the tree can do today." add "Done 2026-09-23: in `HoSpin1FF` the C2 shares the C0 monopole, so the model's own \|F_L\|² minimum IS q₀ (to 1.7e−11 fm⁻¹) — the C2 fill-in cannot reconcile it with 2.83; only Coulomb distortion, not computed, remains between the two."
* **F-4** `docs/OPEN_ITEMS_SOLUTIONS.md` (the RC paragraph ending "…withdrawn — `run_2026-09-03/phase_B_numbers.md` §B1.9"): "(the uncited "\|t\| ≈ 0.31 GeV²" that `PHYSICS_CHANNELS.md` carried until then is withdrawn — …)" → "(the "\|t\| ≈ 0.31 GeV²" that `PHYSICS_CHANNELS.md` carried until then was withdrawn as uncited — §B1.9; it is q = 2.82 fm⁻¹, Li *et al.* 1971's measured minimum, and the UVa FB density's zero 2.694 fm⁻¹ sits lower still: both lie below T11's window, `run_2026-09-23/phase_B2_nuclear.md` §3)"; and after "Q1 is unchanged by it." add "A data-derived C0 shape now exists and lies outside that band (§3.3 of the same record)."
* **F-5** `docs/open_items/run_2026-09-03/phase_B_numbers.md` §B1.9 — dated record, appended note: "*Note 2026-09-23: reason (i) is superseded — 0.31 GeV² ⇔ q = 2.82 fm⁻¹ is Li, Sick, Whitney, Yearian, NPA 162 (1971) 583's q² = 8 fm⁻² minimum (`benchmarking/03_data_nuclear.md` §2), and the UVa FB density puts the C0 zero at 2.694 fm⁻¹; see `run_2026-09-23/phase_B2_nuclear.md` §3. Reasons (ii) and (iii) and the withdrawal as a record of its date stand.*"
* **F-6** (code comment, not in this phase's named sites) `include/lipolgen/rc.hpp`, the `LI6_FF_HO_A_FM` block: "q_0 = 3.1 fm^-1 IS A STARTING GUESS, NOT A SOURCE.  WS98 does not locate the zero;" — keep, and replace the sentence "Fit q_0 over [2.9, 3.3] fm^-1 and pin the REFIT result, with the refit's own uncertainty (T11)." with "Two data-derived numbers sit BELOW [2.9, 3.3]: the UVa FB zero 2.694 and Li71's minimum 2.828 fm^-1 (t3_li6_charge_ff_fb.py, 2026-09-23)." — same line count required (the block is cited).

### Row 6 — RFY

* **R-1** `docs/references/REFERENCES.md` §1 #8 and #9: append "Harness ready and BLOCKED on this document: `validation/benchmarks/t3_li_magnetization_rfy.py` (2026-09-23)."
* **R-2** `docs/benchmarking/00_in_tree_checks.md` row C-3 / `06_critic.md` §4.6 "`HoSpin1FF` anchors F_m(0) on μ (`C-3`) and leaves the slope free" → "… and fixes the slope with the unfitted (q_z, b) = (1.30 fm⁻¹, 1.85 fm): r_mag(⁶Li) = 2.9469 fm (measured 2026-09-23), uncompared until RFY66 / ADNDT 14 is held".

## 9. Final tallies (base `d735abc` + this phase + the concurrent agents' files)

Three other agents worked in the same checkout at the same time (their files:
`tests/test_bench_spin.cpp`, `python/tests/test_bench_chain_rc.py`,
`python/tests/test_bench_spin_deuteron.py`, `include/lipolgen/coherent.hpp`,
`tests/test_tagged.cpp`, `validation/reference/*.json`, …), so the totals
below are the shared tree's; this phase's own share is stated beside each.

| gate | command | baseline (`239ac88`) | measured 2026-09-23, final | this phase's share |
|---|---|---|---|---|
| build | `cmake --build build -j` | — | exit 0 | header comment edit recompiled, no warning added |
| doctest | `build/lipolgen_tests` | 414 / 17 242 629 / 0 skipped | **416 cases / 17 242 712 assertions / 0 skipped / 0 failed** | 0 (no C++ test added; the +2 cases are `tests/test_bench_spin.cpp`) |
| pytest | `python -m pytest python/tests -q` | 1107 passed / 151 skipped | **1134 passed / 153 skipped / 0 failed** (287 s) | **+8 passed, +1 skipped** (`test_bench_nuclear.py`: 8 passed, 1 skipped = row 6 BLOCKED) |
| strict link gate | `python3 validation/check_physics_channels_links.py` | 1240 / 95 / 7 / 0 broken / 24 allow-listed | **1240 / 95 / 7 / 0 broken / 24 allow-listed** | 0 (the `cluster_config.hpp` edit holds its 733-line count) |
| ranges audit | `… --audit-ranges` | 0 REFUSED | **209 range citations, 0 REFUSED** | 0 |
| records | `… --records` | 0 broken | **764 citations in 46 dated run records, 0 broken** | this file adds a record and no `path:line` citation |
| SPDX | `python3 validation/check_spdx_headers.py` | 103/103 | **104/104** | 0 (the +1 is `tests/test_bench_spin.cpp`; `validation/benchmarks/*.py` sits outside the checker's `validation/*.py` glob, and every file of this phase carries the header anyway) |

`git diff include src` for this phase: `include/lipolgen/cluster_config.hpp`,
11 lines replaced by 11, every one a `///` comment line (checked:
`git diff include/lipolgen/cluster_config.hpp | grep '^[-+]' | grep -v '^[-+][-+]' | grep -v '^[-+] *///' | wc -l` → 0).
The baseline doctest had already moved to 17 242 631 assertions before this
phase's header edit (first build of the session, same 414 cases) — a
concurrent agent's, not this phase's.
