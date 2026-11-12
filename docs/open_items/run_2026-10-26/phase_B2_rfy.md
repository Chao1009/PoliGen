<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — phase B2: row 6, Rand–Frosch–Yearian, the measured lithium magnetic form factors

`BENCHMARK_PLAN.md` §4 row 6; `PLAN.md` phase B2 (rewritten after review);
`../../benchmarking/08_channel_priority.md` §4, ingredient **F**, and
`BENCHMARK_PLAN.md` §6, U-10. Base `2d89a64`. Owned and written by this
phase: `validation/benchmarks/t3_li_magnetization_rfy.py` (re-wired),
`validation/benchmarks/data/rfy1966_li_magnetic_ff.csv` (new),
`validation/benchmarks/data/rfy1966_li_magnetization.csv` (new, in a new
format), `python/tests/test_bench_rfy.py` (new), this file, and
`docs/benchmarking/figures/t3_li_magnetization_rfy.png` (written once by the
harness's `--plot`). **No physics default, constant, registry label,
tolerance of another row or reference file moved; nothing was committed.**
The paper stays outside the repository in the session scratch directory;
only numbers are vendored.

## 0. Verdict

| item | verdict |
|---|---|
| row 6 | **BLOCKED → PASS.** RFY Table I's seven measured ⁶Li points F_M²(q²), q² = 0.72–1.92 fm⁻² (q = 0.85–1.39 fm⁻¹), against the shipped `HoSpin1FF` shape F_mag² (q_z = 1.30 fm⁻¹, b = 1.85 fm, "UNFITTED STARTING VALUES"), with one free scale: **χ² = 13.351 / 6, ŝ = 1.419 ± 0.206, p = 0.0378**. The tolerance, p(χ², 6) ≥ 0.01, was fixed by `PLAN.md` phase B2 and written into a hashed design block before the harness computed anything (§2). |
| recorded | s = 1 (both sides normalised to 1 at q = 0, as the generator uses F_m): **χ² = 17.494 / 7, p = 0.0145**. With q² recomputed from E₀ at 180° instead of the printed q²: 13.248 / 6, ŝ = 1.423, p = 0.0393. The q² column check: max \|q²(E₀) − q²(printed)\| = 0.0058 fm⁻² over the 31 rows of Tables I and II (limit 0.006). |
| ⁷Li | **not compared** — `HoSpin1FF.for_ion(li7())` throws (no measured-moment block in `rc.hpp`). Table II's 24 points are vendored and drawn. |
| transcription closures | The paper's own fits re-evaluated at the vendored points with the vendored parameters: ⁶Li **8 / 8 hold**; ⁷Li **4 / 5** by the declared rule (Table VIII (C) a₀ = 1.66: 26.91 against 25.6 + 1) — a recorded FAIL of the closure, not of the row. A rounding-free refit added after it (declared post hoc, §4.2) holds **5 / 5**. |
| radius rule | **retired** (dated note in the harness header, §5). RFY print no magnetization radius; a₀ is vendored verbatim. ADNDT 14's lithium "magnetization radii" are √(5/2)·a₀, the 1p-orbital rms of RFY's fitted well. |

## 1. The document and what it prints

R. E. Rand, R. F. Frosch, M. R. Yearian, *Elastic Electron Scattering from
the Magnetic Multipole Distributions of Li⁶, Li⁷, Be⁹, B¹⁰, B¹¹, and N¹⁴*,
Phys. Rev. **144** (1966) 859–873, doi:10.1103/PhysRev.144.859
(`../../references/REFERENCES.md` §1 #8). The APS PDF: 15 pages (PDF page 1
is journal p. 859 and carries the end of the preceding article, so journal
page = PDF page + 858), 1 303 505 bytes, sha256
`7cbb2d21615d5e500a9c3235a026382e8529fc6c7828fa12a8fa6cfd31d5bbe2`, fetched
2026-10-26 by `curl` through the institution's licence (no client identity
spoofed), kept in the session scratch directory. The PDF was never opened by
the Read tool: pages were rendered with `pdftoppm` at 150, 300 and 450 dpi and
the PNGs read; `pdftotext -layout` only located the tables.

What it prints, read from the page images:

- **Tables I (⁶Li, 7 rows) and II (⁷Li, 24 rows), p. 865**: incident energy
  E₀ (MeV), "(Momentum transfer)² (F⁻²)", "Measured form factor squared" ±
  error. All at 180°. q is the **four-momentum** transfer (p. 860, under
  Eq. (2)); Eq. (9) evaluates it at θ = π with the recoil η = 1 + 2E₀/Mc².
  F_M² is normalised to 1 at q² = 0 through the measured dipole moment μ₀
  (Eq. (9); Table I: μ₀ = 0.822 nm; Table II: μ₀ = 3.256 nm). For ⁶Li
  (J = 1) only M1 contributes; for ⁷Li (J = 3/2) M3 adds and, per Table II's
  caption, the unresolved 0.478-MeV level is inside the measured number.
- **The model** (p. 868): Eq. (17), F_M² = [(⟨j₀⟩ + R₁⟨j₂⟩)² + R₂²⟨j₂⟩²]
  F_N², with the harmonic-well ⟨j₀⟩ = (1 − 2x/3)e⁻ˣ(c.m.), ⟨j₂⟩ = (2/3)x
  e⁻ˣ(c.m.), x = q²a₀²/4, (c.m.) = e^(x/A), F_N = (1 + q²/7.1)⁻¹, and
  "⟨r²⟩ = (5/2)a₀²"; Eq. (18) for ⁷Li with the inelastic level folded in,
  [A⟨j₀⟩² + B⟨j₀⟩⟨j₂⟩ + C⟨j₂⟩²] F_N². The α + d (1s) model (p. 869):
  F_M = F_d exp[−(q²/4q₀²)(1 − 2/A)], F_d a Gaussian of rms radius 1.8 F.
- **Tables VII (⁶Li, p. 869), VIII (⁷Li, p. 870), XIV (p. 873)**: a₀, R₁,
  (R₂, A, B, C, Ω₀/μ₀ for ⁷Li), χ² and degrees of freedom of every fit; the
  coupling-scheme predictions; Table XIV's a₀ from charge and from magnetic
  scattering (⁶Li: 1.78 and "1.6 −2.2"; ⁷Li: 1.72 and 1.70 ± 0.08) with the
  sentence that its errors are "the regions over which the respective
  chi-squared's are within 1 standard deviation of the expected values. (The
  smaller errors in Tables VII–XII are for fixed R₁, R₂.)"

**No magnetization radius is printed anywhere.** The OCR text of the scan
renders ± as `&`, `%`, `a`, `~` or `m` (e.g. "0.017a 0.007"); every number was
taken from the image, and at 300 and 450 dpi no digit is ambiguous. Both CSV
files carry the provenance, the convention column by column and the sha256
of the PDF; the harness re-checks the sha256 of every byte of each on every
read.

| file | rows | first | last | sha256 (every byte) |
|---|---|---|---|---|
| `rfy1966_li_magnetic_ff.csv` | 7 (⁶Li) + 24 (⁷Li) | ⁶Li I.1: 85 MeV, 0.72 fm⁻², 0.10 ± 0.03 | ⁷Li II.24: 200 MeV, 3.87 fm⁻², 0.016 ± 0.001 | `a3ada6b0…4e593aa8` |
| `rfy1966_li_magnetization.csv` | 13 (⁶Li) + 11 (⁷Li) | ⁶Li theory L-S, R₁ = 0 | ⁷Li XIV mag. a₀ = 1.70 ± 0.08 | `8f687426…5bed1e8` |

The second file's format is the task's (nucleus, fit, a0_fm, a0_err_fm,
coupling_note, q_range_fm, table, page), not the 2026-09-23 header's radius
format; `coupling_note` carries the printed parameters as `KEY=VALUE[+-ERR]`
pairs before a ` | ` and the paper's words after it, and `q_range_fm` is
labelled DERIVED (√ of the extreme printed q² of the table the fit used).

## 2. The design, fixed before the comparison

The comparison, the tolerance, the seven recorded items and what had been
seen were written into a `==== DESIGN BLOCK ====` of the harness header and
hashed (`d00b1307a0e397dcad2830513abdcbe21da9fb7eef3d1931e5976a7a0b3a3589`)
before any code computed a number. `design_block_sha256()` re-checks the hash
and reads the three numerical rules back out of the block's sentences (the
tolerance, the closure margin, the q² check limit) against the constants the
code uses, on every run; an edited block exits 2 (tested).

| UTC 2026-10-26 | event |
|---|---|
| 22:53:49 | `rfy1966_li_magnetic_ff.csv` written (Tables I, II) |
| 22:54:51 | `rfy1966_li_magnetization.csv` written (Tables VII, VIII, XIV) |
| **22:57:20** | **design block written and hashed** |
| 23:00:57 | first run of the harness |
| 23:02:36 | post-hoc addendum (a), (b) written outside the block, with the refit diagnostic's code |
| 23:02:41 | second run (the diagnostic); addendum (c), its outcome, appended after |
| 23:07:19 | `--plot`, the figure written once |

The design, in short (the header has the full text): the seven Table I
points against [F_m(t)/F_m(0)]² read through the shipped
`HoSpin1FF.for_ion(li6()).fm` at t = q²(ħc)² with the **printed** q² (checked
against the closed form to 1e−12; measured 5.6e−17); one free scale on the
model, closed form, unconstrained; diagonal, the printed errors, Gaussian;
all seven points (the three above 120 MeV, which the 0.7 % ⁷Li impurity
"almost entirely dominated", are in the paper's own fits); ndf = 6;
**PASS iff p(χ², 6) ≥ 0.01**. Recorded and never in the verdict: (r1) s = 1;
(r2) the pulls; (r3) q² from E₀ and the q² column check; (r4) ⁷Li not
compared; (r5) the transcription closures, rule χ²(table) ≤ χ²(paper) + 1;
(r6) a₀ verbatim with √(5/2)·a₀ labelled as the 1p-orbital rms; (r7) the tree's
slope radius, uncompared.

**What had been seen before the block was written**, and is said so in it:
every row of the five tables (the transcription); the committed report's
tree-side r_mag(⁶Li) = 2.9469 fm; a **by-hand** evaluation of F_mag² at the
seven points with a χ² estimate — about 17.5 for 7 points at s = 1 and about
13 for 6 dof with the free scale (p ≈ 0.04 by eye) — made while reading
Table I, after the tolerance had been fixed by the plan and the task; q² from
E₀ for six rows; √(5/2)·a₀ for the Table VII/VIII rows; and a slope radius
from Eq. (17). No rule was changed after any of it. The harness then measured
13.351 / 6 and 17.494 / 7.

## 3. The numbers (one harness run; window: RFY Table I, 7 points)

| row | E₀ (MeV) | q² printed (fm⁻²) | q (fm⁻¹) | q² from E₀ | F_M² ± σ | tree F_mag² | ŝ·F_mag² | pull |
|---|---|---|---|---|---|---|---|---|
| I.1 | 85 | 0.72 | 0.849 | 0.7203 | 0.100 ± 0.030 | 9.609e−02 | 1.363e−01 | −1.21 |
| I.2 | 100 | 0.99 | 0.995 | 0.9918 | 0.050 ± 0.010 | 3.153e−02 | 4.473e−02 | +0.53 |
| I.3 | 110 | 1.20 | 1.095 | 1.1960 | 0.017 ± 0.007 | 1.078e−02 | 1.530e−02 | +0.24 |
| I.4 | 120 | 1.42 | 1.192 | 1.4185 | 0.010 ± 0.003 | 2.247e−03 | 3.188e−03 | +2.27 |
| I.5 | 127 | 1.59 | 1.261 | 1.5850 | 0.005 ± 0.003 | 2.304e−04 | 3.270e−04 | +1.56 |
| I.6 | 130 | 1.66 | 1.288 | 1.6590 | 0.006 ± 0.003 | 1.840e−05 | 2.610e−05 | +1.99 |
| I.7 | 140 | 1.92 | 1.386 | 1.9175 | 0.001 ± 0.001 | 6.930e−04 | 9.833e−04 | +0.02 |

- **Verdict**: χ² = 13.3509 / 6, ŝ = 1.4188 ± 0.2058, p = 0.03779 ≥ 0.01 →
  **PASS**.
- **Where the χ² sits**: the three points just below the tree's zero
  (q² = 1.42, 1.59, 1.66 fm⁻²; the zero is at q_z² = 1.69 fm⁻²) carry 11.5 of
  the 13.35 (pulls +2.27, +1.56, +1.99): the data stay at 0.005–0.010 where
  the tree's shape falls to its zero. The points at q ≤ 1.1 fm⁻¹ and the last
  one are fine.
- **The scale**: ŝ = 1.42 ± 0.21 is 2.0 of its standard deviations above 1.
  It is the weighted mean of data / shape over the seven points, which is
  1.04 at the first point (q = 0.85 fm⁻¹), 1.59 and 1.58 at the next two and
  4.4 at q² = 1.42 fm⁻² (from the table above). At s = 1 (how the generator
  uses F_m, the normalisation fixed by μ) χ² = 17.494 / 7, p = 0.0145,
  recorded.
- **q² from E₀** (recorded): 13.2481 / 6, ŝ = 1.4227, p = 0.0393 — the
  abscissa is not what decides anything. The column check holds: the largest
  \|Δq²\| is 0.0058 fm⁻² (⁷Li II.17, 147 MeV) against the 0.006 limit;
  E₀ is printed as an integer, so a ±0.5 MeV rounding alone would allow
  ±0.0145 fm⁻² there — the limit was tighter than it needed to be and held
  anyway (header addendum (b)). The recoil is visible in the data: the
  100-MeV rows read 0.99 fm⁻² for ⁶Li and 1.00 for ⁷Li, as η predicts (tested).
- **Runtime**: `run()` 0.38–0.39 s; the script 0.48 s wall.

## 4. The transcription closures

### 4.1 The declared rule (r5): the paper's fits, no free parameter

| fit | χ²(table points) | χ²(paper, incl. other authors' points) | holds |
|---|---|---|---|
| ⁶Li VII (A) a₀ = 2.16, R₁ = 0.53 | 1.436 (7 pts) | 1.7 / 6 dof | yes |
| ⁶Li VII (C) a₀ = 1.70, R₁ = −0.02 | 8.430 | 8.4 / 7 | yes |
| ⁶Li VII (B) a₀ = 1.78, R₁ = 0.10 (Bernheim's a₀) | 6.044 | 6.0 / 7 | yes |
| ⁶Li VII (C) a₀ = 1.90, R₁ = 0.25 | 3.265 | 3.4 / 7 | yes |
| ⁶Li VII (C) a₀ = 2.00, R₁ = 0.36 | 2.017 | 2.3 / 7 | yes |
| ⁶Li VII (B) a₀ = 2.10, R₁ = 0.47 (Burleson–Hofstadter's) | 1.526 | 1.8 / 7 | yes |
| ⁶Li VII (C) a₀ = 2.20, R₁ = 0.57 | 1.465 | 1.7 / 7 | yes |
| ⁶Li α + d (1s), q₀ = 0.39 F⁻¹ | 2.441 | 2.8 / 7 | yes |
| ⁷Li VIII (A) a₀ = 1.70 | 21.814 (24 pts) | 23.3 / 22 | yes |
| ⁷Li VIII (C) a₀ = 1.62 | 29.374 | 30.8 / 23 | yes |
| ⁷Li VIII (C) a₀ = 1.66 | **26.909** | 25.6 / 23 | **no** (by 0.31 beyond the margin) |
| ⁷Li VIII (B) a₀ = 1.72 (Bernheim's) | 21.475 | 23.5 / 23 | yes |
| ⁷Li VIII (C) a₀ = 1.78 | 25.684 | 30.0 / 23 | yes |

The ⁶Li half closes on all eight fits, the α + d formula included as printed
(the exponent's (1 − 2/A) unsquared) — so the seven Table I numbers and this
reading of Eq. (17) reproduce the paper's own fit quality. The `closure[7Li
Table VIII]` sub-row is a **recorded FAIL**: one of five does not hold.

### 4.2 The post-hoc refit (declared after 4.1 was seen; never a verdict)

The page image was re-read at 450 dpi: Table VIII's (C) 1.66 row and the 24
points are as transcribed. The block's premise that rounding the printed
parameters moves χ² "by far less than 1" is not guaranteed when the printed
parameters are not the optimum on the table's own points, so a
rounding-free check was added, its rule written in the header addendum (a)
before it ran: Eq. (18) is linear in (A, B, C), so their weighted
least-squares minimum on the 24 points at the printed a₀ is closed-form, and
since the paper's χ² is a minimum over a superset of the points,
χ²_min(24) ≤ χ²(paper) must hold if the transcription, the reading of
Eq. (18) and the errors are the paper's (exactly for (B), (C); for (A) up to
the rounding of a₀, +0.25).

| ⁷Li fit | χ²_min(24) | χ²(paper) | holds | fitted − printed (A, B, C), in printed errors |
|---|---|---|---|---|
| (A) a₀ = 1.70 | 21.358 | 23.3 (+0.25) | yes | +1.12 / −0.35 / −0.23 |
| (C) 1.62 | 28.816 | 30.8 | yes | −1.03 / +0.41 / −0.02 |
| (C) 1.66 | 24.531 | 25.6 | yes | +0.50 / +0.02 / +1.41 |
| (B) 1.72 | 20.328 | 23.5 | yes | +1.85 / −0.51 / −0.33 |
| (C) 1.78 | 20.195 | 30.0 | yes | +3.78 / −0.99 / −0.82 |

All five hold. At a₀ = 1.66 the 24-point optimum lies 1.41 printed errors
from the printed C — farther than rounding (≤ 0.25 of an error) — which is
what the paper's fit including the two other-author points of Fig. 6
(26 points, 23 dof) would do, and not what a transcription error would do.
The 4.1 FAIL stays recorded as it came out.

## 5. The radius rule, retired — and what ADNDT 14's "magnetization radii" are

The 2026-09-23 header waited for "r_mag rows" from RFY or ADNDT 14 Table V
and planned "the tree's r_mag inside the envelope of all the vendored ⁶Li
rows, each widened by its own quoted error". Retired (dated note in the
harness header) for two reasons:

1. **RFY print no magnetization radius.** The only radius in their model is
   p. 868's ⟨r²⟩ = (5/2)a₀², the rms radius of the 1p-shell orbital of the
   fitted well, which the paper keeps distinct from the slope of F_M at
   q² = 0 (p. 873: "the almost horizontal form factor curve at q² = 0 does not
   imply ⟨r²⟩ = 0, as may be seen by inspection of (17)"). A slope radius of
   Eq. (17) needs R₁ as well as a₀, the c.m. factor and F_N, and a₀ itself is
   known only to 1.6–2.2 F for ⁶Li (p. 869, Table XIV).
2. **The pre-registration condition was not met.** The task allowed a
   slope-radius sub-row only if its formula, the centre-of-mass and
   nucleon-size choices and the RFY error entering it were written in the
   header before anything was derived. A slope radius was evaluated by hand
   before the header existed (§2), so no radius is compared, and no radius is
   recorded as a comparison.

a₀ is vendored verbatim instead. The harness records √(5/2)·a₀ for each a₀
row, labelled as the 1p-orbital rms and never compared, and the values
identify the numbers `../../benchmarking/06_critic.md` §4.6 transcribed from a
third-party copy of ADNDT 14 (1974) Table V (not vendored, not used):

| ADNDT 14 row as 06 §4.6 quotes it | RFY a₀ | √(5/2)·a₀ (harness) |
|---|---|---|
| ⁶Li 3.42(6) fm, "R = 2.16(4)" | Table VII (A) 2.16 ± 0.04 | 3.415 ± 0.063 |
| ⁶Li 2.81 fm, "a₀ fixed from charge scattering" | Table VII (B) 1.78 | 2.814 |
| ⁷Li 2.69(13) fm, "R = 1.70(8)" | Table XIV mag. 1.70 ± 0.08 | 2.688 ± 0.126 |
| ⁷Li 2.72 fm, "a₀ fixed" | Table VIII (B) 1.72 | 2.720 |

All four agree at the printed digit, errors included (the ⁷Li ±0.13 is
Table XIV's ±0.08, not Table VIII's ±0.01). They are **orbital radii of the
fitted well, not magnetization slope radii**; over RFY's own ⁶Li range
a₀ = 1.6–2.2 F they span 2.530–3.479 fm, which brackets the charge rms
(2.54–2.57 fm, de Vries) rather than lying above it. The tree's own slope
radius, 6/q_z² + (3/2)b² → r_mag(⁶Li) = 2.9469 fm, stays a recorded tree-side
number (`config.tree_r_mag_li6_fm_uncompared`), as it was while BLOCKED.

## 6. What the verdict means physically

**For the RC tail's F_m (ingredient F).** LiPolGen puts F_m into POLRAD's
spin-1 elastic tail in both of its t-peak integrands: the unpolarised one
carries (2/3)η F_m² X_t − (2/3)(1 + η) F_m² beside F_c² X_t
(`src/core/rc.cpp:268-284`, `polrad_sigma_el_u`), and in the tensor one
(1 + η + (¾x_A² − η)X_t) F_m² is the only term that does not carry F_q
(`src/core/rc.cpp:286-305`, `polrad_sigma_el_t`). So the magnetic shape is
an input of the tensor elastic tail, i.e. of the A_zz RC band (rank 4;
`../../benchmarking/08_channel_priority.md` §4, ingredient F). Until today
that shape, `include/lipolgen/rc.hpp:651-665`
(`LI6_FF_FM_QZ_FM` = 1.30 fm⁻¹, `LI6_FF_FM_B_FM` = 1.85 fm, "UNFITTED STARTING
VALUES", q_z the midpoint of a design range and b a geometric mean of two
single-landmark solutions), had no lithium datum behind it. Now: over
q = 0.85–1.39 fm⁻¹ — the first lobe and the zero — the only measured ⁶Li
magnetic form factor held does not exclude it at the declared 1 % level,
p = 0.038 with a free normalisation and p = 0.0145 at the moment-fixed one.
The fit is not comfortable: the χ² is carried by the three points just below
q_z, where the data show F_M² = 0.005–0.010 and the tree's shape goes to
zero, and the scale wants +42 ± 21 %. A PASS here
says the unfitted (q_z, b) are not ruled out by seven 1966 points with
30–100 % errors; it does not say they are fitted or right, and nothing about
the second lobe near q = 2 fm⁻¹, which `rc.hpp`'s comment on `F_mag` calls
"exactly where the tail lives". No knob was run, so no observable is priced
by this row.

**For the α + d clustering signature.** RFY's own words are narrower than
`BENCHMARK_PLAN.md` row 6's "the only measured electromagnetic signature of
α+d clustering". The abstract: "The outstanding exception is Li⁷, whose form
factor is accurately described by the extreme independent-particle model.
(L-S coupling is quite incompatible with the data.) The Li⁶ results may also
be interpreted in terms of an 'alpha plus deuteron' model." For ⁷Li the
4-parameter fit's A, B, C = 1.24, 0.23, 1.18 sit beside the
independent-particle 1.29, 0.28, 1.35 and far from L-S's 1.57, −1.33, 2.73
(Table VIII; p. 871 "extremely close to the predictions of the independent
particle model"). For ⁶Li the paper does **not** say the independent-particle
model fails — Table VII lists only L-S (R₁ = 0) and j-j (R₁ = 0.438)
predictions. Its free fit gives R₁ = 0.53 ± 0.03, near j-j, which it calls
"unreasonable" given Lane and Radicati's L-S preference; with Bernheim's
charge-scattering a₀ = 1.78 F the fit gives R₁ = 0.10 ± 0.03, "in closer
agreement with the L-S prediction"; and the α + d (1s) model fits with
χ² = 2.8 / 7. Then (p. 869): "There is no way to choose between this model
(1s) and the L-S-coupled shell model with the present data." This harness
reproduces both: the α + d (1s) closure at 2.441 and the Bernheim-a₀ L-S-like
fit at 6.044 on the seven table points. So the RFY ⁶Li data are **compatible
with** α + d clustering and **do not discriminate** it from the L-S shell
model; and the "magnetization rms larger than the charge rms" argument of
06 §4.6 rests on √(5/2)·a₀ (§5), which over the allowed a₀ range brackets the
charge rms. For the channels that ride on the α–d picture — tagged ⁶Li α
(rank 6), the α–d transfer of b₁(⁶Li) (ingredient D, U-1) — this row neither
supports nor refutes it.

**For ⁷Li** the tree has no magnetic form factor at all
(`HoSpin1FF::for_ion` refuses the ion, and the spin-1 line of POLRAD
Eq. (A.4) it implements is not ⁷Li's spin 3/2). The 24 Table II points are
vendored for the day it has one; note that they include the 0.478-MeV level
(about 25 % of the scattering, p. 871), so an elastic-only ⁷Li F_m would be
compared to Eq. (18)'s elastic part, not to the table.

## 7. What this does NOT validate

F_m outside q = 0.85–1.39 fm⁻¹, in particular the second lobe near 2 fm⁻¹
(Bergstrom, Kowalski, Neuhausen, PRC 25 (1982) 1156, REFERENCES.md #22,
reaches higher and is not held); F_m(0), the measured moment by construction;
the q → 0 slope (the data start at 0.85 fm⁻¹ and the free scale absorbs a
normalisation); C0 and C2 (row 5, the de Vries radius); ⁷Li; any radiative
correction itself (row 7, the POLRAD rows); and whether the α + d picture is
right (above). The points of other authors in Figs. 5–6 (Peterson;
Goldemberg and Torizuka) are not vendored. The errors are taken as printed,
uncorrelated and Gaussian; the paper does not state their composition, and
the three ⁶Li points above 120 MeV carry a large ⁷Li-impurity subtraction.

## 8. Files touched, commands, runtimes

| file | change |
|---|---|
| `validation/benchmarks/t3_li_magnetization_rfy.py` | re-wired: the design block, the two sha256-checked readers, the comparison, the recorded items, the closures and the post-hoc refit, `plot()`, `--plot`; `tree_side()` keeps its 2026-09-23 keys |
| `validation/benchmarks/data/rfy1966_li_magnetic_ff.csv` | new: Tables I, II (31 rows) |
| `validation/benchmarks/data/rfy1966_li_magnetization.csv` | new: Tables VII, VIII, XIV (24 rows), the task's format |
| `python/tests/test_bench_rfy.py` | new: 16 tests |
| `docs/benchmarking/figures/t3_li_magnetization_rfy.png` | new: written once by `--plot` |
| this file | new |

```
cd /home/cpeng/Projects/polli/PoliGen && source env.sh
python3 validation/benchmarks/t3_li_magnetization_rfy.py            # exit 0, PASS, 0.48 s wall
python3 validation/benchmarks/t3_li_magnetization_rfy.py --plot     # the figure, once
python3 -m pytest -p no:cacheprovider python/tests/test_bench_rfy.py \
    "python/tests/test_bench_nuclear.py::test_rfy_tree_side_closed_form" \
    "python/tests/test_bench_nuclear.py::test_rfy_magnetization_benchmark" -q
                                                                    # 18 passed in 1.7 s
python3 validation/benchmarks/make_report.py --only t3_li_magnetization_rfy \
    --out-dir <scratchpad>/report_B2                                # 1 pass; the row renders, both data files hashed
python3 validation/check_spdx_headers.py                            # 131 files, all carry the line
```

The tests: the two file guards (counts, first and last values, one flipped
byte in the body and one in the provenance refused); the note parser; the
design block's hash, its stated tolerances and a moved constant refused; the
verdict and every recorded sub-row pinned (the ⁷Li closure as a FAIL); the
model functions and the 180° kinematics in closed form; the tree side and
rule 3 (F_m/F_m(0) unchanged by `tail_tensor_scale` and M_A); the verdict
following from the numbers (ŝ minimises χ², the scipy-free p agrees with
scipy); the a₀ record against 06 §4.6's four numbers; the REPORT and SUBROW
lines; exit 0 on a run and exit 2 on an altered data file of either kind and
on an edited design block; `plot()` writing a PNG.

## 9. For the integration stage, and open items

The two `test_rfy_*` tests of `python/tests/test_bench_nuclear.py` still pass
against the re-wired harness (`tree_side()` keeps its keys; the BLOCKED skip
can no longer fire) and are superseded by `test_bench_rfy.py`; retiring them,
the README / `BENCHMARK_PLAN.md` / `REFERENCES.md` rows, the stale sentences of
`00_in_tree_checks.md` C-3 and 06 §4.6, and the report regeneration are
listed for the integrator with exact text in this phase's return.

Open: (1) whether the verdict should be the s = 1 comparison, the way the
generator uses F_m (it is recorded: p = 0.0145, also ≥ 0.01) — the plan chose
the shape; (2) Bergstrom–Kowalski–Neuhausen 1982 (REFERENCES.md #22), the
⁶Li magnetic form factor to higher q, is the document that would test the
second lobe and is probably reachable the same way; (3) ADNDT 14 (#9) is no
longer needed for row 6; (4) the licence question for vendoring the numbers
and the filing of the PDF in the author's corpus, as for BONuS; (5) a
candidate for `DECISIONS_PRICED.md`, none taken and nothing priced: whether
F_mag's (q_z, b), which `rc.hpp` assigns to a fit against Wiringa–Schiavilla's
F_T landmarks, should also answer to RFY Table I.
