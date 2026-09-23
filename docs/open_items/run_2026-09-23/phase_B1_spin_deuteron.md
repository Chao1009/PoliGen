# Run 2026-09-23 — Phase B1: benchmark harness group SPIN AND DEUTERON

*BENCHMARK_PLAN.md §4 rows 1, 3, 9, 10. Base commit `d735abc` (tree clean at
start). Written 2026-09-23. Nothing in this phase moved a default, a registry
row, a constant or a reference JSON; every disagreement below is RECORDED with
its distance and priced, never tuned. Every number carries the command that
printed it; all commands run from the repository root after
`source env.sh`.*

## 0. Status

| row | harness | status | the claim, as measured |
|---|---|---|---|
| 1 | `validation/benchmarks/t5_epios_source_modes.py` | **pass** | all eight EPIOS Table II (P_z, P_zz) modes round-trip through `spin1_populations` with residual **0** (tolerance 1e-12), populations are the exact rationals, non-negative; the EST ladder reproduces the ideal P_zz at modes 0, 5, 6 only (5 and 6 in the \|P_z\| = 1 limit, where the tree throws) |
| 3 | `validation/benchmarks/t2_hermes_b1_table2.py` | **fail** (a measurement) | no configuration puts all six bins inside the quadrature error; χ²/6 on b₁ᵈ: 22.26 (convolution, MSTW), 22.67 (convolution, ToyF2), **5.26 (shipped Miller, ×0.5)**, **5.30 (Miller ×1)**, 21.78 (digitized CDKS), 21.80 (b₁ = 0) |
| 9 | `validation/benchmarks/t5_est_identity.py` | **pass** at 2.31e-14 | max \|`spin_temperature_pzz(1, P_z)` − (2 − √(4 − 3P_z²))\| = **1.138e-14** over 19 999 interior points; **the requested 1e-14 does NOT hold** (496 points exceed it); the bisection bound is 2.13e-14. COMPASS ⁶LiD numbers reproduced at the printed digit |
| 10 | `validation/benchmarks/t2_bonus_spectator_shape.py` | **blocked** | the CLAS Physics Database holds BONuS only as F2n/F2p (one measurement); the BONuS p_s tables are in PRC 89 045206's Supplemental Material, HTTP 401 without an APS login |

Tests: `python/tests/test_bench_spin_deuteron.py` (16 pass, 1 skip — row 10,
skipping with the name of the document to download) and
`tests/test_bench_spin.cpp` (2 doctest cases, 81 assertions: rows 1 and 9 on
the C++ side, the EPIOS modes read from the vendored JSON). No row costs more
than 30 s (the slowest, row 3, 4.0–4.4 s over the 2026-09-23 runs), so none is behind `LIPOLGEN_BENCH=1`.

## 1. Row 1 — EPIOS Table II vs the spin bookkeeping

**Source.** EPIOS, arXiv:2510.10794 = PRC 113 (2026) 060501, Table II (printed
p. 18), transcribed from `PolarizedLithiumSim/refs/2510.10794.pdf`
(md5 012cbaa4…) into
`validation/benchmarks/data/epios2026_table2_deuteron_source_modes.json`
with provenance; original data ref. [81] = Chiladze *et al.*, PRST-AB 9
(2006) 050101 (COSY/ANKE, LEP, Nov. 2003). Convention: Cartesian, per ion.
**The P_zz column is ideal**; the caption says "measured vector and tensor
polarizations" and "relative beam intensities" but the table prints neither a
measured tensor column nor an intensity column — only P_z^LEP is measured.

**Measured** (`python3 validation/benchmarks/t5_epios_source_modes.py`):

| mode | (P_z, P_zz) ideal | (p₊, p₀, p₋) | round-trip | EST P_zz at P_z | EST − ideal | EST at P_z^LEP | P_z^LEP/P_z |
|---|---|---|---|---|---|---|---|
| 0 | (0, 0) | (1/3, 1/3, 1/3) | 0 | 0 | 0 | 0 | — |
| 1 | (−2/3, 0) | (0, 1/3, 2/3) | 0 | 0.36701 | **+0.36701** | 0.26263 | 0.858 |
| 2 | (+1/3, +1) | (2/3, 0, 1/3) | 0 | 0.08515 | **−0.91485** | 0.06188 | 0.855 |
| 3 | (−1/3, −1) | (0, 2/3, 1/3) | 0 | 0.08515 | **+1.08515** | 0.06961 | 0.906 |
| 4 | (+1/2, −1/2) | (1/2, 1/2, 0) | 0 | 0.19722 | **+0.69722** | 0.12066 | 0.790 |
| 5 | (−1, +1) | (0, 0, 1) | 0 | tree throws; limit 1 | 0 | 0.49126 | 0.758 |
| 6 | (+1, +1) | (1, 0, 0) | 0 | tree throws; limit 1 | 0 | 0.45180 | 0.731 |
| 7 | (−1/2, −1/2) | (0, 1/2, 1/2) | 0 | 0.19722 | **+0.69722** | 0.13497 | 0.834 |

Round-trip = max over (P_z, P_zz) of the residual by both routes
(`moments_along_axis`, and `rho_from_populations` →
`vector_polarization` / `tensor_polarization`); it is exactly 0 for all
eight. Modes 1–7 sit on the domain edge (min p_m = 0); mode 0 is interior.

**The ladder statement, as measured.** The equal-spin-temperature fill (the
`--pzz-mode ladder` default of `helicity_flip_plan`) reproduces the ideal P_zz
of **3 of 8** source modes, and two of those three are the \|P_z\| = 1 pure
states, where every fill coincides and where the tree's
`spin_temperature_pzz` throws ("maxent populations need |pz| < 1"). At the
five partially polarized modes the ladder misses the source's P_zz by
+0.367, −0.915, +1.085, +0.697, +0.697 (in P_zz ∈ [−2, 1]). The ladder's
P_zz is ≥ 0 at every P_z, so it cannot produce modes 3, 4, 7 (P_zz < 0) at
any P_z. What this measures: the atomic-beam source is not an equal-spin-
temperature system. What it does **not** measure: which fill an EIC lithium
beam will have — no lithium source mode table exists (EPIOS §V D: the ⁶Li/⁷Li
source is under development). The default is unchanged.

The measured/ideal vector ratio is 0.731–0.906 (mode 6 … mode 3), from the
table's own P_z^LEP column.

## 2. Row 9 — the equal-spin-temperature identity

**Measured** (`python3 validation/benchmarks/t5_est_identity.py`;
`build/lipolgen_tests -tc="bench row 9*" -s`):

* J = 1, grid P_z = −1 + 2i/20000, i = 1 … 19 999: max \|diff\| =
  **1.1380e-14** at P_z = −0.5996; **496 of 19 999 points exceed 1e-14**. The
  plan row's "(8e-15)" and 06_critic.md §4.1's table are five points (max
  7.7e-15), not a bound. The residual is `populations_maxent`'s bisection
  (bracket < 1e-13 in β): bound 0.5 × 1e-13 × max\|dP_zz/dβ\| = 0.5e-13 ×
  0.42609 = 2.1305e-14 (maximum of 6 sinh β/(2 cosh β + 1)² at β = 1.113,
  P_z = 0.621). The closed form's own double rounding is ≤ 3.2e-16 (40-digit
  check, scratch). Asserted tolerance: 2.3065e-14. **At P_z = ±1 the tree
  throws**; the closed-form limit is 1.
* J = 3/2 (⁷Li): the ladder is geometric (ratio rel. err ≤ 2.50e-16), obeys
  P_z = B_{3/2}(3/2 · ln t) to 3.89e-16, and its rank-2 moment matches the
  partition-function closed form to 1.10e-14.
* **COMPASS ⁶LiD, numerically.** Koivuniemi *et al.*, SPIN 2004, pp. 796–798,
  fetched 2026-09-23 from
  `wwwcompass.cern.ch/compass/publications/technical/koivuniemi_spin2004.pdf`
  (240 451 bytes, md5 52a3c608…; numbers in
  `validation/benchmarks/data/est_identity_sources.json`, the PDF not
  vendored). From the tree's ladder at their P_d = 0.530 and a common spin
  temperature (β_i ∝ μ_i/J_i): populations (63.55, 25.90, 10.55) % vs printed
  (63.6, 25.9, 10.6) % (max 0.048 %); T_S = 0.8756 mK vs 0.87; P(⁶Li) =
  0.5127 vs 0.513; P(⁷Li) = 0.9238 vs 0.924. Asserted at the printed digit.
* **A convention trap found.** Their Fig. 1 caption defines
  **T = ½(1 − 3p₀) = 11 %**: HALF the Cartesian P_zz = 1 − 3p₀ this tree uses.
  The tree's EST P_zz at 0.530 is 0.22312, i.e. T = 0.11156. Anyone quoting
  "11 % tensor polarization" of the COMPASS deuteron into this tree is off by 2.
* **A transcription error outside the tree.** `PolarizedLithiumSim/refs/
  refs_dict.json`, entry `koivuniemi2004-compass-polarization-buildup`, quotes
  "51 ± 3 %" and "92 ± 4 %"; the rendered caption reads "±51.3 %" and
  "±92.4 %" (the ± is the sign of the polarization). Not edited (read-only
  repository); listed for its owner.

**The "never measured" line, with its evidence.** What this phase can
support is a finding of absence, stated with its searches: in the COMPASS ⁶LiD
target the unsplit 3 kHz NMR line measures p₊ − p₋ only, and the tensor
polarization was *estimated* from EST (Koivuniemi, Fig. 1, read today); the
INSPIRE full-text search of 06_critic.md §4.2 returned no measurement; and
EPIOS §VIII B 4 lists a "Li-6/Li-7 beam polarimeter" only as a *proposed* R&D
item (an atomic-beam target, Li–Li CNI scattering, Breit–Rabi calibration) —
no lithium beam polarimeter of any kind exists, and the section does not
address tensor polarization. The universal "has never been measured by
anyone" is stronger than that; the replacement text in §6 (E4) states the
searches.

**A misattribution found while checking it.** 06_critic.md §4.1 quotes EPIOS
"tensor polarization presents additional challenges, requiring specialized
polarimeter development such as Lamb-shift-based systems or BRP-type" as
evidence that "there is not yet a tensor polarimeter for a lithium beam".
Read in context (`pdftotext -layout` of the local PDF), that sentence is in
§VIII A 2, **"Polarized deuteron ion source"**, and ends "…polarimeters
specifically designed for deuteron ions". It is about deuterons. Correction
E17.

## 3. Row 3 — HERMES b₁ᵈ Table II

**Source.** HERMES, PRL 95 (2005) 242001, Table II, transcribed from
`PolarizedLithiumSim/refs/hep-ex_0506018.pdf` into
`validation/benchmarks/data/hermes2005_b1d_table2.json` (identical to the
refs_dict entry). R1990 from `PolarizedLithiumSim/refs/whitlow1990_R.pdf`
(SLAC-PUB-5284 scan, Eqs. (6)–(7)) into
`validation/benchmarks/data/whitlow1990_r1990.json`. **The scan prints
b₁ = 0.635**; the value used is 0.0635 (R(0.3, 3) = 0.187 against r1998's
0.193; 0.635 would give 0.49). Guarded by `test_r1990_typo_guard`.

**Normalisations (rule 1).** HERMES: per nucleon (F₂ᵈ = (F₂ᵖ + F₂ⁿ)/2).
`DeuteronConvolutionB1`, `CdksB1`: per nucleon, the 0.5 dropped
(`B1_CDKS_TABLE_TO_PER_NUCLEON` = 1). `MillerB1`: Miller's table per deuteron,
0.5 kept (`B1_MILLER_TABLE_TO_PER_NUCLEON` = 0.5, registry row 2). No A-scaling
(A = 2 both sides; the 2/6 inert-α factor is not applied). x on the nucleon
mass both sides.

**Measured** (`python3 validation/benchmarks/t2_hermes_b1_table2.py`), b₁ᵈ per
nucleon in 10⁻²:

| x | Q² | HERMES ± quad. | conv. MSTW | conv. ToyF2 | Miller ×0.5 (shipped) | Miller ×1 | CDKS digitized |
|---|---|---|---|---|---|---|---|
| 0.012 | 0.51 | +11.20 ± 6.17 | +0.0435 | −0.0300 | 5.7147 | 11.4293 | 0.2392 |
| 0.032 | 1.06 | +5.50 ± 3.13 | +0.0046 | −0.0282 | 2.8805 | 5.7609 | 0.0459 |
| 0.063 | 1.65 | +3.82 ± 1.26 | −0.0089 | −0.0203 | 2.0920 | 4.1841 | 0.0024 |
| 0.128 | 2.33 | +0.29 ± 0.69 | −0.0285 | −0.0203 | 0.8896 | 1.7792 | −0.0259 |
| 0.248 | 3.11 | +0.29 ± 0.37 | −0.0592 | −0.0118 | 0.2741 | 0.5482 | −0.0530 |
| 0.452 | 4.69 | −0.38 ± 0.16 | +0.0048 | +0.0219 | −0.2066 | −0.4132 | −0.0033 |

A_zz = −(2/3) b₁/F₁ᵈ with F₁ᵈ rebuilt HERMES's way from MstwSF + R1990
(chi² with R1998 in brackets):

| configuration | χ²(b₁)/6 | bins within | χ²(A_zz)/6 | bins within |
|---|---|---|---|---|
| `DeuteronConvolutionB1` + `MstwSF` (the gate configuration) | 22.264 | 2/6 | 22.009 (22.008) | 2/6 |
| `DeuteronConvolutionB1()` shipped Options (ToyF2) | 22.669 | 2/6 | 22.392 (22.393) | 2/6 |
| `MillerB1` shipped (×0.5) | **5.262** | 4/6 | **5.550** (5.410) | 4/6 |
| Miller ×1 (row 2's other option) | **5.297** | 5/6 | **4.514** (4.587) | 5/6 |
| `CdksB1` digitized (Q² = 2.5) | 21.782 | 2/6 | 21.559 (21.548) | 2/6 |
| b₁ = 0 (baseline, not a model) | 21.799 | 2/6 | 21.586 (21.586) | 2/6 |

**What the numbers say, and no more.**
* The standard convolution (the A = 2 gate's configuration, which passes its
  CDKS Fig. 4 gate) is **indistinguishable from b₁ = 0** against HERMES
  (22.26 vs 21.80): it is two orders of magnitude below the data at
  x < 0.1, which is CDKS's own conclusion. The gate tests the kernel against
  CDKS, not against data.
* Miller's curve beats zero by Δχ² ≈ 16.5 — but Miller tuned P₆q to one
  HERMES point (AUTHOR_DECISIONS §B2(c)), so this is not an independent test.
* **Registry row 2 is not decided by HERMES.** Δχ²(b₁) = +0.035 and
  Δχ²(A_zz) = −1.036 for dropping the 0.5 (R1998: −0.823). At §B2(c)'s own
  test point x = 0.012 the shipped value is 5.71 (−0.89σ) and ×1 is 11.43
  (+0.04σ); at x = 0.452 it is the other way (−0.207, +1.07σ, against −0.413,
  −0.20σ). Six points with these errors cannot tell a factor 2.
* **F₁ price.** F₁ᵈ rebuilt with MSTW + R1990 against the F₁ᵈ HERMES's own
  table implies (−(2/3) b₁/A_zz, 3-digit entries): ratio 1.0256–1.0727 per
  bin (ALLM97/NMC → MSTW LO plus rounding). R1990 vs R1998 moves F₁ by
  (1 + R90)/(1 + R98) = 0.9143 in bin 1 (Q² = 0.51) and by ≤ 0.7 % in bins 2–6.
* Caveats stated in the harness: HERMES's ~1e-3 A_zz normalisation is
  correlated over bins but treated as uncorrelated; bin 1 is below MSTW2008's
  1 GeV² grid floor (PYTHIA's reader extrapolates) and below R1990's SLAC x
  range; the comparison is Born-level against RC-unfolded data.

## 4. Row 10 — the CLAS Physics Database

**Query** (read-only, 2026-09-23): `https://clas.sinp.msu.ru/cgi-bin/jlab/db.cgi`
with `tgt=deuteron`, `quantity=0` (any), `limit=0` → HTTP 200, 320 804 bytes,
**918 measurements in 19 experiment ids**; also run per quantity (56, 57, 58,
93, 94, 173, 13). Catalogue vendored as
`validation/benchmarks/data/clas_db_deuteron_query_2026-09-23.json`.

* **BONuS** (E-03-012, eid 135): **one** measurement, E135M1, the ratio
  F₂ⁿ/F₂ᵖ vs W* (Baillie 2012; 0.65 < Q² < 4.52, 1.1 < W < 2.7). No spectator
  spectrum; no Tkachenko 2014 entry.
* **EG1b** (eid 95: 13 measurements; eid 146: 17): A₁, g₁, g₁/F₁, Γ₁, … —
  inclusive, no spectator.
* **Deeps** (eid 90; Klimenko *et al.*, PRC 73 (2006) 035212,
  nucl-ex/0510032): 251 measurements, of which **115 are F₂ᴺ × P(p_s,
  cos θ_pq)** at p_s = 0.30–0.53 GeV/c, Q² = 1.8 and 2.8 GeV². The only
  deuteron spectator-momentum distribution in the database. Vendored (all 115,
  parsed from the measurement pages) as
  `validation/benchmarks/data/clas_deeps_klimenko2006_f2P.json` — **not
  wired**: it is the high-p_s tail where the paper finds PWIA adequate only at
  cos θ_pq < −0.3. Two label discrepancies kept verbatim: the database's last
  p_s bin is "560 MeV" where the paper's average is 0.53 GeV/c; 55 blocks
  carry database Q² labels 0.18–0.65 that are not the paper's Q² bins.
* also present: eid 20 (e6a, dσ/(dQ² dP_n [dΩ_n]), Egiyan/Griffioen/Strikman).

**The BONuS tables themselves.** Both papers are on disk
(`refs/1110.2770.pdf`, `refs/1402.2477.pdf`); neither body prints a p_s
table. Tkachenko *et al.* ref. [71] puts "tables of numerical results" (the
cos θ_pq spectra, 6 W* × 5 Q² × 4 p_s bins, 70–150 MeV/c, as R_D/S = data over
their own simulation) in the Supplemental Material;
`https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.89.045206`
returned **HTTP 401** (APS login); arXiv:1402.2477 has no ancillary files.
The row is **blocked**; the harness exists, and its test skips naming that
document.

## 5. Files (all new; nothing existing was edited)

- `validation/benchmarks/README.md` — conventions and the `LIPOLGEN_BENCH=1` opt-in convention (created first; no README existed; no test reads the variable yet, since no row needs it)
- `validation/benchmarks/t5_epios_source_modes.py`, `t5_est_identity.py`, `t2_hermes_b1_table2.py`, `t2_bonus_spectator_shape.py`
- `validation/benchmarks/data/epios2026_table2_deuteron_source_modes.json`, `est_identity_sources.json`, `hermes2005_b1d_table2.json`, `whitlow1990_r1990.json`, `clas_db_deuteron_query_2026-09-23.json`, `clas_deeps_klimenko2006_f2P.json`
- `python/tests/test_bench_spin_deuteron.py`
- `tests/test_bench_spin.cpp` (picked up by the existing `file(GLOB … tests/test_*.cpp)`; **CMakeLists.txt not edited**)
- this record

## 6. Doc edits requested (for the Documents stage — none applied here)

Format: **file — site — exact replacement text.**

**E1.** `docs/benchmarking/BENCHMARK_PLAN.md` — §4 table, row 1, last cell
("the FIRST external anchor … the wrong model for an EIC ion beam") — replace with:
> the FIRST external anchor on the spin bookkeeping (today only a polligen self-pin). **WIRED 2026-09-23** (`validation/benchmarks/t5_epios_source_modes.py`, pass): all eight modes round-trip with residual 0 at tolerance 1e-12. Measured: the `ladder` (EST) fill reproduces the source's ideal P_zz at 3 of 8 modes (0, and the two \|P_z\| = 1 pure states where the tree throws) and misses the five partially polarized ones by +0.37 to +1.09 in P_zz; EST cannot produce P_zz < 0 at all. An atomic-beam source is not an EST system; which fill an EIC lithium beam has is not measured (no Li source table exists)

**E2.** same file — §4 row 3, last cell — replace with:
> converts the only measurement of the target observable from a comment into a gate; forces registry row 2 against data. **WIRED 2026-09-23** (`validation/benchmarks/t2_hermes_b1_table2.py`, fail — a measurement, no tuning): χ²/6 on b₁ᵈ = 22.26 (A = 2 convolution, MSTW; indistinguishable from b₁ = 0 at 21.80), 5.26 (shipped Miller ×0.5), 5.30 (Miller ×1). HERMES does not decide registry row 2 (`../open_items/run_2026-09-23/phase_B1_spin_deuteron.md` §3)

**E3.** same file — §4 row 9, last cell ("`spin_temperature_pzz` IS the polarized-target closed form (8e-15); also the honest limitation: P_zz of ⁶Li has never been measured by anyone") — replace with:
> `spin_temperature_pzz` IS the polarized-target closed form: max deviation 1.14e-14 over 19 999 points (8e-15 was five points; the bisection bound is 2.13e-14), **WIRED 2026-09-23** (`validation/benchmarks/t5_est_identity.py`, pass at 2.31e-14), with the COMPASS ⁶LiD populations, T_S and P(⁶Li), P(⁷Li) reproduced at the printed digit. The honest limitation: no measurement of the tensor polarization of ⁶Li was found by any search this project ran (see 06_critic.md §4.1 N-2)

**E4.** `docs/benchmarking/06_critic.md` — §4.1 N-2, the sentences "So in the one material where ⁶Li has ever been polarized, the NMR line carries the **vector** polarization only and **P_zz was never measured**. The tensor polarization of ⁶Li has never been measured by anybody, by any method." — replace with:
> So in the COMPASS ⁶LiD target the NMR line carries the **vector** polarization only, and the tensor polarization was **estimated** from the EST relation, not measured — Fig. 1's caption: "The tensor polarization can be estimated to be T = 1/2(1 − 3p₀) = 11 %" (for the deuteron; note their T is **half** the Cartesian P_zz = 1 − 3p₀ of this tree). No measurement of the tensor polarization of ⁶Li was found by any search run for this survey: the INSPIRE full-text search of §4.2 returned none, and EPIOS §VIII B 4 lists a ⁶Li/⁷Li beam polarimeter only as proposed R&D (no lithium beam polarimeter exists; tensor polarimetry is not addressed there). (Re-read and checked numerically 2026-09-23, `../open_items/run_2026-09-23/phase_B1_spin_deuteron.md` §2.)

**E5.** same file — §4.1 N-2, "`spin_temperature_pzz(1, ·)` **is** the polarized-target community's equal-spin-temperature relation, to machine precision." — replace "to machine precision" with:
> to the bisection's 1e-13-in-β tolerance: max 1.14e-14 over 19 999 points of (−1, 1) (the five points above are a sample, not the bound), and the tree throws at \|P_z\| = 1

**E6.** same file — §4.1 N-1 item 1, "Every mode lands **on the boundary** of the physical domain (min p_m = 0)." — replace with:
> Every polarized mode (1–7) lands **on the boundary** of the physical domain (min p_m = 0); the unpolarized mode 0 is interior (1/3, 1/3, 1/3).

**E7.** same file — §4.1 N-1 item 3, after "only P_z is measured." (end of the first sentence group) — append:
> The caption's "relative beam intensities" are not printed either.

**E8.** same file — §4.7 CLAS Physics Database paragraph, "Not queried in detail here. **Effort: hours to check. Priority: medium.**" — replace with:
> **Queried 2026-09-23** (`../open_items/run_2026-09-23/phase_B1_spin_deuteron.md` §4): for a deuteron target, 918 measurements in 19 experiment ids. BONuS is one measurement (F₂ⁿ/F₂ᵖ vs W*), EG1b is inclusive (g₁, A₁, g₁/F₁), and the only spectator-momentum distribution is Deeps (Klimenko 2006, F₂ᴺ × P(p_s) at 0.30–0.53 GeV/c, 115 tables, vendored unwired). The BONuS p_s tables are in PRC 89 045206's Supplemental Material (HTTP 401 without an APS login).

**E9.** `docs/benchmarking/02_data_nucleon_deuteron.md` — F-9 row, the HEPData cell "**neither in INSPIRE's HEPData index**" — replace with:
> **neither in INSPIRE's HEPData index**; CLAS Physics Database (2026-09-23): BONuS only as F₂ⁿ/F₂ᵖ vs W* (E135M1); the spectator tables are Tkachenko 2014's Supplemental Material (APS login)

and §5.5 table, row "**HEPData**", append:
> The CLAS Physics Database holds the F₂ⁿ/F₂ᵖ ratio only (queried 2026-09-23); Deeps (Klimenko 2006, p_s = 0.30–0.53 GeV/c) is the only deuteron spectator distribution there and is vendored at `validation/benchmarks/data/clas_deeps_klimenko2006_f2P.json`.

**E10.** `docs/references/REFERENCES.md` — §1 (PLEASE DOWNLOAD), add to §1a or §1d:
> - **Tkachenko** *et al.* (CLAS BONuS), PRC **89** (2014) 045206 — its **Supplemental Material** (the "tables of numerical results", ref. [71]: cos θ_pq spectra for 6 W* × 5 Q² × 4 p_s bins) at https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.89.045206 — HTTP 401 without an APS login (2026-09-23). Unblocks BENCHMARK_PLAN §4 row 10 (`validation/benchmarks/t2_bonus_spectator_shape.py`). Drop the files into `validation/benchmarks/data/`.

and in the §4f blockquote ("… Koivuniemi *et al.* 2004 (COMPASS polarization build-up) … None is cited anywhere in the LiPolGen tree …") — replace "Koivuniemi *et al.* 2004 (COMPASS polarization build-up)" by removing it from the no-free-copy list and add:
> Koivuniemi *et al.* 2004 **has** a free copy (`wwwcompass.cern.ch/compass/publications/technical/koivuniemi_spin2004.pdf`, fetched 2026-09-23, md5 52a3c608…) and is now cited by `validation/benchmarks/t5_est_identity.py` (its Fig. 1 numbers are in `validation/benchmarks/data/est_identity_sources.json`).

**E11.** `docs/PHYSICS_CHANNELS.md` — the row "Spin-temperature (maximum-entropy) fill, any J", reference cell "—" — replace with:
> the polarized-target equal-spin-temperature relation: for J = 1, P_zz = 2 − √(4 − 3P_z²), which `spin_temperature_pzz(1, ·)` reproduces to 1.14e-14 (bisection bound 2.13e-14; `validation/benchmarks/t5_est_identity.py`); the COMPASS ⁶LiD common spin temperature, Brillouin with J = 1 (D, ⁶Li) and J = 3/2 (⁷Li) (Koivuniemi *et al.*, SPIN 2004), whose Fig. 1 populations, T_S and P(⁶Li), P(⁷Li) the ladder reproduces at the printed digit; Keller, Crabb, Day, NIM A 981 (2020) 164504 for EST as the equilibrium baseline RF manipulation breaks. Against an atomic-beam SOURCE the EST fill matches 3 of the 8 EPIOS Table II modes (`validation/benchmarks/t5_epios_source_modes.py`)

(The Documents stage decides whether new bibliography keys are needed for the strict gate.)

**E12.** `docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md` — §B2, end of **(c)** (after "…which needs a pion PDF set the tree does not have.") — append:
> **Priced against data 2026-09-23** (`../run_2026-09-23/phase_B1_spin_deuteron.md` §3, `validation/benchmarks/t2_hermes_b1_table2.py`, no tuning): against HERMES Table II (six bins, stat ⊕ syst), keep 0.5 gives χ²(b₁ᵈ) = 5.262/6 (4 bins within 1σ), drop it 5.297/6 (5 within); χ²(A_zzᵈ, F₁ rebuilt with MSTW + R1990) = 5.550 vs 4.514 (R1998: 5.410 vs 4.587); b₁ = 0 gives 21.80. At x = 0.012 the two options read 5.71 and 11.43 (×10⁻²) against HERMES 11.20 ± 6.17. **HERMES does not decide this row** (Δχ² = +0.035 on b₁), and Miller's P₆q was tuned to a HERMES point, so neither χ² is an independent test.

**E13.** `docs/open_items/run_2026-09-03/STATUS.md` — decision-table row 2 — append to its status cell:
> HERMES-priced 2026-09-23: Δχ²(b₁) = +0.035 for dropping the 0.5 — data do not decide it (`../run_2026-09-23/phase_B1_spin_deuteron.md` §3)

**E14.** `docs/benchmarking/BENCHMARK_PLAN.md` — §6 row U-11, "#1 above fixes this" — replace with:
> #1 above (wired 2026-09-23): the bookkeeping now has an external anchor, EPIOS Table II, passing at residual 0; the anchor is for an atomic-beam deuteron source, not for lithium

**E15.** `docs/benchmarking/BENCHMARK_PLAN.md` — §3, "(iii) writes one row of `docs/benchmarking/REPORT.md`" — append after the sentence "The report is regenerated, never hand-edited.":
> (As built 2026-09-23 the harnesses print and return their row — `run()` returns a dict, a `REPORT |` line is printed — and no `REPORT.md` generator exists yet.)

**E17.** `docs/benchmarking/06_critic.md` — §4.1, the paragraph beginning "The same paper (§V D, read here) records …", from "and (p. 3009 of the text layer) that" to the end of the paragraph — replace with:
> and, in §VIII B 4, lists a "Li-6/Li-7 beam polarimeter" (an atomic-beam target for Li–Li CNI scattering with Breit–Rabi calibration) as **proposed** R&D — i.e. **there is not yet a polarimeter of any kind for a lithium beam**, and the tensor case is not addressed. (The sentence "tensor polarization presents additional challenges, requiring specialized polarimeter development such as Lamb-shift-based systems or BRP-type" is in §VIII A 2, the *deuteron* source, and ends "…specifically designed for deuteron ions"; it was misattributed to lithium here until 2026-09-23.) That is the correct citation for why U-11's measured counterpart does not exist for lithium.

**E16 (not a doc; a gate).** `validation/check_spdx_headers.py` — `GROUPS` — add `"validation/benchmarks/*.py"` (and, if wanted, `"python/tests/*.py"`): the new harnesses carry the header but the gate's `validation/*.py` glob is not recursive, so it does not check them. The SPDX count moved 103 → 104 only through `tests/test_bench_spin.cpp`.

## 7. Suite tallies (this phase, measured at the end, shared working tree)

Measured 2026-09-23 after every file above existed, on the SHARED working
tree (three other agents were adding files and editing tracked files
concurrently, so a tally above the baseline is not this phase's alone; this
phase's own share is stated beside each).

| gate | command | measured | baseline (239ac88) | this phase's share |
|---|---|---|---|---|
| doctest | `build/lipolgen_tests` | **416 cases / 17 242 712 assertions / 0 skipped / 0 failed** (6 m 48 s) | 414 / 17 242 629 / 0 | +2 cases, +81 assertions (`tests/test_bench_spin.cpp`) |
| pytest | `python -m pytest python/tests -q` | **1134 passed, 153 skipped** (298.65 s) | 1107 / 151 | +16 passed, +1 skipped (`test_bench_spin_deuteron.py`) |
| strict gate | `python3 validation/check_physics_channels_links.py` | **1240 / 95 / 7 / 0 broken / 24 allow-listed** | same | none |
| ranges | `… --audit-ranges` | **0 REFUSED** (209 citations, 179 blocks) | 0 | none |
| records | `… --records` | **764 citations in 46 records, 0 broken** | 0 broken | +1 record (this file; no `path:line` citations in it) |
| SPDX | `python3 validation/check_spdx_headers.py` | **104/104** | 103/103 | +1 (`tests/test_bench_spin.cpp`; `validation/benchmarks/*.py` carry the header but are outside the gate's globs — E16) |

## 8. Problems

1. The plan's 1e-14 for row 9 is not a property of the code (1.138e-14 measured, 2.13e-14 bound); the test asserts the bound. Tightening would need `populations_maxent`'s `tol`, a code change outside this phase.
2. `spin_temperature_pzz` throws at \|P_z\| = 1 — two of the eight EPIOS modes can only be compared to the closed-form limit.
3. HERMES's F₂ inputs (ALLM97, NMC 1992) are not in the tree; F₁ is rebuilt from MSTW and the difference priced (2.6–7.3 %), not removed.
4. The R1990 preprint scan prints b₁ = 0.635; 0.0635 is used (evidence in the data file).
5. `PolarizedLithiumSim/refs/refs_dict.json`'s Koivuniemi entry misquotes 51.3 % / 92.4 % as "51 ± 3" / "92 ± 4" (read-only repository; for its owner).
6. Row 10 needs a subscription document; the Deeps table is vendored but a harness for it would be a new plan row (its FSI caveat makes it a joint test of n(k) and `GlauberFsiWeight`).
