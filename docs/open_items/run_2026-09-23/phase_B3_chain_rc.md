<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Phase B3 — harness group CHAIN AND RADIATIVE (run 2026-09-23)

`BENCHMARK_PLAN.md` §4 row 7 and §5 items 3, 4, 5. Base commit `d735abc`.
No physics default, no registry row and no reference-JSON number moved.
Every number below carries the command that produced it. Only this file is
written under `docs/`. Every other doc change is in §6 as
"file: site: replacement text" for the Documents stage.

## 0. Verdicts

| item | verdict |
|---|---|
| §4 row 7, DJANGOH Rad/noRad | **Wired, and BLOCKED.** The ePIC DJANGOH samples ran with the **elastic radiative tail switched off** (`IEL2 = IEL31 = IEL32 = IEL33 = 0` in all eight logs). The elastic tail is the **only** O(α) term LiPolGen's `rc_tail` computes. So the two share no term: nothing either side prints can contradict the other, and no tolerance can be defined. Both sides and their distance are printed. Nothing was tuned. |
| §5.3 eSTARlight claim line | The exact sentence "benchmarks `coherent.hpp`'s channel rate" appears **nowhere** in the tree. Its equivalents do: "baseline for the coherent channel", "outside anchor", "validates the unpolarized coherent rate". The code-header site is fixed in `coherent.hpp` (comment only, line-neutral). Eleven doc and code sites are listed in §6.2. |
| §5.4 POLRAD ADGH | Re-downloaded (sha256 `83703668…1ea7b`) and re-driven. **Rows 1 and 3 reproduce to every printed digit. Row 2 does NOT**: σ_u^d = −1.0303e−03, not −1.001e−03 (+2.93 %), and σ_q/σ_u = +0.0623, not +0.0642. That is a stale number, carried at the nine sites this phase found (§6.3); the Documents stage and its fix stage corrected it at sixteen (§6.3, closing note). The five overstatements of the `D-3` entry are in §3.3. |
| §5.5 `b1_default_li6.json` | Labelled a SELF-PIN in the file's new `provenance` field and in a new `validation/README.md` section. Every numeric block is byte-identical to HEAD, and the rtol 1e−12 gates pass. **The file's sha256 changed** (`d7bd8ce6…` → `8373db8e…`). This conflicts with `phase_A_port_gate.md` §4 of this same run (see §7). |

## 1. §4 row 7 — DJANGOH Rad/noRad four-bin table

### 1.1 What the survey recorded, and what was fetched

The four-bin numbers **were already in the tree**, in `docs/benchmarking/01_generators.md` §2.4. They are quoted from the README of
`github.com/eic/InclusiveDjangohSamples`. So the xrootd store was **not** read. The numbers were re-fetched once, on
2026-09-23, at a pinned commit, and vendored:

```
curl -sS https://api.github.com/repos/eic/InclusiveDjangohSamples/commits/main
   -> 9869d9a0d21b9191a1bbec7294b4851feb5ee20a (2025-12-07T21:57:32Z); licence: null
B=https://raw.githubusercontent.com/eic/InclusiveDjangohSamples/9869d9a0d21b9191a1bbec7294b4851feb5ee20a
curl -sS $B/README.md                                   # 698 B, sha256 19445e90…ca132
curl -sS $B/runcards/ep.Rad={0,1}.NC.Q2_1_10.in          # and Rad=1 Q2_1000_10000
curl -sS $B/metadata/djangoh.NC.{Rad,norad}.18x275_Q2_{1_10,10_100,100_1000,1000_10000}.log
curl -sS -o Djangoh_m.pdf https://raw.githubusercontent.com/eic/documents/master/software/general/Djangoh_m.pdf   # 196 416 B
```

The vendored file is `validation/benchmarks/data/djangoh_rad_norad_ep18x275.json`. It holds the four rows, the pinned URL,
the fetch date and the licence (none stated, so only the eight numbers and the run settings are reproduced). It cites the
papers (HERACLES CPC 69 (1992) 155, DJANGO6 CPC 81 (1994) 381, manual MZ-TH/05-15), the run settings, and the size and
sha256 of each log.

| Q² [GeV²] | Rad=1 [µb] | Rad=0 [µb] | ratio | log check |
|---|---|---|---|---|
| 1–10 | 0.6733778912 | 0.6278853416 | 1.072454 | Rad log **truncated** (9 067 B, ends after PYINIT), so this Rad=1 value rests on the README alone; noRad 627.89 nb ✓ |
| 10–100 | 0.07405000239 | 0.06797276736 | 1.089407 | 74.050 / 67.973 nb ✓ |
| 100–1000 | 0.003760005602 | 0.003239636940 | 1.160626 | 3.7600 / 3.2396 nb ✓ |
| 1000–10000 | 0.00008615325533 | 0.00007015936779 | 1.227965 | 0.086153 / 0.070159 nb ✓ |

The README says "the cross sections are not recorded in the log files". **That is false for seven of the eight.** Each of
those seven prints `Cross-section from HERACLES` to 5 significant figures, and each agrees with the README column to every
printed digit (asserted in `python/tests/test_bench_chain_rc.py`).

### 1.2 The finding that decides the row

Every log echoes `INT-OPT-NC` as `INC2 = 1, INC31..INC34 = 18` (Rad) or `0` (noRad), and **`IEL2 = IEL31 = IEL32 = IEL33 = 0`** (four flags, all zero; this line named three until the 2026-09-23 fix stage — the `INT-OPT-NC` echo is `1 18 18 18 18 0 0 0 0`).
The runcards give all nine values: `1 18 18 18 18 0 0 0 0`.

The DJANGOH manual (§3.1.1, `INT-OPT-NC`) defines the IEL flags. IEL2 is elastic ep → ep. IEL31, IEL32 and IEL33 are
"the quasielastic tail ep → epγ" from initial-state, final-state and Compton radiation. **So Rad=1/Rad=0 in this table is
the inelastic O(α) correction with the elastic radiative tail switched off.**

LiPolGen's `rc_tail` is the opposite piece. It is `1 + σ_tail/σ_Born`: the POLRAD elastic-tail peaks plus the Eq. (44)
quasi-elastic tail, and nothing else. There is no lepton-vertex or inelastic shift anywhere in the library, because
`rc.hpp` prices that piece with the tensor band and never applies it.

The two are disjoint, additive pieces of one O(α) correction. That is why the row is BLOCKED rather than PASS: a
must-not-contradict bound between them holds whatever the numbers are.

It unblocks only with a DJANGOH run at `IEL31..IEL33 ≠ 0`, at these cuts. That needs the source (not public,
`01_generators.md` §5.6) or a new ePIC production.

The two survey sentences that motivated the row are therefore overstated, and §6.1 lists them:

- "constrains the magnitude and Q² trend of LiPolGen's unpolarized RC";
- "turns the RC must-not-contradict bound into a Q²-trend comparison".

### 1.3 The harness, and what it measured

`validation/benchmarks/t4_djangoh_rad_noRad.py` (command: `source env.sh && python3 validation/benchmarks/t4_djangoh_rad_noRad.py`,
19.6 s wall).

**Target.** The closest available target is a **proton**: `PROTON()`, `include/lipolgen/beams.hpp:189`, run at beam config 2
(18 × 275, the DJANGOH beams exactly). The deuteron control was not needed.

**What the harness sets.**

- The only settings off their defaults are a zero coherent spin-1 form factor and `qe_kf_gev = 0`. A proton has no coherent
  spin-1 vertex, and `HoSpin1FF::for_ion` refuses it. With `qe_kf_gev = 0`, the Eq. (44) Z = 1, N = 0 slot **is** the free
  proton's elastic radiative tail, with no Pauli factor.
- The plan is `helicity_flip_plan(0.5, 0, 0)`, which is unpolarised.
- The cuts are DJANGOH's: x ≥ 1e−5, y ≥ 1e−4, W² ≥ 9 GeV².
- The y ceiling is **y ≤ 0.985**. That is the shipped `Scenario.y_max` and the top of the tail table's support (`RcModel`
  refuses to extrapolate). DJANGOH's y runs to 1.

**Integration.** The per-bin σ_RC/σ_Born is ∫Born(1+r)/∫Born over (ln Q², y). Each variable uses 4 × 12 and 2 × 6 × 12
composite Gauss–Legendre nodes. Doubling the panels moves t-peak by ≤ 9.0e−7 and polrad-full by ≤ 1.5e−5 (the harness's
`convergence()`). No node falls outside the tail table.

| Q² bin | DJANGOH | LiPolGen `t-peak` (default) | LiPolGen `polrad-full` | t-peak / DJ excess | polrad-full / DJ excess |
|---|---|---|---|---|---|
| 1–10 | 1.072454 | **1.016169** | **1.027096** | 0.2232 | 0.3740 |
| 10–100 | 1.089407 | **1.013172** | **1.015072** | 0.1473 | 0.1686 |
| 100–1000 | 1.160626 | **1.014970** | **1.015138** | 0.0932 | 0.0942 |
| 1000–10000 | 1.227965 | **1.020686** | **1.020757** | 0.0907 | 0.0911 |

**Q² trend.** d ln(ratio − 1)/d ln Q², by least squares over the four bins, is DJANGOH **+0.1748**, t-peak **+0.0377** and
polrad-full **−0.0345**. So DJANGOH's inelastic excess rises with Q², and the proton elastic tail is flat to ±0.04. The
two pieces behave differently, as expected. On this window the proton elastic tail is 9–37 % of DJANGOH's inelastic excess.

**Diagnostics** (same command family, `build_pipeline`):

- `polrad-full` clips 0.4243 % of its tail-table nodes at `tail_max`, all in the y < 0.5 band. `t-peak` clips 0.
- At the shipped 6Li `qe_kf_gev` = 0.169 (Pauli-suppressing a free proton, which is unphysical and reported only so that
  the setting is visible), t-peak reads 1.003361 / 1.003893 / 1.007193 / 1.017332.
- The Born σ per bin, informational only because the SF and y window differ: LiPolGen (shipped `UnpolSfSource.Toy`,
  y ≤ 0.985) gives 0.884691 / 0.0756728 / 0.00372591 / 7.38902e−05 µb. DJANGOH Rad=0 (CTEQ6.1M, y ≤ 1) gives
  0.627885 / 0.0679728 / 0.00323964 / 7.01594e−05 µb. The ratios are 1.409 / 1.113 / 1.150 / 1.053.

**Pytest.** `python/tests/test_bench_chain_rc.py` has 3 tests, which run always in ~19 s (under the 30 s opt-in threshold).
They check that the vendored table and its log cross-check are intact, that the row is BLOCKED with its reason and prints
both sides, and that the numbers above hold at rel 2e−7.

## 2. §5.3 — eSTARlight

The critic's quote ("It is not a like-for-like benchmark of `coherent.hpp`'s channel rate") is a correction. No
tree site carries that exact sentence. What the tree says instead is below: each site either calls eSTARlight a
baseline or anchor for the **coherent channel**, or says it validates the **coherent rate**.

**The 85–97 % figure, measured.** Here it is from `estarlight_li6_q2_floors()` (source `lipolgen._lipolgen`), computed as
1 − σ(Q² > 0.1)/σ(no floor):

| channel | computation | fraction below Q² = 0.1 |
|---|---|---|
| J/ψ | 1 − 1.773/11.971 | **85.2 %** |
| φ | 1 − 30.16/654.344 | **95.4 %** |
| ρ | 1 − 506.4/17823 | **97.2 %** |

The coherent channel generates nothing there: `make_config(channel="coherent").scenario.q2_min` = 0.7.

**Code site, fixed (comment only).** `include/lipolgen/coherent.hpp:131`–133, the `EstarlightLi6Row` comment, now reads
"an EXCLUSIVE-VM rate scale (item 11.1), NOT a check on the coherent channel rate". It points at a new ESTARLIGHT SCOPE
block (`include/lipolgen/coherent.hpp:776`, headed "ESTARLIGHT SCOPE (2026-09-23, BENCHMARK_PLAN.md sec. 5.3"), which carries the numbers above.

The edit is **line-neutral above line 776**. A first version inserted 10 lines at 131 and broke 42 strict-gate citations
(`PHYSICS_CHANNELS.md` cites `coherent.hpp` by line number -- 168, 191, 407 and more), so it was reverted to this form.

The remaining sites are in tests/, validation/ or docs/, which are outside this phase's partition. They are listed in §6.2.

## 3. §5.4 — POLRAD ADGH

### 3.1 Re-download

The command recorded in `polrad_transcription_check.md` §0 **no longer works as written**. Mendeley answers a plain
`curl -sL` with a 395-byte JSON error. It works with a browser user agent:

```
curl -sL -A "Mozilla/5.0" -o adgh_v1_0.gz "https://data.mendeley.com/public-files/datasets/37vgvzgr2w/files/50bc36c4-d146-4487-b697-218850c0a057/file_downloaded"
# -> 97 335 B, sha256 83703668bacfbffe073a63c75374be261b4b36dab8ea8de951762120fba1ea7b (= the API's sha256_hash)
gunzip -c adgh_v1_0.gz > adgh     # 11 593 lines, sha256 4a823438e24b979c0e7366b59412438f20adff83e74c0652f2f3e18197811761
```

### 3.2 Re-drive, with the command

The pieces below are POLRAD's own. The integrator is POLRAD's too; no quadrature was substituted.

- **Verbatim decks.** The `ffdeu` deck (adgh lines 5480–5550) and the `qunc8` deck (lines 7953–8080), extracted with
  `sed -n`.
- **Transcribed.** The `targ_d` bodies of `elu` (adgh line 8893) and `elq` (line 8985).
- **Driver.** It reproduces `conkin`: `s = snuc·amp/amh` (adgh line 747). It also reproduces `apptai`'s `ita = 2`
  branch: `ter = amh/amp`, `eta1`/`eta2`/`yy1`, `siau = ter α³/s · yy1 · ∫elu` and `siaq` likewise (adgh lines
  8608–8643). As in the check's §8, `barn` is omitted.

The target mass is `amt` = 1.87561 (POLRAD's `targ_d` data). Using 1.8756280 instead moves only the fifth digit of `∫elu`.

```
sed -n 5480,5550p adgh > ffdeu.f ; sed -n 7953,8080p adgh > qunc8.f
gfortran -O0 -std=legacy -ffixed-line-length-none -o drv drv.f ffdeu.f qunc8.f   # GNU Fortran 11.4.0
echo 1.87561 | ./drv
```

`drv.f` is reproduced in the Appendix.

| beam, x, y | ∫dη elu | σ_u^d | σ_q^d | σ_q/σ_u | tree (`polrad_transcription_check.md` §8) |
|---|---|---|---|---|---|
| 27.6, 0.050, 0.60 | −4.31996e+03 | **−2.35222e−05** | −2.49227e−06 | **+0.1060** | −4.320e+03 / −2.352e−05 / −2.492e−06 / +0.1060 ✓ |
| 27.6, 0.012, 0.50 | **−2.19502e+05** | **−1.03034e−03** | −6.42349e−05 | **+0.0623** | −2.133e+05 / **−1.001e−03** / −6.429e−05 / **+0.0642** ✗ |
| 11.0, 0.200, 0.50 | −1.33306e+01 | **−1.57002e−07** | +1.82949e−08 | **−0.1165** | −1.333e+01 / −1.570e−07 / +1.829e−08 / −0.1165 ✓ |

**Convergence.** `qunc8` reports fl = 0 and er/|∫| ≤ 6.3e−11 on all rows. Row 2 was also checked independently with a
log-trapezoid in ln η: 2 × 10⁴, 2 × 10⁵ and 2 × 10⁶ nodes give −219502.1306, −219502.1325 and −219502.1325.

**So row 2's σ_u (and the ratio derived from it) is stale.** σ_q reproduces to 0.09 %, which says the kinematics were
right and the σ_u integral was not. The survey's verifier said "the three numbers reproduce". This re-drive does not
confirm that for row 2.

**The in-tree comparison, re-measured.** The tree's HO deuteron stand-in, via `polrad_sigma_el_u` (a C++ one-off linked
against `build/libLiPolGenCore.so`, as in `tests/test_rc.cpp` T8(c′)), gives su/POLRAD = **0.8834 / 0.9828 / 0.3896**
against the tree's rows AS THEY STOOD (stale row 2). Against the re-driven row 2 the middle value is 0.9548, so the `00` D-3 row's 0.983 was NOT current; corrected to 0.955 on 2026-09-23 (re-verification pass).

Against the re-driven row 2 the middle value is **0.9548**. T8(c′)'s 20 % gate passes either way, so correcting the row
moves no test verdict. The tree's st/su is **+0.0658 / +0.0466 / −0.4441**, against POLRAD's +0.1060 / +0.0623 / −0.1165.

### 3.3 The five overstatements in the `D-3` entry (`docs/benchmarking/00_in_tree_checks.md:132`, the POLRAD 2.0 / Akushevich row)

`06_critic.md` §3.1 records only the verdict ("overstates its scope in five places … one stale number"). The verifier's
own list is **not in the tree**. The five below are this phase's reading of the D-3 row against `tests/test_rc.cpp`
T8(c′) and T10.

1. **"compared: … G_C(0)=1, G_M(0)=1.7139610634, G_Q(0)=25.84".** Only G_M(0) is gated against POLRAD's number (T10,
   1e−7). G_Q(0) is gated against the tree's own 25.83 at 1e−3: the HO stand-in gives 25.8294 from Q_d = 0.2859 fm²,
   0.04 % below POLRAD's 25.84, which nothing asserts. G_C(0) = 1 is Z by construction.
2. **"validates … the right power of A".** For the deuteron a wrong power of A is a factor **2**. The 0.25×–4× gate passes
   that at every row. Only the two x < 0.1 rows' 20 % gate catch it.
3. **"… no missing Y₊".** Y₊ = 2.90 (y = 0.6) and 2.50 (y = 0.5). The same applies: only the two low-x rows catch it.
4. **"(any of those moves the magnitude ≥ 6×)".** This is false at A = 2 (power of A 2×, Y₊ 2.5–2.9×). The "6×" is the
   ⁶Li A. The only one of the four that is ≥ 6× here is the dropped minus, and the `su > 0` sign check catches that
   whatever its size.
5. **"The sign of σ_q/σ_u … is gated tightly".** Only the **sign** is gated (`(st/su)·ratio_polrad > 0`). The magnitudes
   differ by ×0.62 / ×0.75 / ×3.81 (tree/POLRAD, re-driven), because the tree uses a different form factor. "Tightly" is
   true of a sign bit, not of the ratio.

The stale number is row 2 (§3.2). The nine sites this phase found are in §6.3; the sixteen actually corrected are in §6.3's closing note.

## 4. §5.5 — `validation/reference/b1_default_li6.json`

A new top-level `provenance` key was added. With `sort_keys=True` it lands after `kernel`. It opens "SELF-PIN, NOT A
REFERENCE OF ANY KIND", names `dump_b1_default_li6.py` and the last commit (`b1071b1`, 2026-09-03), and says what
agreement with the file proves ("unchanged", never "correct"). The file was re-serialised with the dumper's own format
(`json.dumps(indent=1, sort_keys=True) + "\n"`, which round-trips the HEAD file exactly).

**Byte identity, checked:**

- The HEAD file's 2 657 lines equal the new file's 2 658 lines minus the one `provenance` line. The numeric blocks `x`,
  `q2`, `tables` and `kernel` are therefore byte-identical, and `description` and `generator` are unchanged.
- sha256 moved from `d7bd8ce6a01b581d8ff65002f74b315cd4268dc15b280701ce6465002d78cd9e` (= `git show HEAD:…`) to
  `8373db8ed5637ea870717e1b629f3f5c73f2f30dc442b44591b8e0259bd3d777`.

**The rtol gates still pass.** `tests/test_b1_nuclear.cpp` T9 passes (inside the doctest total in §5), and so does
`python/tests/test_b1_model.py::test_default_kernel_is_the_miller_kernel` (inside the pytest total).

`validation/README.md` gained a `## b1_default_li6.json` section saying the same. It is the only paragraph on the file,
because none existed before.

**Caveat.** `validation/dump_polligen_reference.py` rewrites `README.md` wholesale, and `dump_b1_default_li6.py` does not
write the `provenance` key. Either regeneration drops the label; §6.4 lists both script edits.

## 5. Suites and gates (all at the end of this phase; the tree includes other agents' concurrent work)

| check | result | baseline (239ac88) |
|---|---|---|
| `cmake --build build -j` | clean | — |
| `build/lipolgen_tests` | **416 cases / 17 242 712 assertions / 416 passed / 0 skipped** | 414 / 17 242 629 (this phase adds no C++ test; the +2 cases are another agent's `tests/test_bench_spin.cpp`) |
| `python -m pytest python/tests -q` | **1134 passed, 153 skipped** (291.7 s) | 1107 / 151 (this phase adds 3 passed; the rest are other agents' `test_bench_*`) |
| `check_physics_channels_links.py` (strict) | **1240 / 95 / 7 / 0 broken / 24 allow-listed**, exit 0 | same |
| `--audit-ranges` | **0 REFUSED** (209 citations, 179 blocks) | 0 |
| `--records` | **0 broken** (764 citations in 46 records), exit 0 | 0 |
| `check_spdx_headers.py` | **104 / 104** | 103 (+1 not from this phase; `validation/benchmarks/*.py` is not in its globs — the new harness carries the header anyway) |

## 6. Doc edits requested (for the Documents stage)

### 6.1 DJANGOH (row 7)

- `docs/benchmarking/BENCHMARK_PLAN.md`: §4 row 7, "what it buys" cell. Replace
  "turns the RC "must-not-contradict" bound into a Q²-trend comparison; zero dependencies" with:
  "**wired 2026-09-23 as `t4_djangoh_rad_noRad` and BLOCKED**: the ePIC samples ran with the elastic radiative tail OFF
  (IEL2=IEL31=IEL32=IEL33=0), the only O(α) term `rc_tail` computes. The two share no term, so no tolerance exists
  (`../open_items/run_2026-09-23/phase_B3_chain_rc.md` §1). It unblocks with a DJANGOH run at IEL31..33 ≠ 0."
- `docs/benchmarking/01_generators.md`: §2.4 item 2, last sentence ("This is an **obtainable-today table** … without
  running anything."). Replace with: "This table **excludes the elastic radiative tail**: every run has IEL2 = IEL31 =
  IEL32 = IEL33 = 0 (logs at `9869d9a`). Its ratio is the inelastic O(α) correction only, so it constrains nothing
  LiPolGen computes, and LiPolGen's `rc_tail` is the elastic tail only (phase_B3_chain_rc.md §1.2)."
- `docs/benchmarking/01_generators.md`: §6 item 1 ("Constrains the *magnitude and Q² trend* of the unpolarized RC …
  Payoff: turns `C-8`'s single must-not-contradict bound into a trend check."). Replace with: "Wired 2026-09-23 and
  BLOCKED: the table is the inelastic RC with the elastic tail off, and shares no term with `rc_tail`. Payoff after a
  DJANGOH run with IEL31..33 on: a check of the elastic-tail magnitude on a proton."
- `docs/benchmarking/01_generators.md`: the rows for `D-3` and `D-3`/`C-8` in the §4 table. In "Even the published
  Rad=1/Rad=0 table (§2.4) constrains the RC magnitude without running anything", replace with "the published Rad=1/Rad=0
  table (§2.4) has the elastic tail OFF and so constrains no LiPolGen term". In "DJANGOH `Rad=1` vs `Rad=0` ratios (§2.4
  table) … turns a one-point bound into a Q²-trend comparison", replace with "… would, **with the elastic tail on**; the
  published table has it off".
- `docs/benchmarking/05_chain.md`: the §2.3 DJANGOH store row. Replace "The Rad/noRad ratio in (x, Q²) *is* the QED
  radiative-correction factor, produced externally" with "The Rad/noRad ratio is the **inelastic** QED correction
  factor. The production ran with the elastic radiative tail OFF (IEL2 = IEL31..33 = 0, `eic/InclusiveDjangohSamples`
  logs)".
- `eic/InclusiveDjangohSamples` README sentence "The cross sections are not recorded in the log files". This is
  external; note only, do not edit. Seven of eight logs do record them (§1.1).

### 6.2 eSTARlight (§5.3)

Each site below gets the same added clause (or the given replacement):
"**a rate scale for exclusive ρ/φ/J/ψ, not a check on the coherent channel's M_X ≥ 1.2 GeV continuum rate — 85–97 % of
eSTARlight's rate (J/ψ 85.2 %, φ 95.4 %, ρ 97.2 %) is at Q² < 0.1 GeV², where `coherent.hpp` generates nothing**".

- `docs/benchmarking/00_in_tree_checks.md`: row D-4, "what it validates". Replace "an **unpolarized rate and slope
  baseline** for the coherent channel" with "an **unpolarized rate and |t|-slope scale for exclusive VM production on
  ⁶Li** (not the coherent channel's rate: …clause…)".
- `docs/benchmarking/01_generators.md`: §1 table, row "Coherent ⁶Li diffraction", column 3. Replace "the unpolarized
  coherent rate, the |t| slope and the VM ratios" with "the unpolarized **exclusive-VM** rate scale, the |t| slope and
  the VM ratios (not the coherent continuum channel's rate)".
- `docs/benchmarking/01_generators.md`: §2.7 "Overlap". Replace "Already exhausted for the coherent ⁶Li rate" with
  "Already exhausted for the exclusive-VM coherent ⁶Li rate (a scale for exclusive ρ/φ/J/ψ, not for `coherent.hpp`'s
  continuum)".
- `docs/PHYSICS_CHANNELS.md:421`: the eSTARlight row. Replace "to give this scenario an outside anchor." with "to give
  this scenario an outside **rate scale for exclusive vector mesons** — not a check on the coherent channel's rate (…clause…)."
- `docs/PHYSICS_CHANNELS.md:1046`: bibliography [53]. Replace "used for the unpolarized coherent baseline of §9" with
  "used for the unpolarized exclusive-VM rate scale of §9 (not a check on the coherent continuum channel)".
- `docs/OPEN_ITEMS_SOLUTIONS.md:1611`: the heading "11.1 eSTARlight — the unpolarized rate and slope baseline". Replace
  with "11.1 eSTARlight — the unpolarized exclusive-VM rate and slope scale (not a coherent-channel benchmark)".
- `docs/open_items/run_2026-09-03/mantysaari_collaboration_draft.md:67`: replace "An unpolarized coherent-rate baseline
  that IS citable" with "An unpolarized exclusive-VM coherent-rate scale that IS citable".
- `README.md:177` (as of d735abc) and `docs/DEVELOPMENT_PLAN.md:213` (as of d735abc) — applied 2026-09-23 by the Documents stage: replace "the eSTARlight unpolarized (coherent ⁶Li) baseline" with
  "the eSTARlight unpolarized exclusive-VM rate scale".
- `validation/o5_a2_reach.py:17`: replace "the eSTARlight unpolarized coherent rate baseline" with "the eSTARlight
  unpolarized EXCLUSIVE coherent-VM rate scale (not a coherent-channel benchmark)". The O5 use itself, exclusive J/ψ, is
  legitimate.
- `tests/test_coherent.cpp:145` and `:150` (TEST_CASE name "the eSTARlight 6Li baseline"): add "(an exclusive-VM rate
  scale, not a check on the coherent channel's rate)". The TEST_CASE name may stay.

### 6.3 POLRAD ADGH (§5.4)

**The D-3 entry** (`docs/benchmarking/00_in_tree_checks.md:132`, the POLRAD 2.0 / Akushevich row):

- "compared" cell. Replace "G_C(0)=1, G_M(0)=1.7139610634, G_Q(0)=25.84" with "G_M(0) = 1.7139610634 (1e−7); G_Q(0) is
  the tree's 25.83 (1e−3), 0.04 % from POLRAD's 25.84, which is not asserted; G_C(0) = 1 is Z".
- "validates" cell. Replace "that Eq. (38) is transcribed with the right power of A, the right prefactor, no missing Y₊
  and no dropped leading minus (any of those moves the magnitude ≥ 6×)" with "no dropped leading minus (the sign gate).
  A wrong power of A (×2 at A = 2), a missing Y₊ (×2.5–2.9) or a prefactor error > 20 % is caught **only at the two
  x < 0.1 rows**; the 0.25×–4× x = 0.2 row passes all three".
- "tolerance" cell. Replace "The **sign** of σ_q/σ_u … is gated **tightly**" with "only the **sign** of σ_q/σ_u is
  gated; its magnitude differs from POLRAD's by ×0.62 / ×0.75 / ×3.81 (different form factor)".
- Append to the row: "POLRAD's row 2 was re-driven 2026-09-23: σ_u = −1.0303e−03 (not −1.001e−03), σ_q/σ_u = +0.0623
  (not +0.0642); the middle su/POLRAD becomes 0.955."

**The stale number, +0.064 / −1.001e−03 / −2.133e+05 → +0.0623 / −1.0303e−03 / −2.195e+05, at every site:**

- `docs/open_items/run_2026-09-02/polrad_transcription_check.md:533`: the row "27.6 GeV, 0.012, 0.50 | −2.133e+05 |
  **−1.001e−03** | −6.429e−05 | +0.0642". Replace with "−2.195e+05 | **−1.030e−03** | −6.423e−05 | +0.0623", plus a
  dated correction note citing this record.
- `…/polrad_transcription_check.md:568`: "(+0.106, +0.064, −0.117)" → "(+0.106, +0.062, −0.117)".
- `docs/open_items/run_2026-09-02/phase_C_numbers.md:215`: "(+0.106, +0.064, −0.117 at three points)" →
  "(+0.106, +0.062, −0.117 at three points)".
- `docs/open_items/run_2026-09-03/phase_B_numbers.md:1037`: "+0.064" → "+0.062". Also `:1038`: "**68×**" → "**67×**"
  (9.357e−04 against 0.0623 is 66.6).
- `docs/open_items/run_2026-09-03/SUMMARY.md:59-60`: "a factor **68** resp. **≈ 260** below POLRAD's own deuteron
  elastic-tail value +0.064" → "a factor **67** resp. **≈ 255** below … +0.062".
- `include/lipolgen/rc.hpp:940` and `src/core/rc.cpp:1960`: "(check sec. 8: +0.106, +0.064, -0.117 …)" →
  "+0.106, +0.062, -0.117". This is a comment-only code edit, but it was not named for this phase.
- `include/lipolgen/rc.hpp:1576-1581`: "+0.106 / +0.064 / -0.117 … a factor 35 / 68" → "+0.106 / +0.062 / -0.117 …
  a factor 35 / 67". Also "~260" (the vmc-ft edge, 0.0623/2.442e−04 = 255) → "~255". Comment only, not named for this
  phase.
- `tests/test_rc.cpp:1257`: `{27.6, 0.012, 0.50, 1.001e-03 …, 0.0642}` → `1.03034e-03 …, 0.0623`, and the measured
  comment at `tests/test_rc.cpp:1267` ("Measured with this repository's HO stand-in: 0.883, 0.983, 0.390") → "0.883, 0.955, 0.390". Every CHECK still passes: 0.9548 is within the 20 %
  gate, and the sign is unchanged.
- `docs/open_items/run_2026-09-02/polrad_transcription_check.md` §0, the FORTRAN row: the `curl -sL` command now needs
  `-A "Mozilla/5.0"` (§3.1).
  **Closing note, 2026-09-23 (Documents stage, completed by its fix stage) — the sites actually corrected**, found by `grep -rn -E '\+0\.064[^0-9]|1\.001e-0?3|68 ?×|≈ ?260'` over `docs README.md include src tests python validation`: the nine above (the `run_2026-09-03/SUMMARY.md` one now at line 64; `polrad_transcription_check.md` two; `rc.hpp` two; `rc.cpp` one; `tests/test_rc.cpp` one, value and comment; the 2026-09-02 and 2026-09-03 phase records one each), plus `run_2026-09-03/SUMMARY.md` lines 166–168, `run_2026-09-03/STATUS.md` line 33 (its 67 / ≈ 255 by the Documents stage, its "+0.064" by the fix stage), `OPEN_ITEMS_SOLUTIONS.md` §9 (two rows), `README.md` line 54, the D-3 row of `docs/benchmarking/00_in_tree_checks.md`, and `docs/PHYSICS_CHANNELS.md` line 388 (the polarised quasi-elastic row, missed by the Documents stage and corrected by its fix stage) — sixteen. The old values now survive only as dated "was …" annotations beside the new ones, in this record, and as unrelated numbers the grep also hits (a GFMC +0.064 in `design_G_cluster_config.md`, a ×68.5 weight spread).

### 6.4 `b1_default_li6.json` (§5.5): sites that call it a reference

- `validation/README.md` lines 3–8 (intro, generated by `dump_polligen_reference.py`): "Numeric reference tables dumped
  from the Python event generator `polligen`" → add "— except `b1_default_li6.json`, a self-pin of LiPolGen's own C++
  (see its section)". The same text should go in `dump_polligen_reference.py`'s README writer, which must also emit the
  `## b1_default_li6.json` section, or the next re-dump deletes it.
- `validation/dump_b1_default_li6.py` `build()`: add the same `"provenance"` string as the JSON, so that a regeneration
  keeps the label.
- `docs/PHYSICS_CHANNELS.md:208`: "pinned at rtol 1e−12 against `validation/reference/b1_default_li6.json`" →
  "pinned at rtol 1e−12 against the **self-pin** `validation/reference/b1_default_li6.json` (this library's own C++ output
  — a regression guard, not a reference)".
- `docs/USAGE.md`, the `miller` (default) row of the `--b1-model` table: the same replacement.
- `docs/open_items/run_2026-09-03/phase_A_miller_normalisation.md:386`: "the only b₁ reference, `b1_default_li6.json`"
  → "the only b₁ self-pin, `b1_default_li6.json` (LiPolGen's own C++ output, not a reference)".
- `docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md:135-136` ("every rtol-1e−12 reference carrying a tensor entry
  (`validation/reference/{xsec,b1_default_li6,…}`") and the §B6 sentence "**Reference files:** `b1_default_li6`, …" of the same file: add "(the
  first is a self-pin of LiPolGen's own C++, not a polligen reference)". This is prose in the registry file, not a
  registry row; if that prose counts as frozen, leave it.
- `docs/open_items/run_2026-09-03/phase_D_sf_injection.md:385-388` ("All eight files … are compared against kernels"):
  add "(seven polligen references and one self-pin, `b1_default_li6.json`)".
- `docs/benchmarking/BENCHMARK_PLAN.md` §5 item 5: mark "APPLIED 2026-09-23 (provenance field + README section;
  numeric blocks byte-identical; sha256 d7bd8ce6… → 8373db8e…)".
- **Already correct, no edit:** `docs/benchmarking/00_in_tree_checks.md`, row B-1 ("INTERNAL (self-pin)") and the paragraph after that table.
- **Historical sha statements** (true when written; any "now" reading must be updated):
  - `docs/open_items/run_2026-09-06/phase_CW_numbers.md:385` and `:479` (`d7bd8ce6…` "same");
  - `docs/open_items/run_2026-09-03/SUMMARY.md:24` and `docs/open_items/run_2026-09-06/SUMMARY.md:49` ("sha256-identical
    to `a94fd6e`");
  - `docs/open_items/run_2026-09-23/phase_A_port_gate.md:110` (see §7).

## 7. Problems

1. **Cross-agent conflict.** `docs/open_items/run_2026-09-23/phase_A_port_gate.md:110` (this run, another agent) records
   `b1_default_li6.json` as `d7bd8ce6a01b581d…` "identical". After §5.5's metadata label the file's sha256 is `8373db8e…`,
   with every numeric block byte-identical. That table needs a footnote, or its check needs to compare numeric blocks
   rather than whole-file hashes.
2. **Row 7 is a benchmark that cannot be run as planned.** BENCHMARK_PLAN.md's "zero dependencies" is true of the table
   and false of the comparison. Unblocking needs DJANGOH with the elastic tail on (source not public), or an ePIC
   production request. No LiPolGen-side change could unblock it without being a new physics model.
3. **The DJANGOH window differs from LiPolGen's in two places.** The y ceiling is 1 against 0.985 (LiPolGen's tail table
   refuses more). The Born SF is CTEQ6.1M LO against the shipped Toy SF, and it differs by 5–41 % per bin (§1.3). Neither
   affects the verdict, because no term is shared.
4. **The `polrad_transcription_check.md` §0 download recipe has rotted** (§3.1): it needs a user agent.
5. **POLRAD row 2 is stale** (§3.2). The fix is nine sites, three of them comments in include/src and one in tests, none of them in
   this phase's partition, so all are listed in §6.3 and none was applied by this phase (the Documents stage and its fix stage applied them; §6.3, closing note).
7. **One transient doctest failure.** Of four `build/lipolgen_tests` runs at the end of this phase, one reported
   416 cases / 1 failed / 17 242 711 of 17 242 712 assertions. It was not captured, and three runs, before and after it,
   passed 416/416. Other agents were rewriting `validation/reference/tagged.json` and `_manifest.json` in the same tree at
   the time, and this phase touches no C++ test and no code path, so the cause is unidentified and not attributed.
6. The five §5.4 overstatements are this phase's identification. The survey verifier's own list is not in the tree to
   compare against.

## Appendix — `drv.f` (the POLRAD re-drive wrapper; `ffdeu.f`/`qunc8.f` are `sed -n` extracts of `adgh`)

```fortran
      double precision function elu(eta)
      implicit real*8(a-h,o-z)
      common/kdrv/amp2,y,sx
      xa=y/sx
      xx1=xa**2 + 4.*xa*eta - 4.*eta
      xxt=xx1/(2.*eta*xa**2)
      t=4.*amp2*eta
         call ffdeu(t,fc,fm,fq)
         elu=((fc**2+8./9.*fq**2*eta**2+2./3.*fm**2*eta)*xxt
     .         -2./3.*(1.+eta)*fm**2)/eta
      end
      double precision function elq(eta)
      implicit real*8(a-h,o-z)
      common/kdrv/amp2,y,sx
      xa=y/sx
      xx1=xa**2 + 4.*xa*eta - 4.*eta
      xxt=xx1/(2.*eta*xa**2)
      t=4.*amp2*eta
         call ffdeu(t,fc,fm,fq)
         elq=((1.+eta+(.75*xa**2-eta)*xxt)*fm**2
     .        -xxt*xx1/(1.+eta)*fq*(3.*fc+3.*eta*fm+eta*fq)
     .        -2.*eta*xxt*fq*(4.*fc-3.*xa*fm+4./3.*eta*fq))/eta
      end
      program polrdrv
      implicit real*8(a-h,o-z)
      external elu,elq
      common/kdrv/amp2,y,sx
      dimension be(3),xb(3),yb(3)
      data be/27.6d0,27.6d0,11.0d0/, xb/0.05d0,0.012d0,0.20d0/,
     .     yb/0.60d0,0.50d0,0.50d0/
      alfa=.729735d-2
      amh=.938272d0
      aml2=.261112d-6
      read(*,*) amtar
      bo=1d0
      do i=1,3
        snuc=2.*amh*sqrt(be(i)**2+aml2)
        xs=xb(i)
        ys=yb(i)
        y=snuc*xs*ys
        amp=amtar
        amp2=amp**2
        s=snuc*amp/amh
        sx=s*ys
        sqly=dsqrt(sx**2+4.*amp2*y)
        ter=amh/amp
        um=sx-y
        eta1=(um*(sx-sqly)+2.*amp2*y)/(8.*amp2*(um+amp2))
        eta2=sx/4./amp2
        yy1=(1.+(1.-ys)**2)/(1.-ys)
        call qunc8(elu,eta1,eta2,1d-9*bo,1d-7,relu,er,nn2,fl2,3500)
        siau=ter*alfa**3/s*yy1*relu
        call qunc8(elq,eta1,eta2,1d-9*bo,1d-7,relq,er,nn2,fl2,3500)
        siaq=ter*alfa**3/s*yy1*relq
        write(*,'(f6.1,f7.3,f6.2,3es14.5,f10.4,es10.2,f7.2,i6)')
     .    be(i),xs,ys,relu,siau,siaq,siaq/siau,er,fl2,nn2
      enddo
      end
```
