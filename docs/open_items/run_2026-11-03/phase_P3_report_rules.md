<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-03 — phase P3: the report's rules (decisions A5, A6, A7)

Worker P3, 2026-11-06, on `claude/run-2026-10-26` at `df56035` (phases P1 and
P2 committed; this phase uncommitted, the supervising session commits).
Decisions: `DECISIONS.md` A5, A6 and A7 — the SUPERVISOR's choices under the
author's delegation ("Go with your recommendations for the 17 items, just
record the choice and why"), not the author's words. Every number below was
measured in this phase with `source env.sh` on this machine (4 cores, one long
gate at a time), against the committed `docs/benchmarking/REPORT.json` of
`df56035` (copied before any edit); scripts and logs are in the session scratch
`run_1103/p3/` (`pre_*.log`, `post_*.log`, `a5_survey.py` / `.json`,
`t5_numbers.json`, `t5_ab.txt`, `pixdiff.py`, `gate_*.log`, `result.json`).
The phase ran in two stages: the first applied A5 to Deeps alone and listed
four T2 rows for the supervisor; the supervisor decided them in the
supervisor's note of 2026-11-06 under A5 (and amended A6 from §1.3), both in
`DECISIONS.md`, the supervising session's choices under the author's
delegation; the second stage applied that note. This record describes the
final tree.

## 0. What changed, in one paragraph

A5 is written once, as convention 8 of `validation/benchmarks/README.md` — one
line after convention 7, with convention 7's last two lines joined, so no later
README line moves (a first 13-line insertion broke three run-record citations
of README rows 13 and T6; they hold again): a row's headline verdict is
computed on the library's shipped default configuration, the configuration a
declaration names is a recorded sub-row, both numbers stay visible, a harness
states which object its headline is built on; refined by the supervisor's note
of 2026-11-06 under A5 — a row whose object is an opt-in component builds its
headline with that component ON in its shipped form, the sub-rows built on the
default configuration follow the headline — and with ONE named exception, row
16 (`t2_boeglin_nk`, verdict wave inside its hashed block). Applied: **Deeps**
headline AV18 579.12/187 (Hulthén 627.85/187 recorded), FAIL → FAIL;
**BONuS** headline AV18 with `GlauberFsiWeight` at 40 mb, 3/12 and 0/16 (the
declared Hulthén 4/12, 0/16 recorded; stat ⊕ syst, forward rise and PWIA
prescription follow to AV18, the Hulthén ones recorded), FAIL → FAIL;
**deuteron statics** headline AV18, η −1.39σ, PASS, the declared all-waves
rule a recorded sub-row, FAIL — status **FAIL → PASS by the rule**, no η
moved; **HERMES b₁** headline `miller_shipped`, χ²(b₁)/6 5.26 (4/6), A_zz
7.31/6 (3/6) on the shipped F₁ᵈ (CT18NLO through the shipped deuteron kernel,
`r_sigma_lt`), the declared `cdks_conv_mstw` 22.26 / 22.01 recorded unchanged,
FAIL → FAIL; **Boeglin** unchanged (headline Hulthén, AV18 already the sub-row
`nk[av18]` computed outside the block; a dated addendum outside the block).
The report moves 14 / 9 → **15 / 8**, the one status move the rule's. A6:
the harness carried V6 in the TIGHT edge as `PRE_D4_OVER_JM = 2.0` on M₂, not
as `DS16_OVER_TREE` (1 since D4); the edge is now TIGHT = V1 + V2(−), every
TIGHT saturation halves (toy 0.3207 / 0.3808 / 0.4281 → 0.1603 / 0.1904 /
0.2141), the pre-A6 edge recorded as `saturation_tight_v6`, V6's own fields
unchanged, no status moves (the supervisor's amendment of 2026-11-06 under A6
records this). A7: the grid row carries `m2_error_bound = QUAD_REL_TOL × |M₂|`,
errors 3.26e−14 / 4.90e−14 / 5.53e−14 (2.2–3.7e−10 relative, 448× / 302× /
271× inside), the bound chosen AFTER the values were seen; the toy's 1e−15
stays. Findings recorded, not acted on: on the AV18 headline weight two of
BONuS's header CHECKS miss their thresholds (table convergence 8.3e−4 > 5e−4,
bin averaging −3.2 % > 2 %; no block verdict flips); the shipped inclusive
deuteron kernel's own b₁ is 0.32753 × `MillerB1()` (the ⁶Li mapping). No
physics default, constant, tolerance, registry label or reference JSON moved;
`delta_moment_bound_sources.json` and every hashed design block untouched; no
C++ file changed.

## 1. Every number that moved, before → after, with the recipe

### 1.1 A5 — the Deeps headline

Recipe: `python3 validation/benchmarks/t2_deeps_spectator_tail.py --plot`
before (`pre_t2_deeps_spectator_tail.log`) and after
(`post_t2_deeps_spectator_tail.log`); every χ² is the harness's own, nothing
recomputed.

| quantity | before | after |
|---|---|---|
| headline wave (`DEFAULT_WAVE`, the sub-row `in_verdict`) | `hulthen` (`hulthen/lc`) | **`av18`** (`av18/lc`); `DECLARED_WAVE = "hulthen"` new |
| headline χ²/ndf (p) | 627.85/187 (4.93e−49) | **579.12/187 (1.11e−41)** |
| `hulthen/lc` | 627.85/187, FAIL, the verdict | 627.85/187, FAIL, recorded (same rule) |
| `av18/lc` | 579.12/187, FAIL, recorded | 579.12/187, FAIL, the headline |
| row status | FAIL | **FAIL** |
| per-(Q², W*) table printed under the headline | Hulthén's (e.g. Q² = 1.8, W* = 0.94: 32.56/20) | AV18's (36.26/20 there; 5.98/11 and 4.12/4 at W* = 2.4) |
| `hulthen/lc x GlauberFsiWeight` | 422.92/187 | 422.92/187 (kept on `DECLARED_WAVE`) |
| transverse window (data/PWIA, tree Glauber per p_s) | 0.760 … 2.957; 1.681 … 18.936 | identical (kept on `DECLARED_WAVE`) |
| generator / tolerance text | "Hulthen (the declared verdict wave …)" | "AV18 (the shipped wave since D6 …; the headline under A5) … Hulthen (the declaration's deuteron_channel(), recorded) 627.85/187" |
| figure | Hulthén in the `default` slot | AV18 in the `default` slot, Hulthén `alt`; title "AV18 (shipped) LC 579.12/187" |

Documents, in place (line-neutral): `validation/benchmarks/README.md` row 11 and
`docs/benchmarking/BENCHMARK_PLAN.md` row 11 now lead with 579.12/187 for AV18
"the headline since A5 (2026-11-03)", the Hulthén 627.85 recorded.
`docs/benchmarking/02_data_nucleon_deuteron.md` quotes no Deeps headline (its
Deeps lines are about obtainability) and is not edited; there is no
`docs/benchmarking/README.md`.

### 1.2 A5 — every harness against the rule, and the four rows the supervisor's note decided

The shipped default read from `lipolgen.make_config()` on this tree:
`cluster_wave = VmcAV18`, `unpol_sf = Ct18Nlo`, `pol_sf = NnpdfPol`,
`b1_model = Miller`, `fsi = Off` (σ_XN 40 mb when on), `rc = Off`, R =
`r_sigma_lt` (`r_source` unset; the deuteron kernel's f1 equals
F₂/(2x(1 + r_sigma_lt)) at all six HERMES bins); `EMC_BASELINE_DEFAULT = Epps21`.

| row | headline before | after (the rule) | status | action |
|---|---|---|---|---|
| T1 ×4, `t1_r_world` | shipped (D7) | the same | unchanged | none |
| `t2_deeps_spectator_tail` | Hulthén 627.85/187 | AV18 579.12/187 | FAIL → FAIL | applied (§1.1) |
| `t2_bonus_spectator_shape` | Hulthén 40 mb: 4/12, 0/16 | AV18 40 mb: 3/12, 0/16 | FAIL → FAIL | applied (§1.2.1) |
| `t2_deuteron_static` | all three waves (FAIL, Hulthén +17.2σ) | AV18 η −1.39σ | **FAIL → PASS, by the rule** | applied (§1.2.2) |
| `t2_hermes_b1_table2` | `cdks_conv_mstw` 22.26/6 | `miller_shipped` 5.26/6 | FAIL → FAIL | applied (§1.2.3) |
| `t2_boeglin_nk` | Hulthén 21.98/21 (hashed block) | unchanged: the named exception; AV18 11.54/21 (p = 0.951) is the sub-row `nk[av18]`, already computed in `run()` outside the block | PASS → PASS | dated A5 addendum in the docstring outside the block; no number moved |
| T3 ×5 | shipped `HoSpin1FF`; `Epps21Ratio` = the shipped EMC baseline | the same | unchanged | none (recorded: the inclusive kernels' own F₂(⁶Li)/F₂(d) is exactly 1; on it NMC gives 7.82/4, p = 0.099, and 24.86/15, p = 0.052 — PASS too; `a5_survey.json`) |
| T4 ×2, T5 ×5, T6 | shipped / identities / CLI defaults | the same | unchanged | none |

Before the note, the survey measured BONuS on the literal shipped run (`--fsi
off`: weight 1, 12/12 and 16/16, a trivial PASS) — the reason for the note's
refinement.

#### 1.2.1 BONuS (`pre_t2_bonus_spectator_shape.log` → `post_…`)

| quantity | before | after |
|---|---|---|
| `in_verdict` | `{4,5}GeV/hulthen-40mb` | `{4,5}GeV/av18-40mb` (`DEFAULT_CONFIG = "av18-40mb"`, `DECLARED_CONFIG = "hulthen-40mb"`; CONFIGS reordered, AV18 first) |
| headline blocks in band | 4/12, 0/16 | **3/12, 0/16** |
| worst v/δ | 2.62 (4 GeV), 3.44 (5 GeV) | 2.81, 3.56 (same blocks) |
| tree-shape χ² (84 dof) | 226.3 | 221.6 (χ²_flat 739.4 both) |
| headline v by p_s | 0.0344 / 0.0548 / 0.0848 / 0.1406 | 0.0355 / 0.0555 / 0.0871 / 0.1505 |
| Hulthén 40 mb, 20 mb sub-rows | 4/12, 0/16; 8/12, 13/16 | identical, recorded |
| stat ⊕ syst | Hulthén 8/12, 8/16 (19.02, 9.13) | AV18 **7/12, 8/16** (19.31, 8.64); Hulthén's kept as `… hulthen-40mb stat+syst`, identical |
| forward rise ρ_tree by p_s | Hulthén 0.9902 / 0.9838 / 0.9741 / 0.9552 | AV18 0.9904 / 0.9835 / 0.9733 / 0.9521 (same sign 12/28 both); Hulthén's kept as `forward rise [hulthen-40mb, declared]` |
| PWIA prescription max dev / ends / weight ends | Hulthén 0.0227…0.0571 / 1.0457…1.1177 / 1.0710…1.3103 | AV18 0.0151 / 0.0484 / 0.0598 / 0.0559; 1.0298 / 1.0977 / 1.1211 / 1.1118; 1.0735 / 1.1161 / 1.1851 / 1.3376; Hulthén's kept, identical |
| clipped grid fraction | 0.01451 (Hulthén) | 0.01346 (AV18); Hulthén's kept as `fsi_clipped_grid_fraction_declared` |
| status | FAIL | FAIL |
| figure | Hulthén in the default slot | AV18 in the default slot, title on AV18 |

Finding, recorded (the harness docstring's A5 NOTE; `test_bench_bonus.py`):
the header's two CHECKS were made on the Hulthén weight and still hold there;
on the AV18 headline weight the (k_z, k_T)-table convergence |dv| is 3.1e−4 /
1.1e−4 / 5.5e−4 / 8.3e−4 (threshold 5e−4) and bin averaging moves v by −3.2 /
−1.5 / +0.8 / −0.3 % (threshold 2 %). Neither threshold was moved; the tests
keep them on the declared weight and pin the AV18 values as recorded. Applying
either deviation to every one of the 28 headline blocks flips no verdict.

#### 1.2.2 Deuteron statics

| quantity | before | after |
|---|---|---|
| headline | the conjunction over CD-Bonn, AV18, Hulthén | AV18 (`HEADLINE_MODEL`, checked by `_headline()` against `CLUSTER_WAVE_DEFAULT`), η 0.025045, −1.39σ |
| status | FAIL | **PASS** — by the rule, no η, Q_d, tolerance or reference moved |
| sub-rows | the three waves | the three waves (`in_verdict` on AV18 only) + `all waves (the declared rule, recorded)`, FAIL, failing "hulthen +17.2 sigma" |
| REPORT line suffix | "FAIL (recorded: hulthen +17.2 sigma outside; nothing moved)" | "PASS on the shipped wave (A5); the declared all-waves rule FAIL (hulthen eta +17.2 sigma), recorded" |
| figure | Hulthén in the default slot | AV18 in the default slot, "(shipped; headline)" |

#### 1.2.3 HERMES b₁ (row 3)

| quantity | before | after |
|---|---|---|
| headline | `cdks_conv_mstw` (or `cdks_conv_default` without MSTW) | `miller_shipped` (`HEADLINE`), MSTW-independent |
| headline χ²(b₁)/6, bins | 22.264, 2/6 | **5.262, 4/6** |
| `miller_shipped` A_zz χ²/6, bins | 5.550, 4/6 (MSTW F₂, R1990) | **7.308, 3/6** (shipped F₁ᵈ: F₂ from `shipped_inclusive_kernel(deuteron())`, R = `r_sigma_lt`, HERMES's (1 + γ²)); 5.550, 4/6 kept as `chi2_azz_mstw_r1990`; "with R1998" 5.410 unchanged (MSTW) |
| shipped F₁ᵈ per bin (x = 0.012 … 0.452) | — | 12.402 / 4.843 / 2.582 / 1.298 / 0.5638 / 0.1451 (MSTW/R1990: 7.527 / 3.593 / 2.010 / 1.082 / 0.5318 / 0.1655; HERMES-implied 7.044 / 3.427 / 1.929 / 1.018 / 0.4957 / 0.1614) — the low-x bins sit below CT18NLO's Q² grid, frozen as shipped |
| every other configuration's χ² | — | identical (`cdks_conv_mstw` 22.264 / 22.009) |
| status | FAIL | FAIL |

Finding, recorded (the harness's A5 note): the shipped INCLUSIVE deuteron
kernel's own b₁ is 0.32753 × `MillerB1()` (= (2/6) × `LI6_B1_RANK2_TRANSFER`,
the ⁶Li α–d mapping, also with `default_inclusive_kernel`); on it χ²(b₁)/6 =
14.56, 2/6. `miller_shipped` stays the `MillerB1()` the harness declares as
`--b1-model miller`. Not acted on (outside A5; for the supervisor, §5).

### 1.3 A6 — the TIGHT edge, every model (`t5_ab.txt`)

`saturation_tight` before (the committed REPORT.json, = 2|M₂|/B_tight) → after
(|M₂|/B_tight), at Q² = 2 / 5 / 10 GeV². The ratio is 2.0000000000 ± 2.4e−10
for every entry (the report stores ten digits); `saturation_tight_v6` (new)
equals the before column to that precision.

| model | TIGHT before | TIGHT after | status (central / V6) |
|---|---|---|---|
| `toy_1e-2[shipped]` (verdict) | 0.3206737036 / 0.3808336503 / 0.428129935 | **0.1603368518 / 0.1904168252 / 0.2140649675** | PASS (pass / pass), unchanged |
| `moment_A[letter]` | 3.015591746 / 2.932961894 / 2.938096266 | **1.507795873 / 1.466480947 / 1.469048133** | recorded (pass / fail), unchanged |
| `moment_B[letter]` | 3.114028983 / 2.918750585 / 2.822219793 | **1.557014491 / 1.459375293 / 1.411109896** | recorded (pass / fail), unchanged |
| `moment_A[dilution 1]` | 9.046775237 / 8.798885683 / 8.814288797 | **4.523387619 / 4.399442841 / 4.407144398** | recorded (fail / fail), unchanged |
| `moment_B[dilution 1]` | 9.342086949 / 8.756251756 / 8.466659378 | **4.671043474 / 4.378125878 / 4.233329689** | recorded (fail / fail), unchanged |

Unchanged (checked against the committed report): every `saturation`
(toy 0.0759 / 0.0906 / 0.1022), every `saturation_v6` (toy 0.1518 / 0.1812 /
0.2044), `B`, `B_tight` (9.108080e−04 / 7.781218e−04 / 7.013481e−04),
`saturation_loose`, `f1ll_critical(_v6)`, the positivity sub-row, the row
status PASS. The TIGHT edge carries no status, so no verdict could move.

**The fact the A6 entry needs** (for the supervisor's amendment): DECISIONS A6
reads "since D4 that factor is 1 … the edge's value unchanged, the factor being
1". In the harness the TIGHT edge did not carry V6 through `DS16_OVER_TREE`
(= 1 since D4) but through `PRE_D4_OVER_JM = 2.0` (the pre-D4 kernel / the
letter as printed, ×2 on M₂: `saturation_tight = PRE_D4_OVER_JM × |M₂| /
B_tight`). Removing V6 therefore halves every TIGHT value (table above); the
statuses do not move. `validation/benchmarks/data/delta_moment_bound_sources.json`
(sha256-pinned) is not edited: its V6 note stays a record of the sources as
read.

Documents: `validation/benchmarks/README.md` row 13 ("tight edge 0.160 / 0.190 /
0.214 without V6 since A6 … 0.321 / 0.381 / 0.428 with it, recorded") and
`docs/benchmarking/BENCHMARK_PLAN.md` row 13 ("the tight edge kept it until
A6"), in place. In the harness, each line-neutral above the D4 block's cited
lines (`DS16_OVER_TREE` and `PRE_D4_OVER_JM` keep their line numbers): S3's
text, the `B_tight` comment, `saturation_row` (`saturation_tight`,
`saturation_tight_v6`), the bound-side print, the figure's docstring and
whisker caption, the D4 block's comment and `PRE_D4_OVER_JM`'s comment; a dated
A6 / A7 note after `V6_LABEL`. The figure's whiskers now end at the halved
TIGHT edge (content diff, §2).

### 1.4 A7 — the grid row's bound

| Q² [GeV²] | M₂ | error estimate \|fine − coarse\| | relative | bound QUAD_REL_TOL × \|M₂\| | bound / error |
|---|---|---|---|---|---|
| 2 | 1.4604e−04 | 3.26e−14 | 2.23e−10 | 1.46e−11 | 448 |
| 5 | 1.4817e−04 | 4.90e−14 | 3.31e−10 | 1.48e−11 | 302 |
| 10 | 1.5013e−04 | 5.53e−14 | 3.69e−10 | 1.50e−11 | 271 |

New field `m2_error_bound` on every saturation sub-row (the letter's rows, on
the analytic toy F₁, have errors 0–8.7e−19 against bounds 1e−10). The pytest
asserts the grid row's errors within the bound; the 1e−15 absolute pin on the
toy's analytic integrand (`_toy_m2_errors_below(mod, 1e-15)`) stays. The
harness's `integrate()` has raised (exit 2) above the same rule since its
first pass, so the row could not have passed above it; A7 states it as the
row's bound. **The bound was chosen after the values were seen** (D7 measured
3.3e−14 / 4.9e−14 / 5.5e−14; F7, `../run_2026-10-27/phase_I7_backends.md` §13
item 2, asked that this be said); QUAD_REL_TOL itself (1e−7) pre-dates D7 and
was not moved. Status PASS → PASS.

### 1.5 Test pins re-taken (each with a dated `A5 (2026-11-03)` / `A6` / `A7` comment)

| file | pin | before → after |
|---|---|---|
| `test_bench_deeps.py` | the `in_verdict` sub-row; the verdict test's head | `hulthen/lc` → `av18/lc`; `waves["hulthen"]` → `waves["av18"]`; new: the two wave constants, `hulthen/lc` FAIL recorded, AV18's data/model per p_s, the generator's lead |
| `test_bench_bonus.py` | `in_verdict`, heads, `names[:2]`, sub-row count | `hulthen-40mb` → `av18-40mb`; 12 → 16 sub-rows; new pins: AV18 v by p_s, worst v/δ 2.807 / 3.560 and blocks, the three passing AV18 blocks, AV18 stat ⊕ syst (7/12 19.31, 8/16 8.64), AV18 forward ρ_tree and same-sign 12, AV18 prescription; the forward / prescription pins moved to the `[hulthen-40mb, declared]` keys (values unchanged); the two CHECK tests on `DECLARED_CONFIG` (thresholds unchanged) with the AV18 values pinned as recorded |
| `test_bench_static.py` | status, REPORT line | `status` = the all-waves rule → `declared_status`; new: `status == "pass"` on AV18, `in_verdict` = AV18 only, the declared-rule sub-row FAIL; the REPORT line's lead and suffix |
| `test_bench_spin_deuteron.py` | `HERMES_CHI2["miller_shipped"]` A_zz; the no-MSTW headline | 5.550 → 7.308; `cdks_conv_default` → `miller_shipped`; new: headline, 4/6 and 3/6, `chi2_azz_mstw_r1990` 5.550, the kernel's f1 = F₂/(2x(1 + r_sigma_lt)) |
| `test_bench_t5_delta_bound.py` | toy and `moment_A[letter]` `saturation_tight`; the grid row's error | 0.3206737036 … → 0.1603368518 …; 3.015591746 … → 1.507795873 …; new `saturation_tight_v6`; `m2_error_bound` and every error within it (A7) |

The t5 test file's edits are line-neutral (a first, line-adding version broke
six run-record citations). `t2_hermes_b1_table2.py`'s lines 1–160 are
byte-identical to `df56035` (a run record cites its line 160; the A5 note and
`HEADLINE` sit after `build_models()`).

## 2. What was verified and how (A/B)

- PRE / POST runs of every edited harness, diffed (`pre_*.log`, `post_*.log`):
  only the verdict flags, the headline lines and the sub-rows named in §1
  changed; every recorded number of the declared configurations is identical.
- `t5_ab.txt`: TIGHT before / after = 2 to ten digits for all 15 entries;
  `saturation`, `saturation_v6` and the statuses equal the committed report.
- `make_report.py --all --plots --out-dir docs/benchmarking --compare <the
  df56035 REPORT.json copy>`: 15 pass / 8 fail / 0 blocked / 0 error, **6
  changes**: MOVED Deeps, BONuS, deuteron statics and HERMES headlines; STATUS
  `t2_deuteron_static` FAIL → PASS; SUBMOVED `t2_hermes_b1_table2[miller_shipped]`
  (A_zz 5.54993 → 7.30798, within 4 → 3, the two `_mstw_r1990` fields new);
  five new sub-rows reported informational (BONuS AV18 stat ⊕ syst ×2, the
  declared forward and prescription; the deuteron statics' declared rule).
  These are exactly the moves the note named; nothing else moved.
- Figures: 19 regenerated; 14 differ only in the footer stamp (differing rows
  from ≥ 97.2 % of the height, `pixdiff2.txt`) and were restored; five changed
  in content and were looked at: `t2_deeps_spectator_tail.png`,
  `t2_bonus_spectator_shape.png`, `t2_deuteron_static.png`,
  `t2_hermes_b1_table2.png` (headline in the title), `t5_delta_moment_bound.png`.
- The survey (`a5_survey.py`) and the BONuS flip check call the harnesses' own
  functions read-only.

## 3. Files

- `validation/benchmarks/README.md` — convention 8 (one line, line-neutral, with
  the note's refinement and the one exception); rows 3, 10, 11, 13, T2 §5.2, 16
  in place.
- `validation/benchmarks/t2_deeps_spectator_tail.py`, `t2_bonus_spectator_shape.py`,
  `t2_deuteron_static.py`, `t2_hermes_b1_table2.py` — headlines re-pointed,
  dated A5 notes, labels, figures; `t2_boeglin_nk.py` — a dated A5 addendum
  outside the hashed block only (the block's sha256 check passes).
- `validation/benchmarks/t5_delta_moment_bound.py` — A6 / A7 (§1.3–1.4).
- `python/tests/test_bench_deeps.py`, `test_bench_bonus.py`, `test_bench_static.py`,
  `test_bench_spin_deuteron.py`, `test_bench_t5_delta_bound.py` — §1.5.
- `docs/benchmarking/BENCHMARK_PLAN.md` (rows 3, 10, 11, 13, 16, the T2 status
  paragraph), `02_data_nucleon_deuteron.md` (§5.2 and §5.4 rows),
  `08_channel_priority.md` (the "validated today" paragraph), `README.md` (the
  tally line: 15 / 8 since A5) — in place, line-neutral.
- `docs/benchmarking/REPORT.md`, `REPORT.json` and five figures — regenerated.
- `docs/open_items/run_2026-11-03/phase_P3_report_rules.md` (this record),
  `STATUS.md` (the P3 row only).

Not edited: `DECISIONS.md`, the hashed design blocks, `validation/reference/*.json`,
`delta_moment_bound_sources.json`, any C++ file.

## 4. Gates (final tree), tallies and wall times

| gate | result | wall |
|---|---|---|
| C++ build / `./build/lipolgen_tests` | not run: no C++ file changed | — |
| stage 1, flip-only bench subset (before the report was regenerated) | 1 failed / 367 passed / 1 skipped: `test_the_committed_report_is_what_regeneration_gives` (the committed report not yet regenerated; cleared by regeneration, no pin touched) | 7 m 43 s |
| stage 1, `--compare` / `--check` / full pytest / bench / citation / SPDX / doc-link | 14 / 9, 1 change; OK; 1517 / 155 / 0; 368 / 1 / 0; 0 broken; 136; 112 passed | superseded by the rows below |
| `make_report.py --all --plots --out-dir docs/benchmarking --compare <df56035 copy>` | **15 / 8 / 0 blocked / 0 error; 6 changes** (§2) | 2 m 4 s |
| `make_report.py --check` | OK | 1 m 44 s |
| pytest, full | 1517 passed, 155 skipped, 0 failed | 20 m 30 s |
| bench subset (`LIPOLGEN_BENCH=1 -k bench`) | 368 passed, 1 skipped, 0 failed | 6 m 56 s |
| strict citation gate | 1243 / 95 / 7 / **0 broken** / 24 | 3 s |
| `--records-only` | 1432 citations in 85 records: 922 live, 61 into another record, 121 annotated, 129 read by hand, 199 skipped, **0 broken** | <1 s |
| `--audit-ranges` | 209 range citations, 0 REFUSED | 4 s |
| `check_spdx_headers.py` | 136 files, all carry the line | <1 s |
| `pytest python/tests/test_doc_link_gate.py` | 112 passed | 38 s |

Intermediate `--records-only` runs found broken run-record citations three
times (README insertion: 3; line-adding t5 pins: 6; the HERMES A5 note above
its line 160: 1); each edit was redone line-neutrally; none was re-pinned or
annotated. The per-file tests of each edited harness were run after each edit
(Deeps 16, BONuS 26, statics 18, spin-deuteron 24, Boeglin 13, t5 23 passed).

## 5. Open points

1. **The shipped inclusive deuteron kernel's b₁** is 0.32753 × `MillerB1()`
   ((2/6) × `LI6_B1_RANK2_TRANSFER`, the ⁶Li α–d mapping) — measured here, not
   investigated. If it is not intended, it is a physics item outside A5; if
   it is, HERMES's headline would have to be that kernel's b₁ (14.56/6, 2/6,
   FAIL). Recommendation: the supervisor reads `default_inclusive_kernel`'s b₁
   wiring for the deuteron species before the next phase.
2. **BONuS's CHECKS on the AV18 weight** miss their declared thresholds
   (convergence 8.3e−4 > 5e−4 at p_s = 0.135, bin averaging −3.2 % > 2 % at
   0.078); recorded, thresholds not moved, no verdict flips. Recommendation:
   leave recorded; a finer default table is a library change, not a P3 one.
3. **Deeps' Glauber and transverse-window sub-rows** stay on the declared
   Hulthén as records (the supervisor's note, item 5); no edit.
4. **HERMES's low-x bins** sit below CT18NLO's Q² grid, so the shipped F₁ᵈ
   there is the frozen edge (12.40 against HERMES's implied 7.04 at x = 0.012);
   recorded in §1.2.3.
